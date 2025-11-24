"""
Final MEC-Based Haplotype Phasing with Correct Scoring
All fixes applied:
1. No weight normalization (use raw weights)
2. Correct likelihood function (matching user's specification)
3. Proper FFBS counting for target nodes
4. Verbose output showing candidate matrices and scoring
"""

from typing import List, Tuple, Dict, Set, Optional
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix
import itertools


def predict_haplotypes_mec_based(
    nodes: List[str],
    edges: List[Tuple[str, str]],
    samples: Dict[int, Dict[str, str]],
    ploidy: int,
    genotype_path: str,
    M: csr_matrix,
    snp_to_col: Dict[int, int],
    transitions_dict_extra: Dict[str, dict],
    config,
    priority: str = "combined",
    ffbs_weight: float = 0.5,
    likelihood_weight: float = 0.3,
    mec_weight: float = 1.0,
    allow_ffbs_override: bool = True,
    verbose: bool = False
) -> Tuple[pd.DataFrame, np.ndarray]:
    """Predict haplotypes using MEC-based candidate selection."""
    phasing_samples = {nn: samples[t][nn] for t in samples for nn in samples[t]}
    sorted_nodes = sort_nodes(nodes)
    
    genotype_df = pd.read_csv(genotype_path).T
    n_positions = len(genotype_df.columns)
    all_positions = {int(pos) for node in nodes for pos in node.split('-')}
    
    predicted_haplotypes = pd.DataFrame(
        index=[f'haplotype_{p+1}' for p in range(ploidy)],
        columns=genotype_df.columns,
        dtype=float
    )

    # Block tracking
    block_ids = np.full(n_positions, np.nan)
    position_to_block = {}
    current_block_id = 0

    phased_nodes = set()
    phased_positions = set()
    unphased_nodes = set(nodes)
    unphased_positions = {int(pos) for node in nodes for pos in node.split('-')}

    # Initialize
    first_node = sorted_nodes[0]
    initial_positions = sorted({int(pos) for pos in first_node.split('-')})
    initial_positions_idx = [p - 1 for p in initial_positions]
    predicted_haplotypes.loc[:, initial_positions_idx] = str_2_phas_1(
        phasing_samples[first_node], ploidy
    )
    
    for pos in initial_positions:
        block_ids[pos - 1] = current_block_id
        position_to_block[pos] = current_block_id
    
    phased_nodes.add(first_node)
    unphased_nodes.remove(first_node)
    phased_positions.update(initial_positions)
    unphased_positions -= set(initial_positions)

    iteration = 0
    fallback_count = 0
    
    if verbose:
        print("\n" + "="*80)
        print(f"STARTING MEC-BASED PREDICTION (priority={priority})")
        print("="*80)

    while unphased_positions:
        iteration += 1
        
        neighbor_nodes = (
            {n2 for n1, n2 in edges if n1 in phased_nodes and n2 in unphased_nodes} |
            {n1 for n1, n2 in edges if n2 in phased_nodes and n1 in unphased_nodes}
        )

        position_connections = defaultdict(lambda: {"count": 0, "edges": []})
        for n1, n2 in edges:
            if (n1 in phased_nodes and n2 in neighbor_nodes) or \
               (n2 in phased_nodes and n1 in neighbor_nodes):
                for pos in map(int, n1.split('-') + n2.split('-')):
                    if pos in unphased_positions:
                        position_connections[pos]["count"] += 1
                        position_connections[pos]["edges"].append((n1, n2))

        # Handle graph disconnection
        if not position_connections:
            for next_node in sorted_nodes:
                if next_node in unphased_nodes:
                    next_positions = sorted({int(pos) for pos in next_node.split('-')})
                    next_positions_idx = [p - 1 for p in next_positions]
                    next_sample_np = str_2_phas_1(phasing_samples[next_node], ploidy)
                    predicted_haplotypes.loc[:, next_positions_idx] = next_sample_np
                    
                    already_phased_in_node = [p for p in next_positions if p in position_to_block]
                    if already_phased_in_node:
                        inherited_block = position_to_block[already_phased_in_node[0]]
                        for pos in next_positions:
                            if pos not in position_to_block:
                                block_ids[pos - 1] = inherited_block
                                position_to_block[pos] = inherited_block
                    else:
                        current_block_id += 1
                        for pos in next_positions:
                            block_ids[pos - 1] = current_block_id
                            position_to_block[pos] = current_block_id
                    
                    phased_nodes.add(next_node)
                    unphased_nodes.remove(next_node)
                    phased_positions.update(next_positions)
                    unphased_positions -= set(next_positions)
                    break
            continue

        selected_position = max(position_connections, key=lambda p: position_connections[p]["count"])
        relevant_edges = position_connections[selected_position]["edges"]
        
        if verbose:
            print(f"\n{'='*80}")
            print(f"ITERATION {iteration}: Position {selected_position}")

        # Generate candidates
        candidates_info = compute_candidate_phasings_all_transitions(
            selected_position,
            relevant_edges,
            phased_positions,
            transitions_dict_extra,
            phasing_samples,
            ploidy,
            predicted_haplotypes,
            M,
            snp_to_col,
            genotype_df,
            config.error_rate,
            verbose=verbose
        )
        
        if not candidates_info:
            fallback_count += 1
            if verbose:
                print(f"⚠️  FALLBACK (total: {fallback_count})")
            
            # Block assignment
            inherited_block = None
            for node in nodes:
                node_positions = {int(p) for p in node.split('-')}
                if selected_position in node_positions:
                    phased_in_node = node_positions & phased_positions
                    if phased_in_node:
                        for pp in phased_in_node:
                            if pp in position_to_block:
                                inherited_block = position_to_block[pp]
                                break
                        if inherited_block is not None:
                            break
            
            if inherited_block is not None:
                block_ids[selected_position - 1] = inherited_block
                position_to_block[selected_position] = inherited_block
            else:
                current_block_id += 1
                block_ids[selected_position - 1] = current_block_id
                position_to_block[selected_position] = current_block_id
            
            # FFBS fallback
            phased_this = False
            for node in nodes:
                if str(selected_position) in node.split('-') and node in phasing_samples:
                    node_positions = [int(pos) for pos in node.split('-')]
                    node_phasing = str_2_phas_1(phasing_samples[node], ploidy)
                    node_position_index = node_positions.index(selected_position)
                    predicted_haplotypes.loc[:, selected_position - 1] = node_phasing[:, node_position_index]
                    phased_this = True
                    break
            
            if not phased_this:
                g_value = int(genotype_df.iloc[:, selected_position - 1].sum())
                random_phasing = np.zeros(ploidy, dtype=np.int8)
                if 0 < g_value < ploidy:
                    alt_indices = np.random.choice(ploidy, size=g_value, replace=False)
                    random_phasing[alt_indices] = 1
                elif g_value == ploidy:
                    random_phasing = np.ones(ploidy, dtype=np.int8)
                predicted_haplotypes.loc[:, selected_position - 1] = random_phasing
        else:
            # Block assignment
            connected_block = None
            for phased_pos in phased_positions:
                for src, tgt in relevant_edges:
                    src_positions = {int(p) for p in src.split('-')}
                    tgt_positions = {int(p) for p in tgt.split('-')}
                    edge_positions = src_positions | tgt_positions
                    
                    if selected_position in edge_positions and phased_pos in edge_positions:
                        if phased_pos in position_to_block:
                            connected_block = position_to_block[phased_pos]
                            break
                if connected_block is not None:
                    break
            
            if connected_block is not None:
                block_ids[selected_position - 1] = connected_block
                position_to_block[selected_position] = connected_block
            else:
                current_block_id += 1
                block_ids[selected_position - 1] = current_block_id
                position_to_block[selected_position] = current_block_id
            
            # Score candidates
            scored_candidates = score_candidates_with_mec_fixed(
                candidates_info, M, snp_to_col,
                priority=priority, ffbs_weight=ffbs_weight,
                likelihood_weight=likelihood_weight, mec_weight=mec_weight,
                verbose=verbose
            )
            
            # Select best
            if priority == "mec":
                best = min(scored_candidates, key=lambda x: x['score'])
            else:
                best = max(scored_candidates, key=lambda x: x['score'])
            
            if verbose:
                print(f"Best: MEC={best['mec']:.2f}, L={best['likelihood']:.2f}, F={best['ffbs_matches']}, Score={best['score']:.2f}")
            
            # Update
            selected_position_index = best['all_positions'].index(selected_position)
            predicted_haplotypes.loc[:, selected_position - 1] = best['candidate'][:, selected_position_index]

        # Update sets
        for node1, node2 in relevant_edges:
            if node1 in neighbor_nodes:
                phased_nodes.add(node1)
                unphased_nodes.discard(node1)
            if node2 in neighbor_nodes:
                phased_nodes.add(node2)
                unphased_nodes.discard(node2)

        phased_positions.add(selected_position)
        unphased_positions.remove(selected_position)

    # Renumber blocks
    unique_blocks = np.unique(block_ids[~np.isnan(block_ids)])
    if len(unique_blocks) > 0:
        block_mapping = {old_id: new_id for new_id, old_id in enumerate(sorted(unique_blocks))}
        for i in range(len(block_ids)):
            if not np.isnan(block_ids[i]):
                block_ids[i] = block_mapping[block_ids[i]]

    if verbose:
        n_blocks = len(unique_blocks)
        n_unphased = np.sum(np.isnan(block_ids))
        print(f"\n{'='*80}")
        print(f"COMPLETE: {iteration} iterations, {fallback_count} fallbacks")
        print(f"Blocks: {n_blocks}, Unphased: {n_unphased}")
        print("="*80)

    return predicted_haplotypes, block_ids


# ============================================================================
# CANDIDATE GENERATION
# ============================================================================

def compute_candidate_phasings_all_transitions(
    selected_position: int,
    relevant_edges: List[Tuple[str, str]],
    phased_positions: Set[int],
    transitions_dict_extra: Dict[str, dict],
    phasing_samples: Dict[str, str],
    ploidy: int,
    predicted_haplotypes: pd.DataFrame,
    M: csr_matrix,
    snp_to_col: Dict[int, int],
    genotype_df: pd.DataFrame,
    error_rate: float,
    verbose: bool = False
) -> List[Dict]:
    """
    Generate candidates from ALL matched_phasings in transitions_dict_extra.
    Count FFBS matches for all target nodes being extended.
    """
    all_positions = sorted(set(
        int(p) for edge in relevant_edges 
        for node in edge 
        for p in node.split('-')
    ))
    
    candidates_info = []
    seen_canonical_keys = set()
    
    # Identify target nodes (nodes being extended from phased set)
    target_nodes_info = {}
    for src, tgt in relevant_edges:
        # Determine which node is the target (being extended)
        src_phased = all(int(p) in phased_positions for p in src.split('-'))
        tgt_phased = all(int(p) in phased_positions for p in tgt.split('-'))
        
        if src_phased and not tgt_phased:
            target_node = tgt
        elif tgt_phased and not src_phased:
            target_node = src
        elif not src_phased and not tgt_phased:
            # Both partially phased - consider both as targets
            for node in [src, tgt]:
                if node in phasing_samples:
                    node_positions = [int(p) for p in node.split('-')]
                    if all(p in all_positions for p in node_positions):
                        node_cols = [all_positions.index(p) for p in node_positions]
                        node_sample = str_2_phas_1(phasing_samples[node], ploidy)
                        target_nodes_info[node] = {
                            'positions': node_positions,
                            'cols': node_cols,
                            'sample': node_sample
                        }
            continue
        else:
            continue
        
        if target_node in phasing_samples:
            tgt_positions = [int(p) for p in target_node.split('-')]
            if all(p in all_positions for p in tgt_positions):
                tgt_cols = [all_positions.index(p) for p in tgt_positions]
                tgt_sample = str_2_phas_1(phasing_samples[target_node], ploidy)
                target_nodes_info[target_node] = {
                    'positions': tgt_positions,
                    'cols': tgt_cols,
                    'sample': tgt_sample
                }
    
    for src, tgt in relevant_edges:
        edge_key = f"{src}--{tgt}"
        if edge_key not in transitions_dict_extra:
            continue
        
        edge_positions = sorted(set(int(p) for node in [src, tgt] for p in node.split('-')))
        fixed_in_edge = [p for p in edge_positions if p in phased_positions]
        unfixed_in_edge = [p for p in edge_positions if p not in phased_positions]
        
        if selected_position not in unfixed_in_edge:
            continue
        
        # Collect ALL matched_phasings from ALL entries
        all_triples = set()
        for entry_val in transitions_dict_extra[edge_key].values():
            matched_phasings = entry_val.get("matched_phasings", {})
            if not matched_phasings:
                continue
            all_triples.update(matched_phasings.keys())
        
        # Generate candidates from all triples
        for triple_str in all_triples:
            triple_matrix = str_2_phas_1(triple_str, ploidy)
            
            if fixed_in_edge:
                fixed_cols_in_edge = [edge_positions.index(p) for p in fixed_in_edge]
                fixed_vals_expected = predicted_haplotypes.loc[:, [p-1 for p in fixed_in_edge]].values
                valid_perms = _permute_rows_groupwise(triple_matrix, fixed_cols_in_edge, fixed_vals_expected)
            else:
                valid_perms = [triple_matrix]
            
            if not valid_perms:
                continue
            
            for perm_triple in valid_perms:
                candidate = np.full((ploidy, len(all_positions)), -1, dtype=float)
                
                for p in phased_positions:
                    if p in all_positions:
                        col = all_positions.index(p)
                        candidate[:, col] = predicted_haplotypes.loc[:, p - 1].values
                
                for p in unfixed_in_edge:
                    edge_col = edge_positions.index(p)
                    all_col = all_positions.index(p)
                    candidate[:, all_col] = perm_triple[:, edge_col]
                
                if not _verify_genotype_constraint(candidate, all_positions, genotype_df):
                    continue
                
                candidate_key = _canonical_key(candidate)
                if candidate_key in seen_canonical_keys:
                    continue
                seen_canonical_keys.add(candidate_key)
                
                # Count FFBS matches (row-permutation aware)
                ffbs_matches = 0
                for tgt_node, tgt_info in target_nodes_info.items():
                    candidate_tgt = candidate[:, tgt_info['cols']]
                    if _matches_any_permutation(candidate_tgt, tgt_info['sample']):
                        ffbs_matches += 1
                
                # Compute likelihood using user's function
                likelihood = compute_likelihood_from_reads(
                    candidate, all_positions, M, snp_to_col, error_rate
                )
                
                candidates_info.append({
                    'candidate': candidate,
                    'likelihood': likelihood,
                    'ffbs_matches': ffbs_matches,
                    'respects_ffbs': False,
                    'all_positions': all_positions
                })
    
    if verbose:
        print(f"  Generated {len(candidates_info)} unique candidates")
        if len(target_nodes_info) > 0:
            print(f"  Checking FFBS agreement with {len(target_nodes_info)} target nodes: {list(target_nodes_info.keys())}")
    
    return candidates_info


# ============================================================================
# LIKELIHOOD COMPUTATION (User's specification)
# ============================================================================

def compute_likelihood_from_reads(
    phasing: np.ndarray,
    all_positions: List[int],
    M: csr_matrix,
    snp_to_col: Dict[int, int],
    error_rate: float
) -> float:
    """
    Compute likelihood using user's specification:
    For each read covering these positions, compute Σ_k P(read | haplotype_k)
    Sum over all reads.
    """
    K, T = phasing.shape
    
    # Get reads
    cols = [snp_to_col[p - 1] for p in all_positions]
    sub_M = M[:, cols].toarray()
    
    coverage = (sub_M > 0).sum(axis=1)
    covering_reads = np.where(coverage > 0)[0]
    
    if len(covering_reads) == 0:
        return 0.0
    
    observed = (sub_M[covering_reads] == 2).astype(np.int8)
    mask = (sub_M[covering_reads] > 0)
    
    total_likelihood = 0.0
    
    # For each read
    for r in range(len(covering_reads)):
        read = observed[r, :]
        read_mask = mask[r, :]
        
        if not read_mask.any():
            continue
        
        # Use user's formula: Σ_k [(1-ε)^matches × ε^mismatches]
        read_covered = read[read_mask]
        phasing_covered = phasing[:, read_mask]
        
        # For each haplotype k
        prob_sum = 0.0
        for k in range(K):
            hap = phasing_covered[k, :]
            matches = (read_covered == hap).sum()
            mismatches = len(read_covered) - matches
            prob_k = ((1.0 - error_rate) ** matches) * (error_rate ** mismatches)
            prob_sum += prob_k
        
        total_likelihood += prob_sum
    
    return total_likelihood


# ============================================================================
# SCORING
# ============================================================================

def score_candidates_with_mec_fixed(
    candidates_info: List[Dict],
    M: csr_matrix,
    snp_to_col: Dict[int, int],
    priority: str = "combined",
    ffbs_weight: float = 0.5,
    likelihood_weight: float = 0.3,
    mec_weight: float = 1.0,
    verbose: bool = False
) -> List[Dict]:
    """
    Score candidates using RAW weights (NO normalization of weights).
    """
    if not candidates_info:
        return []
    
    # Compute MEC for all candidates
    for cand_dict in candidates_info:
        candidate = cand_dict['candidate']
        all_positions = cand_dict['all_positions']
        mec_score = compute_mec_for_candidate(candidate, all_positions, M, snp_to_col)
        cand_dict['mec'] = mec_score
    
    # Extract metrics
    likelihoods = np.array([c['likelihood'] for c in candidates_info])
    mecs = np.array([c['mec'] for c in candidates_info])
    ffbs_matches = np.array([c['ffbs_matches'] for c in candidates_info])
    
    lik_min, lik_max = likelihoods.min(), likelihoods.max()
    mec_min, mec_max = mecs.min(), mecs.max()
    ffbs_min, ffbs_max = ffbs_matches.min(), ffbs_matches.max()
    
    if verbose:
        print(f"\n  Weights: L={likelihood_weight}, M={mec_weight}, F={ffbs_weight}")
        print(f"  Metric ranges:")
        print(f"    Likelihood: [{lik_min:.2f}, {lik_max:.2f}]")
        print(f"    MEC: [{mec_min:.2f}, {mec_max:.2f}]")
        print(f"    FFBS matches: [{ffbs_min}, {ffbs_max}]")
    
    # Score each candidate
    for i, cand_dict in enumerate(candidates_info):
        # Normalize metrics to [0, 1]
        lik_norm = 0.5 if abs(lik_max - lik_min) < 1e-10 else (cand_dict['likelihood'] - lik_min) / (lik_max - lik_min)
        mec_norm = 0.5 if abs(mec_max - mec_min) < 1e-10 else 1.0 - (cand_dict['mec'] - mec_min) / (mec_max - mec_min)
        ffbs_norm = 0.5 if abs(ffbs_max - ffbs_min) < 1e-10 else (cand_dict['ffbs_matches'] - ffbs_min) / (ffbs_max - ffbs_min)
        
        # WEIGHTED SUM with RAW weights (no normalization of weights!)
        if priority == "mec":
            score = -cand_dict['mec']
        elif priority == "likelihood":
            score = cand_dict['likelihood']
        else:  # combined
            score = (
                likelihood_weight * lik_norm + 
                mec_weight * mec_norm + 
                ffbs_weight * ffbs_norm
            )
        
        cand_dict['score'] = score
        cand_dict['lik_norm'] = lik_norm
        cand_dict['mec_norm'] = mec_norm
        cand_dict['ffbs_norm'] = ffbs_norm
        
        if verbose and i < 5:
            print(f"\n  Cand {i}:")
            print(f"    Matrix:\n{cand_dict['candidate']}")
            print(f"    L={cand_dict['likelihood']:.2f}(norm={lik_norm:.2f}), "
                  f"M={cand_dict['mec']:.1f}(norm={mec_norm:.2f}), "
                  f"F={cand_dict['ffbs_matches']}(norm={ffbs_norm:.2f})")
            print(f"    Score = {likelihood_weight}*{lik_norm:.2f} + {mec_weight}*{mec_norm:.2f} + {ffbs_weight}*{ffbs_norm:.2f} = {score:.3f}")
    
    return candidates_info


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def compute_mec_for_candidate(candidate, all_positions, M, snp_to_col):
    """Compute MEC score."""
    K, n_pos = candidate.shape
    cols = [snp_to_col[p - 1] for p in all_positions]
    sub_M = M[:, cols].toarray()
    coverage = (sub_M > 0).sum(axis=1)
    covering_reads = np.where(coverage > 0)[0]
    
    if len(covering_reads) == 0:
        return 0.0
    
    obs = (sub_M[covering_reads] == 2).astype(np.int8)
    coverage_mask = (sub_M[covering_reads] > 0)
    mec_total = 0.0
    
    for read_idx in range(len(covering_reads)):
        obs_read = obs[read_idx]
        mask = coverage_mask[read_idx]
        if not mask.any():
            continue
        
        mismatches_per_haplotype = []
        for k in range(K):
            hap = candidate[k]
            mismatches = ((obs_read != hap) & mask).sum()
            mismatches_per_haplotype.append(mismatches)
        mec_total += min(mismatches_per_haplotype)
    
    return float(mec_total)


def _verify_genotype_constraint(candidate, all_positions, genotype_df):
    K, T = candidate.shape
    for col_idx, pos in enumerate(all_positions):
        col_values = candidate[:, col_idx]
        if np.any(col_values < 0) or np.any(np.isnan(col_values)):
            continue
        candidate_sum = int(col_values.sum())
        expected_sum = int(genotype_df.iloc[:, pos - 1].sum())
        if candidate_sum != expected_sum:
            return False
    return True


def _canonical_key(matrix):
    sorted_mat = _sort_rows_lexicographically(matrix)
    return tuple(sorted_mat.ravel())


def _sort_rows_lexicographically(matrix):
    row_tuples = [tuple(row) for row in matrix]
    sorted_indices = sorted(range(len(row_tuples)), key=lambda i: row_tuples[i])
    return matrix[sorted_indices, :]


def _matches_any_permutation(mat1, mat2):
    if mat1.shape != mat2.shape:
        return False
    canonical1 = _sort_rows_lexicographically(mat1)
    canonical2 = _sort_rows_lexicographically(mat2)
    return np.allclose(canonical1, canonical2, equal_nan=True)


def _permute_rows_groupwise(H, fixed_cols, fixed_target):
    K, m = H.shape
    if not fixed_cols:
        return [H]
    
    sig_H = [tuple(H[r, fixed_cols].tolist()) for r in range(K)]
    sig_T = [tuple(fixed_target[r, :].tolist()) for r in range(K)]
    grp_H = defaultdict(list)
    grp_T = defaultdict(list)
    for r, s in enumerate(sig_H):
        grp_H[s].append(r)
    for r, s in enumerate(sig_T):
        grp_T[s].append(r)
    
    if set(grp_H.keys()) != set(grp_T.keys()):
        return []
    
    blocks = []
    for s in grp_H.keys():
        H_rows, T_rows = grp_H[s], grp_T[s]
        if len(H_rows) != len(T_rows):
            return []
        blocks.append((H_rows, T_rows))
    
    perms = [[]]
    for H_rows, T_rows in blocks:
        local = []
        for p in itertools.permutations(H_rows):
            local.append(list(zip(T_rows, p)))
        perms = [x + y for x in perms for y in local]
    
    result = []
    for mapping in perms:
        P = np.empty(K, dtype=np.int32)
        for dst, src in mapping:
            P[dst] = src
        result.append(H[P, :])
    return result


def sort_nodes(nodes):
    return sorted(nodes, key=lambda n: int(n.split('-')[0]))


def str_2_phas_1(phasing_str, ploidy):
    arr = np.array([int(c) for c in phasing_str], dtype=np.int8)
    n_positions = len(arr) // ploidy
    return arr.reshape(ploidy, n_positions)