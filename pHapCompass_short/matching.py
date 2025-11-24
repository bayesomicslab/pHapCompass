import numpy as np
import itertools
from collections import defaultdict
from typing import Dict, List, Tuple, Set
from scipy.sparse import csr_matrix
import pandas as pd
from utils import *
import random


def select_best_candidate(candidates, prioritize="probabilities"):
    """
    Select the best candidate based on the chosen priority (counts or probabilities).
    
    Parameters:
        candidates (list): List of tuples (candidate_array, probability, count).
        prioritize (str): Primary criterion to prioritize ("counts" or "probabilities").
    
    Returns:
        tuple: The best candidate (candidate_array, probability, count).
        candidates = [('a', 0.9, 5), ('b',0.9, 5), ('c',0.9, 4), ('d',0.8, 10), ('e',0.9, 10), ('f',0.9, 10)]
    """
    if prioritize == "counts":
        # Sort by count first, then probability
        candidates.sort(key=lambda x: (x[2], x[1]), reverse=True)
    elif prioritize == "probabilities":
        # Sort by probability first, then count
        candidates.sort(key=lambda x: (x[1], x[2]), reverse=True)
    else:
        raise ValueError("Invalid priority. Choose 'counts' or 'probabilities'.")

    # Get the best criteria based on the priority
    best_criteria = candidates[0][1 if prioritize == "probabilities" else 2]
    
    # Filter candidates with the best primary criterion
    filtered_candidates = [c for c in candidates if c[1 if prioritize == "probabilities" else 2] == best_criteria]
    
    # If there are still ties, use the secondary criterion to filter
    if len(filtered_candidates) > 1:
        secondary_criteria = filtered_candidates[0][2 if prioritize == "probabilities" else 1]
        filtered_candidates = [c for c in filtered_candidates if c[2 if prioritize == "probabilities" else 1] == secondary_criteria]
    
    # If there are still ties, pick one randomly
    return random.choice(filtered_candidates)

# ---------- canonicalization & small helpers ----------

def _canonicalize_rows(mat: np.ndarray) -> np.ndarray:
    order = np.lexsort([mat[:, c] for c in range(mat.shape[1]-1, -1, -1)])
    return mat[order, :]

def _rows_equal_up_to_permutation(A: np.ndarray, B: np.ndarray) -> bool:
    """
    Return True iff A and B (K x m) are equal up to a row permutation.
    Uses the same canonicalization scheme as _canonicalize_rows.
    """
    if A.shape != B.shape:
        return False
    return np.array_equal(_canonicalize_rows(A), _canonicalize_rows(B))



def _canonical_key(mat: np.ndarray) -> Tuple[int, ...]:
    return tuple(_canonicalize_rows(mat).ravel().tolist())

def _positions_from_edge(edge: Tuple[str,str]) -> List[int]:
    u, v = edge
    pu = list(map(int, u.split('-')))
    pv = list(map(int, v.split('-')))
    return sorted(set(pu + pv))

def _col_index_map(all_positions: List[int]) -> Dict[int,int]:
    return {pos: idx for idx, pos in enumerate(all_positions)}

def _phas_str_to_mat(phas_str: str, K: int) -> np.ndarray:
    a = np.fromiter((ord(c)-48 for c in phas_str), dtype=np.int8)
    return a.reshape(K, -1)

def _permute_rows_groupwise(H: np.ndarray,
                            fixed_cols: List[int],
                            fixed_target: np.ndarray) -> List[np.ndarray]:
    """Permute rows only within groups sharing the same signature on fixed columns."""
    K, m = H.shape
    if not fixed_cols:
        return [H]
    sig_H = [tuple(H[r, fixed_cols].tolist()) for r in range(K)]
    sig_T = [tuple(fixed_target[r, :].tolist()) for r in range(K)]
    grp_H = defaultdict(list);  grp_T = defaultdict(list)
    for r,s in enumerate(sig_H): grp_H[s].append(r)
    for r,s in enumerate(sig_T): grp_T[s].append(r)
    if set(grp_H.keys()) != set(grp_T.keys()):
        return []
    blocks = []
    for s in grp_H.keys():
        H_rows, T_rows = grp_H[s], grp_T[s]
        if len(H_rows)!=len(T_rows): return []
        blocks.append((H_rows, T_rows))
    perms = [[]]
    for H_rows, T_rows in blocks:
        local = []
        for p in itertools.permutations(H_rows):
            local.append(list(zip(T_rows, p)))
        perms = [x+y for x in perms for y in local]
    out=[]
    for mapping in perms:
        P = np.empty(K, dtype=int)
        for dst, src in mapping:
            P[dst]=src
        out.append(H[P,:])
    return out

# ---------- vectorized likelihood over M without per-read loop ----------

def _submatrix_and_groups(M: csr_matrix,
                          snp_to_col: Dict[int,int],
                          pos_list_1based: List[int],
                          min_cov: int = 1):
    """
    Extract a read×T submatrix for the given positions (T=len(pos_list)),
    turn it into obs bits (0/1) and a coverage mask, group rows by the
    coverage mask (which subset of columns are covered). Returns:
      groups: dict{mask_tuple (bool T) -> np.ndarray rows}
      obs   : np.ndarray (R x T) in {0,1}; rows aligned to concatenation of groups
      idxs  : np.ndarray of row indices (length R) aligned to 'obs'
      cols  : np.ndarray of 0-based columns in M (length T)
    Only keeps reads with coverage >= min_cov on this set.
    """
    T = len(pos_list_1based)
    cols = np.array([snp_to_col[p-1] for p in pos_list_1based], dtype=np.int64)
    sub = M[:, cols]
    # coverage per row on those T columns
    cov = (sub > 0).astype(np.int8)
    cov_counts = np.asarray(cov.sum(axis=1)).ravel()
    keep_rows = np.where(cov_counts >= min_cov)[0]
    if keep_rows.size == 0:
        return {}, np.empty((0,T), dtype=np.int8), keep_rows, cols
    subk = sub[keep_rows, :].toarray()  # (R x T)
    obs = (subk == 2).astype(np.int8)   # ALT->1, REF->0; no-call retained in mask
    mask = (subk > 0)                   # (R x T) boolean
    # group rows by mask pattern
    groups = defaultdict(list)
    for r in range(keep_rows.size):
        groups[tuple(mask[r, :].tolist())].append(r)
    # convert lists to arrays
    for k in list(groups.keys()):
        groups[k] = np.array(groups[k], dtype=np.int64)
    return groups, obs, keep_rows, cols

def _score_mixture_by_groups(obs: np.ndarray,
                             groups: Dict[Tuple[bool,...], np.ndarray],
                             H: np.ndarray,
                             error_rate: float) -> float:
    """
    Sum of mixture likelihoods across all rows in 'obs', grouped by
    coverage mask to avoid per-read loops. H is (K x T).
    """
    if obs.size == 0:
        return 0.0
    total = 0.0
    T = H.shape[1]
    for mask_tuple, ridx in groups.items():
        mask = np.array(mask_tuple, dtype=bool)
        if not mask.any():
            continue
        # restrict to covered columns for these rows
        obs_sub = obs[ridx][:, mask]         # (Rg x Tm)
        H_sub   = H[:, mask]                  # (K  x Tm)
        # mixture likelihood over Tm columns
        Rg = obs_sub.shape[0]
        # vectorized mixture
        # expand: (Rg,1,Tm) vs (1,K,Tm)
        match = (obs_sub[:,None,:] == H_sub[None,:,:])
        m  = match.sum(axis=2, dtype=np.int32)
        Tm = H_sub.shape[1]
        mm = Tm - m
        lk = ((1.0 - error_rate) ** m) * (error_rate ** mm)  # (Rg x K)
        L  = lk.sum(axis=1)                                  # (Rg,)
        total += float(L.sum())
    return total


def predict_haplotypes(
    nodes: List[str],
    edges: List[Tuple[str,str]],
    samples: Dict[int, Dict[str,str]],         # {t: {node: phase_str}}
    ploidy: int,
    genotype_path: str,
    M: csr_matrix,
    snp_to_col: Dict[int,int],
    transitions_dict_extra: Dict[str, dict],
    config,
    priority: str = "probabilities"
) -> pd.DataFrame:
    """
    Phase ONE variant per iteration (Algorithm 1), but score candidates
    from M in vectorized form (no per-read loops).
    """
    # flatten samples once
    phasing_samples = {n: samples[t][n] for t in samples for n in samples[t]}
    sorted_nodes = sort_nodes(nodes)

    genotype_df = pd.read_csv(genotype_path).T
    predicted = pd.DataFrame(index=[f"haplotype_{p+1}" for p in range(ploidy)],
                             columns=genotype_df.columns, dtype=float)

    phased_nodes: Set[str] = set()
    unphased_nodes = set(nodes)
    phased_pos: Set[int]   = set()
    unphased_pos = {int(x) for s in nodes for x in s.split('-')}

    # seed from first node (unchanged)
    first_node = sorted_nodes[0]
    init_pos = sorted(map(int, first_node.split('-')))
    predicted.loc[:, [p-1 for p in init_pos]] = str_2_phas_1(phasing_samples[first_node], ploidy)
    phased_nodes.add(first_node); unphased_nodes.remove(first_node)
    phased_pos.update(init_pos); unphased_pos -= set(init_pos)

    while unphased_pos:
        # neighbor nodes touching phased set
        neighbor_nodes = (
            {v for u,v in edges if u in phased_nodes and v in unphased_nodes} |
            {u for u,v in edges if v in phased_nodes and u in unphased_nodes}
        )

        # count connections per unphased pos
        pos_conn = defaultdict(lambda: {"count": 0, "edges": []})
        for u, v in edges:
            if (u in phased_nodes and v in neighbor_nodes) or (v in phased_nodes and u in neighbor_nodes):
                for pos in map(int, u.split('-') + v.split('-')):
                    if pos in unphased_pos:
                        pos_conn[pos]["count"] += 1
                        pos_conn[pos]["edges"].append((u, v))

        if not pos_conn:
            # fallback: take next unphased node whole
            for nxt in sorted_nodes:
                if nxt in unphased_nodes:
                    p2 = sorted(map(int, nxt.split('-')))
                    predicted.loc[:, [p-1 for p in p2]] = str_2_phas_1(phasing_samples[nxt], ploidy)
                    phased_nodes.add(nxt); unphased_nodes.discard(nxt)
                    phased_pos.update(p2); unphased_pos -= set(p2)
                    break
            continue

        selected_pos = max(pos_conn, key=lambda p: pos_conn[p]["count"])
        relevant_edges = pos_conn[selected_pos]["edges"]

        combos, probs, counts, all_positions = _compute_singlevar_candidates_vectorized(
            selected_pos,
            relevant_edges,
            phased_pos,
            transitions_dict_extra,
            phasing_samples,
            ploidy,
            predicted,
            M,
            snp_to_col,
            config.error_rate
        )

        if combos:
            best = select_best_candidate(list(zip(combos, probs, counts)), prioritize=priority)
            best_mat, _, _ = best
            sel_idx = all_positions.index(selected_pos)
            predicted.loc[:, selected_pos - 1] = best_mat[:, sel_idx]
        else:
            # fallback to sampled column from target node
            for (u, v) in relevant_edges:
                for node in (u, v):
                    if str(selected_pos) in node.split('-'):
                        tpos = list(map(int, node.split('-')))
                        tph = str_2_phas_1(phasing_samples[node], ploidy)
                        col = tpos.index(selected_pos)
                        predicted.loc[:, selected_pos - 1] = tph[:, col]
                        break

        # update sets
        for (u, v) in relevant_edges:
            if u in neighbor_nodes:
                phased_nodes.add(u); unphased_nodes.discard(u)
            if v in neighbor_nodes:
                phased_nodes.add(v); unphased_nodes.discard(v)
        phased_pos.add(selected_pos); unphased_pos.discard(selected_pos)

    return predicted

def _compute_singlevar_candidates_vectorized(
    selected_position: int,
    relevant_edges: List[Tuple[str,str]],
    phased_positions: Set[int],
    transitions_dict_extra: Dict[str, dict],
    samples_brief: Dict[str, str],
    ploidy: int,
    predicted_haplotypes: pd.DataFrame,
    M: csr_matrix,
    snp_to_col: Dict[int,int],
    error_rate: float):
    all_positions = sorted(set(int(p) for e in relevant_edges for n in e for p in n.split('-')))
    pos2col = _col_index_map(all_positions)

    combos: List[np.ndarray] = []
    probs:  List[float] = []
    counts: List[int] = []
    seen: Set[Tuple[int,...]] = set()

    for (src, tgt) in relevant_edges:
        edge_key = f"{src}--{tgt}"
        if edge_key not in transitions_dict_extra:
            continue
        # find (source,target) entry matching sampled phasings
        matched = None
        for key, val in transitions_dict_extra[edge_key].items():
            if val["source_phasing"] == samples_brief[src] and val["target_phasing"] == samples_brief[tgt]:
                matched = val["matched_phasings"]
                break
        if not matched:
            continue

        edge_positions = _positions_from_edge((src, tgt))
        fixed_positions = [pos for pos in edge_positions if pos != selected_position and pos in phased_positions]
        fixed_cols_local = [edge_positions.index(pos) for pos in fixed_positions]
        fixed_vals = predicted_haplotypes.loc[:, [p-1 for p in fixed_positions]].values  # (K x |fixed|)

        # pull reads & group by mask (>=1 covered pos in edge)
        groups, obs, ridx, _ = _submatrix_and_groups(M, snp_to_col, edge_positions, min_cov=1)

        # score each merged triple phasing
        for triple_str in matched.keys():
            H = str_2_phas_1(triple_str, ploidy)   # (K x 3)
            # only row permutations consistent with fixed columns
            permlist = _permute_rows_groupwise(H, fixed_cols_local, fixed_vals) if fixed_cols_local else [H]
            if not permlist:
                continue
            # target node constraint
            tgt_pos = list(map(int, tgt.split('-')))
            tgt_cols = [edge_positions.index(p) for p in tgt_pos]
            tgt_mat = str_2_phas_1(samples_brief[tgt], ploidy)
            tgt_perms = np.array(list(itertools.permutations(tgt_mat)))

            for Hp in permlist:
                # require target matches some permutation
                if not any(np.array_equal(Hp[:, tgt_cols], perm) for perm in tgt_perms):
                    continue

                # total likelihood over reads (vectorized, group by mask)
                total = _score_mixture_by_groups(obs, groups, Hp, error_rate)

                # align to all_positions
                aligned = np.full((ploidy, len(all_positions)), np.nan, dtype=float)
                for pos in all_positions:
                    if pos in phased_positions:
                        aligned[:, pos2col[pos]] = predicted_haplotypes.loc[:, pos-1].values
                for j, pos in enumerate(edge_positions):
                    aligned[:, pos2col[pos]] = Hp[:, j]

                key = _canonical_key(aligned)
                if key in seen:
                    continue
                seen.add(key)
                combos.append(aligned)
                probs.append(total)

    # count FFBS-respected target nodes (unchanged)
    for cand in combos:
        count = 0
        for (_, tgt) in relevant_edges:
            tgt_pos = list(map(int, tgt.split('-')))
            idx = [pos2col[p] for p in tgt_pos]
            tgt_mat = str_2_phas_1(samples_brief[tgt], ploidy)
            tgt_perms = np.array(list(itertools.permutations(tgt_mat)))
            if any(np.array_equal(cand[:, idx], perm) for perm in tgt_perms):
                count += 1
        counts.append(count)

    return combos, probs, counts, all_positions

def predict_haplotypes_multiple_variants(
    nodes: List[str],
    edges: List[Tuple[str,str]],
    samples: Dict[int, Dict[str,str]],
    ploidy: int,
    genotype_path: str,
    M: csr_matrix,
    snp_to_col: Dict[int,int],
    transitions_dict_extra: Dict[str, dict],
    config,
    priority: str = "probabilities"
) -> pd.DataFrame:
    """
    Phase MULTIPLE variants per iteration (your Algorithm 2), but score candidates
    from M in vectorized form (no per-read loops).
    """
    samples_brief = {n: samples[t][n] for t in samples for n in samples[t]}
    sorted_nodes = sort_nodes(nodes)

    genotype_df = pd.read_csv(genotype_path).T
    predicted = pd.DataFrame(index=[f"haplotype_{p+1}" for p in range(ploidy)],
                             columns=genotype_df.columns, dtype=float)

    phased_nodes: Set[str] = set()
    unphased_nodes = set(nodes)
    phased_pos: Set[int]   = set()
    unphased_pos = {int(x) for s in nodes for x in s.split('-')}

    # seed
    if sorted_nodes:
        first = sorted_nodes[0]
        p2 = sorted(map(int, first.split('-')))
        predicted.loc[:, [p-1 for p in p2]] = str_2_phas_1(samples_brief[first], ploidy)
        phased_nodes.add(first); unphased_nodes.discard(first)
        phased_pos.update(p2); unphased_pos -= set(p2)

    while unphased_pos:
        neighbor_nodes = (
            {v for u,v in edges if u in phased_nodes and v in unphased_nodes} |
            {u for u,v in edges if v in phased_nodes and u in unphased_nodes}
        )
        pos_conn = defaultdict(lambda: {"count": 0, "edges": []})
        for u, v in edges:
            if (u in phased_nodes and v in neighbor_nodes) or (v in phased_nodes and u in neighbor_nodes):
                for pos in map(int, u.split('-') + v.split('-')):
                    if pos in unphased_pos:
                        pos_conn[pos]["count"] += 1
                        pos_conn[pos]["edges"].append((u, v))

        if not pos_conn:
            # take next node whole
            nxt = next((n for n in sorted_nodes if n in unphased_nodes), None)
            if nxt:
                p2 = sorted(map(int, nxt.split('-')))
                predicted.loc[:, [p-1 for p in p2]] = str_2_phas_1(samples_brief[nxt], ploidy)
                phased_nodes.add(nxt); unphased_nodes.discard(nxt)
                phased_pos.update(p2); unphased_pos -= set(p2)
            continue

        selected_positions = set(pos_conn.keys())
        relevant = list({e for pos in selected_positions for e in pos_conn[pos]["edges"]})

        # combos, probs, counts, all_positions = _compute_multivar_candidates_vectorized(
        #     selected_positions, relevant, phased_pos, samples_brief,
        #     ploidy, predicted, M, snp_to_col, config.error_rate
        # )

        combos, probs, counts, all_positions = _compute_multivar_candidates_vectorized_beam(
        selected_positions, relevant, phased_pos, samples_brief,
        ploidy, predicted, M, snp_to_col, config.error_rate,
        beam_size=64,              # try 64–128 for K=8
        max_nodes_per_round=6      # try 6–10; start small
        )

        if combos:
            best = select_best_candidate(list(zip(combos, probs, counts)), prioritize=priority)
            best_mat, _, _ = best
            for pos in selected_positions:
                predicted.loc[:, pos - 1] = best_mat[:, all_positions.index(pos)]
        else:
            # fallback: write sampled target-node columns
            assigned = set()
            for (u, v) in relevant:
                for node in (u, v):
                    tpos = list(map(int, node.split('-')))
                    tph = str_2_phas_1(samples_brief[node], ploidy)
                    for pos in selected_positions:
                        if pos in tpos and pos not in assigned:
                            predicted.loc[:, pos - 1] = tph[:, tpos.index(pos)]
                            assigned.add(pos)

        phased_pos.update(selected_positions); unphased_pos -= selected_positions
        for (u, v) in relevant:
            if u in neighbor_nodes: phased_nodes.add(u); unphased_nodes.discard(u)
            if v in neighbor_nodes: phased_nodes.add(v); unphased_nodes.discard(v)

    return predicted

def _compute_multivar_candidates_vectorized(
    selected_positions: Set[int],
    relevant_edges: List[Tuple[str,str]],
    phased_positions: Set[int],
    samples_brief: Dict[str, str],
    ploidy: int,
    predicted_haplotypes: pd.DataFrame,
    M: csr_matrix,
    snp_to_col: Dict[int,int],
    error_rate: float
):
    all_positions = sorted(set(int(p) for e in relevant_edges for n in e for p in n.split('-')))
    pos2col = _col_index_map(all_positions)

    # target nodes = nodes containing any selected position
    target_nodes = sorted({
        n for e in relevant_edges for n in e
        if any(int(p) in selected_positions for p in n.split('-'))
    })
    target_pos = sorted({int(p) for n in target_nodes for p in n.split('-')})

    # build block candidates K x len(target_pos) by stitching FFBS phasings with groupwise permutations
    # start from first
    blocks = [np.full((ploidy, len(target_pos)), np.nan, dtype=float)]
    seen_blocks: Set[Tuple[int,...]] = set()

    first = target_nodes[0]
    c_first = [target_pos.index(int(p)) for p in first.split('-')]
    H_first = str_2_phas_1(samples_brief[first], ploidy)
    b0 = blocks[0].copy(); b0[:, c_first] = H_first
    blocks = [b0]; seen_blocks.add(_canonical_key(b0))

    # expand with remaining target nodes
    for node in target_nodes[1:]:
        cols_local = [target_pos.index(int(p)) for p in node.split('-')]
        H_node    = str_2_phas_1(samples_brief[node], ploidy)
        new_list = []
        for blk in blocks:
            fixed_cols = [c for c in cols_local if not np.isnan(blk[:, c]).all()]
            fixed_vals = blk[:, fixed_cols] if fixed_cols else None
            if fixed_cols:
                perms = _permute_rows_groupwise(H_node, [cols_local.index(c) for c in fixed_cols], fixed_vals)
            else:
                perms = [H_node]
            for P in perms:
                upd = blk.copy()
                upd[:, cols_local] = P
                k = _canonical_key(upd)
                if k not in seen_blocks:
                    seen_blocks.add(k); new_list.append(upd)
        if new_list:
            blocks = new_list

    # align to all_positions and dedup
    combos: List[np.ndarray] = []
    seen_full: Set[Tuple[int,...]] = set()
    for blk in blocks:
        full = np.full((ploidy, len(all_positions)), np.nan, dtype=float)
        for pos in all_positions:
            if pos in phased_positions:
                full[:, pos2col[pos]] = predicted_haplotypes.loc[:, pos-1].values
        for j, pos in enumerate(target_pos):
            full[:, pos2col[pos]] = blk[:, j]
        k = _canonical_key(full)
        if k not in seen_full:
            seen_full.add(k); combos.append(full)

    # score candidates: group rows by mask on all positions and sum mixture likelihoods
    groups, obs, ridx, _ = _submatrix_and_groups(M, snp_to_col, all_positions, min_cov=1)
    probs: List[float] = []
    for cand in combos:
        probs.append(_score_mixture_by_groups(obs, groups, cand, error_rate))

    # FFBS-respect counts
    counts: List[int] = []
    for cand in combos:
        cnt = 0
        for (u, v) in relevant_edges:
            for node in (u, v):
                if node in target_nodes:
                    npos = list(map(int, node.split('-')))
                    idx  = [pos2col[p] for p in npos]
                    tgt  = str_2_phas_1(samples_brief[node], ploidy)
                    tgt_perms = np.array(list(itertools.permutations(tgt)))
                    if any(np.array_equal(cand[:, idx], perm) for perm in tgt_perms):
                        cnt += 1
        counts.append(cnt)

    return combos, probs, counts, all_positions

def _score_groups_with_mask(
    obs_all: np.ndarray,
    groups_all: Dict[Tuple[bool,...], np.ndarray],
    H_full: np.ndarray,               # (K x T_all) over all_positions (target_pos)
    col_mask: np.ndarray,             # (T_all,) boolean: which columns are currently filled
    error_rate: float
) -> float:
    """
    Sum mixture likelihoods over all reads using only columns where col_mask=True.
    Reuses the grouping computed once for all positions.
    """
    if obs_all.size == 0 or not col_mask.any():
        return 0.0
    total = 0.0
    for mask_tuple, ridx in groups_all.items():
        gm = np.array(mask_tuple, dtype=bool)
        use = gm & col_mask
        if not use.any() or ridx.size == 0:
            continue
        obs_sub = obs_all[ridx][:, use]
        H_sub   = H_full[:,   use]
        match = (obs_sub[:, None, :] == H_sub[None, :, :])
        m  = match.sum(axis=2, dtype=np.int32)
        Tm = H_sub.shape[1]
        mm = Tm - m
        lk = ((1.0 - error_rate) ** m) * (error_rate ** mm)
        L  = lk.sum(axis=1)
        total += float(L.sum())
    return total

def _topk(items, scores, k):
    """Keep top-k (items) by scores."""
    if len(items) <= k:
        return items, scores
    idx = np.argpartition(scores, -k)[-k:]
    # sort those top k
    order = idx[np.argsort(scores[idx])[::-1]]
    return [items[i] for i in order], [scores[i] for i in order]


def _compute_multivar_candidates_vectorized_beam(
    selected_positions: Set[int],
    relevant_edges: List[Tuple[str,str]],
    phased_positions: Set[int],
    samples_brief: Dict[str, str],
    ploidy: int,
    predicted_haplotypes: pd.DataFrame,
    M: csr_matrix,
    snp_to_col: Dict[int,int],
    error_rate: float,
    beam_size: int = 64,
    max_nodes_per_round: int = 6
):
    """
    Build K×|selected| candidates with a beam and incremental scoring.
    Returns (combos, final_scores, ffbs_counts, all_positions).
    """
    # ----- choose target nodes (limit per round to keep search feasible) -----
    # nodes touching any selected position, prefer those with higher degree in relevant_edges
    deg = defaultdict(int)
    for u,v in relevant_edges:
        deg[u]+=1; deg[v]+=1
    target_nodes_all = sorted({
        n for e in relevant_edges for n in e
        if any(int(p) in selected_positions for p in n.split('-'))
    }, key=lambda n: -deg[n])
    target_nodes = target_nodes_all[:max_nodes_per_round]

    target_pos = sorted({int(p) for n in target_nodes for p in n.split('-')})
    pos2idx    = {p:i for i,p in enumerate(target_pos)}
    T_all      = len(target_pos)

    # pre-group reads for ALL target_pos once
    groups_all, obs_all, _, _ = _submatrix_and_groups(M, snp_to_col, target_pos, min_cov=1)

    # seed block: place first node without permuting rows (break symmetry)
    first = target_nodes[0]
    c_first = [pos2idx[int(p)] for p in first.split('-')]
    H_first = str_2_phas_1(samples_brief[first], ploidy)

    init_block = np.full((ploidy, T_all), np.nan, dtype=float)
    init_block[:, c_first] = H_first
    col_mask = np.zeros(T_all, dtype=bool); col_mask[c_first] = True

    # initial partial score from columns we just placed
    init_score = _score_groups_with_mask(obs_all, groups_all, init_block, col_mask, error_rate)

    # beam holds (block, score, mask)
    beam_blocks = [init_block]
    beam_scores = [init_score]
    beam_masks  = [col_mask]

    # expand node by node
    for node in target_nodes[1:]:
        cols_local = [pos2idx[int(p)] for p in node.split('-')]
        H_node     = str_2_phas_1(samples_brief[node], ploidy)

        new_blocks = []
        new_scores = []
        new_masks  = []

        for blk, sc, mask in zip(beam_blocks, beam_scores, beam_masks):
            # fixed columns in this node
            fixed_cols = [c for c in cols_local if mask[c]]
            fixed_vals = blk[:, fixed_cols] if fixed_cols else None

            # permute H_node rows only within fixed-col groups
            if fixed_cols:
                # indices in node's 2-cols corresponding to fixed_cols
                local_fixed_idx = [cols_local.index(c) for c in fixed_cols]
                perms = _permute_rows_groupwise(H_node, local_fixed_idx, fixed_vals)
            else:
                perms = [H_node]

            for P in perms:
                upd = blk.copy()
                upd[:, cols_local] = P
                m2 = mask.copy(); m2[cols_local] = True
                # incremental score on new columns only: full mask m2 vs old mask
                inc = _score_groups_with_mask(obs_all, groups_all, upd, m2, error_rate) \
                      - _score_groups_with_mask(obs_all, groups_all, blk, mask, error_rate)
                new_blocks.append(upd)
                new_scores.append(sc + inc)
                new_masks.append(m2)

        # keep top-k
        beam_blocks, beam_scores = _topk(new_blocks, np.array(new_scores), beam_size)
        beam_masks  = [new_masks[i] for i in range(len(new_masks)) if any(np.all(beam_blocks[j]==new_blocks[i]) for j in range(len(beam_blocks)))]

    # align each beam block into all_positions==target_pos (already aligned),
    # canonicalize, dedup, and compute final full score (optional; we already have scores)
    seen = set()
    combos = []
    final_scores = []
    for blk, sc in zip(beam_blocks, beam_scores):
        key = _canonical_key(blk)
        if key in seen: 
            continue
        seen.add(key)
        combos.append(blk)
        final_scores.append(sc)

    # FFBS respect counts (unchanged)
    counts = []
    for cand in combos:
        cnt = 0
        for (u, v) in relevant_edges:
            for node in (u, v):
                if node in target_nodes:
                    npos = list(map(int, node.split('-')))
                    idx  = [pos2idx[p] for p in npos]
                    tgt  = str_2_phas_1(samples_brief[node], ploidy)
                    tgt_perms = np.array(list(itertools.permutations(tgt)))
                    if any(np.array_equal(cand[:, idx], perm) for perm in tgt_perms):
                        cnt += 1
        counts.append(cnt)

    return combos, final_scores, counts, target_pos


