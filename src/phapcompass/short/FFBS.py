from collections import defaultdict, deque, Counter
from typing import Dict, List, Tuple, Optional
import numpy as np
from scipy.sparse import csr_matrix
import math
import random
from itertools import permutations, product
from phapcompass.utils import *


def compute_sampled_path_probability(
    sampled_states: Dict[int, Dict[str, str]],  # Output from sample_states_book
    slices: Dict[int, List[str]],
    edges: List[Tuple[str, str]],
    forward_messages: Dict[int, Dict[str, Dict[str, float]]],
    transitions_dict: Dict[str, np.ndarray]) -> float:
    """
    Compute the probability of a sampled path from FFBS.
    
    The probability is:
    P(path) = P(states at last slice) × ∏ P(parent_state | child_state, forward messages)
    
    Since FFBS samples backward, we compute:
    - At last slice: P(state) ∝ forward[state]
    - At earlier slices: P(state) ∝ forward[state] × ∏ transition[state→child_state]
    """
    
    slice_keys = sorted(slices.keys())
    
    # Build children adjacency
    children = defaultdict(list)
    for u, v in edges:
        children[u].append(v)
    
    log_prob = 0.0
    
    # 1. Probability at last slice (normalize forward messages)
    t_last = slice_keys[-1]
    for node in slices[t_last]:
        sampled_state = sampled_states[t_last][node]
        state_labels = list(forward_messages[t_last][node].keys())
        
        # Get forward messages for all states
        log_alphas = np.array([forward_messages[t_last][node][s] for s in state_labels])
        
        # Normalize (softmax)
        max_log = np.max(log_alphas)
        if np.isfinite(max_log):
            alphas = np.exp(log_alphas - max_log)
            total = alphas.sum()
            if total > 0:
                probs = alphas / total
                sampled_idx = state_labels.index(sampled_state)
                log_prob += np.log(probs[sampled_idx])
    
    # 2. For earlier slices, add backward sampling probabilities
    for t in reversed(slice_keys[:-1]):
        # Find next slice
        t_next = None
        for tt in slice_keys:
            if tt > t:
                t_next = tt
                break
        
        if t_next is None:
            continue
        
        for node in slices[t]:
            sampled_state = sampled_states[t][node]
            state_labels = list(forward_messages[t][node].keys())
            S = len(state_labels)
            
            if S == 0:
                continue
            
            # Get forward messages
            log_alpha = np.array([forward_messages[t][node][s] for s in state_labels])
            
            # Get children in next slice
            child_nodes = [ch for ch in children[node] if ch in slices[t_next]]
            
            # Compute backward probability: forward × ∏ transitions
            add_log = np.zeros(S)
            
            for ch in child_nodes:
                ch_states = list(forward_messages[t_next][ch].keys())
                if len(ch_states) == 0:
                    continue
                
                ch_sampled = sampled_states[t_next][ch]
                try:
                    cj = ch_states.index(ch_sampled)
                except ValueError:
                    continue
                
                # Get transition probabilities for each parent state
                sampled_idx = state_labels.index(sampled_state)
                
                # Get transition matrix
                key = f"{node}--{ch}"
                T = transitions_dict.get(key, None)
                
                if T is not None and T.shape == (len(state_labels), len(ch_states)):
                    # T[sampled_idx, cj] is P(child=cj | parent=sampled_idx)
                    trans_prob = T[sampled_idx, cj]
                    if trans_prob > 0 and np.isfinite(trans_prob):
                        add_log[sampled_idx] += np.log(trans_prob)
                else:
                    # Try reversed key
                    key_rev = f"{ch}--{node}"
                    T_rev = transitions_dict.get(key_rev, None)
                    if T_rev is not None and T_rev.shape == (len(ch_states), len(state_labels)):
                        # Convert P(parent|child) to P(child|parent)
                        col_sum = T_rev.sum(axis=0, keepdims=True).clip(min=1e-12)
                        T_conv = (T_rev / col_sum).T
                        trans_prob = T_conv[sampled_idx, cj]
                        if trans_prob > 0 and np.isfinite(trans_prob):
                            add_log[sampled_idx] += np.log(trans_prob)
            
            # Compute posterior probability (forward × backward)
            log_post = log_alpha + add_log
            
            # Normalize
            max_log = np.max(log_post)
            if np.isfinite(max_log):
                post = np.exp(log_post - max_log)
                total = post.sum()
                if total > 0:
                    probs = post / total
                    sampled_idx = state_labels.index(sampled_state)
                    if probs[sampled_idx] > 0:
                        log_prob += np.log(probs[sampled_idx])
    
    return np.exp(log_prob)


def _rows_covering_all(M, snp_to_col, pos_1based_list):
    """Return (rows, cols) for reads covering ALL 1-based positions."""
    O = (M > 0).astype(np.int8)
    cols = np.array([snp_to_col[p-1] for p in pos_1based_list], dtype=np.int64)
    sub = O[:, cols]
    rows = np.where(np.asarray(sub.getnnz(axis=1)).ravel() == cols.size)[0]
    return rows, cols


def _phas_str_to_mat(phas_str: str, K: int) -> np.ndarray:
    """'001011' -> (K x 2) int array (0/1)."""
    a = np.fromiter((ord(c)-48 for c in phas_str), dtype=np.int8)
    return a.reshape(K, 2)


def _node_label(ci, cj, col_to_snp):
    """0-based columns -> 1-based canonical '{i}-{j}' with i<j."""
    pi = int(col_to_snp[ci]) + 1
    pj = int(col_to_snp[cj]) + 1
    if pi > pj: pi, pj = pj, pi
    return f"{pi}-{pj}"


def _parse_node(lbl: str) -> Tuple[int,int]:
    """'{i}-{j}' -> (i,j) as ints (1-based)."""
    i, j = lbl.split("-")
    i = int(i); j = int(j)
    if i > j: i, j = j, i
    return i, j


def _read_likelihoods_mixture(obs_RxT: np.ndarray, H_KxT: np.ndarray, e: float, divide_by_K: bool=False) -> np.ndarray:
    """
    Per-read mixture likelihood on T loci:
      L_r = sum_h prod_t [(1-e) if obs==H else e].
    """
    R, K = obs_RxT.shape[0], H_KxT.shape[0]
    obs = obs_RxT[:, None, :]  # (R,1,T)
    H   = H_KxT[None, :, :]    # (1,K,T)
    match = (obs == H)
    m  = match.sum(axis=2, dtype=np.int32)  # (R,K)
    T  = H_KxT.shape[1]
    mm = T - m
    lk = ((1.0 - e) ** m) * (e ** mm)      # (R,K)
    L  = lk.sum(axis=1)                    # (R,)
    if divide_by_K and K > 0:
        L = L / K
    return L


def _match_pair_states_to_triples(H_left: np.ndarray, H_right: np.ndarray) -> List[np.ndarray]:
    """
    Given Kx2 H_left = [a,s] and Kx2 H_right = [s,b], enumerate all Kx3 H_[a,s,b]
    by permuting rows of H_right so that the shared column 's' matches:
       H_left[:,1] == H_right_P[:,0].
    """
    K = H_left.shape[0]
    s_left  = H_left[:, 1]
    s_right = H_right[:, 0]
    idx0_L = np.where(s_left  == 0)[0]
    idx1_L = np.where(s_left  == 1)[0]
    idx0_R = np.where(s_right == 0)[0]
    idx1_R = np.where(s_right == 1)[0]
    if len(idx0_L) != len(idx0_R) or len(idx1_L) != len(idx1_R):
        return []
    out = []
    for p0 in permutations(idx0_R):
        for p1 in permutations(idx1_R):
            P = np.empty(K, dtype=np.int64)
            P[idx0_L] = np.array(p0, dtype=np.int64)
            P[idx1_L] = np.array(p1, dtype=np.int64)
            HR = H_right[P, :]
            if not np.array_equal(H_left[:,1], HR[:,0]):
                continue
            H_triple = np.stack([H_left[:,0], H_left[:,1], HR[:,1]], axis=1)  # [a,s,b]
            out.append(H_triple)
    return out


def _node_cols_from_label(node_lbl: str, snp_to_col: Dict[int,int]) -> Tuple[int,int]:
    """'{i}-{j}' (i<j, 1-based) -> (col_i, col_j) 0-based matrix columns."""
    i_str, j_str = node_lbl.split("-")
    i = int(i_str); j = int(j_str)
    if i > j: i, j = j, i
    return snp_to_col[i-1], snp_to_col[j-1]

def _pattern_counts_for_node(
    M: csr_matrix,
    rows: np.ndarray,     # read indices to use
    col_i: int,
    col_j: int
) -> np.ndarray:
    """
    Count patterns among *specific* reads for node columns (col_i, col_j).
    Values in M are {0:no-call, 1:REF, 2:ALT}. We only count rows where both are nonzero.
    Returns counts in order [n00, n01, n10, n11], with 0/1 meaning REF/ALT.
    """
    if rows.size == 0:
        return np.zeros(4, dtype=np.int64)
    # pull sub-columns for those rows
    vi = M[rows, col_i].toarray().ravel()
    vj = M[rows, col_j].toarray().ravel()
    mask = (vi > 0) & (vj > 0)
    if not np.any(mask):
        return np.zeros(4, dtype=np.int64)
    vi = vi[mask]; vj = vj[mask]
    bi = (vi == 2).astype(np.int8)  # ALT->1, REF->0
    bj = (vj == 2).astype(np.int8)
    patt = (bi << 1) | bj           # 0=00, 1=01, 2=10, 3=11
    cnt = np.bincount(patt, minlength=4).astype(np.int64)
    return cnt

def compute_forward_messages_from_M(
    slices: Dict[int, List[str]],
    edges: List[Tuple[str, str]],
    assignment_dict: Dict[str, Dict[str, List[int]]],  # {'states': {node: [read_ids]}, ...}
    emission_dict: Dict[str, Dict[str, Dict[str, float]]],  # node -> phase -> {'00','01','10','11': prob}
    transitions_dict: Dict[str, np.ndarray],           # "u--v" -> (S_u x S_v) P(child|parent)
    M: csr_matrix,
    snp_to_col: Dict[int,int],
    add_uniform_prior: bool = True
) -> Dict[int, Dict[str, Dict[str, float]]]:
    """
    Fast forward messages using M. Exactly reproduces your old semantics:
      alpha_t(node, phase) = sum_reads log P(obs|phase) + logsumexp over parents of alpha_{t-1} + log P(child|parent).

    Inputs:
      - slices: {t: [node_labels]} (layered topological order)
      - edges : [(parent, child), ...]
      - assignment_dict: must contain 'states' mapping node->list[read_idx]
      - emission_dict: node->{phase_str: {'00':p00,'01':p01,'10':p10,'11':p11}}
      - transitions_dict: "parent--child" -> (S_parent x S_child) array
      - M, snp_to_col: sparse matrix + position->column map

    Returns:
      forward_messages[t][node][phase_str] = log alpha.
    """
    # parents adjacency
    parents = defaultdict(list)
    for u, v in edges:
        parents[v].append(u)

    # precompute log emissions per node/phase as a 4-vector (order 00,01,10,11)
    # and the state order (phase names) we'll consistently use
    node_state_order: Dict[str, List[str]] = {}
    node_log_emisvec: Dict[str, Dict[str, np.ndarray]] = {}
    for node, phases in emission_dict.items():
        state_order = list(phases.keys())
        node_state_order[node] = state_order
        node_log_emisvec[node] = {}
        for phase in state_order:
            p = phases[phase]
            # ensure order 00,01,10,11
            log_vec = np.log([
                float(p['00']),
                float(p['01']),
                float(p['10']),
                float(p['11']),
            ])
            node_log_emisvec[node][phase] = log_vec

    # precompute pattern counts per node (from M and assignment_dict['states'][node])
    node_counts4: Dict[str, np.ndarray] = {}
    for node in {n for L in slices.values() for n in L}:
        rows = np.array(assignment_dict['states'].get(node, []), dtype=np.int64)
        ci, cj = _node_cols_from_label(node, snp_to_col)
        cnt4 = _pattern_counts_for_node(M, rows, ci, cj)  # [n00,n01,n10,n11]
        node_counts4[node] = cnt4

    # initialize output structure
    forward_messages: Dict[int, Dict[str, Dict[str, float]]] = {t: {} for t in slices}

    # handle slice order explicitly
    slice_keys = sorted(slices.keys())
    if not slice_keys:
        return forward_messages

    # --- base slice ---
    t0 = slice_keys[0]
    for node in slices[t0]:
        forward_messages[t0][node] = {}
        state_order = node_state_order[node]
        cnt4 = node_counts4[node]
        for phase in state_order:
            log_emit = float(np.dot(cnt4, node_log_emisvec[node][phase]))  # sum over reads via counts
            if add_uniform_prior and len(state_order) > 0:
                log_emit += -math.log(len(state_order))
            forward_messages[t0][node][phase] = log_emit

    # --- recursion ---
    node_to_slice = {n: t for t in slice_keys for n in slices[t]}
    for t in slice_keys[1:]:
        for node in slices[t]:
            forward_messages[t][node] = {}
            state_order = node_state_order[node]
            S_child = len(state_order)
            if S_child == 0:
                continue

            # emission term (same for all parents)
            cnt4 = node_counts4[node]
            log_emit_vec = np.array([np.dot(cnt4, node_log_emisvec[node][ph]) for ph in state_order], dtype=np.float64)
            if add_uniform_prior and S_child > 0:
                log_emit_vec += -math.log(S_child)

            # aggregate over parents
            log_trans_accum = np.full(S_child, -np.inf, dtype=np.float64)
            for pa in parents.get(node, []):
                pt = node_to_slice.get(pa, None)
                if pt is None or pt >= t:
                    continue
                pa_order = node_state_order[pa]
                S_parent = len(pa_order)
                if S_parent == 0:
                    continue

                # alpha_{t-1} for parent, in consistent order
                alpha_pa = np.array([forward_messages[pt][pa][ph] for ph in pa_order], dtype=np.float64)

                # transitions pa->node; if reversed exists, you can add a conversion here if needed
                key = f"{pa}--{node}"
                T = transitions_dict.get(key, None)
                if T is None:
                    # try reversed and convert: P(pa|child) -> P(child|pa) by col-normalize + transpose
                    T_rev = transitions_dict.get(f"{node}--{pa}", None)
                    if T_rev is None:
                        continue
                    col_sum = T_rev.sum(axis=0, keepdims=True).clip(min=1e-12)
                    T = (T_rev / col_sum).T
                if T.shape != (S_parent, S_child):
                    raise ValueError(f"Transition shape mismatch for {key}: {T.shape} vs ({S_parent},{S_child})")

                with np.errstate(divide='ignore'):
                    logT = np.log(np.clip(T, 1e-300, 1.0))  # guard zeros
                tmp = alpha_pa[:, None] + logT            # (S_parent, S_child)
                m = np.max(tmp, axis=0)
                contrib = m + np.log(np.sum(np.exp(tmp - m), axis=0))
                log_trans_accum = np.logaddexp(log_trans_accum, contrib)

            # final
            log_alpha = log_emit_vec + log_trans_accum
            for j, ph in enumerate(state_order):
                forward_messages[t][node][ph] = float(log_alpha[j])

    return forward_messages











def assign_slices_and_interfaces(nodes, edges):
    """
    Assign nodes to slices in a directed acyclic graph (DAG).
    Compute incoming and outgoing interfaces for a graph.
    
    Parameters:
    - nodes: List of nodes in the graph.
    - edges: List of directed edges (tuples) in the graph, where each edge is (source, target).
    
    Returns:
    - slices: List of slices, where each slice is a set of nodes.
    - interfaces: Dictionary containing incoming and outgoing interfaces for each slice.
    """
    # Build adjacency and in-degree and reverse adjacency list
    adjacency_list = defaultdict(list)
    reverse_adjacency_list = defaultdict(list)
    in_degree = {node: 0 for node in nodes}
    for s, t in edges:
        adjacency_list[s].append(t)
        reverse_adjacency_list[t].append(s)
        in_degree[t] += 1


    # Initial slice with zero in-degree nodes
    current_slice = {n for n in nodes if in_degree[n] == 0}
    slice_index = 1
    slices = {slice_index: current_slice}


    while len(current_slice) > 0:
        # print(slice_index, slices)
        slice_index += 1
        successors = set()
        for node in current_slice:
            for nbr in adjacency_list[node]:
                successors.add(nbr)
        current_slice = successors
        slices[slice_index]= current_slice

    slices = {i: sorted(list(slices[i])) for i in slices if len(slices[i]) > 0}


    # Step 3: Compute interfaces
    slice_keys = sorted(slices.keys())  # Ensure slices are processed in order
    interfaces = {"incoming": {}, "outgoing": {}}

    for t in slice_keys:
        current_slice = set(slices[t])
        next_slice = set(slices[t + 1]) if t + 1 in slices else set()
        prev_slice = set(slices[t - 1]) if t - 1 in slices else set()

        # Outgoing interface for slice t
        interfaces["outgoing"][t] = {
            node for node in current_slice if any(nbr in next_slice for nbr in adjacency_list[node])
        }

        # Incoming interface for slice t
        interfaces["incoming"][t] = {
            node for node in current_slice if any(nbr in prev_slice for nbr in reverse_adjacency_list[node])
        }

    interfaces["incoming"] = {i: sorted(list(interfaces["incoming"][i])) for i in slices}
    interfaces["outgoing"] = {i: sorted(list(interfaces["outgoing"][i])) for i in slices}

    return slices, interfaces


def assign_evidence_to_states_and_transitions_from_M(nodes, edges, M: csr_matrix, snp_to_col: Dict[int, int] = None) -> Dict[str, Dict[str, List[int]]]:
    # print(snp_to_col)
    assert snp_to_col is not None, "Pass frag.snp_to_col (maps 0-based positions -> columns)"
    O = (M > 0).astype(np.int8)

    def cols_from_label(lbl: str) -> np.ndarray:
        i_str, j_str = lbl.split("-")
        i = int(i_str); j = int(j_str)   # 1-based positions
        return np.array([snp_to_col[i-1], snp_to_col[j-1]], dtype=np.int64)

    def rows_covering_all(cols: np.ndarray) -> np.ndarray:
        sub = O[:, cols]
        return np.where(np.asarray(sub.getnnz(axis=1)).ravel() == cols.size)[0].astype(np.int64)

    states = {}
    for node in nodes:
        states[node] = rows_covering_all(cols_from_label(node)).tolist()

    trans = {}
    for u, v in edges:
        ui, uj = map(int, u.split("-"))
        vi, vj = map(int, v.split("-"))
        cols_pos = sorted({ui, uj, vi, vj})
        cols = np.array([snp_to_col[p-1] for p in cols_pos], dtype=np.int64)
        trans[f"{u}--{v}"] = rows_covering_all(cols).tolist()

    return {"states": states, "transitions": trans}


def _tcounts_to_phas_strs(K: int, t_counts: np.ndarray) -> List[str]:
    """Convert each t=[t00,t01,t10,t11] (in bank order) into a canonical Kx2 hap-matrix string."""
    out = []
    for t00, t01, t10, t11 in t_counts:
        rows = []
        if t00: rows.append(np.tile([0,0], (int(t00), 1)))
        if t01: rows.append(np.tile([0,1], (int(t01), 1)))
        if t10: rows.append(np.tile([1,0], (int(t10), 1)))
        if t11: rows.append(np.tile([1,1], (int(t11), 1)))
        H = np.vstack(rows) if rows else np.zeros((0,2), dtype=int)
        if H.shape[0] != K:                     # pad/truncate to K rows (degenerate edge cases)
            if H.shape[0] < K:
                H = np.vstack([H, np.tile([0,0], (K - H.shape[0], 1))])
            else:
                H = H[:K]
        out.append(''.join(map(str, H.ravel().tolist())))
    return out


def build_state_names_from_bank(pair_layer, gi, gj, bank, ploidy: int, one_based_labels: bool = True, col_to_snp: List[int] = None) -> Dict[str, List[str]]:
    assert col_to_snp is not None, "Pass frag.col_to_snp"
    names = {}
    for p, (i, j) in enumerate(pair_layer.nodes):
        node_lbl = _node_label_from_cols(int(i), int(j), col_to_snp)
        model = bank[(int(gi[p]), int(gj[p]))]
        names[node_lbl] = _tcounts_to_phas_strs(ploidy, model.t_counts)
    return names


def _node_label_from_cols(col_i: int, col_j: int, col_to_snp: List[int]) -> str:
    """Return 1-based 'i-j' label from 0-based columns."""
    pos_i = int(col_to_snp[col_i]) + 1
    pos_j = int(col_to_snp[col_j]) + 1
    if pos_i > pos_j:
        pos_i, pos_j = pos_j, pos_i
    return f"{pos_i}-{pos_j}"


# def vector_emissions_to_string_keyed(log_emissions: Dict[str, np.ndarray], state_names: Dict[str, List[str]]) -> Dict[str, Dict[str, float]]:
#     """node -> {phas_str: log_emission}"""
#     out = {}
#     for node, vec in log_emissions.items():
#         names = state_names[node]
#         assert len(names) == len(vec), f"state length mismatch for {node}"
#         out[node] = {names[i]: float(vec[i]) for i in range(len(names))}
#     return out


def sample_states_book(slices: Dict[int, List[str]], edges: List[Tuple[str, str]], forward_messages: Dict[int, Dict[str, Dict[str, float]]], transitions_dict: Dict[str, np.ndarray]) -> Dict[int, Dict[str, str]]:
    """
    Robust backward sampler:
    - respects actual slice indices (not assumed to be 1..T without gaps),
    - handles both transition orientations ("u--v" and "v--u"),
    - uses string state labels (from forward_messages keys) consistently,
    - guards against zero/NaN/Inf weights; falls back gracefully.

    Returns: sampled_states[t][node] = selected state label (string).
    """
    # Build slice order and a quick map node->slice
    slice_keys = sorted(slices.keys())
    node_to_slice = {}
    for t in slice_keys:
        for n in slices[t]:
            node_to_slice[n] = t

    # Children adjacency by explicit edges
    children = defaultdict(list)
    for u, v in edges:
        children[u].append(v)

    def _safe_softmax(log_w: np.ndarray) -> np.ndarray:
        """Return a finite probability vector; if degenerate, fallback to uniform; if still degenerate, return zeros."""
        if not np.all(np.isfinite(log_w)):
            # set -inf to very negative; ignore NaN by treating them as -inf
            log_w = np.where(np.isfinite(log_w), log_w, -np.inf)
        m = np.max(log_w)
        if not np.isfinite(m):
            # all -inf -> uniform fallback
            n = log_w.shape[0]
            if n == 0:
                return np.array([])
            return np.full(n, 1.0 / n, dtype=np.float64)
        w = np.exp(log_w - m)
        s = w.sum()
        if not np.isfinite(s) or s <= 0.0:
            # degeneracy -> uniform
            n = log_w.shape[0]
            if n == 0:
                return np.array([])
            return np.full(n, 1.0 / n, dtype=np.float64)
        return w / s


    def _get_transition_row(parent: str, child: str, parent_states: List[str], child_states: List[str], pi: int) -> np.ndarray:
        """
        Return the row of P(child|parent) for parent-state index pi, aligned to child_states order.
        If "parent--child" is missing, try "child--parent" and convert:
           P(parent|child) (S_child x S_parent)  ->  col-normalize, then transpose -> P(child|parent).
        If still missing or degenerate, return None.
        """
        key = f"{parent}--{child}"
        T = transitions_dict.get(key, None)
        if T is not None:
            # expected shape (S_parent x S_child)
            if T.shape == (len(parent_states), len(child_states)):
                row = T[pi, :]
                # guard for negatives / nans
                row = np.clip(row, 0.0, np.inf)
                if row.sum() == 0.0:
                    return None
                return row
            # shape mismatch: try to adapt if transposed accidental
            if T.shape == (len(child_states), len(parent_states)):
                # Looks like P(child|parent) but transposed; transpose
                T2 = T.T
                row = T2[pi, :]
                row = np.clip(row, 0.0, np.inf)
                if row.sum() == 0.0:
                    return None
                return row
            # wrong shape
            return None

        # try reversed
        key_rev = f"{child}--{parent}"
        T_rev = transitions_dict.get(key_rev, None)
        if T_rev is None:
            return None
        # T_rev is P(parent|child) with shape (S_child x S_parent) ideally
        if T_rev.shape == (len(child_states), len(parent_states)):
            # convert: col-normalize then transpose -> (S_parent x S_child)
            col_sum = T_rev.sum(axis=0, keepdims=True).clip(min=1e-12)
            T_conv = (T_rev / col_sum).T
            row = T_conv[pi, :]
            row = np.clip(row, 0.0, np.inf)
            if row.sum() == 0.0:
                return None
            return row

        # maybe T_rev is already (S_parent x S_child) by mistake
        if T_rev.shape == (len(parent_states), len(child_states)):
            row = T_rev[pi, :]
            row = np.clip(row, 0.0, np.inf)
            if row.sum() == 0.0:
                return None
            return row

        return None

    sampled_states: Dict[int, Dict[str, str]] = {}

    # Step 1: sample at the last slice
    t_last = slice_keys[-1]
    sampled_states[t_last] = {}
    for node in slices[t_last]:
        state_labels = list(forward_messages[t_last][node].keys())  # preserve order
        log_alphas = np.array([forward_messages[t_last][node][s] for s in state_labels], dtype=np.float64)
        probs = _safe_softmax(log_alphas)
        # in case of any corner-case numerical drift: re-normalize
        probs = probs / probs.sum() if probs.size and probs.sum() > 0 else probs
        sampled_states[t_last][node] = random.choices(state_labels, weights=probs, k=1)[0]

    # Step 2: go backward through slices
    for t in reversed(slice_keys[:-1]):
        sampled_states[t] = {}
        for node in slices[t]:
            state_labels = list(forward_messages[t][node].keys())
            S = len(state_labels)
            if S == 0:
                sampled_states[t][node] = None
                continue

            log_alpha = np.array([forward_messages[t][node][s] for s in state_labels], dtype=np.float64)

            # only consider children in the *next* slice (as in your original function)
            t_next = None
            # find the very next slice that exists (in case slice indices skip)
            for tt in slice_keys:
                if tt > t:
                    t_next = tt
                    break
            child_nodes = []
            if t_next is not None:
                # filter children that are actually in that next slice
                child_nodes = [ch for ch in children[node] if ch in slices[t_next]]

            # accumulate child constraints
            add_log = np.zeros(S, dtype=np.float64)
            for ch in child_nodes:
                # child's state-space in *the next slice*
                ch_states = list(forward_messages[t_next][ch].keys())
                if len(ch_states) == 0:
                    continue
                ch_sampled = sampled_states[t_next][ch]
                try:
                    cj = ch_states.index(ch_sampled)
                except ValueError:
                    # Shouldn't happen, but guard anyway: skip this child
                    continue

                # Build a vector of child likelihood given each parent state: T[parent_idx, cj]
                contrib = np.full(S, -np.inf, dtype=np.float64)
                for pi in range(S):
                    row = _get_transition_row(node, ch, state_labels, ch_states, pi)
                    if row is None:
                        continue
                    p = float(row[cj]) if np.isfinite(row[cj]) else 0.0
                    # avoid log(0)
                    if p <= 0.0:
                        contrib[pi] = -np.inf
                    else:
                        contrib[pi] = math.log(p)
                # combine: add to running sum (log-space)
                add_log = add_log + contrib  # if contrib has -inf, _safe_softmax will handle later

            # final log-probs for the node
            log_post = log_alpha + add_log
            probs = _safe_softmax(log_post)
            probs = probs / probs.sum() if probs.size and probs.sum() > 0 else probs

            # guard one more time: if still degenerate, fall back to alpha-only
            if probs.size == 0 or not np.all(np.isfinite(probs)) or probs.sum() <= 0.0:
                probs = _safe_softmax(log_alpha)

            # and if that still fails (all -inf), uniform
            if probs.sum() <= 0.0 or not np.isfinite(probs.sum()):
                probs = np.full(S, 1.0 / S, dtype=np.float64)

            sampled_states[t][node] = random.choices(state_labels, weights=probs, k=1)[0]

    return sampled_states


def compute_forward_messages_vectorized(slices: Dict[int, List[str]], edges: List[Tuple[str, str]], emission_dict: Dict[str, object], transitions_dict: Dict[str, np.ndarray], state_names: Optional[Dict[str, List[str]]] = None, add_uniform_prior: bool = True) -> Dict[int, Dict[str, Dict[str, float]]]:
    """
    Log-space forward pass over a layered DAG.

    Inputs
    ------
    slices : {t: [node labels]}
        Layering of nodes (as returned by your assign_slices_and_interfaces).
    edges : list of (u, v)
        Directed edges; transitions are keyed by "u--v".
    emission_dict :
        Either dict[node] -> dict[state_label] -> log_emission   (Mode A)
        or     dict[node] -> np.ndarray(shape=(S,))              (Mode B)
    transitions_dict :
        dict["u--v"] -> np.ndarray of shape (S_u, S_v), row-normalized probabilities.
    state_names :
        Optional in Mode B: dict[node] -> list[str] of state labels (same order as vector).
        If None in Mode B, labels default to 's0','s1',...
        Ignored in Mode A.

    Returns
    -------
    fwd : {t: {node: {state_label: log_alpha}}}
        Same nested structure as your old code, with phasing-string (or provided) labels.
    """
    # Build parents and quick node->slice map
    parents = defaultdict(list)
    for u, v in edges:
        parents[v].append(u)
    slice_keys = sorted(slices.keys())
    node_to_slice = {}
    for t in slice_keys:
        for n in slices[t]:
            node_to_slice[n] = t

    # Normalize emission_dict into Mode A: node -> {state_label: log_emission}
    emA: Dict[str, Dict[str, float]] = {}
    if isinstance(next(iter(emission_dict.values())), dict):
        # Mode A: already label->log emission
        for node, d in emission_dict.items():
            # ensure values are float
            emA[node] = {k: float(v) for k, v in d.items()}
    else:
        # Mode B: vector + (optional) state_names
        for node, vec in emission_dict.items():
            vec = np.asarray(vec, dtype=np.float64)
            if state_names is not None and node in state_names:
                labels = state_names[node]
                assert len(labels) == vec.shape[0], f"state_names length mismatch for {node}"
            else:
                labels = [f"s{i}" for i in range(vec.shape[0])]
            emA[node] = {labels[i]: float(vec[i]) for i in range(vec.shape[0])}

    # Prepare output dict
    fwd: Dict[int, Dict[str, Dict[str, float]]] = {t: {} for t in slice_keys}

    # Base slice
    t0 = slice_keys[0]
    for node in slices[t0]:
        states = list(emA[node].keys())
        log_emit = np.array([emA[node][s] for s in states], dtype=np.float64)
        if add_uniform_prior and len(states) > 0:
            log_emit = log_emit + (-math.log(len(states)))
        fwd[t0][node] = {states[i]: float(log_emit[i]) for i in range(len(states))}

    # Subsequent slices
    for t in slice_keys[1:]:
        for node in slices[t]:
            states = list(emA[node].keys())
            S_node = len(states)
            if S_node == 0:
                fwd[t][node] = {}
                continue

            # accumulate parent contributions in log-space
            log_msg = np.full(S_node, -np.inf, dtype=np.float64)

            # for each parent, pull alphas from its **actual** slice
            for pa in parents.get(node, []):
                pt = node_to_slice.get(pa, None)
                if pt is None or pt >= t:
                    # parent not placed earlier? skip defensively
                    continue

                pa_states = list(emA[pa].keys())
                S_pa = len(pa_states)
                if S_pa == 0:
                    continue

                # alpha vector of the parent
                alpha_pa = np.array([fwd[pt][pa][s] for s in pa_states], dtype=np.float64)

                # transition pa->node
                key = f"{pa}--{node}"
                T = transitions_dict.get(key, None)
                if T is None:
                    # try reverse and convert to P(child|parent)
                    T_rev = transitions_dict.get(f"{node}--{pa}", None)
                    if T_rev is None:
                        continue
                    # column-normalize then transpose -> row-normalized P(child|parent)
                    col_sum = T_rev.sum(axis=0, keepdims=True).clip(min=1e-12)
                    T = (T_rev / col_sum).T

                if T.shape != (S_pa, S_node):
                    raise ValueError(f"Transition shape mismatch for {key}: {T.shape} vs ({S_pa},{S_node})")

                with np.errstate(divide='ignore'):
                    logT = np.log(np.clip(T, 1e-300, 1.0))  # (S_pa, S_node)
                tmp = alpha_pa[:, None] + logT             # (S_pa, S_node)
                m = np.max(tmp, axis=0)
                per_parent = m + np.log(np.sum(np.exp(tmp - m), axis=0))
                log_msg = np.logaddexp(log_msg, per_parent)

            # add local emission + (optional) uniform prior
            log_emit = np.array([emA[node][s] for s in states], dtype=np.float64)
            if add_uniform_prior and S_node > 0:
                log_emit = log_emit + (-math.log(S_node))
            log_alpha = log_emit + log_msg

            fwd[t][node] = {states[i]: float(log_alpha[i]) for i in range(S_node)}

    return fwd


def build_pair_emissions_from_state_names(
    state_names: Dict[str, List[str]],  # node -> list of phasing strings (bank order)
    ploidy: int,
    error_rate: float,
) -> Dict[str, Dict[str, Dict[str, float]]]:

    """
    For each node '{i}-{j}' and each phasing string in state_names[node],
    compute the emission likelihoods for the four 2-bit observation patterns:
      '00', '01', '10', '11'
    using your existing compute_likelihood on the Kx2 phasing matrix.

    Returns:
      emission_dict[node][phasing_str][obs_str] = likelihood (float)
    """
    # local cache keyed by phasing_str -> {'00': val, '01': val, '10': val, '11': val}
    cache: Dict[str, Dict[str, float]] = {}

    # Patterns in a fixed, deterministic order
    patterns = [''.join(bits) for bits in product('01', repeat=2)]

    out: Dict[str, Dict[str, Dict[str, float]]] = {}
    for node, phasings in state_names.items():
        out[node] = {}
        for ph_str in phasings:
            if ph_str in cache:
                out[node][ph_str] = cache[ph_str]
                continue

            H = str_2_phas_1(ph_str, ploidy)        # (K x 2) matrix of 0/1
            # compute per-pattern likelihoods (just like your old code)
            entry = {}
            for pat in patterns:
                obs = np.fromiter((int(c) for c in pat), dtype=int)
                entry[pat] = float(compute_likelihood(obs, H, error_rate))
            out[node][ph_str] = entry
            cache[ph_str] = entry
    return out



def _get_pair_reads(
    M: csr_matrix,
    snp_to_col: Dict[int, int],
    positions: List[int]  # 1-based, sorted [i, j]
) -> Dict:
    """
    Get reads covering both positions.
    Returns observations in {0,1} format (REF=0, ALT=1).
    """
    cols = [snp_to_col[p - 1] for p in positions]
    sub = M[:, cols].toarray()  # (n_reads x 2)
    
    # Only keep reads covering both positions
    complete = (sub > 0).all(axis=1)
    ridx = np.where(complete)[0]
    n_reads = len(ridx)
    
    if n_reads == 0:
        return {
            'n_reads': 0,
            'observations': np.empty((0, 2), dtype=np.int8),
            'ridx': np.empty(0, dtype=np.int64)
        }
    
    # Convert: 1=REF→0, 2=ALT→1
    observations = (sub[ridx] == 2).astype(np.int8)
    
    return {
        'n_reads': n_reads,
        'observations': observations,
        'ridx': ridx
    }

def compute_emissions_from_reads(
    M: csr_matrix,
    pair_layer,
    state_names: Dict[str, List[str]],
    ploidy: int,
    error_rate: float,
    col_to_snp: np.ndarray,
    snp_to_col: Dict[int, int],
    return_extra: bool = True
):
    """
    Compute emission probabilities using actual reads covering each node.
    
    For node '1-2': uses reads covering positions {1, 2}
    For node '2-3': uses reads covering positions {2, 3}
    
    This matches the transition computation: use real data, not parametric model.
    
    Args:
        M: (n_reads x n_snps) sparse matrix {0,1,2}
        pair_layer: graph structure with nodes
        state_names: Dict[node_label, List[phasing_str]]
        ploidy: number of haplotypes
        error_rate: sequencing error rate
        col_to_snp: column → SNP mapping
        snp_to_col: SNP → column mapping
        return_extra: return detailed information
    
    Returns:
        emissions_dict: Dict[node_label, np.ndarray]
            emissions_dict[node][state_idx] = P(reads at node | state)
        emissions_extra: Dict with detailed info (if return_extra=True)
    """
    
    emissions_dict = {}
    emissions_extra = {} if return_extra else None
    
    # Build node labels
    node_labels = []
    for (ci, cj) in pair_layer.nodes.astype(int):
        node_labels.append(_node_label(ci, cj, col_to_snp))
    
    print(f"\nComputing emissions for {len(node_labels)} nodes...")
    print("Using reads covering both positions of each node")
    
    for n_idx, (ci, cj) in enumerate(pair_layer.nodes.astype(int)):
        node_lbl = node_labels[n_idx]
        pos_i, pos_j = _parse_node(node_lbl)
        positions = sorted([pos_i, pos_j])
        
        if n_idx < 3:
            print(f"\nNode {n_idx}: {node_lbl}, positions {positions}")
        
        # Get reads covering both positions
        pair_reads = _get_pair_reads(M, snp_to_col, positions)
        n_reads = pair_reads['n_reads']
        observations = pair_reads['observations']  # (n_reads x 2) in {0,1}
        
        if n_idx < 3:
            print(f"  Reads covering both: {n_reads}")
        
        phasings_str = state_names[node_lbl]
        n_states = len(phasings_str)
        
        emissions = np.zeros(n_states, dtype=np.float64)
        
        if return_extra:
            emissions_extra[node_lbl] = {}
        
        if n_reads == 0:
            # No reads covering this node - uniform (non-informative)
            if n_idx < 3:
                print(f"  → Uniform (no reads)")
            emissions = np.ones(n_states) / n_states
            
            if return_extra:
                for i, phas_str in enumerate(phasings_str):
                    emissions_extra[node_lbl][str(i)] = {
                        'phasing': phas_str,
                        'n_reads': 0,
                        'likelihood': emissions[i],
                        'prior': 'uniform'
                    }
        else:
            # Compute from actual reads
            for i, phas_str in enumerate(phasings_str):
                phas_matrix = _str_to_phasing_matrix(phas_str, ploidy, n_cols=2)
                
                # ADAPTIVE: same as transitions
                if ploidy <= n_reads:
                    # Loop over ploidy (typical)
                    likelihood = 0.0
                    for k in range(ploidy):
                        lk = compute_likelihood(phas_matrix[k,:], observations, error_rate)
                        likelihood += lk
                else:
                    # Loop over reads (high ploidy, few reads)
                    likelihood = 0.0
                    for r in range(n_reads):
                        lk = compute_likelihood(observations[r], phas_matrix, error_rate)
                        likelihood += lk
                
                emissions[i] = likelihood
                
                if return_extra:
                    emissions_extra[node_lbl][str(i)] = {
                        'phasing': phas_str,
                        'n_reads': n_reads,
                        'likelihood': likelihood
                    }
            
            # Normalize to get probabilities
            total = emissions.sum()
            if total > 0 and np.isfinite(total):
                emissions = emissions / total
            else:
                # Fallback to uniform if all zero or non-finite
                emissions = np.ones(n_states) / n_states
            
            if n_idx < 3:
                print(f"  Emission probabilities: {emissions}")
        
        emissions_dict[node_lbl] = emissions
    
    print(f"\nCompleted {len(emissions_dict)} nodes")
    return emissions_dict, emissions_extra if return_extra else {}


def _softmax_safe(log_w: np.ndarray) -> np.ndarray:
    if log_w.size == 0:
        return log_w
    m = np.max(log_w)
    if not np.isfinite(m):
        # all -inf -> uniform
        return np.full(log_w.shape, 1.0/log_w.size, dtype=np.float64)
    w = np.exp(log_w - m)
    s = w.sum()
    if not np.isfinite(s) or s <= 0.0:
        return np.full(log_w.shape, 1.0/log_w.size, dtype=np.float64)
    return w / s


def sample_states_book_classic(
    slices: Dict[int, List[str]],
    edges: List[Tuple[str, str]],
    forward_messages: Dict[int, Dict[str, Dict[str, float]]],
    transitions_dict: Dict[str, np.ndarray],
) -> Dict[int, Dict[str, str]]:
    """
    Backward sampling using ONLY alpha_t and P(child|parent), identical in spirit to your old function:
      1) Sample last slice from alpha_T
      2) For t = T-1..1, for each node:
           log p(state) ∝ alpha_t(state) + sum_{children in next slice} log T_parent->child[state, sampled_child]
    Assumes transitions are keyed 'parent--child' and in the correct shape.

    Returns: sampled_states[t][node] = state_label (string)
    """
    # Precompute child adjacency by next slice only (as in the old code)
    slice_keys = sorted(slices.keys())
    if not slice_keys:
        return {}
    next_children = defaultdict(list)  # (t, node) -> [children in t+1]
    # Build a quick set per slice for membership test
    slice_sets = {t: set(slices[t]) for t in slice_keys}
    # full children
    children = defaultdict(list)
    for u, v in edges:
        children[u].append(v)
    for idx, t in enumerate(slice_keys[:-1]):
        t_next = slice_keys[idx+1]
        next_set = slice_sets[t_next]
        for u in slices[t]:
            # children that are actually in the next slice
            next_children[(t, u)] = [ch for ch in children[u] if ch in next_set]

    # Precompute state orders and index maps for all nodes in all slices
    state_order: Dict[Tuple[int,str], List[str]] = {}
    state_index: Dict[Tuple[int,str], Dict[str,int]] = {}
    for t in slice_keys:
        for node in slices[t]:
            so = list(forward_messages[t][node].keys())
            state_order[(t,node)] = so
            state_index[(t,node)] = {s:i for i,s in enumerate(so)}

    sampled_states: Dict[int, Dict[str, str]] = {}

    # Step 1: sample last slice from alpha_T
    t_last = slice_keys[-1]
    sampled_states[t_last] = {}
    for node in slices[t_last]:
        so = state_order[(t_last,node)]
        log_alpha = np.array([forward_messages[t_last][node][s] for s in so], dtype=np.float64)
        probs = _softmax_safe(log_alpha)
        sampled_states[t_last][node] = random.choices(so, weights=probs, k=1)[0]

    # Step 2: go backward t = ... T-1,..., first
    for idx in range(len(slice_keys)-2, -1, -1):
        t = slice_keys[idx]
        t_next = slice_keys[idx+1]
        sampled_states[t] = {}
        for node in slices[t]:
            so = state_order[(t,node)]
            S = len(so)
            if S == 0:
                sampled_states[t][node] = None
                continue

            # alpha_t for this node (vectorized)
            log_alpha = np.array([forward_messages[t][node][s] for s in so], dtype=np.float64)

            # sum over children in the *next* slice
            add = np.zeros(S, dtype=np.float64)
            for ch in next_children[(t,node)]:
                # child's sampled state index (in next slice's order)
                ch_so = state_order[(t_next,ch)]
                ch_idx = state_index[(t_next,ch)][ sampled_states[t_next][ch] ]
                # transition matrix parent->child must exist in this direction
                T = transitions_dict.get(f"{node}--{ch}", None)
                if T is None:
                    # if you want to support reversed keys add a conversion here; we assume canonical u->v
                    continue
                # add log T[:, ch_idx]
                col = T[:, ch_idx]
                # guard zero
                col = np.clip(col, 1e-300, 1.0)
                add += np.log(col)

            log_post = log_alpha + add
            probs = _softmax_safe(log_post)
            sampled_states[t][node] = random.choices(so, weights=probs, k=1)[0]

    return sampled_states


def _group_counts(vec: np.ndarray) -> Tuple[int,int]:
    """Return (#zeros, #ones) in a 1D 0/1 vector."""
    ones = int(vec.sum())
    zeros = vec.size - ones
    return zeros, ones

def _enumerate_contingencies(z_a: int, o_a: int, z_b: int, o_b: int) -> List[Tuple[int,int,int,int]]:
    """
    Enumerate feasible contingency counts (c00,c01,c10,c11) for pairing two 0/1 multiset vectors
    of the same length m, where left has z_a zeros and o_a ones, right has z_b zeros and o_b ones.
    Constraints:
      c00 + c01 = z_a
      c10 + c11 = o_a
      c00 + c10 = z_b
      c01 + c11 = o_b
      all counts >= 0
    This is tiny (≤ O(m^2)).
    """
    m = z_a + o_a
    out = []
    # choose c00 from feasible range
    c00_min = max(0, z_b - o_a)
    c00_max = min(z_a, z_b)
    for c00 in range(c00_min, c00_max+1):
        c01 = z_a - c00
        c10 = z_b - c00
        c11 = o_a - c10
        if c01 < 0 or c10 < 0 or c11 < 0: 
            continue
        if c01 + c11 != o_b:
            continue
        out.append((c00, c01, c10, c11))
    return out

def _triple_from_contingencies(c0: Tuple[int,int,int,int], c1: Tuple[int,int,int,int]) -> np.ndarray:
    """
    Build K×3 triple matrix [a,s,b] from two group contingency tuples:
      c0 = (c00,c01,c10,c11) for s=0 rows
      c1 = (d00,d01,d10,d11) for s=1 rows
    where c00=#(a=0,b=0), c01=#(a=0,b=1), c10=#(a=1,b=0), c11=#(a=1,b=1)
    and similarly for d.. in the s=1 group.
    The row order inside each group is irrelevant; we canonicalize later anyway.
    """
    rows = []
    c00,c01,c10,c11 = c0
    d00,d01,d10,d11 = c1
    # s=0 rows
    rows += [[0,0,0]] * c00
    rows += [[0,0,1]] * c01
    rows += [[1,0,0]] * c10
    rows += [[1,0,1]] * c11
    # s=1 rows
    rows += [[0,1,0]] * d00
    rows += [[0,1,1]] * d01
    rows += [[1,1,0]] * d10
    rows += [[1,1,1]] * d11
    return np.array(rows, dtype=np.int8)

def build_transitions_readdriven_from_M(
    M,
    pair_layer,
    state_names: Dict[str, List[str]],
    ploidy: int,
    error_rate: float,
    col_to_snp,
    snp_to_col,
    coverage_mode: str = "atleast2",   # "exact3" or "atleast2"
    return_extra: bool = True          # set False to save RAM
):
    """
    Read-driven transitions with counts-based merge (no factorial perms).
    For each edge {x,y}--{u,v} sharing exactly one SNP s:
      - For each (src,target) phasing pair:
          * compute counts in s=0 and s=1 groups for a (from src) and b (from tgt)
          * enumerate small set of feasible contingency 4-tuples per group
          * build each merged triple and score by groups of reads (no per-read loop)
          * sum weights; row-normalize to P(target|source)
    coverage_mode:
      - "exact3"   : use only reads covering {a,s,b}
      - "atleast2" : include (a,s), (s,b), (a,b) as in your latest
    """
    assert coverage_mode in ("exact3", "atleast2")
    transitions = {}
    extra = {} if return_extra else None

    # Precompute per-node Kx2 matrices (bank order)
    node_labels = []
    for (ci,cj) in pair_layer.nodes.astype(int):
        node_labels.append(_node_label(ci, cj, col_to_snp))  # 1-based label i<j
    pair_H = {lbl: [_phas_str_to_mat(s, ploidy) for s in state_names[lbl]] for lbl in node_labels}

    for e_idx, (_u_idx, _v_idx) in enumerate(pair_layer.edges.astype(int)):
        u_lbl = node_labels[int(_u_idx)]
        v_lbl = node_labels[int(_v_idx)]
        u_i, u_j = _parse_node(u_lbl)  # (i<j)
        v_i, v_j = _parse_node(v_lbl)

        U = {u_i, u_j}; V = {v_i, v_j}
        shared = list(U & V); assert len(shared) == 1
        s = shared[0]
        a = list(U - {s})[0]
        b = list(V - {s})[0]

        # We need H_left=[a,s] from u, H_right=[s,b] from v (reorder columns)
        def reorder(Hpairs, first, second, wantL, wantR):
            out=[]
            for H in Hpairs:
                col0 = H[:,0] if first==wantL else H[:,1]
                col1 = H[:,0] if first==wantR else H[:,1]
                out.append(np.stack([col0, col1], axis=1))
            return out

        first_u, second_u = min(u_i,u_j), max(u_i,u_j)
        first_v, second_v = min(v_i,v_j), max(v_i,v_j)
        H_lefts  = reorder(pair_H[u_lbl], first_u, second_u, a, s)   # [a,s]
        H_rights = reorder(pair_H[v_lbl], first_v, second_v, s, b)   # [s,b]
        Su, Sv = len(H_lefts), len(H_rights)

        # Gather reads & group by mask for this triple set
        # Always build triple groups (a,s,b)
        pos_trip = [a, s, b]
        groups3, obs3, rows3, _ = _submatrix_and_groups(M, snp_to_col, pos_trip, min_cov=3 if coverage_mode=="exact3" else 1)

        # Pairs for "atleast2"
        if coverage_mode == "atleast2":
            # avoid double counting triples in pairs
            groups_as = obs_as = None
            groups_sb = obs_sb = None
            groups_ab = obs_ab = None
            # (a,s)
            g_as, o_as, r_as, _ = _submatrix_and_groups(M, snp_to_col, [a,s], min_cov=2)
            if rows3.size and r_as.size:
                r_as = r_as[~np.isin(r_as, rows3)]
                if r_as.size:
                    o_as = M[r_as][:, [snp_to_col[a-1], snp_to_col[s-1]]].toarray()
                    o_as = (o_as == 2).astype(np.int8); g_as = _group_by_mask(o_as> -1)  # trivial mask
            if r_as.size:
                groups_as, obs_as = g_as, o_as
            # (s,b)
            g_sb, o_sb, r_sb, _ = _submatrix_and_groups(M, snp_to_col, [s,b], min_cov=2)
            if rows3.size and r_sb.size:
                r_sb = r_sb[~np.isin(r_sb, rows3)]
                if r_sb.size:
                    o_sb = M[r_sb][:, [snp_to_col[s-1], snp_to_col[b-1]]].toarray()
                    o_sb = (o_sb == 2).astype(np.int8); g_sb = _group_by_mask(o_sb> -1)
            if r_sb.size:
                groups_sb, obs_sb = g_sb, o_sb
            # (a,b)
            g_ab, o_ab, r_ab, _ = _submatrix_and_groups(M, snp_to_col, [a,b], min_cov=2)
            if rows3.size and r_ab.size:
                r_ab = r_ab[~np.isin(r_ab, rows3)]
                if r_ab.size:
                    o_ab = M[r_ab][:, [snp_to_col[a-1], snp_to_col[b-1]]].toarray()
                    o_ab = (o_ab == 2).astype(np.int8); g_ab = _group_by_mask(o_ab> -1)
            if r_ab.size:
                groups_ab, obs_ab = g_ab, o_ab

        # build transitions matrix and (optional) extra
        W = np.zeros((Su, Sv), dtype=np.float64)
        if return_extra:
            extra[f"{u_lbl}--{v_lbl}"] = {}

        for iu, H_as in enumerate(H_lefts):   # [a,s]
            # counts for group s=0 and s=1 on the left side
            s_left = H_as[:,1]
            a_left = H_as[:,0]
            z0_a, o0_a = _group_counts(a_left[s_left==0])
            z1_a, o1_a = _group_counts(a_left[s_left==1])

            for jv, H_sb in enumerate(H_rights):  # [s,b]
                entry_key = f"{iu}-{jv}"
                if return_extra:
                    extra[f"{u_lbl}--{v_lbl}"][entry_key] = {
                        "source_phasing": state_names[u_lbl][iu],
                        "target_phasing": state_names[v_lbl][jv],
                        "matched_phasings": {}
                    }

                s_right = H_sb[:,0]
                b_right = H_sb[:,1]
                z0_b, o0_b = _group_counts(b_right[s_right==0])
                z1_b, o1_b = _group_counts(b_right[s_right==1])

                # enumerate contingency tuples per group
                C0 = _enumerate_contingencies(z0_a, o0_a, z0_b, o0_b)
                C1 = _enumerate_contingencies(z1_a, o1_a, z1_b, o1_b)
                if not C0 or not C1:
                    continue

                total = 0.0
                for c0 in C0:
                    for c1 in C1:
                        H_triple = _triple_from_contingencies(c0, c1)  # (K x 3)

                        # score triples (vectorized by mask groups)
                        if obs3 is not None and len(groups3) > 0:
                            total += _score_mixture_by_groups(obs3, groups3, H_triple, error_rate)

                        if coverage_mode == "atleast2":
                            # add pair contributions by restricting H_triple to observed columns
                            if obs_as is not None and len(groups_as or {})>0:
                                total += _score_mixture_by_groups(obs_as, groups_as, H_triple[:, [0,1]], error_rate)
                            if obs_sb is not None and len(groups_sb or {})>0:
                                total += _score_mixture_by_groups(obs_sb, groups_sb, H_triple[:, [1,2]], error_rate)
                            if obs_ab is not None and len(groups_ab or {})>0:
                                total += _score_mixture_by_groups(obs_ab, groups_ab, H_triple[:, [0,2]], error_rate)

                        if return_extra:
                            tstr = ''.join(map(str, H_triple.ravel().tolist()))
                            extra[f"{u_lbl}--{v_lbl}"][entry_key]["matched_phasings"][tstr] = \
                                extra[f"{u_lbl}--{v_lbl}"][entry_key]["matched_phasings"].get(tstr, 0.0) + 0.0  # set or keep 0 (weights kept in W)

                W[iu, jv] = total

        # row-normalize; uniform for zero rows
        row_sum = W.sum(axis=1, keepdims=True)
        bad = (row_sum.squeeze() <= 0.0) | (~np.isfinite(row_sum.squeeze()))
        if np.any(bad):
            for r in np.where(bad)[0]:
                if Sv > 0:
                    W[r,:] = 1.0 / Sv
        else:
            W = W / row_sum

        transitions[f"{u_lbl}--{v_lbl}"] = W

        # free big temporaries before next edge
        del W, groups3, obs3
        if coverage_mode == "atleast2":
            del groups_as, obs_as, groups_sb, obs_sb, groups_ab, obs_ab

    return transitions, extra if return_extra else {}



def _submatrix_and_groups(
    M: csr_matrix,
    snp_to_col: Dict[int,int],
    pos_list_1based: list,
    min_cov: int = 1,
    rows_filter: Optional[np.ndarray] = None,
):
    """
    Extract a submatrix for positions (1-based) and group rows by coverage mask.

    Returns:
      groups : dict{mask_tuple (bool T) -> np.ndarray of row indices (local in obs)}
      obs    : np.ndarray (R x T) in {0,1} (ALT=1, REF=0)
      ridx   : np.ndarray of original row indices of M (length R) aligned to obs
      cols   : np.ndarray of 0-based column indices in M for the positions (length T)
    Only keeps rows whose coverage on these T columns is >= min_cov.
    If rows_filter is provided (array of row indices), restrict to those rows first.
    """
    T = len(pos_list_1based)
    cols = np.array([snp_to_col[p-1] for p in pos_list_1based], dtype=np.int64)

    # slice rows first if requested
    if rows_filter is not None:
        sub = M[rows_filter][:, cols]
        ridx_global = np.asarray(rows_filter, dtype=np.int64)
    else:
        sub = M[:, cols]
        ridx_global = None

    # coverage mask on these T columns
    cov = (sub > 0).astype(np.int8)
    cov_counts = np.asarray(cov.sum(axis=1)).ravel()

    keep_mask = (cov_counts >= min_cov)
    if not np.any(keep_mask):
        return {}, np.empty((0, T), dtype=np.int8), np.empty((0,), dtype=np.int64), cols

    subk = sub[keep_mask, :].toarray()                # (R x T) {0,1,2}
    obs  = (subk == 2).astype(np.int8)               # ALT=1, REF=0
    mask = (subk > 0)                                 # (R x T) boolean

    if ridx_global is None:
        ridx = np.where(keep_mask)[0].astype(np.int64)
    else:
        ridx = ridx_global[keep_mask]

    # group rows by their coverage pattern
    groups = defaultdict(list)
    for r in range(obs.shape[0]):
        groups[tuple(mask[r, :].tolist())].append(r)
    for k in list(groups.keys()):
        groups[k] = np.array(groups[k], dtype=np.int64)

    return groups, obs, ridx, cols

def _score_mixture_by_groups(
    obs: np.ndarray,
    groups: Dict[Tuple[bool, ...], np.ndarray],
    H: np.ndarray,                # (K x T) haplotype matrix over the same T columns
    error_rate: float
) -> float:
    """
    Sum of mixture likelihoods over all read rows in 'obs', grouped by their coverage mask.
    For a group with mask M, we restrict both obs and H to columns where M is True, then compute:
      L_r = sum_{h=1..K} prod_{t in covered} [ (1-e) if obs==H else e ]
    and accumulate ∑_r L_r across all groups.
    """
    if obs.size == 0 or not groups:
        return 0.0
    total = 0.0
    for mask_tuple, ridx in groups.items():
        mask = np.array(mask_tuple, dtype=bool)
        if not mask.any() or ridx.size == 0:
            continue
        obs_sub = obs[ridx][:, mask]   # (Rg x Tm)
        H_sub   = H[:, mask]           # (K  x Tm)

        # vectorized mixture
        match = (obs_sub[:, None, :] == H_sub[None, :, :])   # (Rg x K x Tm)
        m  = match.sum(axis=2, dtype=np.int32)               # (Rg x K)
        Tm = H_sub.shape[1]
        mm = Tm - m
        lk = ((1.0 - error_rate) ** m) * (error_rate ** mm)  # (Rg x K)
        L  = lk.sum(axis=1)                                  # (Rg,)
        total += float(L.sum())
    return total



def _group_by_mask(mask: np.ndarray) -> Dict[Tuple[bool, ...], np.ndarray]:
    """
    Group row indices by a boolean coverage mask.

    Parameters
    ----------
    mask : (R x T) boolean array
        True where the row covers the column; False = no-call.

    Returns
    -------
    groups : dict
        { mask_tuple (bool T) -> np.ndarray of row indices (local) }
    """
    groups = defaultdict(list)
    R = mask.shape[0]
    for r in range(R):
        groups[tuple(mask[r, :].tolist())].append(r)
    for k in list(groups.keys()):
        groups[k] = np.array(groups[k], dtype=np.int64)
    return groups


############### NEW CODE : ##################



def compute_likelihood(observed, phasing, error_rate):
    """
    Base likelihood function (user provided).
    Tiles first argument to match second argument's rows.
    """
    y = np.tile(observed, (phasing.shape[0], 1))
    diff = y - phasing
    diff[diff != 0] = 1
    comp_diff = 1 - diff
    term1 = diff * error_rate
    term2 = comp_diff * (1 - error_rate)
    terms = term1 + term2
    probs = np.prod(terms, axis=1)
    likelihood = np.sum(probs)
    return likelihood


def compute_transitions_optimized(
    M: csr_matrix,
    pair_layer,
    state_names: Dict[str, List[str]],
    ploidy: int,
    error_rate: float,
    col_to_snp: np.ndarray,
    snp_to_col: Dict[int, int],
    return_extra: bool = True
):
    """
    Optimized transition computation using complete triple reads only.
    
    Adaptive optimization: Chooses to loop over ploidy or reads based on which is smaller.
    - If ploidy ≤ n_reads: loop over ploidy
    - If ploidy > n_reads: loop over reads
    
    This is especially important for high ploidy (6, 8) with sparse triple read coverage.
    """
    
    transitions_dict = {}
    transitions_dict_extra = {}
    
    # Build node labels
    node_labels = []
    for (ci, cj) in pair_layer.nodes.astype(int):
        node_labels.append(_node_label(ci, cj, col_to_snp))
    
    # print(f"\nComputing transitions for {len(pair_layer.edges)} edges...")
    # print(f"Using ONLY complete triple reads")
    # print(f"Adaptive: loops over min(ploidy={ploidy}, n_reads)")
    
    for e_idx, (_u_idx, _v_idx) in enumerate(pair_layer.edges.astype(int)):
        u_lbl = node_labels[int(_u_idx)]
        v_lbl = node_labels[int(_v_idx)]
        u_i, u_j = _parse_node(u_lbl)
        v_i, v_j = _parse_node(v_lbl)
        
        # Find positions
        U = {u_i, u_j}
        V = {v_i, v_j}
        shared = list(U & V)
        assert len(shared) == 1, f"Edge must share exactly 1 SNP"
        
        s = shared[0]
        a = list(U - {s})[0]
        b = list(V - {s})[0]
        poss = sorted([a, s, b])
        
        # if e_idx < 3:
        #     print(f"\nEdge {e_idx}: {u_lbl}--{v_lbl}, positions {poss}")
        
        # Get complete triple reads ONLY
        complete_reads = _get_complete_triple_reads(M, snp_to_col, poss)
        n_reads = complete_reads['n_reads']
        observations = complete_reads['observations']  # (n_reads x 3)
        
        # if e_idx < 3:
        #     print(f"  Complete triple reads: {n_reads}")
        
        source_phasings_str = state_names[u_lbl]
        target_phasings_str = state_names[v_lbl]
        Su, Sv = len(source_phasings_str), len(target_phasings_str)
        
        transitions_mtx = np.zeros((Su, Sv), dtype=np.float64)
        edge_key = f"{u_lbl}--{v_lbl}"
        transitions_dict_extra[edge_key] = {}
        
        if n_reads == 0:
            # Uniform prior
            # if e_idx < 3:
            #     print(f"  → Uniform prior (no complete reads)")
            transitions_mtx = np.ones((Su, Sv)) / Sv
            
            for i in range(Su):
                for j in range(Sv):
                    transitions_dict_extra[edge_key][f"{i}-{j}"] = {
                        'source_phasing': source_phasings_str[i],
                        'target_phasing': target_phasings_str[j],
                        'matched_phasings': {},
                        'prior': 'uniform'
                    }
        else:
            # Compute from data
            for i, source_str in enumerate(source_phasings_str):
                for j, target_str in enumerate(target_phasings_str):
                    entry_key = f"{i}-{j}"
                    transitions_dict_extra[edge_key][entry_key] = {
                        'source_phasing': source_str,
                        'target_phasing': target_str,
                        'matched_phasings': {}
                    }
                    
                    source_phas = _str_to_phasing_matrix(source_str, ploidy)
                    target_phas = _str_to_phasing_matrix(target_str, ploidy)
                    
                    # Find matched phasings
                    matched_phasings = _find_matched_triple_phasings(
                        source_phas, target_phas, u_i, u_j, v_i, v_j, a, s, b, ploidy
                    )
                    
                    # Canonicalize
                    matched_canonical = []
                    for mp in matched_phasings:
                        mp_sorted = mp[np.lexsort([mp[:, col] for col in range(2, -1, -1)])]
                        matched_canonical.append(mp_sorted)
                    
                    matched_strs = list(set([
                        ''.join(map(str, mp.ravel())) for mp in matched_canonical
                    ]))
                    
                    # Sum over matched phasings
                    total_weight = 0.0
                    for phas_str in matched_strs:
                        phas_matrix = _str_to_phasing_matrix(phas_str, ploidy, n_cols=3)
                        
                        # ADAPTIVE: Choose loop strategy based on ploidy vs n_reads
                        # Important for high ploidy (6, 8) with few complete reads
                        this_phas_weight = 0.0
                        if ploidy <= n_reads:
                            # Loop over ploidy (fewer iterations)
                            for k in range(ploidy):
                                lk = compute_likelihood(phas_matrix[k,:], observations, error_rate)
                                this_phas_weight += lk
                        else:
                            # Loop over reads (fewer iterations, happens with high ploidy)
                            for r in range(n_reads):
                                lk = compute_likelihood(observations[r], phas_matrix, error_rate)
                                this_phas_weight += lk
                        
                        transitions_dict_extra[edge_key][entry_key]['matched_phasings'][phas_str] = this_phas_weight
                        total_weight += this_phas_weight
                    
                    transitions_mtx[i, j] = total_weight
            
            # Normalize
            row_sums = transitions_mtx.sum(axis=1, keepdims=True)
            bad_rows = (row_sums.squeeze() <= 0) | (~np.isfinite(row_sums.squeeze()))
            if np.any(bad_rows):
                for r in np.where(bad_rows)[0]:
                    transitions_mtx[r, :] = 1.0 / Sv
            else:
                transitions_mtx = transitions_mtx / row_sums
            
            # if e_idx < 3:
            #     print(f"  Transition matrix:\n{transitions_mtx}")
        
        transitions_dict[edge_key] = transitions_mtx
    
    # print(f"\nCompleted {len(transitions_dict)} edges")
    return transitions_dict, transitions_dict_extra


# ============================================================================
# Helper Functions
# ============================================================================

def _get_complete_triple_reads(
    M: csr_matrix,
    snp_to_col: Dict[int, int],
    positions: List[int]
) -> Dict:
    """Get reads covering ALL positions. Returns {0,1} observations."""
    cols = [snp_to_col[p - 1] for p in positions]
    sub = M[:, cols].toarray()
    complete = (sub > 0).all(axis=1)
    ridx = np.where(complete)[0]
    n_reads = len(ridx)
    
    if n_reads == 0:
        return {
            'n_reads': 0,
            'observations': np.empty((0, len(positions)), dtype=np.int8),
            'ridx': np.empty(0, dtype=np.int64)
        }
    
    observations = (sub[ridx] == 2).astype(np.int8)
    return {'n_reads': n_reads, 'observations': observations, 'ridx': ridx}


def _find_matched_triple_phasings(
    source_phas: np.ndarray, target_phas: np.ndarray,
    u_i: int, u_j: int, v_i: int, v_j: int,
    a: int, s: int, b: int, ploidy: int
) -> List[np.ndarray]:
    """Find all valid triple phasings consistent with source and target."""
    K = ploidy
    matched = []
    
    if u_i == a and u_j == s:
        source_col_a, source_col_s = 0, 1
    elif u_i == s and u_j == a:
        source_col_a, source_col_s = 1, 0
    else:
        raise ValueError(f"Source positions {u_i},{u_j} don't match a={a}, s={s}")
    
    if v_i == s and v_j == b:
        target_col_s, target_col_b = 0, 1
    elif v_i == b and v_j == s:
        target_col_s, target_col_b = 1, 0
    else:
        raise ValueError(f"Target positions {v_i},{v_j} don't match s={s}, b={b}")
    
    source_a = source_phas[:, source_col_a]
    source_s = source_phas[:, source_col_s]
    target_s = target_phas[:, target_col_s]
    target_b = target_phas[:, target_col_b]
    
    for perm in permutations(range(K)):
        perm = list(perm)
        target_s_perm = target_s[perm]
        target_b_perm = target_b[perm]
        
        if np.array_equal(source_s, target_s_perm):
            triple = np.stack([source_a, source_s, target_b_perm], axis=1)
            matched.append(triple)
    
    return matched


def _str_to_phasing_matrix(s: str, K: int, n_cols: Optional[int] = None) -> np.ndarray:
    """Convert phasing string to matrix."""
    arr = np.array([int(x) for x in s], dtype=np.int8)
    if n_cols is None:
        n_cols = len(arr) // K
    return arr.reshape(K, n_cols)


def _node_label(col_i: int, col_j: int, col_to_snp: np.ndarray) -> str:
    """Convert columns to 1-based node label."""
    pos_i = int(col_to_snp[col_i]) + 1
    pos_j = int(col_to_snp[col_j]) + 1
    if pos_i > pos_j:
        pos_i, pos_j = pos_j, pos_i
    return f"{pos_i}-{pos_j}"


def _parse_node(lbl: str) -> Tuple[int, int]:
    """Parse node label."""
    parts = lbl.split('-')
    return int(parts[0]), int(parts[1])



def compute_forward_messages_simple(
    slices: Dict[int, List[str]],
    edges: List[Tuple[str, str]],
    emission_dict: Dict[str, Dict[str, Dict[str, float]]],  # node -> phase -> {'00','01','10','11': prob}
    transitions_dict: Dict[str, np.ndarray],
    M: csr_matrix,
    snp_to_col: Dict[int, int],
    add_uniform_prior: bool = True
) -> Dict[int, Dict[str, Dict[str, float]]]:
    """
    Forward messages using M directly - NO assignment_dict needed!
    
    For each node:
    1. Look at M to find reads covering the node's positions
    2. Count pattern occurrences directly
    3. Compute emission = product over reads
    
    Simple and correct!
    """
    
    # Build parent adjacency
    parents = defaultdict(list)
    for u, v in edges:
        parents[v].append(u)
    
    # Precompute log emission vectors (00, 01, 10, 11 order)
    node_state_order = {}
    node_log_emisvec = {}
    
    for node, phases in emission_dict.items():
        state_order = list(phases.keys())
        node_state_order[node] = state_order
        node_log_emisvec[node] = {}
        
        for phase in state_order:
            p = phases[phase]
            log_vec = np.log([
                float(p['00']),
                float(p['01']),
                float(p['10']),
                float(p['11']),
            ])
            node_log_emisvec[node][phase] = log_vec
    
    # Initialize forward messages
    forward_messages = {t: {} for t in slices}
    slice_keys = sorted(slices.keys())
    
    if not slice_keys:
        return forward_messages
    
    # print("\nComputing forward messages...")
    
    # Base case (t = 1)
    t0 = slice_keys[0]
    # print(f"\nSlice {t0} (base case):")
    
    for node in slices[t0]:
        forward_messages[t0][node] = {}
        state_order = node_state_order[node]
        
        # Get pattern counts DIRECTLY from M
        cnt4 = _get_pattern_counts_from_M(node, M, snp_to_col)
        
        # print(f"  Node {node}: {cnt4.sum():.0f} reads, patterns {cnt4}")
        
        for phase in state_order:
            # Emission: sum of log probs weighted by counts
            log_emit = float(np.dot(cnt4, node_log_emisvec[node][phase]))
            
            # Uniform prior
            if add_uniform_prior and len(state_order) > 0:
                log_emit += -math.log(len(state_order))
            
            forward_messages[t0][node][phase] = log_emit
    
    # Recursion (t = 2, ..., T)
    node_to_slice = {n: t for t in slice_keys for n in slices[t]}
    
    for t in slice_keys[1:]:
        # print(f"\nSlice {t}:")
        
        for node in slices[t]:
            forward_messages[t][node] = {}
            state_order = node_state_order[node]
            S_child = len(state_order)
            
            if S_child == 0:
                continue
            
            # Get pattern counts DIRECTLY from M
            cnt4 = _get_pattern_counts_from_M(node, M, snp_to_col)
            
            # Emission term (same for all parents)
            log_emit_vec = np.array([
                np.dot(cnt4, node_log_emisvec[node][ph]) 
                for ph in state_order
            ], dtype=np.float64)
            
            if add_uniform_prior and S_child > 0:
                log_emit_vec += -math.log(S_child)
            
            # Aggregate over parents
            log_trans_accum = np.full(S_child, -np.inf, dtype=np.float64)
            
            for pa in parents.get(node, []):
                pt = node_to_slice.get(pa, None)
                if pt is None or pt >= t:
                    continue
                
                pa_order = node_state_order[pa]
                S_parent = len(pa_order)
                if S_parent == 0:
                    continue
                
                # Parent forward messages
                alpha_pa = np.array([
                    forward_messages[pt][pa][ph] 
                    for ph in pa_order
                ], dtype=np.float64)
                
                # Transition matrix
                key = f"{pa}--{node}"
                T = transitions_dict.get(key, None)
                
                if T is None:
                    # Try reversed
                    T_rev = transitions_dict.get(f"{node}--{pa}", None)
                    if T_rev is None:
                        continue
                    # Convert P(child|parent) from P(parent|child)
                    col_sum = T_rev.sum(axis=0, keepdims=True).clip(min=1e-12)
                    T = (T_rev / col_sum).T
                
                if T.shape != (S_parent, S_child):
                    raise ValueError(f"Transition shape mismatch for {key}")
                
                # Compute log-sum-exp
                with np.errstate(divide='ignore'):
                    logT = np.log(np.clip(T, 1e-300, 1.0))
                
                tmp = alpha_pa[:, None] + logT  # (S_parent, S_child)
                m = np.max(tmp, axis=0)
                contrib = m + np.log(np.sum(np.exp(tmp - m), axis=0))
                log_trans_accum = np.logaddexp(log_trans_accum, contrib)
            
            # Final forward message
            log_alpha = log_emit_vec + log_trans_accum
            
            for j, ph in enumerate(state_order):
                forward_messages[t][node][ph] = float(log_alpha[j])
            
            # print(f"  Node {node}: {cnt4.sum():.0f} reads")
    
    # print("\nForward messages computed!")
    return forward_messages


def _get_pattern_counts_from_M(
    node: str,
    M: csr_matrix,
    snp_to_col: Dict[int, int]
) -> np.ndarray:
    """
    Get pattern counts [n_00, n_01, n_10, n_11] directly from M.
    
    This is the KEY simplification - no assignment_dict needed!
    """
    # Parse node positions
    positions = [int(x) for x in node.split('-')]
    cols = [snp_to_col[p - 1] for p in positions]
    
    # Extract submatrix for these columns
    sub = M[:, cols].toarray()  # (n_reads x 2)
    
    # Only keep reads covering both positions
    complete = (sub > 0).all(axis=1)
    
    if not complete.any():
        # No reads - return zeros
        return np.zeros(4, dtype=np.int32)
    
    # Convert to {0,1} (REF=0, ALT=1)
    observations = (sub[complete] == 2).astype(np.int8)
    
    # Count patterns
    counts = np.zeros(4, dtype=np.int32)
    for obs in observations:
        pattern_idx = obs[0] * 2 + obs[1]  # 00→0, 01→1, 10→2, 11→3
        counts[pattern_idx] += 1
    
    return counts









