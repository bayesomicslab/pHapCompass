import numpy as np
from collections import defaultdict
import math
from typing import Dict, List, Tuple
from scipy.sparse import csr_matrix


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

def viterbi_decode(
    slices: Dict[int, List[str]],
    edges: List[Tuple[str, str]],
    assignment_dict: Dict[str, Dict[str, List[int]]],
    emission_dict: Dict[str, Dict[str, Dict[str, float]]],
    transitions_dict: Dict[str, np.ndarray],
    M: csr_matrix,
    snp_to_col: Dict[int, int],
    add_uniform_prior: bool = True
) -> Dict[int, Dict[str, str]]:
    """
    Viterbi algorithm: Find the most probable state sequence (MAP estimation).
    
    Uses the SAME emissions and forward computation as FFBS, but:
    - Forward pass: max instead of log-sum-exp
    - Backward pass: deterministic backtracking instead of sampling
    
    Returns:
        best_states[t][node] = best_phasing_string
    """
        
    def _node_cols_from_label(node_lbl: str, snp_to_col: Dict[int,int]) -> Tuple[int,int]:
        """'{i}-{j}' (i<j, 1-based) -> (col_i, col_j) 0-based matrix columns."""
        i_str, j_str = node_lbl.split("-")
        i = int(i_str); j = int(j_str)
        if i > j: i, j = j, i
        return snp_to_col[i-1], snp_to_col[j-1]
    
    # Build parent adjacency
    parents = defaultdict(list)
    for u, v in edges:
        parents[v].append(u)
    
    # Precompute log emissions (same as forward messages)
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
    
    # Precompute pattern counts (same as forward messages)
    node_counts4 = {}
    for node in {n for L in slices.values() for n in L}:
        rows = np.array(assignment_dict['states'].get(node, []), dtype=np.int64)
        ci, cj = _node_cols_from_label(node, snp_to_col)
        cnt4 = _pattern_counts_for_node(M, rows, ci, cj)
        node_counts4[node] = cnt4
    
    # Viterbi tables
    delta = {}  # delta[t][node][state] = max log prob to reach this state
    psi = {}    # psi[t][node][state] = (parent_node, parent_state) that gave max
    
    slice_keys = sorted(slices.keys())
    if not slice_keys:
        return {}
    
    # Initialize
    for t in slice_keys:
        delta[t] = {}
        psi[t] = {}
        for node in slices[t]:
            delta[t][node] = {}
            psi[t][node] = {}
    
    # print("\nRunning Viterbi algorithm...")
    
    # Base case (t = 1)
    t0 = slice_keys[0]
    for node in slices[t0]:
        state_order = node_state_order[node]
        cnt4 = node_counts4[node]
        
        for phase in state_order:
            log_emit = float(np.dot(cnt4, node_log_emisvec[node][phase]))
            if add_uniform_prior and len(state_order) > 0:
                log_emit += -math.log(len(state_order))
            
            delta[t0][node][phase] = log_emit
            psi[t0][node][phase] = None  # No parent
    
    # Recursion (t = 2, ..., T)
    node_to_slice = {n: t for t in slice_keys for n in slices[t]}
    
    for t in slice_keys[1:]:
        for node in slices[t]:
            state_order = node_state_order[node]
            S_child = len(state_order)
            
            if S_child == 0:
                continue
            
            # Emission term (same as forward)
            cnt4 = node_counts4[node]
            log_emit_vec = np.array([
                np.dot(cnt4, node_log_emisvec[node][ph]) 
                for ph in state_order
            ], dtype=np.float64)
            
            if add_uniform_prior and S_child > 0:
                log_emit_vec += -math.log(S_child)
            
            # For each child state, find best parent
            for j, child_state in enumerate(state_order):
                log_emit = log_emit_vec[j]
                
                # Find max over parents
                best_log_prob = -np.inf
                best_parent = None
                best_parent_state = None
                
                for pa in parents.get(node, []):
                    pt = node_to_slice.get(pa, None)
                    if pt is None or pt >= t:
                        continue
                    
                    pa_order = node_state_order[pa]
                    S_parent = len(pa_order)
                    if S_parent == 0:
                        continue
                    
                    # Get transition matrix
                    key = f"{pa}--{node}"
                    T = transitions_dict.get(key, None)
                    
                    if T is None:
                        T_rev = transitions_dict.get(f"{node}--{pa}", None)
                        if T_rev is None:
                            continue
                        col_sum = T_rev.sum(axis=0, keepdims=True).clip(min=1e-12)
                        T = (T_rev / col_sum).T
                    
                    if T.shape != (S_parent, S_child):
                        continue
                    
                    # For each parent state, compute: delta[parent_state] + log(transition)
                    with np.errstate(divide='ignore'):
                        logT = np.log(np.clip(T[:, j], 1e-300, 1.0))  # Column j for child state j
                    
                    for i, pa_state in enumerate(pa_order):
                        prev_delta = delta[pt][pa][pa_state]
                        candidate_prob = prev_delta + logT[i]
                        
                        if candidate_prob > best_log_prob:
                            best_log_prob = candidate_prob
                            best_parent = pa
                            best_parent_state = pa_state
                
                # Update delta and psi
                delta[t][node][child_state] = log_emit + best_log_prob
                psi[t][node][child_state] = (best_parent, best_parent_state)
    
    # print("Viterbi forward pass complete. Backtracking...")
    
    # Backtracking to find best path
    best_states = {}
    
    # Find best state at final slice
    t_final = slice_keys[-1]
    best_states[t_final] = {}
    
    for node in slices[t_final]:
        state_order = node_state_order[node]
        if not state_order:
            best_states[t_final][node] = None
            continue
        
        # Find state with maximum delta
        best_state = max(state_order, key=lambda s: delta[t_final][node][s])
        best_states[t_final][node] = best_state
    
    # Backtrack through previous slices
    for idx in range(len(slice_keys) - 2, -1, -1):
        t = slice_keys[idx]
        t_next = slice_keys[idx + 1]
        best_states[t] = {}
        
        # Build children mapping for this slice
        children_map = defaultdict(list)
        slice_next_set = set(slices[t_next])
        for u, v in edges:
            if u in slices[t] and v in slice_next_set:
                children_map[u].append(v)
        
        for node in slices[t]:
            # Check if this node has children in next slice
            node_children = children_map.get(node, [])
            
            if not node_children:
                # No children - take best state at this node
                state_order = node_state_order[node]
                if state_order:
                    best_state = max(state_order, key=lambda s: delta[t][node][s])
                    best_states[t][node] = best_state
                else:
                    best_states[t][node] = None
                continue
            
            # Follow backpointers from children
            found = False
            for child in node_children:
                if child not in best_states[t_next]:
                    continue
                
                child_best_state = best_states[t_next][child]
                backpointer = psi[t_next][child][child_best_state]
                
                if backpointer and backpointer[0] == node:
                    best_states[t][node] = backpointer[1]
                    found = True
                    break
            
            if not found:
                # Fallback: take best state
                state_order = node_state_order[node]
                if state_order:
                    best_state = max(state_order, key=lambda s: delta[t][node][s])
                    best_states[t][node] = best_state
                else:
                    best_states[t][node] = None
    
    # print("Viterbi decoding complete!")
    
    return best_states