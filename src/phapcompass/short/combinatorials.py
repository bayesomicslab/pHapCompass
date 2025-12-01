from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Tuple, Iterable, Optional, Literal
import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix, coo_matrix


# ---------------------------------------------------------------------
# Utilities: stable logs, Kronecker error channels
# ---------------------------------------------------------------------

def _log_clip(p: NDArray, eps: float = 1e-12) -> NDArray:
    """Safe log with floor."""
    return np.log(np.clip(p, eps, 1.0))


def error_channel_pair(e: float) -> NDArray:
    """
    4x4 error mixing for a SNP-pair under independent flips.
    State order is (00, 01, 10, 11). Flip matrix for one SNP:
        F = [[1-e, e],
             [  e, 1-e]]
    Pair channel = F ⊗ F  (Kronecker).
    """
    F = np.array([[1.0 - e, e],
                  [e, 1.0 - e]], dtype=np.float64)
    return np.kron(F, F)  # (4,4)


def error_channel_triple(e: float) -> NDArray:
    """8x8 error mixing for a 3-SNP block, order (000..111)."""
    F = np.array([[1.0 - e, e],
                  [e, 1.0 - e]], dtype=np.float64)
    return np.kron(F, np.kron(F, F))  # (8,8)


# ---------------------------------------------------------------------
# 1) Pair-node construction and 4-pattern counts from a sparse matrix
# ---------------------------------------------------------------------

@dataclass
class PairLayer:
    """
    Matrix representation of the quotient layer (pair-nodes).

    nodes:  (P x 2) SNP column indices (i<j)
    counts4: (P x 4) [n00, n01, n10, n11] observed counts from reads
    edges:  (E x 2) indices of pair-nodes that share exactly one SNP
    trip_index: (E x 3) SNP triples (i,j,k) aligned with 'edges' where edges connect (i,j) -- (j,k)
    """
    nodes: NDArray[np.int32]
    counts4: NDArray[np.int32]
    edges: NDArray[np.int32]
    trip_index: NDArray[np.int32]     # (E x 3), columns refer to SNP column indices in M


def build_pair_layer(M: csr_matrix, min_cocov: int = 1) -> PairLayer:
    """
    Build standard quotient layer (pair-nodes) using sparse linear algebra only.

    M: (n_reads x n_snps) values in {0,1,2} (0=no-call, 1=REF, 2=ALT)
    """
    assert M.ndim == 2
    n_reads, n_snps = M.shape

    # Indicators
    R = (M == 1).astype(np.int8)
    A = (M == 2).astype(np.int8)
    O = R + A  # any call

    # Co-coverage (sparse)
    C = (O.T @ O).tocsr()
    C.setdiag(0)
    C = C.tocoo()
    keep = C.data >= min_cocov
    ii = C.row[keep]
    jj = C.col[keep]
    mask = ii < jj
    ii, jj = ii[mask], jj[mask]
    pairs = np.stack([ii, jj], axis=1)
    if pairs.size == 0:
        return PairLayer(np.empty((0,2), np.int32),
                         np.empty((0,4), np.int32),
                         np.empty((0,2), np.int32),
                         np.empty((0,3), np.int32))

    # unique sorted
    order = np.lexsort((pairs[:,1], pairs[:,0]))
    pairs = pairs[order]
    uniq = np.ones(len(pairs), dtype=bool)
    uniq[1:] = (pairs[1:] != pairs[:-1]).any(axis=1)
    nodes = pairs[uniq].astype(np.int32)  # (P,2)

    # Vectorized 4-pattern counts:
    RR = (R.T @ R).tocsr()  # both REF
    AA = (A.T @ A).tocsr()  # both ALT
    RA = (R.T @ A).tocsr()  # i=REF, j=ALT
    AR = (A.T @ R).tocsr()  # i=ALT, j=REF

    r = nodes[:,0]
    c = nodes[:,1]
    n00 = RR[r, c].A.ravel()
    n11 = AA[r, c].A.ravel()
    n01 = RA[r, c].A.ravel()
    n10 = AR[r, c].A.ravel()
    counts4 = np.stack([n00, n01, n10, n11], axis=1).astype(np.int32)  # (P,4)

    # Quotient edges via incidence:
    P = nodes.shape[0]
    rows = np.repeat(np.arange(P, dtype=np.int32), 2)
    cols = nodes.ravel()
    data = np.ones(2*P, dtype=np.int8)
    Inc = coo_matrix((data, (rows, cols)), shape=(P, n_snps)).tocsr()
    S = (Inc @ Inc.T).tocsr()
    S.setdiag(0)
    S = S.tocoo()
    sel = (S.data == 1)
    u = S.row[sel]
    v = S.col[sel]
    keep = u < v
    edges = np.stack([u[keep], v[keep]], axis=1).astype(np.int32)

    # Map each edge (u,v) into the corresponding triple (i,j,k) w.r.t. the shared SNP j
    # nodes[u]=(i,j) or (j,i), nodes[v]=(j,k) or (k,j); find the shared column and order as (i,j,k)
    triplets = np.empty((edges.shape[0], 3), dtype=np.int32)
    for idx, (uu, vv) in enumerate(edges):
        a = nodes[uu]
        b = nodes[vv]
        shared = np.intersect1d(a, b)
        # By construction exactly one shared SNP:
        j = int(shared[0])
        i = int(a[0] if a[1] == j else a[1])  # the non-shared in a
        k = int(b[0] if b[1] == j else b[1])  # the non-shared in b
        triplets[idx] = (i, j, k)

    return PairLayer(nodes=nodes, counts4=counts4, edges=edges, trip_index=triplets)


# ---------------------------------------------------------------------
# 2) Emission models: enumerate genotype-consistent phasings, score LL
# ---------------------------------------------------------------------

@dataclass
class NodePhasingModel:
    """
    For a fixed (g_i, g_j, K) genotype pair:
      n_states = number of unique pair phasings (parameterized by x = n11)
      t_counts: (n_states x 4) true haplotype count vector [t00,t01,t10,t11]
      p_obs:    (n_states x 4) observed pair distribution after error channel
      log_p:    (n_states x 4) log p_obs (cached)
      x_vals:   (n_states,)    the n11 values used for these states
    """
    g_i: int
    g_j: int
    K: int
    t_counts: NDArray[np.int16]
    p_obs: NDArray[np.float64]
    log_p: NDArray[np.float64]
    x_vals: NDArray[np.int16]


def _enumerate_pair_states(K: int, g_i: int, g_j: int) -> Tuple[NDArray, NDArray]:
    """
    Enumerate unique 2-SNP phasing count vectors t=[t00,t01,t10,t11] given ploidy K and genotypes g_i,g_j.
    Parameterize by x = t11 in [max(0, g_i+g_j-K), min(g_i,g_j)].
    """
    lo = max(0, g_i + g_j - K)
    hi = min(g_i, g_j)
    x = np.arange(lo, hi + 1, dtype=np.int64)  # n11
    t11 = x
    t10 = (g_i - x).astype(np.int64)
    t01 = (g_j - x).astype(np.int64)
    t00 = (K - g_i - g_j + x).astype(np.int64)
    t = np.stack([t00, t01, t10, t11], axis=1)  # (S,4)
    return x, t


def build_node_phasing_bank(K: int, e: float,
                            genotype_pairs: Iterable[Tuple[int,int]]) -> Dict[Tuple[int,int], NodePhasingModel]:
    """
    Precompute, for every (g_i,g_j) that appears, the phasing states and their (log) observation probabilities.
    """
    E2 = error_channel_pair(e)  # (4,4)
    bank: Dict[Tuple[int,int], NodePhasingModel] = {}
    for (gi, gj) in sorted(set(genotype_pairs)):
        x, t = _enumerate_pair_states(K, gi, gj)        # (S,),(S,4)
        p_true = t.astype(np.float64) / float(K)        # (S,4)
        p_obs = (p_true @ E2.T)                         # (S,4)
        bank[(gi, gj)] = NodePhasingModel(
            g_i=gi, g_j=gj, K=K,
            t_counts=t,
            p_obs=p_obs,
            log_p=_log_clip(p_obs)
        ,   x_vals=x
        )
    return bank


def node_emission_loglik(
    counts4: NDArray[np.int32],
    gi: NDArray[np.int16],
    gj: NDArray[np.int16],
    bank: Dict[Tuple[int,int], NodePhasingModel]
) -> Tuple[List[NDArray[np.float64]], List[NDArray[np.int16]]]:
    """
    For each node p with genotype (gi[p],gj[p]) and observed counts4[p],
    return its vector of log-likelihoods over the states in bank[(gi,gj)],
    and the corresponding x (n11) values.

    Output:
      ll_list[p]: (S_p,) log-likelihoods
      x_list[p]:  (S_p,) the n11 values indexing the states
    """
    P = counts4.shape[0]
    ll_list: List[NDArray[np.float64]] = [None] * P
    x_list:  List[NDArray[np.int16]]   = [None] * P

    # Group nodes by genotype pair to batch-multiply
    key = np.stack([gi, gj], axis=1)
    # pack into tuples for hashing
    key_tuples = [tuple(map(int, k)) for k in key]
    unique_keys = sorted(set(key_tuples))
    for k in unique_keys:
        model = bank[k]
        idx = [i for i, kk in enumerate(key_tuples) if kk == k]
        C = counts4[idx]               # (N_k, 4)
        # ll = sum_k c_k * log p_obs_k  ==> (N_k x 4) @ (4 x S) -> (N_k x S)
        # but model.log_p is (S x 4), so we do: C @ log_p^T
        ll = C @ model.log_p.T         # (N_k, S)
        for j, row in zip(idx, ll):
            ll_list[j] = row
            x_list[j]  = model.x_vals
    return ll_list, x_list


# ---------------------------------------------------------------------
# 3) Triple counts per edge and log-transition matrices (vectorized)
#     Model: i ⟂ k | j  (conditional independence given shared SNP)
# ---------------------------------------------------------------------

@dataclass
class EdgeTripleLayer:
    """
    Triple count tensors and per-edge transition matrices.

    counts8:  (E x 8) observed triple counts for each (i,j,k) from reads
    log_A:    list[E] of (S_ij x S_jk) log transition matrices
    """
    counts8: NDArray[np.int32]
    log_A: List[NDArray[np.float64]]


def _triple_count_mats(M: csr_matrix) -> Tuple[List[csr_matrix], List[csr_matrix]]:
    """
    Return R, A indicator matrices and fast access to column vectors.
    """
    R = (M == 1).astype(np.int8).tocsr()
    A = (M == 2).astype(np.int8).tocsr()
    return R, A


def _counts8_for_j(R: csr_matrix, A: csr_matrix, j: int) -> Tuple[csr_matrix, csr_matrix, csr_matrix, csr_matrix,
                                                                  csr_matrix, csr_matrix, csr_matrix, csr_matrix]:
    """
    For a fixed middle SNP 'j', compute eight (n_snps x n_snps) sparse matrices:
      n000 = R^T @ (R .* R[:,j])   etc.
    where (X .* v) is row-wise multiply by column vector v (mask by reads with X[:,j]=1).

    Ordering of returned matrices matches pattern bits (i,j,k) from 000..111.
    """
    rj = R[:, j]       # (n_reads x 1) sparse col
    aj = A[:, j]

    # Row-wise mask by middle j
    RRj = R.multiply(rj)  # D_rj @ R
    ARj = A.multiply(rj)
    RAj = R.multiply(aj)
    AAj = A.multiply(aj)

    n000 = (R.T @ RRj).tocsr()
    n001 = (R.T @ ARj).tocsr()
    n100 = (A.T @ RRj).tocsr()
    n101 = (A.T @ ARj).tocsr()

    n010 = (R.T @ RAj).tocsr()
    n011 = (R.T @ AAj).tocsr()
    n110 = (A.T @ RAj).tocsr()
    n111 = (A.T @ AAj).tocsr()

    return n000, n001, n010, n011, n100, n101, n110, n111


def build_edge_triples_and_transitions(
    M: csr_matrix,
    layer: PairLayer,
    K: int,
    e: float,
    # Per-node phasing models so we can read (gi,gj) state counts quickly
    node_genotypes: NDArray[np.int16],  # (n_snps,) genotype alt-count per SNP
    node_state_bank: Dict[Tuple[int,int], NodePhasingModel],
    node_state_keys: List[Tuple[int,int]],    # per node p → (g_i,g_j)
    node_state_x:  List[NDArray[np.int16]]    # per node p → x-values for its states (n11)
) -> EdgeTripleLayer:
    """
    For each quotient edge (u,v) joining nodes (i,j) and (j,k):
      1) Compute the 8-pattern observed triple counts c8(i,j,k) from M (vectorized per 'j').
      2) Build a (S_ij x S_jk) log transition matrix using the conditional-independence model i ⟂ k | j:
         - From pair states (t00,t01,t10,t11) on (i,j) and (j,k), derive triple *true* distribution by
           splitting the j=0 and j=1 groups independently:
             K0 = K - g_j,  p_i1|j0 = t10/K0,  p_k1|j0 = m01/K0
             K1 = g_j,      p_i1|j1 = t11/K1,  p_k1|j1 = m11/K1
           → expected triple counts over 8 cells (000..111).
         - Convert to probabilities by /K and apply error channel (8x8).
         - Score c8 with multinomial log-likelihood: sum_u c8[u] * log p_obs[u].
    """
    R, A = _triple_count_mats(M)
    E3 = error_channel_triple(e)

    n_snps = M.shape[1]
    P = layer.nodes.shape[0]
    E = layer.edges.shape[0]

    # Precompute per-middle-SNP the 8 matrices (n_snps x n_snps)
    # We'll fetch (i,k) entry from them for each edge (i,j,k)
    per_j_mats: Dict[int, Tuple[csr_matrix,...]] = {}

    counts8 = np.zeros((E, 8), dtype=np.int32)
    log_A: List[NDArray[np.float64]] = [None] * E

    # Group edges by middle SNP 'j' to avoid recompute
    trip = layer.trip_index  # (E,3)
    # group indices by j
    by_j: Dict[int, List[int]] = {}
    for edx, (_, j, _) in enumerate(trip):
        by_j.setdefault(int(j), []).append(edx)

    for j, eidx in by_j.items():
        mats = per_j_mats.get(j)
        if mats is None:
            mats = _counts8_for_j(R, A, j)
            per_j_mats[j] = mats

        m000, m001, m010, m011, m100, m101, m110, m111 = mats

        # Fill counts8 for all edges with this middle 'j'
        idx = np.array(eidx, dtype=np.int32)
        ik = trip[idx][:, [0, 2]]  # (i,k)
        i_cols = ik[:, 0]
        k_cols = ik[:, 1]

        counts8[idx, 0] = m000[i_cols, k_cols].A.ravel()
        counts8[idx, 1] = m001[i_cols, k_cols].A.ravel()
        counts8[idx, 2] = m010[i_cols, k_cols].A.ravel()
        counts8[idx, 3] = m011[i_cols, k_cols].A.ravel()
        counts8[idx, 4] = m100[i_cols, k_cols].A.ravel()
        counts8[idx, 5] = m101[i_cols, k_cols].A.ravel()
        counts8[idx, 6] = m110[i_cols, k_cols].A.ravel()
        counts8[idx, 7] = m111[i_cols, k_cols].A.ravel()

        # Build transition matrices for these edges
        for ed in idx:
            u, v = layer.edges[ed]      # node indices
            i, j2, k = trip[ed]
            assert j2 == j

            # Pull node models
            gi, gj = node_state_keys[u]
            gj2, gk = node_state_keys[v]  # gj2 == gj (same middle SNP)
            model_ij = node_state_bank[(gi, gj)]
            model_jk = node_state_bank[(gj, gk)]
            x_u = node_state_x[u]   # (S_u,)
            x_v = node_state_x[v]   # (S_v,)

            # Precompute constants
            K0 = float(K - gj)
            K1 = float(gj)
            # For each state, extract t-vecs:
            # t_ij = [t00,t01,t10,t11] for (i,j)
            # t_jk = [u00,u01,u10,u11] for (j,k), note indexing matches (j first)
            T_ij = model_ij.t_counts.astype(np.float64)   # (S_u,4)
            U_jk = model_jk.t_counts.astype(np.float64)   # (S_v,4)

            # Within j=0 group (size K0):
            # p_i1|j0 = t10/K0, p_k1|j0 = u01/K0
            # Within j=1 group (size K1):
            # p_i1|j1 = t11/K1, p_k1|j1 = u11/K1
            # Build expected triple *counts* under i ⟂ k | j:
            # For states u in S_u, v in S_v, we want an (S_u x S_v x 8) tensor.
            # We can do this by outer-products over the j=0 and j=1 blocks.

            # Extract columns:
            t00, t01, t10, t11 = [T_ij[:,k_] for k_ in range(4)]
            u00, u01, u10, u11 = [U_jk[:,k_] for k_ in range(4)]

            # Avoid div by zero
            with np.errstate(divide='ignore', invalid='ignore'):
                pi1_j0 = np.where(K0 > 0, (t10 / K0)[:,None], 0.0)  # (S_u,1)
                pk1_j0 = np.where(K0 > 0, (u01 / K0)[None,:], 0.0)  # (1,S_v)
                pi1_j1 = np.where(K1 > 0, (t11 / K1)[:,None], 0.0)
                pk1_j1 = np.where(K1 > 0, (u11 / K1)[None,:], 0.0)

            # j=0 block expected counts (S_u x S_v):
            # 000: (1-pi)*(1-pk)*K0; 001: (1-pi)*pk*K0; 100: pi*(1-pk)*K0; 101: pi*pk*K0
            one_u = 1.0 - pi1_j0      # (S_u,1)
            one_v = 1.0 - pk1_j0      # (1,S_v)
            c000 = (one_u @ one_v) * K0
            c001 = (one_u @ pk1_j0) * K0
            c100 = (pi1_j0 @ one_v) * K0
            c101 = (pi1_j0 @ pk1_j0) * K0

            # j=1 block expected counts:
            # 010: (1-pi)*(1-pk)*K1; 011: (1-pi)*pk*K1; 110: pi*(1-pk)*K1; 111: pi*pk*K1
            one_u1 = 1.0 - pi1_j1
            one_v1 = 1.0 - pk1_j1
            c010 = (one_u1 @ one_v1) * K1
            c011 = (one_u1 @ pk1_j1) * K1
            c110 = (pi1_j1 @ one_v1) * K1
            c111 = (pi1_j1 @ pk1_j1) * K1

            # Stack into (S_u, S_v, 8)
            # Order: 000,001,010,011,100,101,110,111
            C_true = np.stack([c000, c001, c010, c011, c100, c101, c110, c111], axis=2)  # (Su,Sv,8)

            # Convert to probabilities by /K and apply error channel
            P_true = C_true / float(K)             # (Su,Sv,8)
            # (8,) <- (8,8) @ (8,)
            # Do mixing by right-multiplying with E3^T: p_obs = p_true @ E3^T
            P_obs = P_true @ E3.T                  # (Su,Sv,8)
            logP = _log_clip(P_obs)                # (Su,Sv,8)

            # Score observed counts for this edge
            c8 = counts8[ed].astype(np.float64)    # (8,)
            # ll(u,v) = sum_u c8[u] * logP(u,v,u)
            # Broadcast multiply-sum over last axis
            LL = (logP * c8[None, None, :]).sum(axis=2)   # (Su,Sv)
            log_A[ed] = LL

    return EdgeTripleLayer(counts8=counts8, log_A=log_A)


def adapt_emissions_from_ll_nodes(
    pair_layer,
    ll_nodes,
    gi, gj,
    bank,
    one_based_labels=True,
    col_to_snp: List[int] = None,   # <<< NEW (required)
):
    assert col_to_snp is not None, "Pass frag.col_to_snp"
    labels = [
        _node_label_from_cols(int(i), int(j), col_to_snp)
        for (i, j) in pair_layer.nodes
    ]

    # emissions: node -> (S_p,) log-likelihood vector (same order as ll_nodes)
    log_emissions = {labels[p]: ll_nodes[p] for p in range(len(labels))}

    node_key = list(zip(gi.tolist(), gj.tolist()))
    state_t = {labels[p]: bank[k].t_counts for p, k in enumerate(node_key)}
    state_x = {labels[p]: bank[k].x_vals   for p, k in enumerate(node_key)}
    return log_emissions, state_t, state_x, labels


def adapt_transitions_from_edge_layer(
    pair_layer,
    edge_layer,
    col_to_snp: List[int],  # <<< REQUIRED
):
    """
    Transitions keyed by TOPOLOGICAL 'u--v':
      For edge e with triple (i,j,k) in column indices:
        u = label(i,j)  using positions (1-based)
        v = label(j,k)  using positions (1-based)
      The direction is from smaller "other" to larger "other" around the shared SNP.
    """
    transitions_prob = {}
    transitions_log_raw = {}

    def label_from_cols(ci, cj):
        return _node_label_from_cols(int(ci), int(cj), col_to_snp)

    for e_idx, (_u_idx, _v_idx) in enumerate(pair_layer.edges):
        i_col, j_col, k_col = map(int, pair_layer.trip_index[e_idx])

        # source and target labels in topological order
        u_label = label_from_cols(i_col, j_col)  # (i,j)
        v_label = label_from_cols(j_col, k_col)  # (j,k)

        # convert raw log-scores to row-normalized probabilities
        A_log = edge_layer.log_A[e_idx]  # (S_u x S_v) in bank order
        max_row = np.max(A_log, axis=1, keepdims=True)
        A_exp = np.exp(A_log - max_row)
        row_sum = A_exp.sum(axis=1, keepdims=True)
        with np.errstate(invalid='ignore', divide='ignore'):
            A_prob = np.divide(A_exp, row_sum, where=row_sum > 0)

        zero_rows = (row_sum.squeeze() == 0)
        if np.any(zero_rows):
            Su, Sv = A_prob.shape
            A_prob[zero_rows, :] = 1.0 / max(Sv, 1)

        key = f"{u_label}--{v_label}"
        transitions_prob[key] = A_prob
        transitions_log_raw[key] = A_log

    return transitions_prob, transitions_log_raw



# ──────────────────────────────────────────────────────────────────────────────
# Fragment adjacency from M  (reads×SNPs; entries 0/1/2 → indicator >0)
# Optional: min co-coverage threshold and column-gap window
# ──────────────────────────────────────────────────────────────────────────────

def _fragment_adjacency_from_M(M: csr_matrix,
                               min_cocov: int = 1,
                               max_col_gap: Optional[int] = None) -> csr_matrix:
    O = (M > 0).astype(np.int8)
    C = (O.T @ O).tocsr()
    C.setdiag(0)
    C.eliminate_zeros()

    if max_col_gap is not None:
        rows, cols = C.nonzero()
        mask = (np.abs(rows - cols) <= int(max_col_gap))
        C = coo_matrix((C.data[mask], (rows[mask], cols[mask])), shape=C.shape).tocsr()

    if min_cocov > 1:
        C.data = (C.data >= min_cocov).astype(np.int8)
        C.eliminate_zeros()
    else:
        C.data[:] = 1

    A = ((C + C.T) > 0).astype(np.int8).tocsr()
    A.setdiag(0)
    A.eliminate_zeros()
    return A

# ──────────────────────────────────────────────────────────────────────────────
# MCS → PEO and chordality check
# ──────────────────────────────────────────────────────────────────────────────

def _mcs_peo(A: csr_matrix) -> Tuple[np.ndarray, np.ndarray]:
    A = A.tocsr()
    n = A.shape[0]
    labels = np.zeros(n, dtype=np.int32)
    unnum  = np.ones(n, dtype=bool)
    order  = np.empty(n, dtype=np.int32)

    for t in range(n - 1, -1, -1):
        cand = np.where(unnum)[0]
        v = cand[np.argmax(labels[cand])]
        order[t] = v
        unnum[v] = False
        nbrs = A.indices[A.indptr[v]:A.indptr[v+1]]
        nbrs = nbrs[unnum[nbrs]]
        labels[nbrs] += 1

    pos = np.empty(n, dtype=np.int32)
    pos[order] = np.arange(n, dtype=np.int32)
    return order, pos

def _is_chordal_with_peo(A: csr_matrix, order: np.ndarray, pos: np.ndarray) -> bool:
    """
    Check: for every vertex v, its later neighbors form a clique.
    Equivalent to 'order' being a PEO → chordal graph.
    """
    A = A.tocsr()
    for v in order:
        nbrs = A.indices[A.indptr[v]:A.indptr[v+1]]
        later = nbrs[pos[nbrs] > pos[v]]
        if later.size <= 1:
            continue
        # all pairs among 'later' must be adjacent
        # we can check via degrees into the 'later' set:
        later_set = set(later.tolist())
        for u in later:
            u_nbrs = A.indices[A.indptr[u]:A.indptr[u+1]]
            # count neighbors of u within 'later'
            if len(later_set.intersection(u_nbrs.tolist())) < len(later) - 1:
                return False
    return True

# ──────────────────────────────────────────────────────────────────────────────
# Maximal cliques for chordal graphs (PEO)
# ──────────────────────────────────────────────────────────────────────────────

def _maximal_cliques_from_peo(A: csr_matrix,
                              order: np.ndarray,
                              pos: np.ndarray) -> List[np.ndarray]:
    A = A.tocsr()
    n = A.shape[0]
    candidates = []
    for v in order:
        nbrs = A.indices[A.indptr[v]:A.indptr[v+1]]
        later = nbrs[pos[nbrs] > pos[v]]
        C = np.sort(np.concatenate(([v], later)).astype(np.int32))
        candidates.append(C)

    # dedup identical
    uniq_map = {}
    for C in candidates:
        key = C.tobytes()
        if key not in uniq_map:
            uniq_map[key] = C
    uniq = list(uniq_map.values())
    if len(uniq) <= 1:
        return uniq

    # prune subsets using sparse incidence
    c = len(uniq)
    sizes = np.array([len(C) for C in uniq], dtype=np.int32)
    rows = np.concatenate(uniq)
    cols = np.concatenate([np.full(len(C), j, dtype=np.int32) for j, C in enumerate(uniq)])
    data = np.ones_like(rows, dtype=np.int8)
    B = coo_matrix((data, (rows, cols)), shape=(n, c)).tocsr()

    W = (B.T @ B).tocoo()
    drop = np.zeros(c, dtype=bool)
    for u, v, w in zip(W.row, W.col, W.data):
        if u == v: continue
        if w == sizes[u] and sizes[v] >= sizes[u]:
            drop[u] = True
    cliques = [uniq[i] for i in range(c) if not drop[i]]
    cliques.sort(key=lambda arr: (int(arr.min()), -len(arr)))
    return cliques

# ──────────────────────────────────────────────────────────────────────────────
# General maximal cliques (Bron–Kerbosch with Tomita pivoting, bitsets)
# Works for ANY graph. Output-sensitive; fast on sparse graphs.
# ──────────────────────────────────────────────────────────────────────────────


def _maximal_cliques_general(A: csr_matrix) -> List[np.ndarray]:
    n, adj_bits = _neighbors_bitsets(A)
    R = 0
    P = int((1 << n) - 1)   # ← ensure Python int (works for any n)
    X = 0
    out: List[np.ndarray] = []
    _bron_kerbosch_pivot(adj_bits, R, P, X, out)
    out = [np.sort(c).astype(np.int32) for c in out]
    out.sort(key=lambda arr: (int(arr.min()) if arr.size else -1, -len(arr)))
    return out

# ──────────────────────────────────────────────────────────────────────────────
# Public API: find_maximal_cliques (auto chordal fast-path + general fallback)
# ──────────────────────────────────────────────────────────────────────────────

def find_maximal_cliques(M: csr_matrix,
                         min_cocov: int = 1,
                         max_col_gap: Optional[int] = None) -> List[np.ndarray]:
    """
    Return ALL maximal cliques of the fragment graph derived from M.
    If the graph is chordal (by MCS test), use linear-time PEO method.
    Otherwise, fall back to Bron–Kerbosch with pivoting (bitsets).
    """
    A = _fragment_adjacency_from_M(M, min_cocov=min_cocov, max_col_gap=max_col_gap)
    if A.shape[0] == 0:
        return []
    order, pos = _mcs_peo(A)
    if _is_chordal_with_peo(A, order, pos):
        return _maximal_cliques_from_peo(A, order, pos)
    else:
        return _maximal_cliques_general(A)

# ──────────────────────────────────────────────────────────────────────────────
# Junction graph (clique graph) + maximum-weight spanning tree (junction tree)
# ──────────────────────────────────────────────────────────────────────────────


def _neighbors_bitsets(A: csr_matrix) -> Tuple[int, List[int]]:
    """
    Encode each vertex's neighborhood as a Python int bitset.
    Python ints are arbitrary precision; bit ops are fast in C.
    """
    A = A.tocsr()
    n = A.shape[0]
    bits: List[int] = [0] * n
    for v in range(n):
        mask = 0
        start, end = A.indptr[v], A.indptr[v + 1]
        for u in A.indices[start:end]:
            mask |= (1 << int(u))
        bits[v] = int(mask)  # ensure pure Python int
    return n, bits


def _bron_kerbosch_pivot(adj_bits: List[int], R: int, P: int, X: int,
                         out: List[np.ndarray]) -> None:
    """
    Bron–Kerbosch with Tomita pivoting on Python-int bitsets.
    R, P, X are *plain Python ints* used as bit masks over vertices [0..n-1].
    """
    # Make absolutely sure we're using Python ints (not numpy scalars)
    R = int(R); P = int(P); X = int(X)

    if P == 0 and X == 0:
        # report clique R
        # iterate set bits in R
        verts = []
        r = R
        while r:
            lb = r & -r                # lowbit (power of two)
            i  = (lb.bit_length() - 1) # index of lowbit
            verts.append(i)
            r ^= lb
        out.append(np.array(sorted(verts), dtype=np.int32))
        return

    # choose pivot u ∈ P ∪ X maximizing |P ∩ N(u)|
    U = int(P | X)
    best_u = -1
    best_score = -1
    u_bits = U
    while u_bits:
        lb = u_bits & -u_bits
        u  = lb.bit_length() - 1
        u_bits ^= lb
        score = (P & adj_bits[u]).bit_count()
        if score > best_score:
            best_score = score
            best_u = u

    # iterate over P \ N(u)
    targets = int(P & (~adj_bits[best_u])) if best_u != -1 else int(P)
    while targets:
        lb = targets & -targets
        v  = lb.bit_length() - 1
        targets ^= lb

        _bron_kerbosch_pivot(
            adj_bits,
            int(R | (1 << v)),
            int(P & adj_bits[v]),
            int(X & adj_bits[v]),
            out
        )
        P = int(P & ~(1 << v))
        X = int(X |  (1 << v))


# Build B (n_snps x C) where B[i,c] = 1 iff SNP i is in clique c
def _clique_incidence(n_snps: int, cliques: List[np.ndarray]) -> csr_matrix:
    if not cliques:
        return csr_matrix((n_snps, 0), dtype=np.int8)
    rows = np.concatenate(cliques)
    cols = np.concatenate([np.full(len(C), j, dtype=np.int32) for j, C in enumerate(cliques)])
    data = np.ones_like(rows, dtype=np.int8)
    return coo_matrix((data, (rows, cols)), shape=(n_snps, len(cliques))).tocsr()

def build_clique_graph(cliques: List[np.ndarray],
                       n_snps: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Junction (clique) graph via one sparse product.
    Returns:
      edges_uv : (E x 2) clique indices u<v
      weights  : (E,)    |intersection(C_u, C_v)|
    """
    C = len(cliques)
    if C <= 1:
        return np.empty((0,2), np.int32), np.empty((0,), np.int32)

    B = _clique_incidence(n_snps, cliques)  # (n_snps x C)
    W = (B.T @ B).tocsr()                   # (C x C) intersections
    W.setdiag(0)
    W = W.tocoo()

    mask = (W.row < W.col) & (W.data > 0)
    edges_uv = np.stack([W.row[mask], W.col[mask]], axis=1).astype(np.int32)
    weights  = W.data[mask].astype(np.int32)
    return edges_uv, weights

# Simple union-find for Kruskal
class _UF:
    __slots__=("p","r")
    def __init__(self, n): self.p=np.arange(n,dtype=np.int32); self.r=np.zeros(n,dtype=np.int8)
    def find(self,x):
        while self.p[x]!=x:
            self.p[x]=self.p[self.p[x]]
            x=self.p[x]
        return x
    def union(self,a,b):
        ra, rb = self.find(a), self.find(b)
        if ra==rb: return False
        if self.r[ra]<self.r[rb]: self.p[ra]=rb
        elif self.r[ra]>self.r[rb]: self.p[rb]=ra
        else: self.p[rb]=ra; self.r[ra]+=1
        return True

def maximum_weight_spanning_forest(C: int,
                                   edges_uv: np.ndarray,
                                   weights: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Kruskal (descending weights). Returns a maximum-weight spanning forest:
      mst_edges : (K x 2) edge list (K = C - #components)
      mst_w     : (K,)    weights
    """
    if C <= 1 or edges_uv.size == 0:
        return np.empty((0,2), np.int32), np.empty((0,), np.int32)

    order = np.argsort(-weights, kind="stable")
    uf = _UF(C)
    sel = []
    wsel = []
    for k in order:
        u, v = int(edges_uv[k,0]), int(edges_uv[k,1])
        if uf.union(u, v):
            sel.append([u, v]); wsel.append(int(weights[k]))
            # optional early stop when single component: if len(sel) == C-1: break
    mst_edges = np.array(sel, dtype=np.int32) if sel else np.empty((0,2), np.int32)
    mst_w     = np.array(wsel, dtype=np.int32) if sel else np.empty((0,), np.int32)
    return mst_edges, mst_w



def make_junction_graph(cliques: List[np.ndarray],
                        n_snps: int,
                        tree_mode: Literal["none","forest","auto"] = "auto",
                        weight_mode: Literal["overlap","read_support","hybrid"] = "overlap",
                        alpha: float = 1.0,   # weight for overlap size
                        beta: float  = 1.0,   # weight for read support
                        M: Optional[csr_matrix] = None
                        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Build the clique (junction) graph and an optional maximum-weight spanning forest.

    weight_mode:
      - "overlap":    w(u,v) = |C_u ∩ C_v|   (default; theoretical JT on chordal graphs)
      - "read_support": w(u,v) = #reads covering ALL SNPs in C_u ∪ C_v  (needs M)
      - "hybrid":     w = alpha*|∩| + beta*#reads  (needs M)

    Returns:
      cg_edges  : (E x 2) clique-graph edges (u<v)
      cg_w      : (E,)   weights (int if pure modes; float if hybrid)
      jt_edges  : (K x 2) MWST forest edges (subset of cg_edges; order not guaranteed)
      jt_w      : (K,)    weights on forest edges
    """
    C = len(cliques)
    if C <= 1:
        return (np.empty((0,2), np.int32), np.empty((0,), np.int32),
                np.empty((0,2), np.int32), np.empty((0,), np.int32))

    # overlap edges/weights (|C_u ∩ C_v|)
    # re-use your existing build_clique_graph(cliques, n_snps) if you already have it
    def _overlap_edges_weights():
        B = _clique_incidence(n_snps, cliques)
        W = (B.T @ B).tocsr()
        W.setdiag(0)
        W = W.tocoo()
        mask = (W.row < W.col) & (W.data > 0)
        return (np.stack([W.row[mask], W.col[mask]], axis=1).astype(np.int32),
                W.data[mask].astype(np.int32))

    if weight_mode == "overlap":
        cg_edges, cg_w = _overlap_edges_weights()

    elif weight_mode == "read_support":
        if M is None:
            raise ValueError("weight_mode='read_support' requires M.")
        cg_edges, cg_w = _build_clique_graph_read_support(M, cliques)

    elif weight_mode == "hybrid":
        if M is None:
            raise ValueError("weight_mode='hybrid' requires M.")
        edges_ov, w_ov = _overlap_edges_weights()
        edges_rs, w_rs = _build_clique_graph_read_support(M, cliques)

        # merge on edge keys (u,v)
        # put both into a dict keyed by (u,v)
        d = {}
        for (u,v), w in zip(edges_ov, w_ov):
            d[(int(u), int(v))] = [float(w), 0.0]
        for (u,v), w in zip(edges_rs, w_rs):
            key = (int(u), int(v))
            if key in d:
                d[key][1] = float(w)
            else:
                d[key] = [0.0, float(w)]    # edge has no overlap (rare), still allowed

        cg_edges = np.array(list(d.keys()), dtype=np.int32)
        cg_w     = np.array([alpha*ov + beta*rs for (ov, rs) in d.values()], dtype=np.float64)

    else:
        raise ValueError(f"Unknown weight_mode '{weight_mode}'")

    # Max-weight spanning forest (Kruskal). Reuse your existing maximum_weight_spanning_forest
    jt_edges, jt_w = maximum_weight_spanning_forest(C, cg_edges, cg_w)

    return cg_edges, cg_w, jt_edges, jt_w


def _clique_incidence_reads(M: csr_matrix, cliques: List[np.ndarray]) -> csr_matrix:
    """
    Build R (n_reads x C) where R[r,c] = 1 iff read r has nonzero calls at ALL SNPs in clique c.
    Vectorized via submatrix nnz counts per row.
    """
    O = (M > 0).astype(np.int8).tocsr()   # indicator
    n_reads, _ = O.shape
    C = len(cliques)
    if C == 0:
        return csr_matrix((n_reads, 0), dtype=np.int8)

    rows_all = []
    cols_all = []
    data_all = []

    for j, cols in enumerate(cliques):
        if len(cols) == 0:
            continue
        sub = O[:, cols]                               # (n_reads x |C_j|)
        # rows that cover all SNPs in this clique
        covered = np.where(np.asarray(sub.getnnz(axis=1)).ravel() == len(cols))[0]
        if covered.size:
            rows_all.append(covered.astype(np.int32))
            cols_all.append(np.full(covered.size, j, dtype=np.int32))
            data_all.append(np.ones(covered.size, dtype=np.int8))

    if not rows_all:
        return csr_matrix((n_reads, C), dtype=np.int8)

    rows = np.concatenate(rows_all)
    cols = np.concatenate(cols_all)
    data = np.concatenate(data_all)
    R = coo_matrix((data, (rows, cols)), shape=(n_reads, C)).tocsr()
    return R


def _build_clique_graph_read_support(M: csr_matrix,
                                     cliques: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute weights for the clique graph edges as the number of reads covering ALL SNPs
    in BOTH cliques (i.e., support for the union). Vectorized via R^T R.
    Returns edges (u<v) and weights (int).
    """
    C = len(cliques)
    if C <= 1:
        return np.empty((0,2), np.int32), np.empty((0,), np.int32)

    # R[r,c] = 1 if read r fully covers clique c
    R = _clique_incidence_reads(M, cliques)     # (n_reads x C)
    S = (R.T @ R).tocsr()                       # (C x C): read-support co-coverage of cliques
    S.setdiag(0)
    S = S.tocoo()
    mask = (S.row < S.col) & (S.data > 0)
    edges = np.stack([S.row[mask], S.col[mask]], axis=1).astype(np.int32)
    w_read = S.data[mask].astype(np.int32)
    return edges, w_read


def _node_label_from_cols(col_i: int, col_j: int, col_to_snp: List[int]) -> str:
    """Return 1-based 'i-j' label from 0-based column indices (enforce i<j)."""
    pi = int(col_to_snp[col_i]) + 1  # 1-based
    pj = int(col_to_snp[col_j]) + 1  # 1-based
    if pi > pj:
        pi, pj = pj, pi
    return f"{pi}-{pj}"




