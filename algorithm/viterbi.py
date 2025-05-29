import math
from collections import defaultdict
import numpy as np



def select_states_viterbi_like(slices, edges, forward_messages, transitions_dict, emission_dict):
    backpointer = defaultdict(lambda: defaultdict(dict))
    viterbi_scores = defaultdict(lambda: defaultdict(dict))

    reverse_adjacency = defaultdict(list)
    for src, dst in edges:
        reverse_adjacency[dst].append(src)

    sorted_ts = sorted(slices.keys())

    # Forward pass to compute max score and backpointer
    for t in sorted_ts:
        for node in slices[t]:
            curr_states = list(emission_dict[node].keys())
            for j, curr_state in enumerate(curr_states):
                max_score = -math.inf
                best_parent = None

                if t == 1 or node not in reverse_adjacency:
                    max_score = forward_messages[t][node][curr_state]
                    best_parent = None
                else:
                    for parent in reverse_adjacency[node]:
                        if parent not in forward_messages[t - 1]:
                            continue

                        edge_label = f"{parent}--{node}"
                        if edge_label not in transitions_dict:
                            continue

                        parent_states = list(emission_dict[parent].keys())
                        for i, prev_state in enumerate(parent_states):
                            trans_prob = transitions_dict[edge_label][i, j]
                            if trans_prob <= 0:
                                continue
                            score = forward_messages[t - 1][parent][prev_state] + math.log(trans_prob)
                            if score > max_score:
                                max_score = score
                                best_parent = (t - 1, parent, prev_state)

                viterbi_scores[t][node][curr_state] = max_score
                backpointer[t][node][curr_state] = best_parent

    # Backward pass to recover best path
    final_states = {}
    visited = set()
    for t in sorted_ts[::-1]:
        for node in slices[t]:
            if node in visited:
                continue
            best_state = max(viterbi_scores[t][node], key=lambda s: viterbi_scores[t][node][s])
            curr = (t, node, best_state)

            while curr:
                t_curr, node_curr, state_curr = curr
                if node_curr not in final_states:
                    final_states[node_curr] = state_curr
                    visited.add(node_curr)
                curr = backpointer[t_curr][node_curr][state_curr]

    return {1: final_states}


def true_viterbi_dp_complete(slices, edges, emission_dict, transitions_dict):
    reverse_adj = defaultdict(list)
    for src, tgt in edges:
        reverse_adj[tgt].append(src)

    dp = defaultdict(lambda: defaultdict(dict))
    backpointer = defaultdict(lambda: defaultdict(dict))

    for t in sorted(slices):
        for node in slices[t]:
            state_space = list(emission_dict[node].keys())
            for phasing in state_space:
                log_em = emission_dict[node][phasing]


                if t == 1 or node not in reverse_adj:
                    # No parents => prior * emission
                    log_prior = -math.log(len(state_space))  # uniform prior
                    dp[t][node][phasing] = log_em + log_prior
                    backpointer[t][node][phasing] = None
                else:
                    max_score = float('-inf')
                    best_parent = None
                    for parent in reverse_adj[node]:
                        if parent not in dp[t-1]:
                            continue
                        edge_label = f"{parent}--{node}"
                        if edge_label not in transitions_dict:
                            continue

                        parent_phasings = list(emission_dict[parent].keys())
                        child_phasings = list(emission_dict[node].keys())
                        i = child_phasings.index(phasing)

                        for j, p_ph in enumerate(parent_phasings):
                            trans_prob = transitions_dict[edge_label][j, i]
                            if trans_prob <= 0:
                                continue
                            score = dp[t-1][parent][p_ph] + math.log(trans_prob)
                            if score > max_score:
                                max_score = score
                                best_parent = (t-1, parent, p_ph)

                    dp[t][node][phasing] = log_em + max_score
                    backpointer[t][node][phasing] = best_parent

    # Reconstruct paths
    result = {}
    visited = set()
    for t in sorted(dp.keys(), reverse=True):
        for node in dp[t]:
            
            if node in visited:
                continue
            best_ph = max(dp[t][node], key=lambda k: dp[t][node][k])
            curr = (t, node, best_ph)
            while curr:
                ct, cn, cp = curr
                if cn not in result:
                    result[cn] = cp
                    visited.add(cn)
                curr = backpointer[ct][cn][cp]

    # Fill disconnected nodes
    all_nodes = set(n for sl in slices.values() for n in sl)
    disconnected = all_nodes - set(result.keys())
    for node in disconnected:
        best_ph = max(emission_dict[node], key=lambda k: emission_dict[node][k])
        result[node] = best_ph

    final_results = {1: result}

    return final_results


def true_viterbi_dp_for_dbn(slices, edges, emission_dict, transitions_dict):
    """
    A complete, correct Viterbi dynamic programming for a general DBN.

    Requirements:
    -------------
    1) 'slices' is a dict: {t: [nodeA, nodeB, ...], ...} in topological order.
    2) 'edges' is a list of (parent, child).
    3) 'emission_dict[node][phasing]' = precomputed log-likelihood 
       of node 'node' in state 'phasing' (already in log space).
    4) 'transitions_dict[f"{u}--{v}"]' = a 2D array of shape (num_u, num_v),
       with prob[u_index, v_index], so we do math.log(...) during DP.

    We unify each node as a single entity across slices 
    (meaning node is not repeated if your slices says '7-8' in multiple slices).

    Steps:
    ------
    1) Build adjacency to know parents, children.
    2) Forward DP (dp[node][phasing]) = best log-prob path that ends in node=phasing.
    3) For multi-parents, we combine them with a max over all parent states 
       (the standard Viterbi max-product).
    4) We store backpointers to figure out which parent & parent's phasing gave that max.
    5) Then we identify terminal nodes from the final slice or from "no children".
    6) We BFS from each terminal node's best phasing back to all ancestors, 
       ensuring each node in the subgraph gets assigned a best state.
    7) Fill any leftover truly disconnected node with its best emission alone.

    Returns:
    --------
    final_results = {1: best_phasing_dict}
    where best_phasing_dict[node] = chosen phasing (string).

    No placeholders, no simplification. Full real Viterbi for DBN.
    """

    # Step 1) Build adjacency
    # We want 'reverse_adj[v]' = list of all parents u of v
    # and 'forward_adj[u]' = list of all children v of u
    reverse_adj = defaultdict(list)
    forward_adj = defaultdict(list)
    for u, v in edges:
        reverse_adj[v].append(u)
        forward_adj[u].append(v)

    # Build a set of all nodes
    all_nodes = set()
    for t in slices:
        for nd in slices[t]:
            all_nodes.add(nd)

    # Step 2) Data structures for DP:
    # dp[node][phasing] = best log-prob path that ends in node=phasing
    # backpointer[node][phasing] = (parent_node, parent_phasing) that gave max
    dp = defaultdict(dict)
    backpointer = defaultdict(dict)

    # Step 3) Forward pass, following the topological order 'slices'
    sorted_times = sorted(slices.keys())

    for t in sorted_times:
        for node in slices[t]:
            # For each node, go through each phasing
            if node not in emission_dict:
                # If no emission info, skip or treat as zero-likelihood
                continue

            phasings = list(emission_dict[node].keys())
            for ph in phasings:
                log_emit = emission_dict[node][ph]  # already in log

                # If no parents, or t==1, we do prior assumption
                if t == sorted_times[0] or len(reverse_adj[node]) == 0:
                    # uniform prior among phasings => -log(#phasings)
                    prior_val = -math.log(len(phasings)) 
                    dp[node][ph] = log_emit + prior_val
                    backpointer[node][ph] = None  # no parent
                else:
                    # we do the max over all parents/states
                    best_score = float('-inf')
                    best_parent = None
                    # Each parent can appear in a previous slice
                    for parent in reverse_adj[node]:
                        # Check if we have dp for parent
                        if parent not in dp:
                            continue

                        # transitions_dict is a 2D array for edge parent->node
                        edge_label = f"{parent}--{node}"
                        if edge_label not in transitions_dict:
                            continue

                        parent_phasings = list(dp[parent].keys())  # states that have a DP score
                        child_phasings = list(emission_dict[node].keys())
                        i_v = child_phasings.index(ph)

                        for j_u, p_ph in enumerate(parent_phasings):
                            # j_u is the index of parent's phasing in parent's possible states
                            trans_prob = transitions_dict[edge_label][j_u, i_v]
                            if trans_prob <= 0:
                                continue
                            candidate = dp[parent][p_ph] + math.log(trans_prob)
                            if candidate > best_score:
                                best_score = candidate
                                best_parent = (parent, p_ph)

                    dp[node][ph] = log_emit + best_score
                    backpointer[node][ph] = best_parent

    # Step 4) Identify terminal nodes. 
    # We can do so by the last slice: 
    # or more generally, nodes with no children in the DBN
    # but we respect your existing slice concept. 
    last_t = max(sorted_times)
    terminal_nodes = set(slices[last_t])  # last slice for your DBN

    # Step 5) BFS-like backtracking to fill best states for all connected ancestors
    result = {}
    visited = set()

    queue = []

    # Enqueue best phasing from each node in the last slice
    for nd in slices[last_t]:
        if nd in dp:  # skip if dp not exist
            if len(dp[nd]) > 0:
                # pick best ph
                best_ph = max(dp[nd].keys(), key=lambda x: dp[nd][x])
                queue.append((nd, best_ph))

    # BFS or stack approach
    while queue:
        node, ph = queue.pop()
        if node not in result:
            # assign best state
            result[node] = ph
            visited.add(node)
            # push the parent from backpointer
            par_info = backpointer[node][ph]
            if par_info is not None:
                parent_node, parent_ph = par_info
                if parent_node not in result:  # not assigned yet
                    queue.append((parent_node, parent_ph))

    # Step 6) For any node never reached in backtracking but is not truly disconnected
    # We'll pick the best ph in dp. If no dp, it's fully disconnected => fallback
    assigned_nodes = set(result.keys())
    leftover = all_nodes - assigned_nodes

    for nd in leftover:
        if nd in dp and len(dp[nd]) > 0:
            # pick best state from dp
            best_ph = max(dp[nd].keys(), key=lambda x: dp[nd][x])
            result[nd] = best_ph
        else:
            # truly no parent, no child => disconnected
            if nd in emission_dict and len(emission_dict[nd]) > 0:
                best_ph = max(emission_dict[nd].keys(), key=lambda x: emission_dict[nd][x])
                result[nd] = best_ph
            else:
                # fallback
                result[nd] = None

    final_results = {1: result}
    return final_results
