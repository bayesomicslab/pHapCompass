
def compute_log_emission(node, phasing, emission_dict):
    """
    For a given node and candidate phasing, combine the sub-likelihood factors
    (which are not in log space) stored in emission_dict[node][phasing] by multiplying them,
    and return the logarithm of the product.
    
    For example, if:
       emission_dict[node][phasing] = {'00': 1.003996, '01': 0.002004, '10': 5.990004, '11': 1.003996}
    then the final emission probability = (1.003996)*(0.002004)*(5.990004)*(1.003996),
    and this function returns math.log(product).
    
    (If any sub-value is <= 0, return -infinity.)
    """
    sub_dict = emission_dict[node][phasing]
    product_val = 1.0
    for key, val in sub_dict.items():
        if val <= 0:
            return -inf
        product_val *= val
    if product_val <= 0:
        return -inf
    return math.log(product_val)

def handle_disconnected_nodes(slices, viterbi_scores, viterbi_back, emission_dict):
    """
    For each node (except in the first slice) if all candidate phasings have a score of -inf,
    then treat that node as a fresh start by setting its score to (0 + log(emission)).
    This allows disconnected subgraphs to be decoded separately.
    """
    sorted_times = sorted(slices.keys())
    for idx, t in enumerate(sorted_times):
        if idx == 0:
            continue  # first slice already initialized
        if t not in viterbi_scores:
            continue
        for node in slices[t]:
            if node not in viterbi_scores[t]:
                continue
            all_inf = True
            for phasing, val in viterbi_scores[t][node].items():
                if val > -inf:
                    all_inf = False
                    break
            if all_inf:
                # Reinitialize this node as a fresh start.
                for phasing in viterbi_scores[t][node].keys():
                    log_emit = compute_log_emission(node, phasing, emission_dict)
                    viterbi_scores[t][node][phasing] = 0.0 + log_emit
                    viterbi_back[t][node][phasing] = None
    return viterbi_scores, viterbi_back

def compute_forward_messages_viterbi(slices, edges, assignment_dict, emission_dict, transitions_dict, frag_path):
    """
    Computes forward filtering messages in log space using the same recursion as in FFBS.
    (This is essentially the same as your FFBS forward pass.)
    
    Parameters:
      slices: dict mapping slice index to list of nodes.
      edges: list of (parent, child) tuples.
      assignment_dict: dictionary (if needed for observations; here assumed not to be used because
                       emission_dict already contains the final likelihoods).
      emission_dict: as defined above.
      transitions_dict: dict where for each edge "parent--child" we have a NumPy array of shape (m, n)
                        giving the transition probabilities from parent's m states to child's n states.
      frag_path: path to fragment file (not used here because emission_dict is final).
      
    Returns:
      forward_messages: dict mapping slice index -> node -> phasing -> log probability.
    """
    # Build reverse adjacency (not strictly needed here but we follow the original structure)
    reverse_adjacency_list = defaultdict(list)
    for s, t in edges:
        reverse_adjacency_list[t].append(s)
    
    forward_messages = {}
    # Initialize: set all scores to -inf.
    for t in slices:
        forward_messages[t] = {}
        for node in slices[t]:
            forward_messages[t][node] = {}
            for phase in emission_dict[node].keys():
                forward_messages[t][node][phase] = -inf

    # Base case: first slice.
    first_t = min(slices.keys())
    for node in slices[first_t]:
        state_space = list(emission_dict[node].keys())
        for phase in state_space:
            log_prior = 0.0  # Uniform prior
            log_emit = compute_log_emission(node, phase, emission_dict)
            forward_messages[first_t][node][phase] = log_prior + log_emit

    # Recursion: for each subsequent slice, use logaddexp to sum contributions.
    sorted_times = sorted(slices.keys())
    for idx in range(1, len(sorted_times)):
        t = sorted_times[idx]
        prev_t = sorted_times[idx - 1]
        for node in slices[t]:
            state_space = list(emission_dict[node].keys())
            for i, phase in enumerate(state_space):
                log_emit = compute_log_emission(node, phase, emission_dict)
                log_transition_sum = -inf
                # Only consider parents that appear in the previous slice.
                parents = list(set(reverse_adjacency_list[node]).intersection(set(slices[prev_t])))
                for p_node in parents:
                    parent_order = list(emission_dict[p_node].keys())
                    for j, p_phase in enumerate(parent_order):
                        parent_message = forward_messages[prev_t][p_node].get(p_phase, -inf)
                        if parent_message == -inf:
                            continue
                        edge_label = f"{p_node}--{node}"
                        try:
                            trans_prob = transitions_dict[edge_label][j, i]
                        except Exception:
                            trans_prob = 0.0
                        if trans_prob <= 0:
                            continue
                        candidate = parent_message + math.log(trans_prob)
                        log_transition_sum = np.logaddexp(log_transition_sum, candidate)
                forward_messages[t][node][phase] = log_emit + log_transition_sum
    return forward_messages

def max_decoding(forward_messages, slices):
    """
    Given forward_messages in log space (computed as in FFBS) and the slices,
    perform maximum decoding (i.e. pick, for each node, the candidate state with highest score).
    
    Returns:
      decoded: dict mapping slice index -> dict of node -> best phasing.
    """
    decoded = {}
    for t in slices:
        decoded[t] = {}
        for node in slices[t]:
            best_phase = None
            best_val = -inf
            for phase, val in forward_messages[t][node].items():
                if val > best_val:
                    best_val = val
                    best_phase = phase
            decoded[t][node] = best_phase
    return decoded

def run_viterbi(slices, edges, assignment_dict, emission_dict, transitions_dict, data_path):
        # So we now compute forward messages:
    forward_messages = compute_forward_messages_viterbi(slices, edges, assignment_dict, emission_dict, transitions_dict, data_path)
    
    # And then decode by taking the maximum for each node.
    best_assignment = max_decoding(forward_messages, slices)
    return best_assignment