import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Categorical
import itertools
import numpy as np
from torch_model_utils import *


def find_matchings(nodes_part1, nodes_part2):
    # Sort both parts and remember the original indices.
    sorted_part1 = sorted(enumerate(nodes_part1), key=lambda x: x[1])
    sorted_part2 = sorted(enumerate(nodes_part2), key=lambda x: x[1])
    
    # Split nodes by type and collect their original indices.
    def split_by_type(sorted_nodes):
        grouped = {}
        for idx, t in sorted_nodes:
            if t not in grouped:
                grouped[t] = []
            grouped[t].append(idx)
        return grouped
    
    grouped_part1 = split_by_type(sorted_part1)
    grouped_part2 = split_by_type(sorted_part2)
    if grouped_part1.keys() != grouped_part2.keys():
        return []
    if any([len((grouped_part1[i])) != len((grouped_part2[i])) for i in grouped_part1.keys()]):
        return []
    # Start with a single empty matching.
    matchings = [[]]
    for node_type, indices1 in grouped_part1.items():
        indices2 = grouped_part2[node_type]
        
        # For each current matching, extend it with all possible permutations for the current type.
        new_matchings = []
        for perm in itertools.permutations(indices2, len(indices2)):
            for current_matching in matchings:
                # Add new matching to the results only if it doesn't conflict with the current matching.
                if all((i1, i2) not in current_matching for i1, i2 in zip(indices1, perm)):
                    new_matchings.append(current_matching + list(zip(indices1, perm)))
        matchings = new_matchings
    
    return matchings


def str_2_phas_1(phasing, ploidy):
    return np.array([int(p) for p in [*phasing]]).reshape(ploidy, -1)


def phas_2_str(phas):
    return ''.join([str(ph) for ph in list(np.ravel(phas))])


def find_phasings_matches(ff, sf, common_ff, common_sf):
    templates = []
    all_local = find_matchings(list(ff[:, -1]), list(sf[:, 0]))
    for al in all_local:
        ff_ordering = [ii[0] for ii in al]
        sf_ordering = [ii[1] for ii in al]
        assert any(ff[ff_ordering, common_ff] == sf[sf_ordering, common_sf])
        temp = np.hstack([ff[ff_ordering, :], sf[sf_ordering, 1:]])
        byte_set = {a.tobytes() for a in templates}
        if temp.tobytes() not in byte_set:
            templates.append(temp)
    return templates


def find_common_element_and_index(node1, node2):
    # Split the node names into components (e.g., '1-2' -> ['1', '2'])
    node1_parts = node1.split('-')
    node2_parts = node2.split('-')
    
    # Find the common element between node1 and node2
    common_element = None
    for part in node1_parts:
        if part in node2_parts:
            common_element = part
            break
    
    if common_element is None:
        raise ValueError(f"No common element found between {node1} and {node2}")
    
    # Find the index of the common element in both nodes
    common_ff = node1_parts.index(common_element)
    common_sf = node2_parts.index(common_element)
    
    return common_ff, common_sf


def permute_rows(a):
    # Generate all permutations of the rows
    perms = list(itertools.permutations(a))
    
    # Convert each permutation into a NumPy array and return as a list of matrices
    permuted_matrices = [np.array(p) for p in perms]
    
    return permuted_matrices


def compute_likelihood(observed, phasing, error_rate):
    """This likelihood computation assumes the length of observation is the same as the length of phasing"""
    y = np.tile(observed, (phasing.shape[0], 1))
    diff = y - phasing
    diff[diff != 0] = 1
    comp_diff = 1 - diff
    term1 = diff * error_rate
    term2 = comp_diff * (1 - error_rate)
    terms = term1 + term2
    probs = np.prod(terms, axis=1)
    likelihood = np.mean(probs)
    return likelihood


def compute_likelihood_generalized_plus(observed, phasing, obs_pos, phas_pos, error_rate):
    """This likelihood computation can accept different length observed and phasing, but the length of obs_pos and
    phas_pos should be the same. The likelihood is computed on the provided indices on both vectors"""
    new_phasing = phasing[:, phas_pos]
    new_observed = observed[obs_pos]
    y = np.tile(new_observed, (new_phasing.shape[0], 1))
    diff = y - new_phasing
    diff[diff != 0] = 1
    comp_diff = 1 - diff
    term1 = diff * error_rate
    term2 = comp_diff * (1 - error_rate)
    terms = term1 + term2
    probs = np.prod(terms, axis=1)
    likelihood = np.mean(probs)
    return likelihood


def generate_binary_combinations(length):
    return [list(x) for x in itertools.product([0, 1], repeat=length)]


def generate_all_possible_emissions(unique_node_lengths):
    all_emissions = []
    for length in unique_node_lengths:
        all_emissions.extend(generate_binary_combinations(length))
    return all_emissions


def generate_hmm_with_weights_and_emissions(qg, error_rate):
    state_names = []  # To store the names of states in hmmg
    state_index_map = {}  # Map state names to their index in the transition matrix
    emission_index_map = {}  # To map emissions to their index for the matrix
    unique_node_lengths = set()  # Track the unique lengths of node names
    sorted_nodes = sorted(qg.nodes(data=True), key=lambda x: x[0])
    # Step 1: Generate the state names and collect unique node lengths
    for node, node_data in sorted_nodes:
        for state_key in node_data['weight'].keys():
            state_name = f"{node}-{state_key}"
            state_names.append(state_name)
            node_length = len(node.split('-'))  # Get the length of the node part (excluding state_key)
            unique_node_lengths.add(node_length)  # Track the unique lengths

    # Convert unique_node_lengths to a sorted list (just in case)
    unique_node_lengths = sorted(unique_node_lengths)

    # Step 2: Generate all possible emissions (0-1 combinations) based on the actual unique node lengths
    all_emissions = generate_all_possible_emissions(unique_node_lengths)
    
    # Map each emission to a unique index
    for i, emission in enumerate(all_emissions):
        emission_index_map[tuple(emission)] = i  # Map emissions to indices
    
    num_emissions = len(all_emissions)
    num_states = len(state_names)
    
    # Step 3: Create the emission probability matrix (number of states x number of emissions), initialized with zeros
    emission_prob_matrix = np.zeros((num_states, num_emissions))
    
    # Create an index for each state
    for i, state_name in enumerate(state_names):
        state_index_map[state_name] = i

    # Step 4: Initialize transition matrix
    transition_matrix = np.zeros((num_states, num_states))  # Start with a zero matrix

    # Step 5: Define transitions between neighboring nodes with weights
    for node1, node2 in qg.edges():
        # Get the states of node1
        for state_key1 in qg.nodes[node1]['weight'].keys():
            state1 = f"{node1}-{state_key1}"
            state1_idx = state_index_map[state1]
            # Get the states of node2
            for state_key2 in qg.nodes[node2]['weight'].keys():
                state2 = f"{node2}-{state_key2}"
                state2_idx = state_index_map[state2]
                
                # Step 6: Compute the transition weight based on phasings
                ff = str_2_phas_1(state_key1, 3)
                sf = str_2_phas_1(state_key2, 3)
                
                common_ff, common_sf = find_common_element_and_index(node1, node2)
                
                phasings = find_phasings_matches(ff, sf, common_ff, common_sf)
                all_phasings_str = []
                for phas in phasings:
                    all_phasings = permute_rows(phas)
                    phasings_str = [phas_2_str(phasss) for phasss in all_phasings]
                    all_phasings_str += phasings_str
                # print('matched phasings', phasings_str)

                # 4. Calculate the transition weight by summing the values of matching keys
                transition_weight = 0
                for phasing_key in all_phasings_str:
                    if phasing_key in qg[node1][node2]['weight']:
                        transition_weight += qg[node1][node2]['weight'][phasing_key]
                
                transition_matrix[state1_idx][state2_idx] = transition_weight

    # Step 7: Normalize the transition matrix row-wise, skipping zero-sum rows
    for i in range(num_states):
        row_sum = transition_matrix[i].sum()
        if row_sum > 0:
            transition_matrix[i] /= row_sum
        else:
            transition_matrix[i] = np.zeros(num_states)

    # Step 8: Calculate emissions for each state and populate the emission probability matrix
    for state in state_names:
        node_part, state_key = state.rsplit('-', 1)  # Split the state into node part and key
        node_elements = node_part.split('-')  # For example, '1-2' -> ['1', '2']
        node_length = len(node_elements)  # Determine the number of elements in the node

        # Generate all possible combinations of 0 and 1 based on the node length
        possible_emissions = generate_binary_combinations(node_length)

        # Compute the likelihood for each emission using compute_likelihood
        phasing = str_2_phas_1(state_key, 3)  # Compute phasing for the state key
        state_idx = state_index_map[state]  # Get the index for the state in the matrix

        # Set emission probabilities for combinations with the same length as node_elements
        for emission in possible_emissions:
            likelihood = compute_likelihood(np.array(emission), phasing, error_rate)
            emission_idx = emission_index_map[tuple(emission)]  # Get the index for the emission
            emission_prob_matrix[state_idx][emission_idx] = likelihood

    return state_names, transition_matrix, emission_prob_matrix, emission_index_map


def generate_hmm_with_weights_and_emissions2(qg, error_rate):
    state_names = []  # To store the names of states in hmmg
    state_index_map = {}  # Map state names to their index in the transition matrix
    emission_index_map = {}  # To map emissions to their index for the matrix
    unique_node_lengths = set()  # Track the unique lengths of node names
    sorted_nodes = sorted(qg.nodes(data=True), key=lambda x: x[0])
    # Step 1: Generate the state names and collect unique node lengths
    for node, node_data in sorted_nodes:
        for state_key in node_data['weight'].keys():
            state_name = f"{node}-{state_key}"
            state_names.append(state_name)
            node_length = len(node.split('-'))  # Get the length of the node part (excluding state_key)
            unique_node_lengths.add(node_length)  # Track the unique lengths

    # Convert unique_node_lengths to a sorted list (just in case)
    unique_node_lengths = sorted(unique_node_lengths)

    # Step 2: Generate all possible emissions (0-1 combinations) based on the actual unique node lengths
    all_emissions = generate_all_possible_emissions(unique_node_lengths)
    
    # Map each emission to a unique index
    for i, emission in enumerate(all_emissions):
        emission_index_map[tuple(emission)] = i  # Map emissions to indices
    
    num_emissions = len(all_emissions)
    num_states = len(state_names)
    
    # Step 3: Create the emission probability matrix (number of states x number of emissions), initialized with zeros
    emission_prob_matrix = np.zeros((num_states, num_emissions))
    
    # Create an index for each state
    for i, state_name in enumerate(state_names):
        state_index_map[state_name] = i

    # Step 4: Initialize transition matrix
    transition_matrix = np.zeros((num_states, num_states))  # Start with a zero matrix

    sorted_tuple_list = [sorted(iii, key=lambda x: int(x.split('-')[0])) for iii in qg.edges()]


    # Step 5: Define transitions between neighboring nodes with weights
    for node1, node2 in sorted_tuple_list:
        # Get the states of node1
        for state_key1 in qg.nodes[node1]['weight'].keys():
            state1 = f"{node1}-{state_key1}"
            state1_idx = state_index_map[state1]
            # Get the states of node2
            for state_key2 in qg.nodes[node2]['weight'].keys():
                state2 = f"{node2}-{state_key2}"
                state2_idx = state_index_map[state2]
                # print(state1, state2)
                # Step 4: Compute the transition weight based on phasings
                # 1. Convert the state keys to phasings
                ff = str_2_phas_1(state_key1, 3)
                sf = str_2_phas_1(state_key2, 3)
                
                # 2. Find the common index (e.g., common_ff and common_sf)
                common_ff, common_sf = find_common_element_and_index(node1, node2)
                # print('common indices:', common_ff, common_sf)

                # 3. Get the matching phasings
                phasings = find_phasings_matches(ff, sf, common_ff, common_sf)
                all_phasings_str = []
                for phas in phasings:
                    all_phasings = permute_rows(phas)
                    # all_phasings = set(tuple(arr) for arr in permute_rows(phas))
                    # all_phasings = {tuple(arr) for arr in permute_rows(phas)}

                    phasings_str = [phas_2_str(phasss) for phasss in all_phasings]
                    for pe in phasings_str:
                        if pe not in all_phasings_str:
                            all_phasings_str.append(pe)

                # 4. Calculate the transition weight by summing the values of matching keys
                transition_weight = 0
                for phasing_key in all_phasings_str:
                    if phasing_key in qg[node1][node2]['weight']:
                        transition_weight += qg[node1][node2]['weight'][phasing_key]
                    # else:

                # print('edge weights:', qg[node1][node2]['weight'])
                # print('transition prob.', transition_weight)
                # 5. Set the transition weight in the matrix
                transition_matrix[state1_idx][state2_idx] = transition_weight



    # Step 7: Normalize the transition matrix row-wise, skipping zero-sum rows
    for i in range(num_states):
        row_sum = transition_matrix[i].sum()
        if row_sum > 0:
            transition_matrix[i] /= row_sum
        else:
            transition_matrix[i] = np.zeros(num_states)

    # Step 8: Calculate emissions for each state and populate the emission probability matrix
    for state in state_names:
        node_part, state_key = state.rsplit('-', 1)  # Split the state into node part and key
        node_elements = node_part.split('-')  # For example, '1-2' -> ['1', '2']
        node_length = len(node_elements)  # Determine the number of elements in the node

        # Generate all possible combinations of 0 and 1 based on the node length
        possible_emissions = generate_binary_combinations(node_length)

        # Compute the likelihood for each emission using compute_likelihood
        phasing = str_2_phas_1(state_key, 3)  # Compute phasing for the state key
        state_idx = state_index_map[state]  # Get the index for the state in the matrix

        # Set emission probabilities for combinations with the same length as node_elements
        for emission in possible_emissions:
            likelihood = compute_likelihood(np.array(emission), phasing, error_rate)
            emission_idx = emission_index_map[tuple(emission)]  # Get the index for the emission
            emission_prob_matrix[state_idx][emission_idx] = likelihood

    return state_names, transition_matrix, emission_prob_matrix, emission_index_map


def sort_tuples_by_number(tuples):
    sorted_tuples = []
    for t in tuples:
        # Extract the first number from each string in the tuple
        num1 = int(t[0].split('-')[0])
        num2 = int(t[1].split('-')[0])
        
        # Sort the tuple based on the extracted numbers
        if num1 <= num2:
            sorted_tuples.append(t)
        else:
            sorted_tuples.append((t[1], t[0]))
    
    return sorted_tuples


def generate_hmm_with_weights(qg, error_rate=0.001):
    state_names = []  # To store the names of states in hmmg
    state_index_map = {}  # Map state names to their index in the transition matrix
    emissions_map = {}  # To store emissions for each state
    unique_node_lengths = set()  # Track the unique lengths of node names

    # Step 1: Generate the state names
    for node, node_data in qg.nodes(data=True):
        for state_key in node_data['weight'].keys():
            state_name = f"{node}-{state_key}"
            state_names.append(state_name)
    
    # Create an index for each state
    for i, state_name in enumerate(state_names):
        state_index_map[state_name] = i

    # Step 2: Initialize transition matrix
    num_states = len(state_names)
    transition_matrix = np.zeros((num_states, num_states))  # Start with a zero matrix

    # Step 3: Define transitions between neighboring nodes with weights
    for node1, node2 in qg.edges():
        # Get the states of node1
        for state_key1 in qg.nodes[node1]['weight'].keys():
            state1 = f"{node1}-{state_key1}"
            state1_idx = state_index_map[state1]
            # Get the states of node2
            for state_key2 in qg.nodes[node2]['weight'].keys():
                state2 = f"{node2}-{state_key2}"
                state2_idx = state_index_map[state2]
                # print(state1, state2)
                # Step 4: Compute the transition weight based on phasings
                # 1. Convert the state keys to phasings
                ff = str_2_phas_1(state_key1, 3)
                sf = str_2_phas_1(state_key2, 3)
                
                # 2. Find the common index (e.g., common_ff and common_sf)
                common_ff, common_sf = find_common_element_and_index(node1, node2)
                # print('common indices:', common_ff, common_sf)

                # 3. Get the matching phasings
                phasings = find_phasings_matches(ff, sf, common_ff, common_sf)
                all_phasings_str = []
                for phas in phasings:
                    all_phasings = permute_rows(phas)
                    # all_phasings = set(tuple(arr) for arr in permute_rows(phas))
                    # all_phasings = {tuple(arr) for arr in permute_rows(phas)}

                    phasings_str = [phas_2_str(phasss) for phasss in all_phasings]
                    for pe in phasings_str:
                        if pe not in all_phasings_str:
                            all_phasings_str.append(pe)

                # 4. Calculate the transition weight by summing the values of matching keys
                transition_weight = 0
                for phasing_key in all_phasings_str:
                    if phasing_key in qg[node1][node2]['weight']:
                        transition_weight += qg[node1][node2]['weight'][phasing_key]
                    # else:

                # print('edge weights:', qg[node1][node2]['weight'])
                # print('transition prob.', transition_weight)
                # 5. Set the transition weight in the matrix
                transition_matrix[state1_idx][state2_idx] = transition_weight
    


    # transition_matrix = transition_matrix / transition_matrix.sum(axis=1, keepdims=True)
    # Step 5: Normalize the transition matrix row-wise, skipping zero-sum rows
    for i in range(num_states):
        row_sum = transition_matrix[i].sum()
        if row_sum > 0:
            transition_matrix[i] /= row_sum
        else:
            # Optionally, set a uniform distribution for zero-sum rows
            transition_matrix[i] = np.zeros(num_states)  # or set to uniform probabilities
            # Uncomment the following line if you prefer uniform probabilities instead of zeros
            # transition_matrix[i] = np.full(num_states, 1/num_states)

    # Step 6: Calculate emissions for each state
    for state in state_names:
        node_part, state_key = state.rsplit('-', 1)  # Split the state into node part and key
        node_elements = node_part.split('-')  # For example, '1-2' -> ['1', '2']
        node_length = len(node_elements)  # Determine the number of elements in the node

        # Step 7: Generate all possible combinations of 0 and 1 based on node length
        possible_emissions = generate_binary_combinations(node_length)

        # Step 8: Compute the likelihood for each emission using compute_likelihood
        phasing = str_2_phas_1(state_key, 3)  # Compute phasing for the state key
        emission_probs = []
        for emission in possible_emissions:
            likelihood = compute_likelihood(np.array(emission), phasing, error_rate)
            emission_probs.append((emission, likelihood))

        # Step 9: Handle length mismatches (set emissions for longer lengths to zero)
        if node_length < len(state_key):
            extra_combinations = generate_binary_combinations(len(state_key))
            for extra_emission in extra_combinations:
                if extra_emission not in possible_emissions:
                    emission_probs.append((extra_emission, 0))  # Set to zero

        # Store emissions and their probabilities for this state
        emissions_map[state] = emission_probs


    return state_names, transition_matrix, emissions_map


def combine_potentials(self, transition, emission):
    """Mix the transition and emission scores

    Args:
      transition: 
      emission_potentials: torch.Tensor(float), 
        size=[batch, max_len, num_state]

    Returns:
      scores: size=[batch, len, num_state, num_state]
        scores := log phi(batch, x_t, y_{t-1}, y_t)
    """
    batch_size = emission.size(0)
    seq_len = emission.size(1)
    num_state = emission.size(2)

    # scores[batch, t, C, C] = log_potentials(t, from y_{t-1}, to y_t)
    if(len(transition.size()) == 2):
      log_potentials = transition.view(1, 1, num_state, num_state)\
        .expand(batch_size, seq_len, num_state, num_state) + \
        emission.view(batch_size, seq_len, 1, num_state)\
        .expand(batch_size, seq_len, num_state, num_state)
    else: 
      log_potentials = transition + \
        emission.view(batch_size, seq_len, 1, num_state)\
        .expand(batch_size, seq_len, num_state, num_state)
    return log_potentials


def rsample(self, transition_potentials, emission_potentials, seq_lens, 
    log_potentials=None, tau=1.0, return_prob=False, z_st=False):
    """Reparameterized CRF sampling, a Gumbelized version of the 
    Forward-Filtering Backward-Sampling algorithm

    TODO: an autograd based implementation 
    requires to redefine the backward function over a relaxed-sampling semiring
    
    Args:
      emission_potentials: type=torch.tensor(float), 
        size=[batch, max_len, num_state]
      seq_lens: type=torch.tensor(int), size=[batch]
      tau: type=float, anneal strength

    Returns
      sample: size=[batch, max_len]
      relaxed_sample: size=[batch, max_len, num_state]
      sample_log_prob: size=[batch]
    """
    if(transition_potentials is not None):
      transition_potentials, emission_potentials = self.normalize(
        transition_potentials, emission_potentials)

    # Algo 2 line 1
    if(log_potentials is None):
      log_potentials = self.combine_potentials(
        transition_potentials, emission_potentials)
    alpha, log_Z = self.forward_sum(
      None, emission_potentials, seq_lens, log_potentials) 

    batch_size = emission_potentials.size(0)
    max_len = emission_potentials.size(1)
    num_state = emission_potentials.size(2)
    device = emission_potentials.device

    # Backward sampling start
    # The sampling still goes backward, but for simple implementation we
    # reverse the sequence, so in the code it still goes from 1 to T 
    relaxed_sample_rev = torch.zeros(batch_size, max_len, num_state).to(device)
    sample_prob = torch.zeros(batch_size, max_len).to(device)
    sample_rev = torch.zeros(batch_size, max_len).type(torch.long).to(device)
    alpha_rev = reverse_sequence(alpha, seq_lens).to(device)
    log_potentials_rev = reverse_sequence(log_potentials, seq_lens).to(device)
    
    # Algo 2 line 3, log space
    # w.shape=[batch, num_state]
    w = alpha_rev[:, 0, :].clone()
    w -= log_Z.view(batch_size, -1)
    p = w.exp()
    # switching regularization for longer chunk, not mentioned in the paper
    # so do no need to care. In the future this will be updated with posterior
    # regularization
    # if(return_switching): 
    #   switching = 0.
    
    # Algo 2 line 4
    relaxed_sample_rev[:, 0] = reparameterize_gumbel(w, tau)
    # Algo 2 line 5
    sample_rev[:, 0] = relaxed_sample_rev[:, 0].argmax(dim=-1)
    sample_prob[:, 0] = tmu.batch_index_select(p, sample_rev[:, 0]).flatten()
    mask = tmu.length_to_mask(seq_lens, max_len).type(torch.float)
    prev_p = p
    for i in range(1, max_len):
      # y_after_to_current[j, k] = log_potentials(y_{t - 1} = k, y_t = j, x_t)
      # size=[batch, num_state, num_state]
      y_after_to_current = log_potentials_rev[:, i-1].transpose(1, 2)
      # w.size=[batch, num_state]
      w = tmu.batch_index_select(y_after_to_current, sample_rev[:, i-1])
      w_base = tmu.batch_index_select(alpha_rev[:, i-1], sample_rev[:, i-1])
      # Algo 2 line 7, log space
      w = w + alpha_rev[:, i] - w_base.view(batch_size, 1)
      p = F.softmax(w, dim=-1) # p correspond to pi in the paper
      # if(return_switching):
      #   switching += (tmu.js_divergence(p, prev_p) * mask[:, i]).sum()
      prev_p = p
      # Algo 2 line 8
      relaxed_sample_rev[:, i] = tmu.reparameterize_gumbel(w, tau)
      # Algo 2 line 9
      sample_rev[:, i] = relaxed_sample_rev[:, i].argmax(dim=-1)
      sample_prob[:, i] = tmu.batch_index_select(p, sample_rev[:, i]).flatten()

    # Reverse the sequence back
    sample = tmu.reverse_sequence(sample_rev, seq_lens)
    relaxed_sample = tmu.reverse_sequence(relaxed_sample_rev, seq_lens)
    if(z_st):
      sample_size = relaxed_sample.size(2)
      sample_one_hot = tmu.ind_to_one_hot(sample.view(-1), sample_size)
      sample_one_hot = sample_one_hot.view(batch_size, -1, sample_size).float()
      relaxed_sample = (sample_one_hot - relaxed_sample).detach() + relaxed_sample

    sample_prob = tmu.reverse_sequence(sample_prob, seq_lens)
    sample_prob = sample_prob.masked_fill(mask == 0, 1.)
    sample_log_prob_stepwise = (sample_prob + 1e-10).log()
    sample_log_prob = sample_log_prob_stepwise.sum(dim=1)

    ret = [sample, relaxed_sample]
    # if(return_switching): 
    #   switching /= (mask.sum(dim=-1) - 1).sum()
    #   ret.append(switching)
    if(return_prob):
      ret.extend([sample_log_prob, sample_log_prob_stepwise])
    return ret


class HMM(nn.Module):
    """Implementation of Hidden Markov Model with forward algorithm and 
    backward sampling."""
    
    def __init__(self, potential_normalization='none', potential_scale=1.0):
        super(HMM, self).__init__()
        self.potential_normalization = potential_normalization
        self.potential_scale = potential_scale

    def normalize(self, transition, emission):
        if self.potential_normalization == 'minmax':
            transition = (transition - transition.mean()) / (transition.max() - transition.min())
            emission = (emission - emission.mean(-1, keepdim=True)) / (emission.max(-1, keepdim=True).values - emission.min(-1, keepdim=True).values)
        elif self.potential_normalization == 'zscore':
            transition = (transition - transition.mean()) / transition.std()
            emission = (emission - emission.mean(-1, keepdim=True)) / emission.std(-1, keepdim=True)
        return self.potential_scale * transition, self.potential_scale * emission

    def combine_potentials(self, transition, emission):
        """For HMM, the emission depends only on the current state."""
        batch_size = emission.size(0)
        num_state = emission.size(2)
        log_potentials = transition.view(1, num_state, num_state).expand(batch_size, num_state, num_state)
        return log_potentials


    def forward_sum_phasing(self, transition_potentials, emission_potentials, seq_lens):
        # max phasings for nodes 
        # n nodes 
        # positions are the number of nodes (# layers in NN)
        # max number of phasings for the nodes (e.g. k + m)
        num_state_node = 4
        positions = 7
        alpha = torch.zeros(positions, num_state_node).to(emission_potentials.device)
        alpha[0, :] = emission_potentials[0, 0, :]
        for t in range(1, positions):
            # Compute alpha_t for each time step
            alpha[t, :] = torch.logsumexp(alpha[t - 1, :].unsqueeze(2) + transition_potentials, dim=1) + emission_potentials[:, t, :]


    def forward_sum(self, transition_potentials, emission_potentials, seq_lens):
        """Forward algorithm for HMM to compute log probability."""
        transition_potentials, emission_potentials = self.normalize(transition_potentials, emission_potentials)

        batch_size = emission_potentials.size(0)
        seq_len = emission_potentials.size(1)
        num_state = emission_potentials.size(2)

        alpha = torch.zeros(batch_size, seq_len, num_state).to(emission_potentials.device)
        alpha[:, 0, :] = emission_potentials[:, 0, :]

        for t in range(1, seq_len):
            # Compute alpha_t for each time step
            alpha[:, t, :] = torch.logsumexp(alpha[:, t - 1, :].unsqueeze(2) + transition_potentials, dim=1) + emission_potentials[:, t, :]

        log_Z = torch.logsumexp(alpha[:, -1, :], dim=-1)
        return alpha, log_Z

    def rsample(self, transition_potentials, emission_potentials, seq_lens):
        """Backward sampling for HMM."""
        alpha, log_Z = self.forward_sum(transition_potentials, emission_potentials, seq_lens)

        batch_size = emission_potentials.size(0)
        seq_len = emission_potentials.size(1)
        num_state = emission_potentials.size(2)

        sample = torch.zeros(batch_size, seq_len).long().to(emission_potentials.device)

        # Backward sampling (start with last state)
        p_T = torch.softmax(alpha[:, -1, :], dim=-1)
        sample[:, -1] = Categorical(p_T).sample()

        for t in reversed(range(1, seq_len)):
            transition_potential_t = transition_potentials[sample[:, t]].view(batch_size, num_state)
            p_t = torch.softmax(alpha[:, t - 1, :] + transition_potential_t, dim=-1)
            sample[:, t - 1] = Categorical(p_t).sample()

        return sample


# Initialize the HMM model
hmm = HMM(potential_normalization='none', potential_scale=1.0)

# Example settings: 3 states, sequence length 5, and batch size 2
num_states = 2
seq_len = 10
batch_size = 1

# # Random transition potentials: size=[num_states, num_states]
# transition_potentials = torch.randn(num_states, num_states)

# # Random emission potentials: size=[batch_size, seq_len, num_states]
# emission_potentials = torch.randn(batch_size, seq_len, num_states)

# # Random sequence lengths (between 1 and seq_len)
# seq_lens = torch.randint(1, seq_len + 1, (batch_size,))



state_names, transition_matrix, emission_prob_matrix, emission_index_map = generate_hmm_with_weights_and_emissions2(qg, error_rate)
state_names_dict = {i: state for i, state in enumerate(state_names)}

transition_potentials = torch.tensor(transition_matrix, dtype=torch.float32)

emission_prob_matrix = emission_prob_matrix.T
emission_prob_matrix_new = np.expand_dims(emission_prob_matrix, axis=0)


emission_potentials = torch.tensor(emission_prob_matrix_new, dtype=torch.float32)


# emission_prob_tensor = torch.tensor(emission_prob_matrix).unsqueeze(0)






seq_lens = 10

print("Sampled sequences from FFBS:")
# Generate samples using FFBS (rsample function)
for rr in range(10):
    sample = hmm.rsample(transition_potentials, emission_potentials, seq_lens)
    state_to_name_func = np.vectorize(lambda x: state_names_dict[x])
    state_names_mapped = state_to_name_func(sample)

    # print("Sampled sequence from FFBS:")
    print(sample)
# samples = hmm.rsample(transition_potentials, emission_potentials, seq_lens)

# Print the sampled sequences
# print("Sampled sequences from FFBS:")
# print(samples)

# import numpy as np

# def forward_sampling(transition_probs, initial_probs, T):
#     """
#     Perform forward sampling to generate a sequence of hidden states.
    
#     Parameters:
#     - transition_probs: (N, N) transition matrix, where N is the number of hidden states.
#     - initial_probs: (N,) initial state distribution.
#     - T: Length of the sequence to sample.
    
#     Returns:
#     - hidden_states: (T,) array of sampled hidden states.
#     """
#     N = transition_probs.shape[0]
#     hidden_states = np.zeros(T, dtype=int)

#     # Step 1: Sample the first hidden state from the initial distribution
#     hidden_states[0] = np.random.choice(N, p=initial_probs)

#     # Step 2: Sample the remaining states based on the transition probabilities
#     for t in range(1, T):
#         previous_state = hidden_states[t - 1]
#         hidden_states[t] = np.random.choice(N, p=transition_probs[previous_state])
    
#     return hidden_states


# def ffbs_hmm_sampling(transition_probs, initial_probs, T):
#     """
#     Perform FFBS to sample a sequence of hidden states from an HMM.
    
#     Parameters:
#     - transition_probs: (N, N) transition matrix.
#     - initial_probs: (N,) initial state distribution.
#     - T: Length of the hidden state sequence to sample.
    
#     Returns:
#     - sampled_states: (T,) array of sampled hidden states.
#     """
#     # Perform forward sampling to generate a hidden state sequence
#     sampled_states = forward_sampling(transition_probs, initial_probs, T)
    
#     return sampled_states


# Example usage
if __name__ == "__main__":
    # HMM parameters
    N = 3  # Number of hidden states
    T = 10  # Length of hidden state sequence


    initial_probs = np.array([0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    # Run FFBS sampling
    state_names, transition_matrix, emission_prob_matrix, emission_index_map = generate_hmm_with_weights_and_emissions2(qg, error_rate)
    sampled_states = ffbs_hmm_sampling(transition_matrix, initial_probs, T)

    print("Sampled hidden states:", sampled_states)
