import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Categorical
import itertools
import numpy as np
import os
import argparse
import sys
import networkx as nx
from data.input_handler import InputHandler
from data.configuration import Configuration
from algorithm.haplotype_assembly import HaplotypeAssembly
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from models.factor_graph import Factorgraph
from algorithm.chordal_contraction import chordal_contraction_cycle_base, chordal_contraction
from utils.utils import *
from algorithm.inference import *


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


def select_top_k_keys(node_states_dict, k, m):
    sorted_keys = sorted(node_states_dict, key=node_states_dict.get, reverse=True)
    top_k_keys = sorted_keys[:k]

    # Step 2: Remove top k keys from the dictionary
    remaining_keys = sorted_keys[k:]

    # Step 3: Randomly choose m keys from the remaining keys
    random_m_keys = random.sample(remaining_keys, m) if len(remaining_keys) >= m else remaining_keys
    result_keys = top_k_keys + random_m_keys
    sel_mk_dict = {key: node_states_dict[key] for key in result_keys}
    return sel_mk_dict


def generate_hmm_with_weights_and_emissions_random(qg, error_rate, k, m):
    state_names = []  # To store the names of states in hmmg
    state_index_map = {}  # Map state names to their index in the transition matrix
    emission_index_map = {}  # To map emissions to their index for the matrix
    unique_node_lengths = set()  # Track the unique lengths of node names
    
    sorted_nodes = sorted(qg.nodes(data=True), key=lambda x: x[0])

    sel_states_keys = {}

    # Step 1: Generate the state names and collect unique node lengths
    for node, node_data in sorted_nodes:
        node_states_dict = node_data['weight']
        if len(node_states_dict) > k + m:
            print(f"Warning: Node {node} has more than {k+m} states. Selecting top {k+m} states.")
            sel_mk_dict = select_top_k_keys(node_states_dict, k, m)
        else:
            sel_mk_dict = node_states_dict  # If k is greater than the number of states, use all states
        sel_states_keys[node] = sel_mk_dict
        for state_key in sel_mk_dict.keys():
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
        # for state_key1 in qg.nodes[node1]['weight'].keys():
        for state_key1 in sel_states_keys[node1].keys():
            state1 = f"{node1}-{state_key1}"
            state1_idx = state_index_map[state1]
            # Get the states of node2
            # for state_key2 in qg.nodes[node2]['weight'].keys():
            for state_key2 in sel_states_keys[node2].keys():
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



class Args:
    def __init__(self):
        self.vcf_path = 'example/62_ID0.vcf'
        self.data_path = 'example/test.txt'
        self.bam_path = 'example/example.bam'
        self.genotype_path = 'example/genotype.txt'
        self.ploidy = 3
        self.error_rate = 0.001
        self.epsilon = 0.0001
        self.output_path = 'output'
        self.root_dir = 'D:/UCONN/HaplOrbit'
        self.alleles = [0, 1]



def forward_filtering_no_observations(A, B, pi, T):
    """
    Perform forward filtering in FFBS without explicit observations, handling zero sums.

    A : transition matrix (NxN)
    B : emission matrix (NxM), where M is the number of possible emissions
    pi : initial state distribution (N)
    T : length of the sequence to sample
    
    Returns:
    alpha : forward probabilities (T x N)
    """
    N = A.shape[0]  # Number of hidden states
    
    # Alpha array to store forward probabilities
    alpha = np.zeros((T, N))
    
    # Step 1: Initialize the forward probabilities with the initial state distribution (pi)
    alpha[0, :] = pi  # No observation, so we just start with the initial probabilities
    
    # Step 2: Forward recursion step (no conditioning on emissions, just transitions)
    for t in range(1, T):
        for j in range(N):
            alpha[t, j] = np.sum(alpha[t-1, :] * A[:, j])
        
        # Check if sum of alpha is zero and handle it
        alpha_sum = np.sum(alpha[t, :])
        if alpha_sum == 0:
            alpha[t, :] = 0  # Keep it as zeros if all transition probabilities are zero
        else:
            # Normalize to avoid underflow and ensure valid probabilities
            alpha[t, :] /= alpha_sum
    
    return alpha


def backward_sampling_no_observations(alpha, A, B, T):
    """
    Perform backward sampling in FFBS without explicit observations.

    alpha : forward probabilities (T x N)
    A : transition matrix (NxN)
    B : emission matrix (NxM), where M is the number of possible emissions
    T : length of the sequence to sample
    
    Returns:
    sampled_states : list of sampled hidden states
    sampled_emissions : list of sampled emissions corresponding to the hidden states
    """
    N = A.shape[0]  # Number of hidden states
    M = B.shape[1]  # Number of possible emissions
    
    # Initialize arrays to store the sampled states and emissions
    sampled_states = np.zeros(T, dtype=int)
    sampled_emissions = np.zeros(T, dtype=int)
    
    if np.sum(alpha[T-1, :]) == 0:
        probs = np.ones(alpha[T-1, :].shape[0])/alpha[T-1, :].shape[0]
    else: 
        probs = alpha[T-1, :] / np.sum(alpha[T-1, :])
    # Step 1: Sample the final hidden state from the forward probabilities at time T-1
    sampled_states[T-1] = np.random.choice(N, p=probs)
    
    # Step 2: Backward sampling of the hidden states
    for t in range(T-2, -1, -1):
        current_state = sampled_states[t+1]
        prob = alpha[t, :] * A[:, current_state]
        if np.sum(prob) == 0:
            prob = np.ones(prob.shape[0])/prob.shape[0]  # Keep it as zeros if all transition probabilities are zero
        else:  # Normalize probabilities
            prob /= np.sum(prob)  # Normalize probabilities
        sampled_states[t] = np.random.choice(N, p=prob)
    
    # Step 3: Sample the emissions based on the sampled hidden states
    for t in range(T):
        sampled_emissions[t] = np.random.choice(M, p=B[sampled_states[t], :])
    
    return sampled_states, sampled_emissions


def ffbs_no_observations(A, B, pi, T):
    """
    Full FFBS algorithm that combines forward filtering and backward sampling, without explicit observations.
    
    A : transition matrix (NxN)
    B : emission matrix (NxM)
    pi : initial state distribution (N)
    T : length of the sequence to sample
    
    Returns:
    sampled_states : list of sampled hidden states
    sampled_emissions : list of sampled emissions
    """
    # Step 1: Forward Filtering
    alpha = forward_filtering_no_observations(A, B, pi, T)
    
    # Step 2: Backward Sampling
    sampled_states, sampled_emissions = backward_sampling_no_observations(alpha, A, B, T)
    
    return sampled_states, sampled_emissions


# Create the mock args object
args = Args()


# Initialize classes with parsed arguments
input_handler = InputHandler(args)

config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)

# fragment_model = FragmentGraph(args.data_path, args.genotype_path, args.ploidy, input_handler.alleles)
# frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)

plot_graph(frag_graph)
print('Fragment Graph constructed.')

quotient_g = QuotientGraph(frag_graph).construct(fragment_list, input_handler, config)
plot_graph(quotient_g)
print('Quotient Graph constructed.')

# qg = chordal_contraction_cycle_base(quotient_g, fragment_list, input_handler, config)
qg = chordal_contraction(quotient_g, fragment_list, input_handler, config)
plot_graph(qg)    
print('Chordal Graph constructed.')



error_rate = 0.001
k = 2
m = 2
# state_names, transition_matrix, emission_prob_matrix, emission_index_map = generate_hmm_with_weights_and_emissions(qg, error_rate)
state_names, transition_matrix, emission_prob_matrix, emission_index_map = generate_hmm_with_weights_and_emissions_random(qg, error_rate, k, m)
pi = np.array([1/2 , 1/2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
T = 10
ffbs_no_observations(transition_matrix, emission_prob_matrix, pi, T)


