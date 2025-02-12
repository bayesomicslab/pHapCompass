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


def str_2_phas_1(phasing, ploidy):
    return np.array([int(p) for p in [*phasing]]).reshape(ploidy, -1)


def phas_2_str(phas):
    return ''.join([str(ph) for ph in list(np.ravel(phas))])


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




