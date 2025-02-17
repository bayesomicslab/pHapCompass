import os
import pandas as pd
import random
import sys
import numpy as np
# sys.path.append('/home/mok23003/BML/HaplOrbit')
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from collections import defaultdict, deque
from collections import Counter, defaultdict
import itertools
from utils.utils import *


def extract_positions(node_label):
    """
    Extracts positions from a node label of the form 'int-int'.
    """
    return set(map(int, node_label.split('-')))


def find_most_connected_positions(nodes, edges, phased_nodes):
    """
    Finds positions that are not part of the set but have the most edges to nodes in the set.
    
    Parameters:
    - nodes: list of node labels (e.g., ['1-2', '2-3', ...]).
    - edges: list of tuples representing edges between node labels (e.g., [('1-2', '2-3'), ...]).
    - node_set: a set of node labels to compare against.
    
    Returns:
    - A set of positions with the most connections to the input node_set.
    """
    
    
    # Convert node labels to positions in the node_set
    node_set_positions = set()
    for node in phased_nodes:
        node_set_positions.update(extract_positions(node))

    # Dictionary to count connections for positions outside the set
    position_counter = Counter()

    # Create an adjacency list for quick lookup of neighbors
    adjacency_list = defaultdict(list)
    for u, v in edges:
        adjacency_list[u].append(v)
        adjacency_list[v].append(u)

    # Iterate over all nodes in the node set
    for node in phased_nodes:
        for neighbor in adjacency_list[node]:
            # Skip neighbors that are already in the node_set
            if neighbor in phased_nodes:
                continue

            # Get positions of the neighbor node
            neighbor_positions = extract_positions(neighbor)

            # Count positions that are outside the set but connected to it
            new_positions = neighbor_positions - node_set_positions
            position_counter.update(new_positions)

    # Find positions with the highest number of connections
    max_connections_value = np.max(list(position_counter.values()))
    most_connected_positions = {pos for pos, count in position_counter.items() if count == max_connections_value}
    # Find nodes that contain the most connected positions and are not in the node_set
    most_connected_nodes = {
        node for node in nodes
        if extract_positions(node).intersection(most_connected_positions) and node not in phased_nodes
    }
    return most_connected_nodes
    

def select_max_prob_key(probabilities):
    """
    Selects the key with the highest probability from a dictionary.
    If multiple keys have the same max probability, selects one randomly.
    
    Parameters:
        probabilities (dict): A dictionary where keys are options and values are probabilities.

    Returns:
        The key with the highest probability.
    """
    if not probabilities:
        raise ValueError("The probabilities dictionary is empty.")

    # Get the maximum probability value
    # max_prob = max(probabilities.values())
    max_prob = np.max(list(probabilities.values()))

    # Find all keys with the maximum probability
    max_keys = [key for key, prob in probabilities.items() if prob == max_prob]

    # Randomly select a key if there are ties
    return random.choice(max_keys)


def permute_rows(a):
    # Generate all permutations of the rows
    perms = list(itertools.permutations(a))
    
    # Convert each permutation into a NumPy array and return as a list of matrices
    permuted_matrices = [np.array(p) for p in perms]
    
    return permuted_matrices

# @profile
def hstack_with_order(source_matrix, target_matrix, source_label, target_label):
    """
    Combines columns from source and target matrices based on the given labels.

    Parameters:
    - source_matrix: np.ndarray, source matrix with columns mapped by source_label
    - target_matrix: np.ndarray, target matrix with columns mapped by target_label
    - source_label: str, label with numbers associated with source_matrix columns
    - target_label: str, label with numbers associated with target_matrix columns

    Returns:
    - np.ndarray, combined matrix with columns arranged based on sorted unique labels
    """

    # Parse the labels and extract unique, sorted numbers
    source_keys = list(map(int, source_label.split('-')))
    target_keys = list(map(int, target_label.split('-')))
    all_keys = sorted(set(source_keys + target_keys))

    # Create dictionaries for source and target matrix columns
    source_dict = {key: source_matrix[:, i] for i, key in enumerate(source_keys)}
    target_dict = {key: target_matrix[:, i] for i, key in enumerate(target_keys)}

    # Merge the two dictionaries: prefer source_dict values if keys overlap
    combined_dict = {**target_dict, **source_dict}

    # Retrieve the columns based on the sorted keys
    combined_matrix = np.column_stack([combined_dict[key] for key in all_keys])

    return combined_matrix

# @profile
def find_phasings_matches(ff, sf, common_ff, common_sf, source_label, target_label):
    # Correct
    templates = []
    # byte_set = {a.tobytes() for a in templates}
    all_local = find_matchings(list(ff[:, common_ff]), list(sf[:, common_sf]))
    for al in all_local:
        ff_ordering = [ii[0] for ii in al]
        sf_ordering = [ii[1] for ii in al]
        assert any(ff[ff_ordering, common_ff] == sf[sf_ordering, common_sf])
        ordered_ff = ff[ff_ordering, :]
        ordered_sf = sf[sf_ordering, :]
        temp = hstack_with_order(ordered_ff, ordered_sf, source_label, target_label)

        # temp = np.hstack([ff[ff_ordering, :], sf[sf_ordering, 1:]])
        byte_set = {a.tobytes() for a in templates}
        if temp.tobytes() not in byte_set:
            templates.append(temp)
    return templates


def find_phasings_matches_generalized(ff, sf, common_ff, common_sf, source_label, target_label):
    # Correct
    templates = []
    # byte_set = {a.tobytes() for a in templates}

    all_local = find_matchings_generalized(list(ff[:, common_ff]), list(sf[:, common_sf]))
    # all_local = find_matchings(list(ff[:, common_ff]), list(sf[:, common_sf]))
    for al in all_local:
        ff_ordering = [ii[0] for ii in al]
        sf_ordering = [ii[1] for ii in al]
        # assert any(ff[ff_ordering, common_ff] == sf[sf_ordering, common_sf])
        ordered_ff = ff[ff_ordering, :]
        ordered_sf = sf[sf_ordering, :]
        temp = hstack_with_order(ordered_ff, ordered_sf, source_label, target_label)

        # temp = np.hstack([ff[ff_ordering, :], sf[sf_ordering, 1:]])
        byte_set = {a.tobytes() for a in templates}
        if temp.tobytes() not in byte_set:
            templates.append(temp)
    return templates


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


def generate_binary_combinations(length):
    return [list(x) for x in itertools.product([0, 1], repeat=length)]


def generate_all_possible_emissions(unique_node_lengths):
    all_emissions = []
    for length in unique_node_lengths:
        all_emissions.extend(generate_binary_combinations(length))
    return all_emissions


def transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model, config):
    transitions_dict = {}
    transitions_dict_extra = {}
    for edge in edges_map_quotient.keys():
        transitions_dict_extra[edge] = {}
        source = edges_map_quotient[edge][0]
        target = edges_map_quotient[edge][1]
        source_weights = quotient_g.graph.vertex_properties["v_weights"][source]['weight']
        target_weights = quotient_g.graph.vertex_properties["v_weights"][target]['weight']
        source_label = quotient_g.graph.vertex_properties["v_label"][source]
        target_label = quotient_g.graph.vertex_properties["v_label"][target]
        common_ff, common_sf = find_common_element_and_index(source_label, target_label)
        source_phasings = list(source_weights.keys())
        target_phasings = list(target_weights.keys())
        # transitions_dict = {'source': source_phasings, 'target': target_phasings}
        transitions_mtx = np.zeros((len(source_phasings), len(target_phasings)))

        for i, ffstr in enumerate(source_phasings):
            for j, sfstr in enumerate(target_phasings):
                transitions_dict_extra[edge][str(i) + '-' + str(j)] = {}
                transitions_dict_extra[edge][str(i) + '-' + str(j)]['source_phasing'] = ffstr
                transitions_dict_extra[edge][str(i) + '-' + str(j)]['target_phasing'] = sfstr
                transitions_dict_extra[edge][str(i) + '-' + str(j)]['matched_phasings'] = {}
                matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, ploidy), str_2_phas_1(sfstr, ploidy), common_ff, common_sf, source_label, target_label)
                sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)
                
                matched_phasings_str = list(set([phas_2_str(pm) for pm in sorted_phasings]))
                # print(i, ffstr, j, sfstr)
                # print('matched phasings:', matched_phasings_str, len(matched_phasings_str))
                # if len(matched_phasings_str) > 1:
                #     print('More than one matching phasing')
                #     # stop
                poss = sorted(list(set([int(ss) for ss in source_label.split('-')] + [int(tt) for tt in target_label.split('-')])))
                match_reads = get_matching_reads_for_positions([int(i) for i in poss], fragment_model.fragment_list)
                wei = 0
                for phas in matched_phasings_str:
                    this_phas_weight = 0
                    for indc, this_po, obs in match_reads:
                        this_phas_read_weight = compute_likelihood_generalized_plus(np.array(obs), str_2_phas_1(phas, ploidy), indc, list(range(len(indc))), 
                                                                   config.error_rate)
                        wei += this_phas_read_weight
                        this_phas_weight += this_phas_read_weight
                    transitions_dict_extra[edge][str(i) + '-' + str(j)]['matched_phasings'][phas] = this_phas_weight
                transitions_mtx[i, j] = wei
        transitions_mtx = transitions_mtx / transitions_mtx.sum(axis=1, keepdims=True)
        transitions_dict[edge] = transitions_mtx
    return transitions_dict, transitions_dict_extra


def emissions(ploidy, quotient_g, quotient_g_v_label_reversed, error_rate):
    emission_dict = {}
    # Calculate emissions for each state and populate the emission probability matrix
    for state in quotient_g_v_label_reversed.keys():
        emission_dict[state] = {}

        node_elements = state.split('-')  # For example, '1-2' -> ['1', '2']
        node_length = len(node_elements)  # Determine the number of elements in the node

        # Generate all possible combinations of 0 and 1 based on the node length
        possible_emissions = generate_binary_combinations(node_length)
        v = quotient_g_v_label_reversed[state]
        phasings = quotient_g.graph.vertex_properties["v_weights"][v]['weight'].keys()
        for phasing in phasings:
            emission_dict[state][phasing] = {}
            phasing_np = str_2_phas_1(phasing, ploidy)  # Compute phasing for the state key
            for emission in possible_emissions:
                likelihood = compute_likelihood(np.array(emission), phasing_np, error_rate)
                emission_dict[state][phasing][''.join([str(e) for e in emission])] = likelihood

    return emission_dict


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


def assign_evidence_to_states_and_transitions(nodes, edges, input_file):
    """
    Assigns evidence (line indices) from the input file to states (nodes) and transitions (edges).
    
    Parameters:
    - nodes: List of states (nodes) in the graph.
    - edges: List of transitions (edges) in the graph.
    - input_file: Path to the input file.
    
    Returns:
    - assignment_dict: Dictionary containing line indices assigned to states and transitions.
    """
    # Initialize the assignment dictionary
    assignment_dict = {
        'states': {node: [] for node in nodes},
        'transitions': {f"{s}--{t}": [] for s, t in edges}
    }

    # Open and process the input file line by line
    with open(input_file, 'r') as file:
        for line_idx, line in enumerate(file): 
            line = line.strip()
            parts = line.split()

            # Parse the number of parts in the read
            num_parts = int(parts[0])
            current_positions = []

            # Process each part of the read
            for i in range(num_parts):
                start_pos = int(parts[2 + i * 2])  # Starting position
                read = parts[3 + i * 2]  # Read sequence
                length = len(read)  # Length of the read
                current_positions.extend(range(start_pos, start_pos + length))  # Covered positions

            # Assign to states
            for node in nodes:
                node_positions = list(map(int, node.split('-')))
                if all(pos in current_positions for pos in node_positions):  # Check explicit positions
                    assignment_dict['states'][node].append(line_idx)

            # Assign to transitions
            for s, t in edges:
                s_positions = list(map(int, s.split('-')))
                t_positions = list(map(int, t.split('-')))
                all_positions = s_positions + t_positions
                if all(pos in current_positions for pos in all_positions):  # Check explicit positions
                    assignment_dict['transitions'][f"{s}--{t}"].append(line_idx)

    return assignment_dict


def extract_observations_for_state(state, indices, observation_file):
    """
    Extracts observations associated with a specific state from the observation file.
    
    Parameters:
    - state (str): The state (e.g., '1-2'), representing genomic positions covered by the node.
    - indices (list): Line indices from the observation file relevant to the state.
    - observation_file (str): Path to the observation file.
    
    Returns:
    - List of observations (strings) associated with the state (same length as indices).
    """
    # Parse positions covered by the state
    state_positions = list(map(int, state.split('-')))
    
    observations = []
    with open(observation_file, 'r') as file:
        lines = file.readlines()
        
        for idx in indices:
            # Read the specific line
            line = lines[idx].strip()
            parts = line.split()
            
            # Parse the read structure
            num_parts = int(parts[0])  # Number of parts in the read
            
            # For each part of the read, check if it overlaps with the state
            observation = ""
            for part_idx in range(num_parts):
                start_pos = int(parts[2 + part_idx * 2])  # Start position of the current part
                read = parts[3 + part_idx * 2]  # Read sequence
                
                # Extract the relevant substring for the state
                for pos in state_positions:
                    if start_pos <= pos < start_pos + len(read):
                        observation += read[pos - start_pos]
            
            observations.append(observation)
    
    return observations


# def compute_forward_backward_with_interfaces(slices, interfaces, observations, emission_dict, transitions_dict):
#     """
#     Compute forward and backward messages for FFBS in a DBN with arbitrary state spaces using interfaces.
    
#     Parameters:
#     - slices: Dictionary of slices {t: list of nodes in slice t}.
#     - interfaces: Dictionary of interfaces {"outgoing": {t: list of nodes}, "incoming": {t: list of nodes}}.
#     - observations: Dictionary of observations {t: list of observations relevant to slice t}.
#     - emission_dict: Hierarchical dictionary of emission probabilities for each node and its state space.
#     - transitions_dict: Hierarchical dictionary of transition probabilities for each edge.
    
#     Returns:
#     - forward_messages: Forward messages {t: {node: array of forward probabilities for each state-space value}}.
#     - backward_messages: Backward messages {t: {node: array of backward probabilities for each state-space value}}.
#     """
#     adjacency_list = defaultdict(list)
#     reverse_adjacency_list = defaultdict(list)
#     for s, t in edges:
#         adjacency_list[s].append(t)
#         reverse_adjacency_list[t].append(s)

#     forward_messages = {}
#     backward_messages = {}
    
#     # Initialize forward messages
#     for t in slices:
#         forward_messages[t] = {}
#         for node in slices[t]:  # Use outgoing interface
#             forward_messages[t][node] = {}
#             state_space = list(emission_dict[node].keys())
#             for phase in state_space:
#                 forward_messages[t][node][phase] = 0
#                 # forward_messages[t][node] = {phase: 0}
    
#     # Initialize backward messages
#     for t in slices:
#         backward_messages[t] = {}
#         for node in slices[t]:  # Use outgoing interface
#             backward_messages[t][node] = {}
#             state_space = list(emission_dict[node].keys())
#             for phase in state_space:
#                 backward_messages[t][node][phase] = 0

#     # Compute forward messages
#     # Base case
#     t = 1  
#     for node in slices[t]:
#         observations = extract_observations_for_state(node, assignment_dict['states'][node], frag_path)  # Compute for outgoing interface
#         state_space = list(emission_dict[node].keys())
#         for phase in state_space:
#             forward_messages[t][node][phase] += -1 * (np.sum([np.log(emission_dict[node][phase][r]) for r in observations]) + np.log(1/len(state_space)))

#     # Recursion
#     for t in range(2, len(slices) + 1):
#         for node in slices[t]:
#             observations = extract_observations_for_state(node, assignment_dict['states'][node], frag_path)  # Compute for outgoing interface
#             state_space = list(emission_dict[node].keys())
#             parents = reverse_adjacency_list[node] # interfaces["outgoing"][t-1]
#             parents_in_prev_t = list(set(slices[t-1]).intersection(set(parents)))
            
#             for i, phase in enumerate(state_space):
#                 log_likelihood = np.sum([np.log(emission_dict[node][phase][r]) for r in observations])
#                 transition_sum = 0
#                 for prev_node in parents_in_prev_t:
#                     prev_state_space = list(emission_dict[prev_node].keys())
#                     for j, prev_phase in enumerate(prev_state_space):

#                         transition_prob = transitions_dict[f"{prev_node}--{node}"][j, i]
#                         prev_meassage = forward_messages[t - 1][prev_node][prev_phase]
#                         # transition_prob = transitions_dict.get(f"{prev_node}--{node}", np.zeros((len(prev_state_space), len(state_space))))[j, i]
#                         # print(node, phase, prev_node, prev_phase, 'addition:', prev_meassage * transition_prob)
#                         transition_sum += prev_meassage * transition_prob
#                 forward_messages[t][node][phase] += -1 * log_likelihood * np.log(transition_sum)
    
#     # Compute for outgoing interface
#     # Base case
#     t = len(slices)
#     for node in slices[t]:  
#         state_space = list(emission_dict[node].keys())
#         for phase in state_space:  
#             backward_messages[t][node][phase] += 1

#     # Recursion
#     for t in range(len(slices) - 1, 0, -1):
#         # print(t)
#         for node in slices[t]:  # Compute for outgoing interface
#             state_space = list(emission_dict[node].keys())
#             children = adjacency_list[node]  # interfaces["incoming"][t+1]
#             children_in_next_t = list(set(slices[t+1]).intersection(set(children)))
#             # if t == len(slices):  # Base case
#             #     backward_messages[t][node] = np.ones(len(state_space))
#             # else:  # Recursion
#             # backward_sum = 0
#             for i, phase in enumerate(state_space):
#                 backward_sum = 0
#                 for next_node in children_in_next_t:
#                     backward_sum_in = 0
#                     observations = extract_observations_for_state(next_node, assignment_dict['states'][next_node], frag_path)
#                     next_state_space = list(emission_dict[next_node].keys())
#                     for j, next_phase in enumerate(next_state_space):
#                         transition_prob = transitions_dict[f"{node}--{next_node}"][i, j]
#                         bm = backward_messages[t + 1][next_node][next_phase]
#                         for r in observations:
                            
#                             # log_likelihood = np.sum([np.log(emission_dict[next_node][next_phase][r]) for r in observations])
#                             likelihood = emission_dict[next_node][next_phase][r]
#                             backward_sum_in += bm * likelihood * transition_prob
#                             print(node, phase, next_node, next_phase, r, 'addition:', bm * likelihood * transition_prob)

#                             # backward_sum_in += np.log(bm * likelihood * transition_prob)
#                         backward_sum += np.log(backward_sum_in)
#                 backward_messages[t][node][phase] += backward_sum
    
#     return forward_messages, backward_messages


def compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, frag_path):
    """
    Compute forward messages in log space for FFBS in a DBN with arbitrary state spaces.
    
    Parameters:
    - slices: Dictionary of slices {t: list of nodes in slice t}.
    - edges: List of directed edges (source, target).
    - assignment_dict: Dictionary mapping states to line indices of observations.
    - emission_dict: Hierarchical dictionary of emission probabilities for each node and its state space.
    - transitions_dict: Dictionary of transition probabilities (log-space).
    - frag_path: Path to the observation file.
    
    Returns:
    - forward_messages: Forward messages {t: {node: {state: log probability}}}.
    """
    # Build adjacency lists for reverse lookups
    reverse_adjacency_list = defaultdict(list)
    for s, t in edges:
        reverse_adjacency_list[t].append(s)
    
    forward_messages = {}

    # Initialize forward messages
    for t in slices:
        forward_messages[t] = {}
        for node in slices[t]:  # Nodes in the current slice
            forward_messages[t][node] = {}
            state_space = list(emission_dict[node].keys())
            for phase in state_space:
                forward_messages[t][node][phase] = -np.inf  # Initialize in log space to negative infinity

    # Compute forward messages
    # Base case (t = 1)
    t = 1
    for node in slices[t]:
        # Extract observations for the current state
        observations = extract_observations_for_state(node, assignment_dict['states'][node], frag_path)
        state_space = list(emission_dict[node].keys())
        for phase in state_space:
            # Compute log P(r | Z_1^{(i)}) for all observations
            log_emission = np.sum([np.log(emission_dict[node][phase][r]) for r in observations])
            # Uniform prior: log P(I_1)
            log_prior = np.log(1 / len(state_space))
            forward_messages[t][node][phase] = log_emission + log_prior

    # Recursion (t = 2, ..., T)
    for t in range(2, len(slices) + 1):
        for node in slices[t]:
            # Extract observations for the current state
            observations = extract_observations_for_state(node, assignment_dict['states'][node], frag_path)
            state_space = list(emission_dict[node].keys())
            parents = reverse_adjacency_list[node]  # Get parents of the current node
            parents_in_prev_t = list(set(slices[t - 1]).intersection(set(parents)))  # Filter parents in slice t-1
            
            for i, phase in enumerate(state_space):
                # Compute log P(r | Z_t^{(i)}) for all observations
                log_likelihood = np.sum([np.log(emission_dict[node][phase][r]) for r in observations])

                # Compute log-sum-exp over all transitions from parents
                log_transition_sum = -np.inf  # Initialize for log-sum-exp
                for prev_node in parents_in_prev_t:
                    prev_state_space = list(emission_dict[prev_node].keys())
                    for j, prev_phase in enumerate(prev_state_space):
                        transition_prob = transitions_dict[f"{prev_node}--{node}"][j, i]
                        prev_message = forward_messages[t - 1][prev_node][prev_phase]
                        log_transition_sum = np.logaddexp(
                            log_transition_sum,
                            prev_message + np.log(transition_prob)
                        )

                # Update forward message
                forward_messages[t][node][phase] = log_likelihood + log_transition_sum

    return forward_messages


def compute_backward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, frag_path):
    """
    Compute backward messages in log space for FFBS in a DBN with arbitrary state spaces.
    
    Parameters:
    - slices: Dictionary of slices {t: list of nodes in slice t}.
    - edges: List of directed edges (source, target).
    - assignment_dict: Dictionary mapping states to line indices of observations.
    - emission_dict: Hierarchical dictionary of emission probabilities for each node and its state space.
    - transitions_dict: Dictionary of transition probabilities (log-space).
    - frag_path: Path to the observation file.
    
    Returns:
    - backward_messages: Backward messages {t: {node: {state: log probability}}}.
    """
    # Build adjacency list for forward traversal
    adjacency_list = defaultdict(list)
    for s, t in edges:
        adjacency_list[s].append(t)

    backward_messages = {}

    # Initialize backward messages
    for t in slices:
        backward_messages[t] = {}
        for node in slices[t]:
            backward_messages[t][node] = {}
            state_space = list(emission_dict[node].keys())
            for phase in state_space:
                backward_messages[t][node][phase] = -np.inf  # Initialize to negative infinity in log space

    # Base case (t = T)
    t = len(slices)
    for node in slices[t]:
        state_space = list(emission_dict[node].keys())
        for phase in state_space:
            backward_messages[t][node][phase] = 0  # log(1) = 0

    # Recursion (t = T-1, ..., 1)
    for t in range(len(slices) - 1, 0, -1):
        for node in slices[t]:
            state_space = list(emission_dict[node].keys())
            children = adjacency_list[node]
            children_in_next_t = list(set(slices[t + 1]).intersection(set(children)))
            if len(children_in_next_t) == 0:
                for i, phase in enumerate(state_space):
                    backward_messages[t][node][phase] = 0
                continue

            for i, phase in enumerate(state_space):
                log_backward_sum = -np.inf  # Initialize for log-sum-exp

                for next_node in children_in_next_t:
                    next_state_space = list(emission_dict[next_node].keys())
                    observations = extract_observations_for_state(next_node, assignment_dict['states'][next_node], frag_path)

                    # Compute the inner summation in probability space
                    prob_sum = 0
                    for j, next_phase in enumerate(next_state_space):
                        # Multiply in probability space
                        transition_prob = transitions_dict[f"{node}--{next_node}"][i, j]
                        bm = np.exp(backward_messages[t + 1][next_node][next_phase])  # Convert back to probability space

                        # Compute emission probabilities for all observations
                        emission_prob = np.prod([emission_dict[next_node][next_phase][r] for r in observations])

                        # Add to the summation in probability space
                        prob_sum += bm * transition_prob * emission_prob

                    # Combine using log after summing in probability space
                    log_backward_sum = np.logaddexp(log_backward_sum, np.log(prob_sum))

                # Update backward message for the current phase
                backward_messages[t][node][phase] = log_backward_sum

    return backward_messages


def sample_states_no_resample_optimized(slices, edges, forward_messages, backward_messages, transitions_dict):
    """
    Sample states using the FFBS algorithm, avoiding resampling nodes already sampled in earlier slices,
    and leveraging precomputed forward and backward messages for efficiency.

    Parameters:
    - slices: Dictionary of slices {t: list of nodes in slice t}.
    - edges: List of directed edges (source, target).
    - forward_messages: Precomputed forward messages {t: {node: {state: log probability}}}.
    - backward_messages: Precomputed backward messages {t: {node: {state: log probability}}}.
    - transitions_dict: Dictionary of transition probabilities (log-space).

    Returns:
    - sampled_states: Dictionary of sampled states for each node in each slice.
    """
    sampled_states = {}
    already_sampled_nodes = set()

    # Step 1: Sample Q_T
    t = len(slices)
    sampled_states[t] = {}
    for node in slices[t]:
        if node in already_sampled_nodes:
            continue
        state_space = list(forward_messages[t][node].keys())
        log_probs = np.array([
            forward_messages[t][node][state] + backward_messages[t][node][state]
            for state in state_space
        ])
        probs = np.exp(log_probs - np.max(log_probs))  # Stabilize for numerical issues
        probs /= np.sum(probs)  # Normalize
        sampled_states[t][node] = random.choices(state_space, weights=probs, k=1)[0]
        already_sampled_nodes.add(node)

    # Step 2: Sample Q_t for t = T-1, ..., 1
    for t in range(len(slices) - 1, 0, -1):
        # print(t)
        sampled_states[t] = {}
        for node in slices[t]:
            if node in already_sampled_nodes:
                sampled_states[t][node] = sampled_states[t + 1][node]  # Use the sampled state from the next slice
                # print('comes here:', t, node)
                continue
            state_space = list(forward_messages[t][node].keys())
            children = [edge[1] for edge in edges if edge[0] == node]
            children_in_next_t = list(set(slices[t + 1]).intersection(set(children)))

            log_probs = []
            for state in state_space:
                log_alpha = forward_messages[t][node][state]
                log_child_contrib = 0

                # Compute contribution from children
                for child in children_in_next_t:
                    if child in already_sampled_nodes:
                        # Use the sampled state of the child
                        sampled_child_state = sampled_states[t + 1][child]
                        transition_prob = transitions_dict[f"{node}--{child}"][
                            state_space.index(state), list(forward_messages[t + 1][child].keys()).index(sampled_child_state)
                        ]
                        log_child_contrib += np.log(transition_prob)
                    else:
                        # Use backward message of the child
                        next_state_space = list(forward_messages[t + 1][child].keys())
                        child_sum = -np.inf  # Log-sum-exp initialization

                        for child_state in next_state_space:
                            transition_prob = transitions_dict[f"{node}--{child}"][
                                state_space.index(state), next_state_space.index(child_state)
                            ]
                            beta = backward_messages[t + 1][child][child_state]
                            child_sum = np.logaddexp(child_sum, beta + np.log(transition_prob))

                        log_child_contrib += child_sum

                log_probs.append(log_alpha + log_child_contrib)

            probs = np.exp(log_probs - np.max(log_probs))  # Stabilize for numerical issues
            probs /= np.sum(probs)  # Normalize
            sampled_states[t][node] = random.choices(state_space, weights=probs, k=1)[0]
            already_sampled_nodes.add(node)

    return sampled_states


def sample_states_ground_truth(slices, nodes, genotype_path):
    true_haplotypes = pd.read_csv(genotype_path).T
    sorted_nodes = sort_nodes(nodes)
    sampled_states = {}

    # **Step 1:** Sample Z_T from the final forward distribution (only uses alpha_T).
    t = len(slices)
    sampled_states[t] = {}
    for node in sorted_nodes:
        positions = [int(i)-1 for i in node.split('-')]
        true_phasing = true_haplotypes.loc[:, positions].values
        true_phasing_str = phas_2_str(true_phasing)
        sampled_states[t][node] = true_phasing_str
    return sampled_states


def sample_states_book(slices, edges, forward_messages, transitions_dict):
    """
    Sample states using forward messages and transition probabilities without backward messages.
    This matches the recursive sampling process described in the theory.
    
    Parameters:
    - slices: Dictionary of slices {t: list of nodes in slice t}.
    - edges: List of directed edges (source, target).
    - forward_messages: Precomputed forward messages {t: {node: {state: log probability}}}.
    - transitions_dict: Dictionary of transition probabilities.

    Returns:
    - sampled_states: Dictionary of sampled states for each node in each slice.
    """
    sampled_states = {}

    # **Step 1:** Sample Z_T from the final forward distribution (only uses alpha_T).
    t = len(slices)
    sampled_states[t] = {}
    for node in slices[t]:
        state_space = list(forward_messages[t][node].keys())
        log_probs = np.array([forward_messages[t][node][state] for state in state_space])
        probs = np.exp(log_probs - np.max(log_probs))  # Prevent numerical underflow
        probs /= np.sum(probs)  # Normalize
        sampled_states[t][node] = random.choices(state_space, weights=probs, k=1)[0]

    # **Step 2:** Recursively sample Z_t for t = T-1, ..., 1 based only on alpha_t and transitions
    for t in range(len(slices) - 1, 0, -1):
        sampled_states[t] = {}
        for node in slices[t]:
            state_space = list(forward_messages[t][node].keys())
            children_in_next_t = [child for _, child in edges if _ == node and child in slices[t + 1]]

            log_probs = []
            for state in state_space:
                log_alpha = forward_messages[t][node][state]
                log_transition_sum = 0

                # Add contributions from the sampled states of children
                for child in children_in_next_t:
                    sampled_child_state = sampled_states[t + 1][child]
                    transition_prob = transitions_dict[f"{node}--{child}"][
                        state_space.index(state), list(forward_messages[t + 1][child].keys()).index(sampled_child_state)
                    ]
                    log_transition_sum += np.log(transition_prob)

                # Final log probability of state
                log_probs.append(log_alpha + log_transition_sum)

            probs = np.exp(log_probs - np.max(log_probs))  # Prevent numerical underflow
            probs /= np.sum(probs)  # Normalize
            sampled_states[t][node] = random.choices(state_space, weights=probs, k=1)[0]

    return sampled_states


def sample_states(slices, edges, forward_messages, transitions_dict):
    """
    Sample states using forward messages and transition probabilities without backward messages.
    This matches the recursive sampling process described in the theory.
    """
    sampled_states = {}

    # **Step 1:** Sample Z_T from the final forward distribution (only uses alpha_T).
    t = len(slices)
    sampled_states[t] = {}
    for node in slices[t]:
        state_space = list(forward_messages[t][node].keys())
        log_probs = np.array([forward_messages[t][node][state] for state in state_space])
        probs = np.exp(log_probs - np.max(log_probs))  # Prevent numerical underflow
        probs /= np.sum(probs)  # Normalize
        sampled_states[t][node] = random.choices(state_space, weights=probs, k=1)[0]

    # **Step 2:** Recursively sample Z_t for t = T-1, ..., 1
    for t in range(len(slices) - 1, 0, -1):
        sampled_states[t] = {}
        for node in slices[t]:
            state_space = list(forward_messages[t][node].keys())
            child_node = [child for _, child in edges if _ == node and child in slices[t + 1]][0]  # Single child node
            sampled_child_state = sampled_states[t + 1][child_node]

            log_probs = []
            for state in state_space:
                log_alpha = forward_messages[t][node][state]
                transition_prob = transitions_dict[f"{node}--{child_node}"][
                    state_space.index(state), list(forward_messages[t + 1][child_node].keys()).index(sampled_child_state)
                ]
                log_probs.append(log_alpha + np.log(transition_prob))

            probs = np.exp(log_probs - np.max(log_probs))  # Prevent numerical underflow
            probs /= np.sum(probs)  # Normalize
            sampled_states[t][node] = random.choices(state_space, weights=probs, k=1)[0]

    return sampled_states


def get_combinations_with_probabilities(phased_nodes, relevant_edges, transitions_dict_extra, samples_brief, ploidy, predicted_haplotypes, config, fragment_model):
    """
    Collapse information from all phased nodes and their outgoing edges, 
    returning all valid phasing combinations and their probabilities.

    Parameters:
        phased_nodes (set): Nodes already phased.
        relevant_edges (list): Relevant edges for the selected position.
        transitions_dict_extra (dict): Contains phasing information for edges.
        samples_brief (dict): Phasing samples for nodes.
        ploidy (int): Ploidy level.   
        predicted_haplotypes (pd.DataFrame): Current phased haplotypes.
        config: Configuration object (e.g., error rate).

    Returns:
        list: Valid combinations of phasings.
        list: Probabilities for the combinations.
    """
    all_combinations = {}  # Store valid combinations
    probabilities = {}     # Store probabilities for each combination

    # Iterate over all phased nodes
    for source_node in phased_nodes:
        # Get edges originating from the source node
        all_combinations[source_node] = []
        probabilities[source_node] = []
        combination_set = set()

        outgoing_edges = [edge for edge in relevant_edges if edge[0] == source_node or edge[1] == source_node]
        
        for edge in outgoing_edges:
            source, target = edge if edge[0] == source_node else edge[::-1]
            # source, target = edge
            edge_label = f"{source}--{target}"

            this_edge_positions = sorted(set(int(pos) for part in [source, target] for pos in part.split('-')))

            match_reads = get_matching_reads_for_positions(this_edge_positions, fragment_model.fragment_list)

            if edge_label not in transitions_dict_extra:
                continue

            # Retrieve matched phasings for the edge
            edge_phasings = transitions_dict_extra[edge_label]
            matched_phasings = None
            for key, value in edge_phasings.items():
                if (
                    value["source_phasing"] == samples_brief[source] and
                    value["target_phasing"] == samples_brief[target]
                ):
                    matched_phasings = value["matched_phasings"]
                    break

            if not matched_phasings:
                continue

            # Get fixed positions and their values
            fixed_positions = [int(pos) for pos in source.split("-")]
            fixed_positions_indices = [int(pos) - 1 for pos in fixed_positions]
            fixed_values = predicted_haplotypes.loc[:, fixed_positions_indices].values

            # Get target node positions and phasing
            target_positions = [int(pos) for pos in target.split("-")]
            target_phasing = str_2_phas_1(samples_brief[target], ploidy)
            target_phasing_permutations = np.array(list(itertools.permutations(target_phasing)))

            target_positions_indices = [this_edge_positions.index(pos) for pos in target_positions if pos in this_edge_positions]


            # Iterate over all matched phasings for the edge
            for matchp in matched_phasings.keys():
                matched_phasings_np = str_2_phas_1(matchp, ploidy)

                for permuted_key in itertools.permutations(matched_phasings_np):
                    permuted_key_np = np.array(permuted_key)
                    # print(permuted_key_np)
                    # Check fixed positions match
                    if not np.array_equal(fixed_values, permuted_key_np[:, fixed_positions_indices]):
                        continue

                    # Check if rows match a permutation of the target node phasing
                    if not any(
                        np.array_equal(permuted_key_np[:, target_positions_indices], perm)
                        for perm in target_phasing_permutations):
                        continue
                    
                    # Flatten the array and convert to tuple for hashable representation
                    flattened_key = tuple(permuted_key_np.flatten())

                    # Skip if the combination is already processed
                    if flattened_key in combination_set:
                        continue

                    # Compute likelihood for the valid permutation
                    this_phas_weight = 0
                    for indc, this_po, obs in match_reads:
                        # Find shared indices between observed positions and the current phasing array
                        shared_indices = [this_edge_positions.index(pos) for pos in this_po if pos in this_edge_positions]

                        # Compute the likelihood for the shared positions
                        this_phas_read_weight = compute_likelihood_generalized_plus(
                            observed=np.array(obs),
                            phasing=permuted_key_np,
                            obs_pos=indc,            # Observed indices in the read
                            phas_pos=shared_indices, # Indices in the phasing array
                            error_rate=config.error_rate
                        )
                        this_phas_weight += this_phas_read_weight

                    combination_set.add(flattened_key)  # Add to set for future lookup

                    # Store the valid combination and its likelihood
                    all_combinations[source_node].append(permuted_key_np)
                    probabilities[source_node].append(this_phas_weight)

    return all_combinations, probabilities


def compute_candidate_phasings_target_based(
    selected_position, relevant_edges, phased_positions, transitions_dict_extra,
    samples_brief, ploidy, predicted_haplotypes, config, fragment_model
):
    """
    Compute candidate phasings and probabilities for nodes connected to the selected position.

    Parameters:
        selected_position (int): The position to be phased.
        relevant_edges (list): Edges connected to the selected position.
        phased_positions (set): Positions already phased.
        transitions_dict_extra (dict): Contains phasing information for edges.
        samples_brief (dict): Phasing samples for nodes.
        ploidy (int): Ploidy level.
        predicted_haplotypes (pd.DataFrame): Current phased haplotypes.
        config: Configuration object (e.g., error rate).
        fragment_model: Contains fragment information.

    Returns:
        dict: Candidate phasings for each node (keys are nodes, values are arrays of phasings).
        dict: Probabilities for each phasing (keys are nodes, values are lists of probabilities).
    """
    candidate_phasings = {}
    probabilities = {}

    # Find nodes containing the selected position
    target_nodes = {edge[1] for edge in relevant_edges if str(selected_position) in edge[1].split('-')}

    # Collect all positions involved in the relevant edges
    all_positions = sorted(set(
        int(pos) for edge in relevant_edges for node in edge for pos in node.split('-')
    ))

    for target_node in target_nodes:
        # Initialize candidate lists for the target node
        candidate_phasings[target_node] = []
        probabilities[target_node] = []

        # Reset combination_set for this target node
        combination_set = set()

        # Get all edges connected to the target node
        incoming_edges = [edge for edge in relevant_edges if edge[1] == target_node or edge[0] == target_node]

        for edge in incoming_edges:
            source, target = edge if edge[1] == target_node else edge[::-1]
            edge_label = f"{source}--{target}"

            if edge_label not in transitions_dict_extra:
                continue

            edge_phasings = transitions_dict_extra[edge_label]
            matched_phasings = None
            for key, value in edge_phasings.items():
                if (
                    value["source_phasing"] == samples_brief[source] and
                    value["target_phasing"] == samples_brief[target]
                ):
                    matched_phasings = value["matched_phasings"]
                    break

            if not matched_phasings:
                continue

            # Compute positions and values specific to this edge
            this_edge_positions = sorted(set(int(pos) for part in [source, target] for pos in part.split('-')))
            fixed_positions = [pos for pos in this_edge_positions if pos != selected_position and pos in phased_positions]
            fixed_positions_indices = [this_edge_positions.index(pos) for pos in fixed_positions]
            fixed_values = predicted_haplotypes.loc[:, [pos - 1 for pos in fixed_positions]].values

            # Retrieve matching reads for all positions in these edges
            match_reads = get_matching_reads_for_positions(this_edge_positions, fragment_model.fragment_list)

            # Iterate over matched phasings for this edge
            for matchp in matched_phasings.keys():
                matched_phasings_np = str_2_phas_1(matchp, ploidy)

                for permuted_key in itertools.permutations(matched_phasings_np):
                    permuted_key_np = np.array(permuted_key)

                    # Ensure fixed positions match for this edge
                    if not np.array_equal(fixed_values, permuted_key_np[:, fixed_positions_indices]):
                        continue

                    # Ensure at least one permutation of the target phasing matches
                    target_positions_indices = [
                        all_positions.index(pos) for pos in [int(i) for i in target.split('-')]
                    ]
                    target_phasing = str_2_phas_1(samples_brief[target], ploidy)
                    target_phasing_permutations = np.array(list(itertools.permutations(target_phasing)))
                    if not any(
                        np.array_equal(permuted_key_np[:, target_positions_indices], perm)
                        for perm in target_phasing_permutations
                    ):
                        continue

                    # Create a full array with NaN for all positions
                    aligned_key_np = np.full((ploidy, len(all_positions)), np.nan)
                    this_edge_indices = [all_positions.index(pos) for pos in this_edge_positions]
                    aligned_key_np[:, this_edge_indices] = permuted_key_np

                    # Flatten and deduplicate
                    flattened_key = tuple(aligned_key_np.flatten())
                    if flattened_key in combination_set:
                        continue

                    # Compute likelihood
                    this_phas_weight = 0
                    for indc, this_po, obs in match_reads:
                        shared_indices = [all_positions.index(pos) for pos in this_po if pos in all_positions]
                        this_phas_read_weight = compute_likelihood_generalized_plus(
                            observed=np.array(obs),
                            phasing=aligned_key_np,
                            obs_pos=indc,
                            phas_pos=shared_indices,
                            error_rate=config.error_rate
                        )
                        this_phas_weight += this_phas_read_weight

                    # Store valid phasing and likelihood
                    candidate_phasings[target_node].append(aligned_key_np)
                    probabilities[target_node].append(this_phas_weight)
                    combination_set.add(flattened_key)

    return candidate_phasings, probabilities


def compute_candidate_phasings0(
    selected_position, relevant_edges, phased_positions, transitions_dict_extra,
    samples_brief, ploidy, predicted_haplotypes, config, fragment_model
):
    """
    Compute candidate phasings and probabilities for nodes connected to the selected position.

    Parameters:
        selected_position (int): The position to be phased.
        relevant_edges (list): Edges connected to the selected position.
        phased_positions (set): Positions already phased.
        transitions_dict_extra (dict): Contains phasing information for edges.
        samples_brief (dict): Phasing samples for nodes.
        ploidy (int): Ploidy level.
        predicted_haplotypes (pd.DataFrame): Current phased haplotypes.
        config: Configuration object (e.g., error rate).
        fragment_model: Contains fragment information.

    Returns:
        list: All candidate phasings (NumPy arrays of uniform size).
        list: Probabilities associated with each candidate phasing.
    """
    all_combinations = []
    all_probabilities = []
    candidate_counts = []
    # Collect all positions involved in the relevant edges
    all_positions = sorted(set(
        int(pos) for edge in relevant_edges for node in edge for pos in node.split('-')
    ))

    # Reset combination_set for deduplication across all edges
    combination_set = set()

    for edge in relevant_edges:
        source, target = edge
        edge_label = f"{source}--{target}"

        if edge_label not in transitions_dict_extra:
            continue

        edge_phasings = transitions_dict_extra[edge_label]
        matched_phasings = None
        for key, value in edge_phasings.items():
            if (
                value["source_phasing"] == samples_brief[source] and
                value["target_phasing"] == samples_brief[target]
            ):
                matched_phasings = value["matched_phasings"]
                break

        if not matched_phasings:
            continue

        # Compute positions and values specific to this edge
        this_edge_positions = sorted(set(int(pos) for part in [source, target] for pos in part.split('-')))
        fixed_positions = [pos for pos in this_edge_positions if pos != selected_position and pos in phased_positions]
        fixed_positions_indices = [this_edge_positions.index(pos) for pos in fixed_positions]
        fixed_values = predicted_haplotypes.loc[:, [pos - 1 for pos in fixed_positions]].values

        # Precompute target phasing permutations (fixed for the edge)
        target_positions_indices = [
            this_edge_positions.index(pos) for pos in [int(i) for i in target.split('-')] if pos in this_edge_positions
        ]
        target_phasing = str_2_phas_1(samples_brief[target], ploidy)
        target_phasing_permutations = np.array(list(itertools.permutations(target_phasing)))

        # Retrieve matching reads for all positions in these edges
        match_reads = get_matching_reads_for_positions(this_edge_positions, fragment_model.fragment_list)

        # Iterate over matched phasings for this edge
        for matchp in matched_phasings.keys():
            matched_phasings_np = str_2_phas_1(matchp, ploidy)

            for permuted_key in itertools.permutations(matched_phasings_np):
                permuted_key_np = np.array(permuted_key)

                # Ensure fixed positions match for this edge
                if not np.array_equal(fixed_values, permuted_key_np[:, fixed_positions_indices]):
                    continue

                # Ensure at least one permutation of the target phasing matches
                if not any(
                    np.array_equal(permuted_key_np[:, target_positions_indices], perm)
                    for perm in target_phasing_permutations
                ):
                    continue

                # Compute likelihood
                this_phas_weight = 0
                for indc, this_po, obs in match_reads:
                    shared_indices = [this_edge_positions.index(pos) for pos in this_po if pos in this_edge_positions]
                    this_phas_read_weight = compute_likelihood_generalized_plus(
                        observed=np.array(obs),
                        phasing=permuted_key_np,
                        obs_pos=indc,
                        phas_pos=shared_indices,
                        error_rate=config.error_rate
                    )
                    this_phas_weight += this_phas_read_weight

                # Create a full array with NaN for all positions
                aligned_key_np = np.full((ploidy, len(all_positions)), np.nan)

                # Fill positions in the current edge
                this_edge_indices = [all_positions.index(pos) for pos in this_edge_positions]
                aligned_key_np[:, this_edge_indices] = permuted_key_np

                # Fill remaining columns with values from predicted_haplotypes
                remaining_positions = [pos for pos in all_positions if pos not in this_edge_positions]
                for pos in remaining_positions:
                    if pos in phased_positions:
                        col_index = all_positions.index(pos)
                        aligned_key_np[:, col_index] = predicted_haplotypes.loc[:, pos - 1].values

                # Flatten and deduplicate
                flattened_key = tuple(aligned_key_np.flatten())
                if flattened_key in combination_set:
                    continue

                # Store valid phasing and likelihood
                all_combinations.append(aligned_key_np)
                all_probabilities.append(this_phas_weight)
                combination_set.add(flattened_key)

    # Count how many target nodes match each candidate
    for candidate in all_combinations:
        count = 0
        for edge in relevant_edges:
            target = edge[1]
            target_positions = [int(pos) for pos in target.split('-')]
            target_indices = [all_positions.index(pos) for pos in target_positions]
            target_phasing = str_2_phas_1(samples_brief[target], ploidy)
            target_phasing_permutations = np.array(list(itertools.permutations(target_phasing)))
            if any(
                np.array_equal(candidate[:, target_indices], perm)
                for perm in target_phasing_permutations
            ):
                count += 1
        candidate_counts.append(count)

    return all_combinations, all_probabilities, candidate_counts, all_positions


def compute_candidate_phasings(
    selected_position,            # Not really used anymore for multi-position
    relevant_edges,
    phased_positions,
    transitions_dict_extra,
    samples_brief,
    ploidy,
    predicted_haplotypes,
    config,
    fragment_model
):
    """
    Compute candidate phasings for *all newly introduced positions* within 'relevant_edges',
    using an iterative "node-by-node" approach. We maintain a list of partial solutions,
    each a (ploidy x len(all_positions)) array. We fill columns for new nodes one by one.
    
    At the end, each candidate solution is fully assigned (no NaN columns) for the subgraph.
    Then we compute the total likelihood and return them.
    
    Parameters:
        selected_position (int): Old param, not used here. Kept for signature compatibility.
        relevant_edges (list): List of edges describing the subgraph of newly introduced nodes
                               plus edges to phased nodes.
        phased_positions (set): Positions already phased. We won't alter those columns.
        transitions_dict_extra (dict): Weighted phasing constraints for edges. 
                                       { "node1--node2": {...} }
        samples_brief (dict): Map node -> string. The "base" phasing for each node.
        ploidy (int): Ploidy level => # rows in the final alignment array.
        predicted_haplotypes (pd.DataFrame): Already phased haplotypes table, shape:
                         (#haplotypes = ploidy) x (#positions).
        config: Contains e.g. error_rate.
        fragment_model: For retrieving matching reads.

    Returns:
        all_combinations (list of np.ndarray):
            Each is shape (ploidy, len(all_positions)) with no NaNs.
        all_probabilities (list of float):
            Sums of read-likelihood for each final candidate.
        candidate_counts (list of int):
            How many edges' target phasing is matched by each candidate (old measure).
        all_positions (list of int):
            The sorted union of positions in this subgraph.
    """

    # -----------------------------
    # 1) Identify subgraph nodes & positions
    # -----------------------------
    subgraph_nodes = set()
    for (n1, n2) in relevant_edges:
        subgraph_nodes.add(n1)
        subgraph_nodes.add(n2)

    all_positions = sorted({
        int(pos)
        for (n1, n2) in relevant_edges
        for node in (n1, n2)
        for pos in node.split('-')
    })

    # Map each node -> sorted list of integer positions
    node_positions_map = {}
    for nd in subgraph_nodes:
        nd_positions = sorted(int(p) for p in nd.split('-'))
        node_positions_map[nd] = nd_positions

    # Distinguish which subgraph nodes are "already phased" vs "new"
    # A node is fully phased if all its positions in phased_positions
    def node_is_fully_phased(nd):
        return all(pos in phased_positions for pos in node_positions_map[nd])

    # We'll pick an ordering of subgraph nodes. E.g. sorted by node name:
    subgraph_nodes_sorted = sorted(subgraph_nodes)

    # For quick "edge label" check
    def get_edge_label(a, b):
        # In your data, it might be "a--b" only if (a<b) or so. But let's be
        # consistent with how transitions_dict_extra was built. We'll check both:
        if f"{a}--{b}" in transitions_dict_extra:
            return f"{a}--{b}"
        if f"{b}--{a}" in transitions_dict_extra:
            return f"{b}--{a}"
        return None

    # -----------------------------
    # 2) Prepare "initial" partial solution(s)
    # We store them as (ploidy x len(all_positions)) arrays. The columns for
    # already phased positions are filled from predicted_haplotypes, the rest are NaN.
    # We'll keep a list "candidate_solutions" of these arrays.
    # Start with exactly one partial solution (with only phased columns filled).
    # -----------------------------
    init_array = np.full((ploidy, len(all_positions)), np.nan, dtype=float)

    # Fill from predicted_haplotypes for columns that are in phased_positions
    for col_i, pos in enumerate(all_positions):
        if pos in phased_positions:
            # Copy from predicted_haplotypes
            init_array[:, col_i] = predicted_haplotypes.loc[:, pos - 1].values

    candidate_solutions = [init_array]

    # -----------------------------
    # 3) Helper to get "all permutations" for a node's base phasing
    #    This duplicates old logic where we permute the rows of str_2_phas_1(...)
    # -----------------------------
    def get_node_all_permutations(node_str):
        """Return list of 2D arrays shape (ploidy, #positions_in_node) for all row permutations."""
        base_np = str_2_phas_1(node_str, ploidy)
        out = []
        for perm in itertools.permutations(base_np):
            arr_ = np.array(perm)
            out.append(arr_)
        return out

    # -----------------------------
    # 4) Helper: For an "edge" constraint, we want to see if the node1 partial columns
    #    are consistent with node2's candidate block. We can also check transitions_dict_extra
    #    for matched phasings. This is somewhat simpler if we only keep solutions that
    #    match the base "samples_brief[node]" anyway, as in older code.
    # -----------------------------
    def is_compatible_edge(nd1, nd2, full_sol):
        """
        Return True if nd1<->nd2 is consistent under transitions_dict_extra + the
        current partial columns for nd1 and nd2 in full_sol.
        We'll check:
         1) Each node's assigned sub-block vs samples_brief[node]
         2) The edge must exist => find matched phasings => see if there's a row assignment that matches.
        A simpler approach is "If the node array in full_sol is the same row arrangement as samples_brief[node], then check transitions_dict_extra edge if it is present".
        But that can get tricky quickly...
        
        For now, we replicate old logic: the "source_phasing" must be samples_brief[nd1], "target_phasing" must be samples_brief[nd2]. Then check if there's a 'matched_phasings'.
        We'll see if the row pattern in full_sol for nd1,nd2 matches one of those 'matched_phasings'.
        """
        edge_label = get_edge_label(nd1, nd2)
        if edge_label is None:
            return True  # No constraints => OK

        if edge_label not in transitions_dict_extra:
            return True  # No constraints in the dictionary => OK

        # If transitions_dict_extra[edge_label] expects e.g. "source_phasing": samples_brief[nd1], "target_phasing": ...
        # we see if there's a matched_phasings
        edge_phasings = transitions_dict_extra[edge_label]
        matched_phasings = None
        for k, v in edge_phasings.items():
            if (v["source_phasing"] == samples_brief[nd1]
                and v["target_phasing"] == samples_brief[nd2]):
                matched_phasings = v["matched_phasings"]
                break
        if not matched_phasings:
            # Means we don't have a recognized pairing => can't be consistent
            return False

        # Now we see if the actual row arrangement in full_sol for nd1, nd2 matches one of matched_phasings
        # We'll gather the sub-block for nd1 from full_sol
        nd1_pos = node_positions_map[nd1]
        nd2_pos = node_positions_map[nd2]

        # shape: (ploidy, #pos_in_nd1)
        nd1_sub = np.zeros((ploidy, len(nd1_pos)), dtype=int)
        for j, p in enumerate(nd1_pos):
            col_idx = all_positions.index(p)
            nd1_sub[:, j] = full_sol[:, col_idx]

        # shape: (ploidy, #pos_in_nd2)
        nd2_sub = np.zeros((ploidy, len(nd2_pos)), dtype=int)
        for j, p in enumerate(nd2_pos):
            col_idx = all_positions.index(p)
            nd2_sub[:, j] = full_sol[:, col_idx]

        # We'll check if nd1_sub is a row permutation of str_2_phas_1(samples_brief[nd1]), etc.
        # Then also if there's a key in matched_phasings that matches nd1_sub => nd2_sub arrangement.
        # But old code used permutations for the sub-block. It's a bit big. We'll do a simpler approach: 
        # We'll invert the logic: if old code was used, the final arrangement is indeed 1 of the permutations that was accepted. So we only confirm that some "mp_key" => row permutation => matches nd2_sub. That can be done, but let's do a direct approach:

        # We'll flatten nd1_sub, nd2_sub. We see if there's an mp_key that matches nd2_sub after we interpret mp_key with row perms, etc.
        # This is fairly elaborate. Possibly we can do a "loose check" or skip. 
        # For brevity, let's do a minimal check: If the node arrays align with the base sample strings => proceed. If not => fail.

        arr_from_brief_nd1 = str_2_phas_1(samples_brief[nd1], ploidy)
        if not any(np.array_equal(nd1_sub, np.array(perm)) for perm in itertools.permutations(arr_from_brief_nd1)):
            return False

        arr_from_brief_nd2 = str_2_phas_1(samples_brief[nd2], ploidy)
        if not any(np.array_equal(nd2_sub, np.array(perm)) for perm in itertools.permutations(arr_from_brief_nd2)):
            return False

        # If we got here => they match the base phasing. We'll skip checking "matched_phasings" as it is quite big. Let's do a direct approach:
        # Actually let's do a minimal "the presence of matched_phasings means it's possible, so we return True".
        return True

    # -----------------------------
    # 5) The main iterative approach:
    #    We'll gather the "new" subgraph nodes that are not fully phased, in sorted order, then do
    #    a node-by-node extension of candidate_solutions.
    # -----------------------------
    new_nodes = [nd for nd in subgraph_nodes_sorted if not node_is_fully_phased(nd)]
    # The already-phased nodes are essentially locked in candidate_solutions (since we filled columns from predicted_haplotypes).

    def integrate_node_into_candidates(nd, current_candidates):
        """
        For each array in current_candidates, attempt to fill columns for node `nd`
        with all row permutations of samples_brief[nd] that are consistent w.r.t edges
        linking nd to the subgraph nodes that are already assigned in that partial solution.
        We'll produce a new candidate list in the end.
        """
        node_pos = node_positions_map[nd]
        # base permutations for node
        node_candidates = get_node_all_permutations(samples_brief[nd])

        new_candidates = []
        for cand_arr in current_candidates:
            # We only fill columns for node_pos if they're currently NaN
            # If they're not NaN => skip or check? Let's forcibly override if they're all NaN or possibly partial. We'll do a check:
            to_fill_cols = [all_positions.index(p) for p in node_pos]
            # We'll find which subgraph nodes are already integrated => any node that is fully phased or has been processed in new_nodes up to now => 
            # Actually we don't store that info. We'll do a simpler approach: just after we fill nd, we check consistency with all subgraph nodes that appear in relevant_edges with nd. 
            for node_phasing in node_candidates:
                # Make a copy of cand_arr
                arr_copy = np.copy(cand_arr)
                # Fill in the columns for nd
                for col_i, pos_i in enumerate(to_fill_cols):
                    arr_copy[:, pos_i] = node_phasing[:, col_i]
                # Now check edge constraints for every edge (nd, x) in subgraph
                # Actually we only need to check subgraph nodes that are either fully phased or already integrated, but let's do a simple approach: check all. If they're not assigned => some columns are NaN => might pass trivially or fail. We'll see.
                # We'll do a for e in relevant_edges: if nd is in e, let's see other node in e => is_compatible_edge
                # but is_compatible_edge expects no NaNs. If the other node is also unfilled => no constraints. We'll do a quick "if not np.isnan" check or node_is_fully_phased?
                # We'll build a function to see if the other node is assigned in arr_copy
                # Actually let's do a small function:
                def node_assigned_in_array(nodeX, arrX):
                    nodeX_pos = node_positions_map[nodeX]
                    col_idxs = [all_positions.index(p) for p in nodeX_pos]
                    # if none is nan => assigned
                    return not np.isnan(arrX[:, col_idxs]).any()

                # Now let's test edges:
                consistent_flag = True
                for (a,b) in relevant_edges:
                    if nd in (a,b):
                        other = b if (a==nd) else a
                        # if other is assigned in arr_copy => check is_compatible_edge
                        if node_assigned_in_array(other, arr_copy):
                            if not is_compatible_edge(nd, other, arr_copy):
                                consistent_flag = False
                                break
                if consistent_flag:
                    new_candidates.append(arr_copy)

        return new_candidates

    candidate_solutions_current = candidate_solutions  # the partial solutions from the start (just phased cols)

    for nd in new_nodes:
        candidate_solutions_current = integrate_node_into_candidates(nd, candidate_solutions_current)
        # at this point, node nd is assigned in all solutions in candidate_solutions_current

    # candidate_solutions_current now has fully assigned columns for all new nodes
    # but we still have to compute the read-likelihood and do the "candidate_counts" measure.

    all_combinations = []
    all_probabilities = []
    candidate_counts = []

    # We'll define a function to compute full-likelihood from read coverage
    def compute_solution_likelihood(arr_):
        # sum over all reads that cover these all_positions
        match_reads = get_matching_reads_for_positions(all_positions, fragment_model.fragment_list)
        tot_like = 0.0
        for (indc, these_pos, obs) in match_reads:
            shared_indices = [all_positions.index(p) for p in these_pos if p in all_positions]
            ll = compute_likelihood_generalized_plus(
                observed=np.array(obs),
                phasing=arr_,
                obs_pos=indc,
                phas_pos=shared_indices,
                error_rate=config.error_rate
            )
            tot_like += ll
        return tot_like

    for sol_arr in candidate_solutions_current:
        # ensure no columns are NaN
        if np.isnan(sol_arr).any():
            # skip partial solution
            continue

        # compute solution likelihood
        sol_like = compute_solution_likelihood(sol_arr)
        all_combinations.append(sol_arr)
        all_probabilities.append(sol_like)

        # replicate old "candidate_count" measure
        cc = 0
        for edge in relevant_edges:
            target = edge[1]
            target_positions = [int(pos) for pos in target.split('-')]
            target_indices = [all_positions.index(x) for x in target_positions]
            target_phasing = str_2_phas_1(samples_brief[target], ploidy)
            target_permutations = np.array(list(itertools.permutations(target_phasing)))
            if any(np.array_equal(sol_arr[:, target_indices], perm) for perm in target_permutations):
                cc += 1
        candidate_counts.append(cc)

    return all_combinations, all_probabilities, candidate_counts, all_positions



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


def sample_states_book_multiple_times(slices, edges, forward_messages, transitions_dict, n=10):
    """
    Samples multiple times and computes the consensus (mode) for each node.

    Parameters:
    - slices: Dictionary of slices {t: list of nodes in slice t}.
    - edges: List of directed edges (source, target).
    - forward_messages: Precomputed forward messages {t: {node: {state: log probability}}}.
    - transitions_dict: Dictionary of transition probabilities.
    - n: Number of times to sample.

    Returns:
    - consensus_samples: Dictionary of consensus samples after multiple runs.
    """
    # Store multiple sampled dictionaries
    all_samples = [sample_states_book(slices, edges, forward_messages, transitions_dict) for _ in range(n)]

    # Compute consensus sample
    consensus_samples = {}

    # Iterate over slices (5, 4, 3, ...)
    for t in slices.keys():
        consensus_samples[t] = {}

        # Iterate over nodes in the current slice
        for node in slices[t]:
            sampled_values = [all_samples[i][t][node] for i in range(n)]  # Collect samples across all runs
            most_common_sample = Counter(sampled_values).most_common(1)[0][0]  # Get mode
            consensus_samples[t][node] = most_common_sample  # Store consensus

    return consensus_samples


def predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities"):
    """
    Predict haplotypes iteratively by phasing one variant at a time using candidate sampling and selection.
    """
    # Step 1: Initialization
    phasing_samples = {nn: samples[t][nn] for t in samples for nn in samples[t]}
    sorted_nodes = sort_nodes(nodes)
    genotype_df = pd.read_csv(genotype_path).T
    predicted_haplotypes = pd.DataFrame(index=[f'haplotype_{p+1}' for p in range(ploidy)], columns=genotype_df.columns)

    # Initialize phased and unphased sets
    phased_nodes = set()
    phased_positions = set()
    unphased_nodes = set(nodes)
    unphased_positions = {int(pos) for node in nodes for pos in node.split('-')}

    # Start with the first node
    first_node = sorted_nodes[0]
    initial_positions = sorted({int(pos) for pos in first_node.split('-')})
    initial_positions = [p - 1 for p in initial_positions]

    predicted_haplotypes.loc[:, initial_positions] = str_2_phas_1(phasing_samples[first_node], ploidy)
    phased_nodes.add(first_node)
    unphased_nodes.remove(first_node)
    phased_positions.update(map(int, first_node.split('-')))
    unphased_positions -= set(map(int, first_node.split('-')))

    # Step 2: Iterative Phasing
    while unphased_positions:
        # Find neighbors and relevant edges
        neighbor_nodes = {n2 for n1, n2 in edges if n1 in phased_nodes and n2 in unphased_nodes} | \
                         {n1 for n1, n2 in edges if n2 in phased_nodes and n1 in unphased_nodes}

        position_connections = defaultdict(lambda: {"count": 0, "edges": []})
        for n1, n2 in edges:
            if (n1 in phased_nodes and n2 in neighbor_nodes) or (n2 in phased_nodes and n1 in neighbor_nodes):
                for pos in map(int, n1.split('-') + n2.split('-')):
                    if pos in unphased_positions:
                        position_connections[pos]["count"] += 1
                        position_connections[pos]["edges"].append((n1, n2))

        # Check if there are no connections
        if not position_connections:
            # print(f"No connections found for unphased positions: {unphased_positions}")
            # Select the first unphased node from the sorted list
            for next_node in sorted_nodes:
                if next_node in unphased_nodes:
                    next_positions = sorted({int(pos) for pos in next_node.split('-')})
                    next_positions = [p - 1 for p in next_positions]

                    # Extract phasing for the new node
                    next_node_sample = phasing_samples[next_node]
                    next_sample_np = str_2_phas_1(next_node_sample, ploidy)

                    # Update predicted haplotypes
                    predicted_haplotypes.loc[:, next_positions] = next_sample_np
                    phased_nodes.add(next_node)
                    unphased_nodes.remove(next_node)
                    phased_positions.update(map(int, next_node.split('-')))
                    unphased_positions -= set(map(int, next_node.split('-')))
                    break
            continue

        # Select the position with the most connections
        selected_position = max(position_connections, key=lambda p: position_connections[p]["count"])
        relevant_edges = position_connections[selected_position]["edges"]

        # Compute candidates
        all_combinations, all_probabilities, candidate_counts, all_positions = compute_candidate_phasings(
            selected_position, relevant_edges, phased_positions, transitions_dict_extra,
            phasing_samples, ploidy, predicted_haplotypes, config, fragment_model)

        # Select and update
        if all_combinations:
            candidates = list(zip(all_combinations, all_probabilities, candidate_counts))
            best_candidate = select_best_candidate(candidates, prioritize=priority)
            best_combination, best_probability, best_count = best_candidate

            selected_position_index = all_positions.index(selected_position)
            predicted_haplotypes.loc[:, selected_position - 1] = best_combination[:, selected_position_index]
        else:
            # Fallback: Use target phasing
            for edge in relevant_edges:
                target_node = edge[1]
                if str(selected_position) in target_node.split('-'):
                    target_positions = [int(pos) for pos in target_node.split('-')]
                    target_phasing = str_2_phas_1(phasing_samples[target_node], ploidy)
                    target_position_index = target_positions.index(selected_position)
                    predicted_haplotypes.loc[:, selected_position - 1] = target_phasing[:, target_position_index]
                    break

        # Update phased sets
        for edge in relevant_edges:
            node1, node2 = edge
            if node1 in neighbor_nodes:
                phased_nodes.add(node1)
                unphased_nodes.discard(node1)
            if node2 in neighbor_nodes:
                phased_nodes.add(node2)
                unphased_nodes.discard(node2)

        phased_positions.add(selected_position)
        unphased_positions.remove(selected_position)

    return predicted_haplotypes


def predict_haplotypes_multiple_variants(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities"):
    """
    Algorithm 1: Iterative selection of positions to be phased.

    This function finds unphased positions connected to the current phased set and
    moves them into the phased sets. 

    Parameters:
        nodes (list): All nodes in the graph (e.g., ['1-2','2-3',...]).
        edges (list): List of edges (tuples) connecting nodes.
        ploidy (int): Ploidy level.
        samples (dict): A map node->string_phasing, for reference only here.
        sorted_nodes (list): Sorted node labels for fallback when no connections exist.

    Returns:
        phased_nodes (set),
        unphased_nodes (set),
        phased_positions (set),
        unphased_positions (set)
    """

        # Step 1: Initialization
    phasing_samples = {nn: samples[t][nn] for t in samples for nn in samples[t]}
    sorted_nodes = sort_nodes(nodes)
    genotype_df = pd.read_csv(genotype_path).T
    predicted_haplotypes = pd.DataFrame(index=[f'haplotype_{p+1}' for p in range(ploidy)], columns=genotype_df.columns)

    phased_nodes = set()
    unphased_nodes = set(nodes)
    phased_positions = set()
    unphased_positions = {int(pos) for node in nodes for pos in node.split('-')}

    # Start by picking a first node from sorted_nodes, same as original code
    if sorted_nodes:
        first_node = sorted_nodes[0]
        phased_nodes.add(first_node)
        unphased_nodes.discard(first_node)

        # Example: for a node '3-4', we turn them into integers {3,4} then mark them phased
        first_positions = {int(pos) for pos in first_node.split('-')}
        phased_positions.update(first_positions)
        unphased_positions -= first_positions

    # -------------------------------------------------------------------------
    # MAIN WHILE LOOP: keep going until all positions are phased
    # -------------------------------------------------------------------------
    while unphased_positions:

        # Identify neighbor_nodes: nodes that are unphased but connected to at least one phased node
        neighbor_nodes = (
            {n2 for n1, n2 in edges if n1 in phased_nodes and n2 in unphased_nodes} |
            {n1 for n1, n2 in edges if n2 in phased_nodes and n1 in unphased_nodes}
        )

        # Build a dictionary to see which unphased positions connect to the phased set
        position_connections = defaultdict(lambda: {"count": 0, "edges": []})
        for n1, n2 in edges:
            # If this edge connects a phased node with a neighbor node...
            if (n1 in phased_nodes and n2 in neighbor_nodes) or (n2 in phased_nodes and n1 in neighbor_nodes):
                # ...then each position in these two nodes is relevant
                for pos in map(int, n1.split('-') + n2.split('-')):
                    if pos in unphased_positions:
                        position_connections[pos]["count"] += 1
                        position_connections[pos]["edges"].append((n1, n2))

        # If there are no connections, pick the next unphased node from sorted list as fallback
        if not position_connections:
            for next_node in sorted_nodes:
                if next_node in unphased_nodes:
                    # Mark its positions as phased
                    next_positions = {int(pos) for pos in next_node.split('-')}
                    phased_positions.update(next_positions)
                    unphased_positions -= next_positions

                    # Mark node as phased
                    phased_nodes.add(next_node)
                    unphased_nodes.discard(next_node)
                    break
            continue

        # ---------------------------------------------------------------------
        # NEW LOGIC: gather *all* unphased positions that appear in position_connections
        # ---------------------------------------------------------------------
        selected_positions = set(position_connections.keys())

        # You could also do something like picking a subset, or picking only
        # positions from a specific node, etc. But here we gather them all.
        # relevant_edges is the union of edges from all these positions
        all_relevant_edges = set()
        for pos in selected_positions:
            for e in position_connections[pos]["edges"]:
                all_relevant_edges.add(e)
        relevant_edges = list(all_relevant_edges)
        print(f"Selected positions: {selected_positions}", f"Relevant edges: {relevant_edges}")
        # ---------------------------------------------------------------------
        # TODO: (Algorithm 2 will actually phase these selected_positions.)
        #
        # e.g., you might do:
        # phased_haplotypes_for_selected = do_actual_phasing(selected_positions, relevant_edges, ...)
        #
        # For now we just skip it and keep a placeholder:
        # ---------------------------------------------------------------------
        # TODO: call the second algorithm with (selected_positions, relevant_edges, phased_positions, ...)
        #       to do the actual phasing
        # ---------------------------------------------------------------------

        # Mark these positions as phased
        phased_positions.update(selected_positions)
        unphased_positions -= selected_positions

        # Mark neighbor_nodes as phased (so we won't pick them again)
        for edge in relevant_edges:
            node1, node2 = edge
            if node1 in neighbor_nodes:
                phased_nodes.add(node1)
                unphased_nodes.discard(node1)
            if node2 in neighbor_nodes:
                phased_nodes.add(node2)
                unphased_nodes.discard(node2)

    # Return updated sets so we know what ended up phased/unphased
    return phased_nodes, unphased_nodes, phased_positions, unphased_positions


def predict_haplotypes_multiple_variants2(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities"):
    """
    Algorithm 1: Iterative selection of positions to be phased.

    This function finds unphased positions connected to the current phased set and
    moves them into the phased sets. It does NOT include the actual phasing logic (Algorithm 2).

    Parameters:
        nodes (list): All nodes in the graph (e.g., ['1-2','2-3',...]).
        edges (list): List of edges (tuples) connecting nodes.
        samples (dict): Original structure with sample phasings. (node->string_phasing)
        ploidy (int): Ploidy level.
        genotype_path (str): CSV file path for genotype info.
        fragment_model: Contains fragment information (unused here, but part of function signature).
        transitions_dict_extra (dict): Extra data about transitions (unused here, but part of function signature).
        config: Configuration object (e.g., error rate) (unused in selection, but part of signature).
        priority (str): Priority mode, unused in selection but included for consistency.

    Returns:
        phased_nodes (set)       : The set of nodes that ended up phased.
        unphased_nodes (set)     : The set of nodes still unphased (should be empty if done).
        phased_positions (set)   : The set of positions that ended up phased.
        unphased_positions (set) : The set of positions still unphased.
    """
    # -------------------------------------------------------------------------
    # Step 1: Initialization (same as old function to keep consistency)
    # -------------------------------------------------------------------------
    phasing_samples = {nn: samples[t][nn] for t in samples for nn in samples[t]}
    sorted_nodes = sort_nodes(nodes)

    genotype_df = pd.read_csv(genotype_path).T
    predicted_haplotypes = pd.DataFrame(
        index=[f'haplotype_{p+1}' for p in range(ploidy)],
        columns=genotype_df.columns
    )

    # Sets to track phased/unphased nodes & positions
    phased_nodes = set()
    unphased_nodes = set(nodes)
    phased_positions = set()
    unphased_positions = {int(pos) for node in nodes for pos in node.split('-')}

    # Pick the first node from sorted_nodes as a starting point
    if sorted_nodes:
        first_node = sorted_nodes[0]
        phased_nodes.add(first_node)
        unphased_nodes.discard(first_node)

        first_positions = {int(pos) for pos in first_node.split('-')}
        phased_positions.update(first_positions)
        unphased_positions -= first_positions

    # -------------------------------------------------------------------------
    # MAIN WHILE LOOP: keep going until all positions are phased
    # -------------------------------------------------------------------------
    while unphased_positions:

        # Identify neighbor_nodes: those unphased but connected to any phased node
        neighbor_nodes = (
            {n2 for n1, n2 in edges if n1 in phased_nodes and n2 in unphased_nodes} |
            {n1 for n1, n2 in edges if n2 in phased_nodes and n1 in unphased_nodes}
        )

        # Build dictionary: which unphased positions connect to the phased set via edges
        position_connections = defaultdict(lambda: {"count": 0, "edges": []})
        for n1, n2 in edges:
            # If edge connects a phased node with one in neighbor_nodes, track them
            if (n1 in phased_nodes and n2 in neighbor_nodes) or (n2 in phased_nodes and n1 in neighbor_nodes):
                # Each position in these two nodes is "connected"
                for pos in map(int, n1.split('-') + n2.split('-')):
                    if pos in unphased_positions:
                        position_connections[pos]["count"] += 1
                        position_connections[pos]["edges"].append((n1, n2))

        # If there are no connections at all, fallback: pick next unphased node from sorted_nodes
        if not position_connections:
            for next_node in sorted_nodes:
                if next_node in unphased_nodes:
                    next_positions = {int(pos) for pos in next_node.split('-')}
                    phased_positions.update(next_positions)
                    unphased_positions -= next_positions

                    phased_nodes.add(next_node)
                    unphased_nodes.discard(next_node)
                    break
            continue

        # Gather *all* unphased positions that appear in position_connections
        selected_positions = set(position_connections.keys())

        # Build a set of "relevant_edges": union of edges from these positions ...
        all_relevant_edges = set()
        for pos in selected_positions:
            for e in position_connections[pos]["edges"]:
                all_relevant_edges.add(e)
        

        # ... plus also edges among neighbor_nodes themselves (if multiple neighbor_nodes exist)
        for e in edges:
            n1, n2 = e
            if n1 in neighbor_nodes and n2 in neighbor_nodes:
                # Even if this edge doesn't directly connect to a phased node, we add it
                # because it belongs to the subgraph of neighbor_nodes.
                all_relevant_edges.add(e)

        relevant_edges = list(all_relevant_edges)
        print(f"Selected positions: {selected_positions}", f"Relevant edges: {relevant_edges}")
        # ---------------------------------------------------------------------
        # TODO (Algorithm 2): Do the actual phasing of these selected positions.
        # e.g. call your "phase_positions(selected_positions, relevant_edges, ...)"
        # This code is left out for clarity.
        # ---------------------------------------------------------------------
        # compute_candidate_phasings(selected_positions, relevant_edges, phased_positions, transitions_dict_extra, phasing_samples, ploidy, predicted_haplotypes, config, fragment_model)

        # Mark these positions as phased
        phased_positions.update(selected_positions)
        unphased_positions -= selected_positions

        # Mark neighbor_nodes as phased, so we won't pick them again in future loops
        for edge in relevant_edges:
            node1, node2 = edge
            if node1 in neighbor_nodes:
                phased_nodes.add(node1)
                unphased_nodes.discard(node1)
            if node2 in neighbor_nodes:
                phased_nodes.add(node2)
                unphased_nodes.discard(node2)

    # Return updated sets at the end
    return phased_nodes, unphased_nodes, phased_positions, unphased_positions



def build_children_dict(edges):
    """
    Builds a dictionary where keys are nodes and values are lists of their child nodes.

    Parameters:
    - edges: list of tuples representing edges between node labels (e.g., [('1-2', '2-3'), ...]).

    Returns:
    - A dictionary with nodes as keys and lists of child nodes as values.
    """
    children_dict = defaultdict(list)

    for parent, child in edges:
        children_dict[parent].append(child)

    return dict(children_dict)


def topological_sort_and_get_parents(nodes, edges):
    """
    Perform topological sorting on a DAG and identify parents for each node.
    
    Parameters:
    - nodes: List of nodes in the graph.
    - edges: List of directed edges (source, target).
    
    Returns:
    - sorted_nodes: List of nodes in topological order.
    - parent_dict: Dictionary where the key is a node and the value is a list of its parents.
    """
    # Build adjacency list and in-degree count
    adjacency_list = defaultdict(list)
    in_degree = {node: 0 for node in nodes}
    
    for s, t in edges:
        adjacency_list[s].append(t)
        in_degree[t] += 1
    
    # Initialize queue with zero in-degree nodes
    zero_in_degree_queue = deque([node for node in nodes if in_degree[node] == 0])
    sorted_nodes = []  # Topological order
    parent_dict = {node: [] for node in nodes}  # Dictionary to store parents of each node
    
    while zero_in_degree_queue:
        node = zero_in_degree_queue.popleft()
        sorted_nodes.append(node)
        
        for neighbor in adjacency_list[node]:
            parent_dict[neighbor].append(node)  # Add current node as a parent of its neighbor
            in_degree[neighbor] -= 1
            if in_degree[neighbor] == 0:
                zero_in_degree_queue.append(neighbor)
    
    if len(sorted_nodes) != len(nodes):
        raise ValueError("The graph is not a DAG (contains cycles).")
    
    return sorted_nodes, parent_dict


def emissions_v2(ploidy, quotient_g, quotient_g_v_label_reversed, error_rate):
    """
    The only difference in version 2 is that we are using quotient graph directly:
    quotient_g.graph ==> quotient_g
    """
    emission_dict = {}
    # Calculate emissions for each state and populate the emission probability matrix
    for state in quotient_g_v_label_reversed.keys():
        emission_dict[state] = {}

        node_elements = state.split('-')  # For example, '1-2' -> ['1', '2']
        node_length = len(node_elements)  # Determine the number of elements in the node

        # Generate all possible combinations of 0 and 1 based on the node length
        possible_emissions = generate_binary_combinations(node_length)
        v = quotient_g_v_label_reversed[state]
        phasings = quotient_g.vertex_properties["v_weights"][v]['weight'].keys()
        for phasing in phasings:
            emission_dict[state][phasing] = {}
            phasing_np = str_2_phas_1(phasing, ploidy)  # Compute phasing for the state key
            for emission in possible_emissions:
                likelihood = compute_likelihood(np.array(emission), phasing_np, error_rate)
                emission_dict[state][phasing][''.join([str(e) for e in emission])] = likelihood

    return emission_dict


def transition_matrices_v2(quotient_g, edges_map_quotient, ploidy, config, fragment_model):
    """
    The only difference in version 2 is that we are using quotient graph directly:
    quotient_g.graph ==> quotient_g
    """
    transitions_dict = {}
    transitions_dict_extra = {}
    for edge in edges_map_quotient.keys():
        transitions_dict_extra[edge] = {}
        source = edges_map_quotient[edge][0]
        target = edges_map_quotient[edge][1]
        source_weights = quotient_g.vertex_properties["v_weights"][source]['weight']
        target_weights = quotient_g.vertex_properties["v_weights"][target]['weight']
        source_label = quotient_g.vertex_properties["v_label"][source]
        target_label = quotient_g.vertex_properties["v_label"][target]
        common_ff, common_sf = find_common_element_and_index(source_label, target_label)
        source_phasings = list(source_weights.keys())
        target_phasings = list(target_weights.keys())
        # transitions_dict = {'source': source_phasings, 'target': target_phasings}
        transitions_mtx = np.zeros((len(source_phasings), len(target_phasings)))

        for i, ffstr in enumerate(source_phasings):
            for j, sfstr in enumerate(target_phasings):
                transitions_dict_extra[edge][str(i) + '-' + str(j)] = {}
                transitions_dict_extra[edge][str(i) + '-' + str(j)]['source_phasing'] = ffstr
                transitions_dict_extra[edge][str(i) + '-' + str(j)]['target_phasing'] = sfstr
                transitions_dict_extra[edge][str(i) + '-' + str(j)]['matched_phasings'] = {}
                matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, ploidy), str_2_phas_1(sfstr, ploidy), common_ff, common_sf, source_label, target_label)
                sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)
                
                matched_phasings_str = list(set([phas_2_str(pm) for pm in sorted_phasings]))
                # print(i, ffstr, j, sfstr)
                # print('matched phasings:', matched_phasings_str, len(matched_phasings_str))
                # if len(matched_phasings_str) > 1:
                #     print('More than one matching phasing')
                #     # stop
                poss = sorted(list(set([int(ss) for ss in source_label.split('-')] + [int(tt) for tt in target_label.split('-')])))
                match_reads = get_matching_reads_for_positions([int(i) for i in poss], fragment_model.fragment_list)
                wei = 0
                for phas in matched_phasings_str:
                    this_phas_weight = 0
                    for indc, this_po, obs in match_reads:
                        this_phas_read_weight = compute_likelihood_generalized_plus(np.array(obs), str_2_phas_1(phas, ploidy), indc, list(range(len(indc))), 
                                                                   config.error_rate)
                        wei += this_phas_read_weight
                        this_phas_weight += this_phas_read_weight
                    transitions_dict_extra[edge][str(i) + '-' + str(j)]['matched_phasings'][phas] = this_phas_weight
                transitions_mtx[i, j] = wei

        transitions_mtx = transitions_mtx / transitions_mtx.sum(axis=1, keepdims=True)
        transitions_dict[edge] = transitions_mtx
    return transitions_dict, transitions_dict_extra

