import os
import argparse
import sys
import networkx as nx
# import torch
from data.input_handler import InputHandler
from data.configuration import Configuration
from algorithm.haplotype_assembly import HaplotypeAssembly
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from models.factor_graph import Factorgraph
from algorithm.chordal_contraction import chordal_contraction_cycle_base, chordal_contraction
from utils.utils import *
from algorithm.inference import *
from algorithm.chordal_contraction import *
import graph_tool.all as gt
from collections import defaultdict, deque
from evaluation.evaluation import compute_vector_error_rate, calculate_accuracy, calculate_mismatch_error, mec
from collections import Counter
from collections import Counter, defaultdict

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


def transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model):
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


def compute_forward_backward_with_interfaces(slices, interfaces, observations, emission_dict, transitions_dict):
    """
    Compute forward and backward messages for FFBS in a DBN with arbitrary state spaces using interfaces.
    
    Parameters:
    - slices: Dictionary of slices {t: list of nodes in slice t}.
    - interfaces: Dictionary of interfaces {"outgoing": {t: list of nodes}, "incoming": {t: list of nodes}}.
    - observations: Dictionary of observations {t: list of observations relevant to slice t}.
    - emission_dict: Hierarchical dictionary of emission probabilities for each node and its state space.
    - transitions_dict: Hierarchical dictionary of transition probabilities for each edge.
    
    Returns:
    - forward_messages: Forward messages {t: {node: array of forward probabilities for each state-space value}}.
    - backward_messages: Backward messages {t: {node: array of backward probabilities for each state-space value}}.
    """
    adjacency_list = defaultdict(list)
    reverse_adjacency_list = defaultdict(list)
    for s, t in edges:
        adjacency_list[s].append(t)
        reverse_adjacency_list[t].append(s)

    forward_messages = {}
    backward_messages = {}
    
    # Initialize forward messages
    for t in slices:
        forward_messages[t] = {}
        for node in slices[t]:  # Use outgoing interface
            forward_messages[t][node] = {}
            state_space = list(emission_dict[node].keys())
            for phase in state_space:
                forward_messages[t][node][phase] = 0
                # forward_messages[t][node] = {phase: 0}
    
    # Initialize backward messages
    for t in slices:
        backward_messages[t] = {}
        for node in slices[t]:  # Use outgoing interface
            backward_messages[t][node] = {}
            state_space = list(emission_dict[node].keys())
            for phase in state_space:
                backward_messages[t][node][phase] = 0

    # Compute forward messages
    # Base case
    t = 1  
    for node in slices[t]:
        observations = extract_observations_for_state(node, assignment_dict['states'][node], frag_path)  # Compute for outgoing interface
        state_space = list(emission_dict[node].keys())
        for phase in state_space:
            forward_messages[t][node][phase] += -1 * (np.sum([np.log(emission_dict[node][phase][r]) for r in observations]) + np.log(1/len(state_space)))

    # Recursion
    for t in range(2, len(slices) + 1):
        for node in slices[t]:
            observations = extract_observations_for_state(node, assignment_dict['states'][node], frag_path)  # Compute for outgoing interface
            state_space = list(emission_dict[node].keys())
            parents = reverse_adjacency_list[node] # interfaces["outgoing"][t-1]
            parents_in_prev_t = list(set(slices[t-1]).intersection(set(parents)))
            
            for i, phase in enumerate(state_space):
                log_likelihood = np.sum([np.log(emission_dict[node][phase][r]) for r in observations])
                transition_sum = 0
                for prev_node in parents_in_prev_t:
                    prev_state_space = list(emission_dict[prev_node].keys())
                    for j, prev_phase in enumerate(prev_state_space):

                        transition_prob = transitions_dict[f"{prev_node}--{node}"][j, i]
                        prev_meassage = forward_messages[t - 1][prev_node][prev_phase]
                        # transition_prob = transitions_dict.get(f"{prev_node}--{node}", np.zeros((len(prev_state_space), len(state_space))))[j, i]
                        # print(node, phase, prev_node, prev_phase, 'addition:', prev_meassage * transition_prob)
                        transition_sum += prev_meassage * transition_prob
                forward_messages[t][node][phase] += -1 * log_likelihood * np.log(transition_sum)
    
    # Compute for outgoing interface
    # Base case
    t = len(slices)
    for node in slices[t]:  
        state_space = list(emission_dict[node].keys())
        for phase in state_space:  
            backward_messages[t][node][phase] += 1

    # Recursion
    for t in range(len(slices) - 1, 0, -1):
        # print(t)
        for node in slices[t]:  # Compute for outgoing interface
            state_space = list(emission_dict[node].keys())
            children = adjacency_list[node]  # interfaces["incoming"][t+1]
            children_in_next_t = list(set(slices[t+1]).intersection(set(children)))
            # if t == len(slices):  # Base case
            #     backward_messages[t][node] = np.ones(len(state_space))
            # else:  # Recursion
            # backward_sum = 0
            for i, phase in enumerate(state_space):
                backward_sum = 0
                for next_node in children_in_next_t:
                    backward_sum_in = 0
                    observations = extract_observations_for_state(next_node, assignment_dict['states'][next_node], frag_path)
                    next_state_space = list(emission_dict[next_node].keys())
                    for j, next_phase in enumerate(next_state_space):
                        transition_prob = transitions_dict[f"{node}--{next_node}"][i, j]
                        bm = backward_messages[t + 1][next_node][next_phase]
                        for r in observations:
                            
                            # log_likelihood = np.sum([np.log(emission_dict[next_node][next_phase][r]) for r in observations])
                            likelihood = emission_dict[next_node][next_phase][r]
                            backward_sum_in += bm * likelihood * transition_prob
                            print(node, phase, next_node, next_phase, r, 'addition:', bm * likelihood * transition_prob)

                            # backward_sum_in += np.log(bm * likelihood * transition_prob)
                        backward_sum += np.log(backward_sum_in)
                backward_messages[t][node][phase] += backward_sum
    
    return forward_messages, backward_messages


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


def predict_haplotypes(nodes, edges, samples, ploidy, transitions_dict, genotype_path, fragment_model, transitions_dict_extra, config):
    samples_brief = {}
    for t in samples.keys():
        for nn in samples[t].keys():
            if nn not in samples_brief.keys():
                samples_brief[nn] = samples[t][nn]
                # samples_brief.update({nn: samples[t][nn]})

    # sorted_nodes = sort_nodes(nodes)
    genotype_df = pd.read_csv(genotype_path).T
    # positions = genotype_df.columns + 1 
    predicted_haplotypes = pd.DataFrame(index=['haplotype_'+str(p+1) for p in range(ploidy)], columns=genotype_df.columns)

    # edges_names = list(transitions_dict.keys())

    # Initial sets
    phased_nodes = set()   # Set to store phased nodes
    phased_positions = set()  # Set to store phased positions
    unphased_nodes = set(nodes)  # Initially all nodes are unphased
    unphased_positions = set(int(pos) for node in nodes for pos in node.split('-'))  # Positions

    # Step 1: Start with the first node
    first_node = nodes[0]

    # fill the poss - 1 of the df
    poss = sorted(list(set([int(ss) for ss in first_node.split('-')]))) # only this node
    poss = [p - 1 for p in poss]

    first_node_sample = samples_brief[first_node]
    first_sample_np = str_2_phas_1(first_node_sample, ploidy)

    predicted_haplotypes.loc[:, poss] = first_sample_np
    

    phased_nodes.add(first_node)
    unphased_nodes.remove(first_node)
    pos1, pos2 = map(int, first_node.split('-'))
    phased_positions.update([pos1, pos2])
    unphased_positions -= {pos1, pos2}

    while unphased_positions:
        # print('unphased_positions:', unphased_positions)
        # Step 2: Find neighbor nodes connected to phased nodes
        neighbor_nodes = set()
        for node1, node2 in edges:
            if node1 in phased_nodes and node2 in unphased_nodes:
                neighbor_nodes.add(node2)
            elif node2 in phased_nodes and node1 in unphased_nodes:
                neighbor_nodes.add(node1)

        # Step 3: Count connections and store associated edges for unphased positions
        position_connections = defaultdict(lambda: {"count": 0, "edges": []})
        for node1, node2 in edges:
            if (node1 in phased_nodes and node2 in neighbor_nodes) or (node2 in phased_nodes and node1 in neighbor_nodes):
                for pos in map(int, node1.split('-') + node2.split('-')):
                    if pos in unphased_positions:
                        position_connections[pos]["count"] += 1
                        position_connections[pos]["edges"].append((node1, node2))

        # Step 4: Select the unphased position with the highest connections
        selected_position = max(position_connections, key=lambda p: position_connections[p]["count"])
        relevant_edges = position_connections[selected_position]["edges"]
        # print('selected_position:', selected_position, '\nrelevant_edges:', relevant_edges)

        #TODO: fill the selected_position - 1 of the df with correct order of rows
        all_poss = []
        for edge in relevant_edges:
            source, target = edge
            all_poss.extend([int(ss) for ss in source.split('-')] + [int(tt) for tt in target.split('-')])
        all_poss = sorted(list(set(all_poss))) 
        
        fixed_positions = [p for p in all_poss if p in phased_positions]
        fixed_positions_df = [p - 1 for p in fixed_positions]
        # existing_indices = [all_poss_id.index(p) for p in fixed_positions_df]
        fixed_values = predicted_haplotypes.loc[:, fixed_positions_df].values
        # sel_pos_index = all_poss_id.index(selected_position - 1)
        fixed_positions_ids = [all_poss.index(p) for p in fixed_positions]
        unphased_position_id = all_poss.index(selected_position)
        match_reads = get_matching_reads_for_positions([int(i) for i in all_poss], fragment_model.fragment_list)

        matched_phasings_str = {}
        matched_phasings_all_conn = {}
        # temp_np = np.empty((ploidy, len(all_poss)))
        # temp_np = np.full((ploidy, len(all_poss)), np.nan)
        # min_pos = all_poss[0]
        # max_pos = all_poss[-1]
        
        edges_as_sets = [set(edge[0].split('-') + edge[1].split('-')) for edge in relevant_edges]
        all_temp_np = []
        # Find shared positions between edges
        # shared_positions = [p for p in [int(i) for i in list(set.intersection(*edges_as_sets))] if p in phased_positions]

        if len(relevant_edges) == 1:
            edge = relevant_edges[0]
            source, target = edge
            edge_label = source + '--' + target
            matched_phasings_all_conn[edge_label] = []
            this_edge_positions = sorted(set(int(x) for part in edge_label.split('--') for x in part.split('-')))
            this_edge_positions_df = [p - 1 for p in this_edge_positions]

            fixed_positions_edge = [p for p in this_edge_positions if p in all_poss and p in phased_positions] # oonayee ke dar all_poss hastan va phased hastan 
            fixed_positions_edge_df = [p - 1 for p in fixed_positions_edge]

            fixed_ids_this_poss = [this_edge_positions.index(p) for p in fixed_positions_edge] # index of the fixed positions in the edge
            this_edge_ind_in_all_poss = [all_poss.index(p) for p in this_edge_positions if p in phased_positions]
            fixed_values_edge = predicted_haplotypes.loc[:, fixed_positions_edge_df].values

            edge_phasings = transitions_dict_extra[edge_label]
            for key, value in edge_phasings.items():
                if value['source_phasing'] == samples_brief[source] and value['target_phasing'] == samples_brief[target]:
                    matched_phasings = value['matched_phasings']
                    break
            if matched_phasings:
                # matched_phasings_all_conn[edge_label] = matched_phasings
                # sampled_key = select_max_prob_key(probabilities)
                for match in matched_phasings.keys():
                    matched_phasings_np = str_2_phas_1(match, ploidy)
                    # print(fixed_values_edge)
                    correct_permuted = None
                    for permuted_key in itertools.permutations(matched_phasings_np):
                        permuted_key_np = np.array(permuted_key)
                        # print(permuted_key_np[:, fixed_ids_this_poss], '\n', fixed_values_edge)
                        # print('-----------------------')
                        if np.array_equal(fixed_values_edge, permuted_key_np[:, fixed_ids_this_poss]):
                            # print('oooooooooooooo')
                            correct_permuted = permuted_key_np
                            break
                    temp_np = np.full((ploidy, len(all_poss)), np.nan)
                    temp_np[:, this_edge_ind_in_all_poss] = fixed_values_edge
                    # this_edge_temp_ids =  [p - min_pos for p in this_edge_positions]
                    this_edge_temp_ids =  [all_poss.index(p) for p in this_edge_positions]
                    temp_np[:, this_edge_temp_ids] = correct_permuted
                    all_temp_np.append(temp_np)
                    matched_phasings_all_conn[edge_label].append(temp_np)          
            
            merged_nps = []
            probs = []
            for np1 in all_temp_np:
                non_nan_ind = np.where(~np.isnan(np1[0]))[0]
                this_phas_weight = 0
                for indc, this_po, obs in match_reads:
                    rd_poss_ids_in_all_poss = [all_poss.index(p) for p in this_po] # index in temp_np for reads
                    shared_indices = sorted(set(rd_poss_ids_in_all_poss).intersection(set(non_nan_ind)))
                    # [this_po.index(p) for p in shared_indices]
                    # rd_poss_ids_in_this_edge = [this_edge_positions.index(p) for p in this_po]
                    this_phas_read_weight = compute_likelihood_generalized_plus(np.array(obs), np1, indc, shared_indices, config.error_rate)
                    # wei += this_phas_read_weight
                    this_phas_weight += this_phas_read_weight
                merged_nps.append(np1)
                probs.append(this_phas_weight)
            # max_np = merged_nps[probs.index(max(probs))]
            max_prob = max(probs)
            max_indices = [i for i, p in enumerate(probs) if p == max_prob]
            chosen_index = random.choice(max_indices)
            max_np = merged_nps[chosen_index]

            correct_permuted = None
            for permuted_key in itertools.permutations(max_np):
                permuted_key_np = np.array(permuted_key)
                if np.array_equal(fixed_values, permuted_key_np[:, fixed_positions_ids]):
                    correct_permuted = permuted_key_np
                    break
            correct_permuted_int = correct_permuted.astype(int)
            predicted_haplotypes.loc[:, selected_position - 1] = correct_permuted_int[:, unphased_position_id]


        if len(relevant_edges) > 1:

            for edge in relevant_edges:
                # print('------------- edge:', edge)
                source, target = edge
                target_phasing = samples_brief[target]
                edge_label = source + '--' + target
                matched_phasings_all_conn[edge_label] = []
                this_edge_positions = sorted(set(int(x) for part in edge_label.split('--') for x in part.split('-')))
                this_edge_positions_df = [p - 1 for p in this_edge_positions]

                fixed_positions_edge = [p for p in this_edge_positions if p in all_poss and p in phased_positions] # oonayee ke dar all_poss hastan va phased hastan 
                fixed_positions_edge_df = [p - 1 for p in fixed_positions_edge]

                fixed_ids_this_poss = [this_edge_positions.index(p) for p in fixed_positions_edge] # index of the fixed positions in the edge
                this_edge_ind_in_all_poss = [all_poss.index(p) for p in this_edge_positions if p in phased_positions]
                fixed_values_edge = predicted_haplotypes.loc[:, fixed_positions_edge_df].values

                edge_phasings = transitions_dict_extra[edge_label]
                for key, value in edge_phasings.items():
                    if value['source_phasing'] == samples_brief[source] and value['target_phasing'] == samples_brief[target]:
                        matched_phasings = value['matched_phasings']
                        break
                if matched_phasings:
                    # matched_phasings_all_conn[edge_label] = matched_phasings
                    # sampled_key = select_max_prob_key(probabilities)
                    for match in matched_phasings.keys():
                        matched_phasings_np = str_2_phas_1(match, ploidy)
                        # print(fixed_values_edge)
                        correct_permuted = None
                        # target_positions = [int(pos) for pos in target.split('-')]  # Extract positions of the target node
                        # target_positions_indices = [all_poss.index(pos) for pos in target_positions]
                        for permuted_key in itertools.permutations(matched_phasings_np):
                            permuted_key_np = np.array(permuted_key)
                            # print(permuted_key_np[:, fixed_ids_this_poss], '\n', fixed_values_edge)
                            # print('-----------------------')
                            if np.array_equal(fixed_values_edge, permuted_key_np[:, fixed_ids_this_poss]):
                                # print('oooooooooooooo')
                                correct_permuted = permuted_key_np
                                break
                        temp_np = np.full((ploidy, len(all_poss)), np.nan)
                        temp_np[:, this_edge_ind_in_all_poss] = fixed_values_edge
                        # this_edge_temp_ids =  [p - min_pos for p in this_edge_positions]
                        this_edge_temp_ids =  [all_poss.index(p) for p in this_edge_positions]
                        temp_np[:, this_edge_temp_ids] = correct_permuted
                        all_temp_np.append(temp_np)
                        matched_phasings_all_conn[edge_label].append(temp_np)

            # # candidate_phasings = list(itertools.product(*matched_phasings_all_conn.values()))
            # candidate_phasings = list(itertools.combinations(matched_phasings_all_conn.values(), 2))

            # Generate combinations of keys
            key_combinations = itertools.combinations(matched_phasings_all_conn.keys(), 2)

            # Generate combinations of NumPy arrays where each comes from a different key
            candidate_phasings = []
            for key1, key2 in key_combinations:
                arrays1 = matched_phasings_all_conn[key1]
                arrays2 = matched_phasings_all_conn[key2]
                # Cartesian product of arrays from the two keys
                candidate_phasings.extend(itertools.product(arrays1, arrays2))


            merged_nps = []
            probs = []
            # counter = 0
            for np1, np2 in candidate_phasings:
                if np.array_equal(np1,np2, equal_nan=True) and np.isnan(np1).any():
                    # print('equal')
                    continue
                # print(counter)
                # counter += 1
                merge_list = match_np(np1, np2)
                for this_merge in merge_list:
                    # print(this_merge)
                    valid_columns_merge = np.all(~np.isnan(this_merge), axis=0)  # Boolean mask of valid columns
                    non_nan_ind_merge = np.where(valid_columns_merge)[0]  
                    this_phas_weight = 0
                    for indc, this_po, obs in match_reads:
                        
                        rd_poss_ids_in_all_poss = [all_poss.index(p) for p in this_po if p in all_poss] # index in temp_np for reads
                        shared_indices = sorted(set(rd_poss_ids_in_all_poss).intersection(set(non_nan_ind_merge)))
                        # this_rd_indc = [p for p in indc if p in shared_indices]
                        this_rd_indc = [this_po.index(all_poss[p]) for p in shared_indices]
                        # this_rd_indc = [this_po.index(p) for p in shared_indices]
                        # [this_po.index(p) for p in shared_indices]
                        # rd_poss_ids_in_this_edge = [this_edge_positions.index(p) for p in this_po]
                        this_phas_read_weight = compute_likelihood_generalized_plus(np.array(obs), this_merge, this_rd_indc, shared_indices, config.error_rate)
                        # wei += this_phas_read_weight
                        this_phas_weight += this_phas_read_weight
                    merged_nps.append(this_merge)
                    probs.append(this_phas_weight)
            
            # # max_np = merged_nps[probs.index(max(probs))]
            # max_prob = max(probs)
            # max_indices = [i for i, p in enumerate(probs) if p == max_prob]
            # chosen_index = random.choice(max_indices)
            # max_np = merged_nps[chosen_index]


            # Sort both lists together based on probs
            sorted_pairs = sorted(zip(probs, merged_nps), key=lambda x: x[0], reverse=True)
            # Separate the sorted probabilities and numpy arrays if needed
            sorted_probs, sorted_nps = zip(*sorted_pairs)
            # Convert back to lists if desired
            sorted_probs = list(sorted_probs)
            sorted_nps = list(sorted_nps)


            # valid_columns_max = np.all(~np.isnan(max_np), axis=0)  # Boolean mask of valid columns
            # non_nan_ind_max = np.where(valid_columns_max)[0]  

            # shared_idxs = sorted(set(non_nan_ind_max).intersection(set(fixed_positions_ids)))

            correct_permuted = None
            for max_np in sorted_nps:
                # print(max_np)
                valid_columns_max = np.all(~np.isnan(max_np), axis=0)  # Boolean mask of valid columns
                non_nan_ind_max = np.where(valid_columns_max)[0]  

                shared_idxs = sorted(set(non_nan_ind_max).intersection(set(fixed_positions_ids)))
                shared_idx_in_df_fixed = [all_poss[p] - 1 for p in shared_idxs]
                shared_fixed = predicted_haplotypes.loc[:, shared_idx_in_df_fixed].values
                for permuted_key in itertools.permutations(max_np):
                    
                    permuted_key_np = np.array(permuted_key)
                    # print(fixed_values[:, shared_idxs], '\n', permuted_key_np[:, shared_idxs])
                    # print('---------------------------------')
                    if np.array_equal(shared_fixed, permuted_key_np[:, shared_idxs]):
                        correct_permuted = permuted_key_np
                        break
                if correct_permuted is not None:
                    break

            correct_permuted_int = correct_permuted.astype(int)
            predicted_haplotypes.loc[:, selected_position - 1] = correct_permuted_int[:, unphased_position_id]

        # Step 5: Add nodes containing the selected position to phased nodes
        for edge in relevant_edges:
            node1, node2 = edge
            if node1 in neighbor_nodes:
                phased_nodes.add(node1)
                unphased_nodes.discard(node1)
            if node2 in neighbor_nodes:
                phased_nodes.add(node2)
                unphased_nodes.discard(node2)

        # Update phased positions and remove from unphased
        phased_positions.add(selected_position)
        unphased_positions.remove(selected_position)

        # # Output the selected position and relevant edges for this iteration
        # print(f"Selected position: {selected_position}")
        # print(f"Relevant edges: {relevant_edges}")

    print("Phasing complete.")
    return predicted_haplotypes



def predict_haplotypes1(nodes, edges, samples, ploidy, transitions_dict, genotype_path, fragment_model, transitions_dict_extra, config):
    samples_brief = {}
    for t in samples.keys():
        for nn in samples[t].keys():
            if nn not in samples_brief.keys():
                samples_brief[nn] = samples[t][nn]
                # samples_brief.update({nn: samples[t][nn]})

    # sorted_nodes = sort_nodes(nodes)
    genotype_df = pd.read_csv(genotype_path).T
    # positions = genotype_df.columns + 1 
    predicted_haplotypes = pd.DataFrame(index=['haplotype_'+str(p+1) for p in range(ploidy)], columns=genotype_df.columns)

    # edges_names = list(transitions_dict.keys())

    # Initial sets
    phased_nodes = set()   # Set to store phased nodes
    phased_positions = set()  # Set to store phased positions
    unphased_nodes = set(nodes)  # Initially all nodes are unphased
    unphased_positions = set(int(pos) for node in nodes for pos in node.split('-'))  # Positions

    # Step 1: Start with the first node
    first_node = nodes[0]

    # fill the poss - 1 of the df
    poss = sorted(list(set([int(ss) for ss in first_node.split('-')]))) # only this node
    poss = [p - 1 for p in poss]

    first_node_sample = samples_brief[first_node]
    first_sample_np = str_2_phas_1(first_node_sample, ploidy)

    predicted_haplotypes.loc[:, poss] = first_sample_np
    

    phased_nodes.add(first_node)
    unphased_nodes.remove(first_node)
    pos1, pos2 = map(int, first_node.split('-'))
    phased_positions.update([pos1, pos2])
    unphased_positions -= {pos1, pos2}

    while unphased_positions:
        # print('unphased_positions:', unphased_positions)
        # Step 2: Find neighbor nodes connected to phased nodes
        neighbor_nodes = set()
        for node1, node2 in edges:
            if node1 in phased_nodes and node2 in unphased_nodes:
                neighbor_nodes.add(node2)
            elif node2 in phased_nodes and node1 in unphased_nodes:
                neighbor_nodes.add(node1)

        # Step 3: Count connections and store associated edges for unphased positions
        position_connections = defaultdict(lambda: {"count": 0, "edges": []})
        for node1, node2 in edges:
            if (node1 in phased_nodes and node2 in neighbor_nodes) or (node2 in phased_nodes and node1 in neighbor_nodes):
                for pos in map(int, node1.split('-') + node2.split('-')):
                    if pos in unphased_positions:
                        position_connections[pos]["count"] += 1
                        position_connections[pos]["edges"].append((node1, node2))

        # Step 4: Select the unphased position with the highest connections
        selected_position = max(position_connections, key=lambda p: position_connections[p]["count"])
        relevant_edges = position_connections[selected_position]["edges"]
        # print('selected_position:', selected_position, '\nrelevant_edges:', relevant_edges)

        #TODO: fill the selected_position - 1 of the df with correct order of rows
        all_poss = []
        for edge in relevant_edges:
            source, target = edge
            all_poss.extend([int(ss) for ss in source.split('-')] + [int(tt) for tt in target.split('-')])
        all_poss = sorted(list(set(all_poss))) 
        
        fixed_positions = [p for p in all_poss if p in phased_positions]
        fixed_positions_df = [p - 1 for p in fixed_positions]
        # existing_indices = [all_poss_id.index(p) for p in fixed_positions_df]
        fixed_values = predicted_haplotypes.loc[:, fixed_positions_df].values
        # sel_pos_index = all_poss_id.index(selected_position - 1)
        fixed_positions_ids = [all_poss.index(p) for p in fixed_positions]
        unphased_position_id = all_poss.index(selected_position)
        match_reads = get_matching_reads_for_positions([int(i) for i in all_poss], fragment_model.fragment_list)

        matched_phasings_str = {}
        matched_phasings_all_conn = {}
        # temp_np = np.empty((ploidy, len(all_poss)))
        # temp_np = np.full((ploidy, len(all_poss)), np.nan)
        # min_pos = all_poss[0]
        # max_pos = all_poss[-1]
        
        edges_as_sets = [set(edge[0].split('-') + edge[1].split('-')) for edge in relevant_edges]
        all_temp_np = []
        # Find shared positions between edges
        # shared_positions = [p for p in [int(i) for i in list(set.intersection(*edges_as_sets))] if p in phased_positions]

        if len(relevant_edges) == 1:
            edge = relevant_edges[0]
            source, target = edge
            edge_label = source + '--' + target
            matched_phasings_all_conn[edge_label] = []
            this_edge_positions = sorted(set(int(x) for part in edge_label.split('--') for x in part.split('-')))
            this_edge_positions_df = [p - 1 for p in this_edge_positions]

            fixed_positions_edge = [p for p in this_edge_positions if p in all_poss and p in phased_positions] # oonayee ke dar all_poss hastan va phased hastan 
            fixed_positions_edge_df = [p - 1 for p in fixed_positions_edge]

            fixed_ids_this_poss = [this_edge_positions.index(p) for p in fixed_positions_edge] # index of the fixed positions in the edge
            this_edge_ind_in_all_poss = [all_poss.index(p) for p in this_edge_positions if p in phased_positions]
            fixed_values_edge = predicted_haplotypes.loc[:, fixed_positions_edge_df].values

            edge_phasings = transitions_dict_extra[edge_label]
            for key, value in edge_phasings.items():
                if value['source_phasing'] == samples_brief[source] and value['target_phasing'] == samples_brief[target]:
                    matched_phasings = value['matched_phasings']
                    break
            if matched_phasings:
                # matched_phasings_all_conn[edge_label] = matched_phasings
                # sampled_key = select_max_prob_key(probabilities)
                for match in matched_phasings.keys():
                    matched_phasings_np = str_2_phas_1(match, ploidy)
                    # print(fixed_values_edge)
                    correct_permuted = None
                    for permuted_key in itertools.permutations(matched_phasings_np):
                        permuted_key_np = np.array(permuted_key)
                        # print(permuted_key_np[:, fixed_ids_this_poss], '\n', fixed_values_edge)
                        # print('-----------------------')
                        if np.array_equal(fixed_values_edge, permuted_key_np[:, fixed_ids_this_poss]):
                            # print('oooooooooooooo')
                            correct_permuted = permuted_key_np
                            break
                    temp_np = np.full((ploidy, len(all_poss)), np.nan)
                    temp_np[:, this_edge_ind_in_all_poss] = fixed_values_edge
                    # this_edge_temp_ids =  [p - min_pos for p in this_edge_positions]
                    this_edge_temp_ids =  [all_poss.index(p) for p in this_edge_positions]
                    temp_np[:, this_edge_temp_ids] = correct_permuted
                    all_temp_np.append(temp_np)
                    matched_phasings_all_conn[edge_label].append(temp_np)          
            
            merged_nps = []
            probs = []
            for np1 in all_temp_np:
                non_nan_ind = np.where(~np.isnan(np1[0]))[0]
                this_phas_weight = 0
                for indc, this_po, obs in match_reads:
                    rd_poss_ids_in_all_poss = [all_poss.index(p) for p in this_po] # index in temp_np for reads
                    shared_indices = sorted(set(rd_poss_ids_in_all_poss).intersection(set(non_nan_ind)))
                    # [this_po.index(p) for p in shared_indices]
                    # rd_poss_ids_in_this_edge = [this_edge_positions.index(p) for p in this_po]
                    this_phas_read_weight = compute_likelihood_generalized_plus(np.array(obs), np1, indc, shared_indices, config.error_rate)
                    # wei += this_phas_read_weight
                    this_phas_weight += this_phas_read_weight
                merged_nps.append(np1)
                probs.append(this_phas_weight)
            # max_np = merged_nps[probs.index(max(probs))]
            max_prob = max(probs)
            max_indices = [i for i, p in enumerate(probs) if p == max_prob]
            chosen_index = random.choice(max_indices)
            max_np = merged_nps[chosen_index]

            correct_permuted = None
            for permuted_key in itertools.permutations(max_np):
                permuted_key_np = np.array(permuted_key)
                if np.array_equal(fixed_values, permuted_key_np[:, fixed_positions_ids]):
                    correct_permuted = permuted_key_np
                    break
            correct_permuted_int = correct_permuted.astype(int)
            predicted_haplotypes.loc[:, selected_position - 1] = correct_permuted_int[:, unphased_position_id]


        if len(relevant_edges) > 1:

            for edge in relevant_edges:
                # print('------------- edge:', edge)
                source, target = edge
                target_phasing = samples_brief[target]
                edge_label = source + '--' + target
                matched_phasings_all_conn[edge_label] = []
                this_edge_positions = sorted(set(int(x) for part in edge_label.split('--') for x in part.split('-')))
                this_edge_positions_df = [p - 1 for p in this_edge_positions]

                fixed_positions_edge = [p for p in this_edge_positions if p in all_poss and p in phased_positions] # oonayee ke dar all_poss hastan va phased hastan 
                fixed_positions_edge_df = [p - 1 for p in fixed_positions_edge]

                fixed_ids_this_poss = [this_edge_positions.index(p) for p in fixed_positions_edge] # index of the fixed positions in the edge
                this_edge_ind_in_all_poss = [all_poss.index(p) for p in this_edge_positions if p in phased_positions]
                fixed_values_edge = predicted_haplotypes.loc[:, fixed_positions_edge_df].values

                edge_phasings = transitions_dict_extra[edge_label]
                for key, value in edge_phasings.items():
                    if value['source_phasing'] == samples_brief[source] and value['target_phasing'] == samples_brief[target]:
                        matched_phasings = value['matched_phasings']
                        break
                if matched_phasings:
                    # matched_phasings_all_conn[edge_label] = matched_phasings
                    # sampled_key = select_max_prob_key(probabilities)
                    for match in matched_phasings.keys():
                        matched_phasings_np = str_2_phas_1(match, ploidy)
                        # print(fixed_values_edge)
                        correct_permuted = None
                        target_positions = [int(pos) for pos in target.split('-')]  # Extract positions of the target node
                        target_positions_indices = [this_edge_positions.index(pos) for pos in target_positions]
                        target_phasing_permutations = np.array(list(itertools.permutations(str_2_phas_1(target_phasing, ploidy))))

                        # for permuted_key in itertools.permutations(matched_phasings_np):
                        #     permuted_key_np = np.array(permuted_key)
                        #     # print(permuted_key_np[:, fixed_ids_this_poss], '\n', fixed_values_edge)
                        #     # print('-----------------------')
                        #     if np.array_equal(fixed_values_edge, permuted_key_np[:, fixed_ids_this_poss]):
                        #         # print('oooooooooooooo')
                        #         correct_permuted = permuted_key_np
                        #         break
                        for permuted_key in itertools.permutations(matched_phasings_np):
                            permuted_key_np = np.array(permuted_key)

                            # Check fixed positions first
                            if not np.array_equal(fixed_values_edge, permuted_key_np[:, fixed_ids_this_poss]):
                                continue

                            # Check if any permutation of target_phasing matches target positions
                            if not any(
                                np.array_equal(permuted_key_np[:, target_positions_indices], perm)
                                for perm in target_phasing_permutations):
                                continue

                            # If both conditions are true, break with the valid permutation
                            correct_permuted = permuted_key_np
                            break

                        temp_np = np.full((ploidy, len(all_poss)), np.nan)
                        temp_np[:, this_edge_ind_in_all_poss] = fixed_values_edge
                        # this_edge_temp_ids =  [p - min_pos for p in this_edge_positions]
                        this_edge_temp_ids =  [all_poss.index(p) for p in this_edge_positions]
                        temp_np[:, this_edge_temp_ids] = correct_permuted
                        all_temp_np.append(temp_np)
                        matched_phasings_all_conn[edge_label].append(temp_np)

            # # candidate_phasings = list(itertools.product(*matched_phasings_all_conn.values()))
            # candidate_phasings = list(itertools.combinations(matched_phasings_all_conn.values(), 2))

            # Generate combinations of keys
            key_combinations = itertools.combinations(matched_phasings_all_conn.keys(), 2)

            # Generate combinations of NumPy arrays where each comes from a different key
            candidate_phasings = []
            for key1, key2 in key_combinations:
                arrays1 = matched_phasings_all_conn[key1]
                arrays2 = matched_phasings_all_conn[key2]
                # Cartesian product of arrays from the two keys
                candidate_phasings.extend(itertools.product(arrays1, arrays2))


            merged_nps = []
            probs = []
            # counter = 0
            for np1, np2 in candidate_phasings:
                if np.array_equal(np1, np2, equal_nan=True) and np.isnan(np1).any():
                    # print('equal')
                    continue
                # print(counter)
                # counter += 1
                merge_list = match_np(np1, np2)
                for this_merge in merge_list:
                    # print(this_merge)
                    valid_columns_merge = np.all(~np.isnan(this_merge), axis=0)  # Boolean mask of valid columns
                    non_nan_ind_merge = np.where(valid_columns_merge)[0]  
                    this_phas_weight = 0
                    for indc, this_po, obs in match_reads:
                        
                        rd_poss_ids_in_all_poss = [all_poss.index(p) for p in this_po if p in all_poss] # index in temp_np for reads
                        shared_indices = sorted(set(rd_poss_ids_in_all_poss).intersection(set(non_nan_ind_merge)))
                        # this_rd_indc = [p for p in indc if p in shared_indices]
                        this_rd_indc = [this_po.index(all_poss[p]) for p in shared_indices]
                        # this_rd_indc = [this_po.index(p) for p in shared_indices]
                        # [this_po.index(p) for p in shared_indices]
                        # rd_poss_ids_in_this_edge = [this_edge_positions.index(p) for p in this_po]
                        this_phas_read_weight = compute_likelihood_generalized_plus(np.array(obs), this_merge, this_rd_indc, shared_indices, config.error_rate)
                        # wei += this_phas_read_weight
                        this_phas_weight += this_phas_read_weight
                    merged_nps.append(this_merge)
                    probs.append(this_phas_weight)
            
            # # max_np = merged_nps[probs.index(max(probs))]
            # max_prob = max(probs)
            # max_indices = [i for i, p in enumerate(probs) if p == max_prob]
            # chosen_index = random.choice(max_indices)
            # max_np = merged_nps[chosen_index]


            # Sort both lists together based on probs
            sorted_pairs = sorted(zip(probs, merged_nps), key=lambda x: x[0], reverse=True)
            # Separate the sorted probabilities and numpy arrays if needed
            sorted_probs, sorted_nps = zip(*sorted_pairs)
            # Convert back to lists if desired
            sorted_probs = list(sorted_probs)
            sorted_nps = list(sorted_nps)


            # valid_columns_max = np.all(~np.isnan(max_np), axis=0)  # Boolean mask of valid columns
            # non_nan_ind_max = np.where(valid_columns_max)[0]  

            # shared_idxs = sorted(set(non_nan_ind_max).intersection(set(fixed_positions_ids)))

            correct_permuted = None
            for max_np in sorted_nps:
                # print(max_np)
                valid_columns_max = np.all(~np.isnan(max_np), axis=0)  # Boolean mask of valid columns
                non_nan_ind_max = np.where(valid_columns_max)[0]  

                shared_idxs = sorted(set(non_nan_ind_max).intersection(set(fixed_positions_ids)))
                shared_idx_in_df_fixed = [all_poss[p] - 1 for p in shared_idxs]
                shared_fixed = predicted_haplotypes.loc[:, shared_idx_in_df_fixed].values
                for permuted_key in itertools.permutations(max_np):
                    
                    permuted_key_np = np.array(permuted_key)
                    # print(fixed_values[:, shared_idxs], '\n', permuted_key_np[:, shared_idxs])
                    # print('---------------------------------')
                    if np.array_equal(shared_fixed, permuted_key_np[:, shared_idxs]):
                        correct_permuted = permuted_key_np
                        break
                if correct_permuted is not None:
                    break

            correct_permuted_int = correct_permuted.astype(int)
            predicted_haplotypes.loc[:, selected_position - 1] = correct_permuted_int[:, unphased_position_id]

        # Step 5: Add nodes containing the selected position to phased nodes
        for edge in relevant_edges:
            node1, node2 = edge
            if node1 in neighbor_nodes:
                phased_nodes.add(node1)
                unphased_nodes.discard(node1)
            if node2 in neighbor_nodes:
                phased_nodes.add(node2)
                unphased_nodes.discard(node2)

        # Update phased positions and remove from unphased
        phased_positions.add(selected_position)
        unphased_positions.remove(selected_position)

        # # Output the selected position and relevant edges for this iteration
        # print(f"Selected position: {selected_position}")
        # print(f"Relevant edges: {relevant_edges}")

    print("Phasing complete.")
    return predicted_haplotypes



def match_np(np1,np2):
    np1_cols = np.all(~np.isnan(np1), axis=0)
    np2_cols = np.all(~np.isnan(np2), axis=0)
    indices1 = np.where(np1_cols)[0]
    indices2 = np.where(np2_cols)[0]
    np1_nan_cols = np.isnan(np1).all(axis=0)
    np2_nan_cols = np.isnan(np2).all(axis=0)
    # np1_nan_indices = np.where(np1_nan_cols)[0]
    np2_nan_indices = np.where(np2_nan_cols)[0]

    premutations = []
    fixed_positions = sorted(list(set(indices1).intersection(set(indices2))))
    fixed_values1 = np1[:, fixed_positions]
    correct_permuted = None
    for permuted_key in itertools.permutations(np2):
        permuted_key_np = np.array(permuted_key)
        if np.array_equal(fixed_values1, permuted_key_np[:, fixed_positions]):
            correct_permuted = permuted_key_np
            correct_permuted[:, np2_nan_indices] = np1[:, np2_nan_indices]
            premutations.append(correct_permuted)
    return premutations



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

# # Run an example
# if __name__ == '__main__':
def main():
    # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/test.frag'
    # frag_path = '/labs/Aguiar/pHapCompass/test/test2.frag'
    frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_3/cov_10/frag/00.frag'
    ploidy= 3
    # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/haplotypes.csv'
    # genotype_path = '/labs/Aguiar/pHapCompass/test/haplotypes.csv'
    genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_3/haplotypes.csv'

    class Args:
        def __init__(self):
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = frag_path
            self.bam_path = 'example/example.bam'
            self.genotype_path = genotype_path
            self.ploidy = 3
            self.error_rate = 0.001
            self.epsilon = 0.0001
            self.output_path = 'output'
            self.root_dir = 'D:/UCONN/HaplOrbit'
            self.alleles = [0, 1]

    # Create the mock args object
    args = Args()

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)

    fragment_model.construct(input_handler, config)
    print('Fragment Graph constructed.')

    e_labels = fragment_model.graph.edge_properties["e_label"]
    v_labels = fragment_model.graph.vertex_properties["v_label"]
    gt.graph_draw(fragment_model.graph, output_size=(500, 500), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=14,  
    edge_font_size=12)

    fragment_model_v_label_reversed = fragment_model.v_label_reversed

    edges_map_fragment = {}
    for k in fragment_model.e_label_reversed.keys():
        edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]


    # create quotient graph
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct(input_handler, config)

    e_labels_q = quotient_g.graph.edge_properties["e_label"]
    v_labels_q = quotient_g.graph.vertex_properties["v_label"]
    gt.graph_draw(quotient_g.graph, output_size=(500, 500), vertex_text=v_labels_q, edge_text=e_labels_q, vertex_font_size=14,  
    edge_font_size=12)

    quotient_g_v_label_reversed = quotient_g.v_label_reversed

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model)
    emission_dict = emissions(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, interfaces =  assign_slices_and_interfaces(nodes, edges)

    assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, frag_path)

    forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, frag_path)

    backward_messages = compute_backward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, frag_path)

    # samples = sample_states_no_resample_optimized(slices, edges, forward_messages, backward_messages, transitions_dict)
    # samples = sample_states(slices, edges, forward_messages, transitions_dict)
    samples = sample_states_book(slices, edges, forward_messages, transitions_dict)
    # samples_brief = {}
    # for t in samples.keys():
    #     for nn in samples[t].keys():
    #         if nn not in samples_brief.keys():
    #             samples_brief[nn] = samples[t][nn]

    predicted_haplotypes = predict_haplotypes1(nodes, edges, samples, ploidy, transitions_dict, genotype_path, fragment_model, transitions_dict_extra, config)
    print('Predicted Haplotypes:\n', predicted_haplotypes)
    print('\nTrue Haplotypes:\n', pd.read_csv(genotype_path).T) 

    predicted_haplotypes_np = predicted_haplotypes.to_numpy()
    true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()

    predicted_haplotypes_np = predicted_haplotypes_np[:,0: -1]
    true_haplotypes = true_haplotypes[:,0: -1]

    vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_np, true_haplotypes)
    accuracy, best_permutation = calculate_accuracy(predicted_haplotypes_np, true_haplotypes)
    mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_np, true_haplotypes)
    mec_ = mec(predicted_haplotypes_np, fragment_model.fragment_list)
    print('Vector Error Rate:', vector_error_rate, '\nVector Error:', vector_error, '\nAccuracy:', accuracy, '\nMismatch Error:', mismatch_error, '\nMEC:', mec_)
