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


def predict_haplotypes(samples, transitions_dict, transitions_dict_extra, nodes, genotype_path, ploidy):
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

    edges_names = list(transitions_dict.keys())
    sorted_edges = sort_edges(edges_names)

    for edge in sorted_edges:
        source, target = edge.split('--')
        # common_ff, common_sf = find_common_element_and_index(source, target)
        source_sample = samples_brief[source]
        target_sample = samples_brief[target]
        edge_phasings = transitions_dict_extra[edge]
        for key, value in edge_phasings.items():
                # print(key, 'source', value['source_phasing'], source_sample, 'target', value['target_phasing'], target_sample)
                if value['source_phasing'] == source_sample and value['target_phasing'] == target_sample:
                    matched_phasings = value['matched_phasings']
                    break
        if matched_phasings:
            total = sum(matched_phasings.values())
            probabilities = {key: value / total for key, value in matched_phasings.items()}
            keys = list(probabilities.keys())
            probs = list(probabilities.values())
            # sampled_key = random.choices(keys, weights=probs, k=1)[0]

            # Sample 100 times and pick the most frequent sample
            samples_100 = random.choices(keys, weights=probs, k=100)
            sampled_key = Counter(samples_100).most_common(1)[0][0]

            sampled_key_np = str_2_phas_1(sampled_key, ploidy)
            poss = sorted(list(set([int(ss) for ss in source.split('-')] + [int(tt) for tt in target.split('-')])))
            poss = [p - 1 for p in poss]
            # print('Edge:', edge, 'Matched phasings:', matched_phasings, 'Positions:', poss)
            if predicted_haplotypes.loc[:, poss].isna().any().any():
                # Extract the current values for the already filled positions
                existing_cols = [pos for pos in poss if not predicted_haplotypes.loc[:, pos].isna().all()]
                existing_values = predicted_haplotypes.loc[:, existing_cols].values
                
                # Extract corresponding positions in sampled_key_np
                existing_indices = [poss.index(pos) for pos in existing_cols]
                
                # Generate permutations and find the matching one
                sampled_permuted = None
                for permuted_key in itertools.permutations(sampled_key_np):
                    permuted_key_np = np.array(permuted_key)
                    if np.array_equal(permuted_key_np[:, existing_indices], existing_values):
                        sampled_permuted = permuted_key_np
                        break

                if sampled_permuted is None:
                    raise ValueError("No valid permutation found for the given positions.")
                
                # Update the dataframe with the permuted sampled key
                predicted_haplotypes.loc[:, poss] = sampled_permuted

                # predicted_haplotypes.loc[:, poss] = sampled_key_np
                print('Edge:', edge, 'Poss', poss, '\nSampled key:\n', sampled_key_np, '\nSampled permuted:\n', sampled_permuted, 
                      '\nPredicted haplotypes:\n', predicted_haplotypes, '\nkeys:', keys, 'probs:', probs)
                # print('nan detected.')
            # else:
            #     print('These positions were already phased.', poss)
        else:
            print('No matched phasings found for the edge:', edge)

    return predicted_haplotypes

# # Run an example
# if __name__ == '__main__':

    frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/test.frag'
    # frag_path = '/labs/Aguiar/pHapCompass/test/test2.frag'
    ploidy= 3
    genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/haplotypes.csv'
    # genotype_path = '/labs/Aguiar/pHapCompass/test/haplotypes.csv'

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

    samples = sample_states_no_resample_optimized(slices, edges, forward_messages, backward_messages, transitions_dict)
    # samples_brief = {}
    # for t in samples.keys():
    #     for nn in samples[t].keys():
    #         if nn not in samples_brief.keys():
    #             samples_brief[nn] = samples[t][nn]

    predicted_haplotypes = predict_haplotypes(samples, transitions_dict, transitions_dict_extra, nodes, genotype_path, ploidy)

    print('Predicted Haplotypes:\n', predicted_haplotypes)
    print('\nTrue Haplotypes:\n', pd.read_csv(genotype_path).T) 

    predicted_haplotypes_np = predicted_haplotypes.to_numpy()
    true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()

    vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_np, true_haplotypes)
    accuracy, best_permutation = calculate_accuracy(predicted_haplotypes_np, true_haplotypes)
    mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_np, true_haplotypes)
    mec_ = mec(predicted_haplotypes_np, fragment_model.fragment_list)


