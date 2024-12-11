import os
import argparse
import sys
import networkx as nx
import torch
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


def transition_matrices(quotient_g, edges_map_quotient):
    transitions_dict = {}
    for edge in edges_map_quotient.keys():
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
                
                matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, ploidy), str_2_phas_1(sfstr, ploidy), common_ff, common_sf, source_label, target_label)
                sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)
                
                matched_phasings_str = list(set([phas_2_str(pm) for pm in sorted_phasings]))
                poss = sorted(list(set([int(ss) for ss in source_label.split('-')] + [int(tt) for tt in target_label.split('-')])))
                match_reads = get_matching_reads_for_positions([int(i) for i in poss], fragment_model.fragment_list)
                wei = 0
                for phas in matched_phasings_str:
                    for indc, this_po, obs in match_reads:
                        wei += compute_likelihood_generalized_plus(np.array(obs), str_2_phas_1(phas, ploidy), indc, list(range(len(indc))), 
                                                                   config.error_rate)
                transitions_mtx[i, j] = wei

        transitions_mtx = transitions_mtx / transitions_mtx.sum(axis=1, keepdims=True)
        transitions_dict[edge] = transitions_mtx
    return transitions_dict


def emissions(ploidy, quotient_g_v_label_reversed, error_rate):
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



def assign_slices(nodes, edges):
    """
    Assign nodes to slices in a directed acyclic graph (DAG).
    
    Parameters:
    - nodes: List of nodes in the graph.
    - edges: List of directed edges (tuples) in the graph, where each edge is (source, target).
    
    Returns:
    - slices: List of slices, where each slice is a set of nodes.
    """
    # Step 1: Compute in-degrees for all nodes
    in_degree = {node: 0 for node in nodes}
    for source, target in edges:
        in_degree[target] += 1

    # Step 2: Initialize queue with nodes having zero in-degree
    queue = deque([node for node in nodes if in_degree[node] == 0])

    # Step 3: Assign nodes to slices
    slices = []

    while queue:
        current_slice = set()
        for _ in range(len(queue)):
            node = queue.popleft()
            current_slice.add(node)
            # Reduce in-degree of neighbors and add to queue if in-degree becomes zero
            for source, target in edges:
                if source == node:
                    in_degree[target] -= 1
                    if in_degree[target] == 0:
                        queue.append(target)
        slices.append(current_slice)

    return slices


frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/test.frag'
ploidy= 3
genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/haplotypes.csv'

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

fragment_model.construct2(input_handler, config)
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
quotient_g.construct3(input_handler, config)

e_labels_q = quotient_g.graph.edge_properties["e_label"]
v_labels_q = quotient_g.graph.vertex_properties["v_label"]
gt.graph_draw(quotient_g.graph, output_size=(500, 500), vertex_text=v_labels_q, edge_text=e_labels_q, vertex_font_size=14,  
edge_font_size=12)


quotient_g_v_label_reversed = quotient_g.v_label_reversed

edges_map_quotient = {}
for k in quotient_g.e_label_reversed.keys():
    edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]


transitions_dict = transition_matrices(quotient_g, edges_map_quotient)
emission_dict = emissions(ploidy, quotient_g_v_label_reversed, config.error_rate)

nodes = list(emission_dict.keys())
edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

slices = assign_slices(nodes, edges)
print(slices)


