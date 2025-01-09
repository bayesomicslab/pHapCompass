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


def transition_matrices(quotient_g, edges_map_quotient):
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


def chordal_contraction_graph_tool(quotient_graph, config, fragment_model):
    # more efficient version of the chordal contraction
    new_graph = quotient_graph.copy()
    e_weights = new_graph.edge_properties["e_weights"]
    new_graph.clear_filters()
    e_entropy = new_graph.new_edge_property("double")

    # Loop over edges and assign entropy from the e_weights property
    for e in new_graph.edges():
        e_entropy[e] = e_weights[e]['entropy']

    new_graph.ep['e_entropy'] = e_entropy

    chordless_cycles = get_chordless_cycles(new_graph)

    to_be_removed_nodes = []

    for cyc_id, cyc in enumerate(chordless_cycles):
        
        print(cyc_id)
        
        edges = [new_graph.edge(cyc[-1], cyc[0])]
        for i in range(len(cyc) - 1):
            edges += [new_graph.edge(cyc[i], cyc[i+1])]
        edges = [x for x in edges if x is not None]
        while len(edges) > 3:
            min_edge = min(edges, key=lambda e: new_graph.ep['e_entropy'][e])
            source_label = new_graph.vp['v_label'][min_edge.source()]
            target_label = new_graph.vp['v_label'][min_edge.target()]
            # new node positions
            poss = sorted(set([int(nn) for nn in source_label.split('-')] + [int(nn) for nn in target_label.split('-')]))
            # new vertex properties:
            new_vertex_name = '-'.join([str(nnn) for nnn in poss])
            vertex_weights = new_graph.ep['e_weights'][min_edge]
            
            new_graph.vertex_properties["v_weights"][min_edge.source()] = vertex_weights
            new_graph.vertex_properties["v_label"][min_edge.source()] = new_vertex_name

            source_nbrs = [n for n in min_edge.source().all_neighbors() if n != min_edge.target()]
            target_nbrs = [n for n in min_edge.target().all_neighbors() if n != min_edge.source()]
            common_nbrs = set(source_nbrs).intersection(set(target_nbrs))

            for n in common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]
                # e_poss = sorted(set([int(nn) for nn in v_label.split('-')] + poss))
                # print(len(e_poss))
                # new_edge_name = '-'.join([str(nnn) for nnn in e_poss])
                sorted_labels = sort_nodes([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)
                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)

                e1 = new_graph.edge(min_edge.source(), n)
                e2 = new_graph.edge(min_edge.target(), n)
                
                new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
                new_graph.remove_edge(e2)

            for n in set(source_nbrs)-common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]

                sorted_labels = sort_nodes([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)
                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)


                e1 = new_graph.edge(min_edge.source(), n)
                # e2 = new_graph.edge(min_edge.target(), n)
                new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
                # new_graph.edge_properties["e_weights"][e2]

            for n in set(target_nbrs)-common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]
                sorted_labels = sort_nodes([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)

                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)

                e2 = new_graph.edge(min_edge.target(), n)
                new_graph.remove_edge(e2)
                e1 = new_graph.add_edge(min_edge.source(), n)
                new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
            
            # to_be_removed_nodes += [min_edge.target()]
            to_be_removed_nodes.append(min_edge.target())
            new_graph.remove_edge(min_edge)
            edges.remove(min_edge)

    new_graph.remove_vertex(to_be_removed_nodes)
    return new_graph
    
# @profile
def get_chordless_cycles(subgraph_copy):
    """
    Check if the given graph is chordal using Lexicographical Breadth-First Search (Lex-BFS).

    Parameters:
    
    graph: A graph_tool Graph object.

    Returns:
    
    True if the graph is chordal, False otherwise."""
    chordless_cycles = []
    labels = {}
    for v in subgraph_copy.vertices():
        labels[v] = []

    ordering = []
    unnumbered = set(subgraph_copy.vertices())

    n = subgraph_copy.num_vertices()
    for i in range(n, 0, -1):
        # print(i)
        # Select the unnumbered vertex with the largest label lexicographically
        max_label_vertex = None
        max_label = None
        for v in unnumbered:
            label = labels[v]
            if max_label is None or label > max_label:
                max_label = label
                max_label_vertex = v

        v = max_label_vertex
        ordering.append(v)
        unnumbered.remove(v)

        # Update labels of unnumbered neighbors
        for u in v.all_neighbors():
            if u in unnumbered:
                labels[u].append(i)

    # Reverse the ordering to get the perfect elimination ordering (PEO)
    ordering.reverse()
    position = {v: idx for idx, v in enumerate(ordering)}
    for idx, v in enumerate(ordering):
        # Check if the neighbors of v that appear later in the PEO form a clique
        
        nbrs = [u for u in v.all_neighbors() if position[u] > idx]
        if not nbrs:
            continue
        # Find the neighbor with the smallest position (earliest in PEO)
        w = min(nbrs, key=lambda u: position[u])
        # print('min neighbor:', w, 'neighbors:', nbrs)
        for u in nbrs:
            if u != w and not subgraph_copy.edge(u, w) and not subgraph_copy.edge(w, u):
                # The neighbors do not form a clique
                # print('ordering:', v, 'neighbor', u, 'min neighbor',w)
                # gt.shortest_path(subgraph_copy, u, w)
                v_filter = subgraph_copy.new_vertex_property("bool", val=True)
                v_filter[v] = False
                subgraph_copy.set_vertex_filter(v_filter)
                # Exclude the specified node by setting its filter to False
                # v_filter[subgraph_copy.vertex(v)] = False
                # subgraph_2 = gt.GraphView(subgraph_copy, vfilt=v_filter)
    
                # Find the shortest path in the filtered graph
                path_v = gt.shortest_path(subgraph_copy, source=w, target=u)[0]
                # path_edges = path_edges[::-1]
                path_v = path_v + [v]
                # print(path_v)
                chordless_cycles.append(path_v)
                # subgraph_copy.set_vertex_filter(None)

                # Reset the vertex filter after the computation
                subgraph_copy.set_vertex_filter(None)

    return chordless_cycles


if __name__ == '__main__':

    frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/test.frag'
    # frag_path = '/labs/Aguiar/pHapCompass/test/test.frag'
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


    chordal_graph = chordal_contraction_graph_tool(quotient_g.graph, config, fragment_model)

    v_labels_c = chordal_graph.vertex_properties["v_label"]
    e_labels_c = chordal_graph.edge_properties["e_label"]
    gt.graph_draw(chordal_graph, output_size=(500, 500), vertex_text=v_labels_c, edge_text=e_labels_c, vertex_font_size=14,  
    edge_font_size=12)

    

