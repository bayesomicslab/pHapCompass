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

def normalize_potentials(clique_graph):
    """
    Normalize the potentials in v_weights so that they sum to 1.
    """
    for v in clique_graph.vertices():
        data = clique_graph.vp['v_weights'][v]
        total_weight = sum(data['weight'].values())
        
        if total_weight > 0:
            normalized_weights = {k: v / total_weight for k, v in data['weight'].items()}
            # Update the vertex property after normalizing
            data['weight'] = normalized_weights
        else:
            # If total weight is zero, assign uniform distribution
            num_phasings = len(data['weight'])
            if num_phasings > 0:
                uniform_weight = 1.0 / num_phasings
                data['weight'] = {k: uniform_weight for k in data['weight'].keys()}
            else:
                print(f"Warning: No phasings for node {v}, skipping normalization.")



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


def chordal_contraction_fragment_graph(frag_graph):
    # frag_graph = fragment_model.graph
    # more efficient version of the chordal contraction
    new_graph = frag_graph.copy()
    # print("Vertex Properties:", list(new_graph.vp.keys()))
    # print("Edge Properties:", list(new_graph.ep.keys()))

    # e_weights = new_graph.edge_properties["e_weights"]
    frag_e_labels = new_graph.edge_properties["e_label"]
    # frag_v_labels = new_graph.vertex_properties["v_label"]

    new_graph.clear_filters()
    # e_entropy = new_graph.new_edge_property("double")
    n_e_reads = new_graph.new_edge_property("int")
    # n_v_reads = new_graph.new_vertex_property("int")

    fragment_list_positions = fragment_model.fragment_list[::2] # [1::2] # [::2]

    for e in new_graph.edges():
        # print(e)
        poss = [int(i) for i in frag_e_labels[e].split('-')]
        relev_reads = [r for r in fragment_list_positions if all(poss) in r]
        n_e_reads[e] = len(relev_reads)
    
    new_graph.ep['n_e_reads'] = n_e_reads

    # # Loop over edges and assign entropy from the e_weights property
    # for e in new_graph.edges():
    #     e_entropy[e] = e_weights[e]['entropy']

    # new_graph.ep['e_entropy'] = e_entropy

    chordless_cycles = get_chordless_cycles(new_graph)

    to_be_removed_nodes = []

    for cyc_id, cyc in enumerate(chordless_cycles):
        # stop
        # print(cyc_id)
        
        
        edges = [new_graph.edge(cyc[-1], cyc[0])]
        for i in range(len(cyc) - 1):
            edges += [new_graph.edge(cyc[i], cyc[i+1])]
        edges = [x for x in edges if x is not None]
        while len(edges) > 3:
            min_edge = min(edges, key=lambda e: new_graph.ep['n_e_reads'][e])
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
                # (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                # first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                # second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                # final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)

                e1 = new_graph.edge(min_edge.source(), n)
                e2 = new_graph.edge(min_edge.target(), n)
                
                # new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name

                n_poss = sorted(set([int(nn) for nn in v_label.split('-')] + poss))
                relev_reads = [r for r in fragment_list_positions if all(n_poss) in r]
                new_graph.edge_properties['n_e_reads'][e1] = len(relev_reads)

                # new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
                new_graph.remove_edge(e2)

            for n in set(source_nbrs)-common_nbrs:
                # stop
                v_label = new_graph.vertex_properties["v_label"][n]

                sorted_labels = sort_nodes([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)
                # (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                # first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                # second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                # final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)

                e1 = new_graph.edge(min_edge.source(), n)
                # e2 = new_graph.edge(min_edge.target(), n)
                # new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                # new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
                # new_graph.edge_properties["e_weights"][e2]
                n_poss = sorted(set([int(nn) for nn in v_label.split('-')] + poss))
                relev_reads = [r for r in fragment_list_positions if all(n_poss) in r]
                new_graph.edge_properties['n_e_reads'][e1] = len(relev_reads)



            for n in set(target_nbrs)-common_nbrs:
                # stop
                v_label = new_graph.vertex_properties["v_label"][n]
                sorted_labels = sort_nodes([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)

                # (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                # first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                # second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                # final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)
                
                e2 = new_graph.edge(min_edge.target(), n)
                flag = False
                if e2 in edges:
                    edges.remove(e2)
                    flag = True
                new_graph.remove_edge(e2)
                e1 = new_graph.add_edge(min_edge.source(), n)
                if flag:
                    edges.append(e1)
                # new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                # new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
                n_poss = sorted(set([int(nn) for nn in v_label.split('-')] + poss))
                relev_reads = [r for r in fragment_list_positions if all(n_poss) in r]
                new_graph.edge_properties['n_e_reads'][e1] = len(relev_reads)
            
            # to_be_removed_nodes += [min_edge.target()]
            to_be_removed_nodes.append(min_edge.target())
            new_graph.remove_edge(min_edge)
            edges.remove(min_edge)
    # print('Nodes to be removed:', to_be_removed_nodes)
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



# Helper function for log-sum-exp to avoid numerical issues
def logsumexp(arr):
    max_val = np.max(arr)
    return max_val + np.log(np.sum(np.exp(arr - max_val)))



# Step 1: Message Passing for Exact Inference
def send_message(clique_graph, clique_potentials, src, tgt):
    """
    Send message from source clique (src) to target clique (tgt) in log space.
    """
    separator_vars = set(clique_potentials[src].keys()) & set(clique_potentials[tgt].keys())
    
    # Log-space message initialization
    message = {sep_state: -np.inf for sep_state in separator_vars}  # Log(0) = -inf

    for sep_state in separator_vars:
        # Sum over states in src clique not in separator
        relevant_phasings = [
            (state, prob) for state, prob in clique_potentials[src].items() if sep_state in state
        ]
        log_probs = np.array([prob for state, prob in relevant_phasings])
        message[sep_state] = logsumexp(log_probs)

    return message

# Step 2: Upward Pass (Collect Evidence)
def collect_evidence(clique_graph, clique_potentials, root_node):
    """
    Collect evidence by passing messages upward from the leaves to the root.
    """
    messages = {}
    stack = [(root_node, None)]  # Stack for DFS traversal
    visited = set()

    # DFS traversal to visit leaves first
    while stack:
        current, parent = stack.pop()
        if current in visited:
            continue

        visited.add(current)
        for neighbor in current.all_neighbors():
            if neighbor not in visited:
                stack.append((neighbor, current))

        # Send message from current to parent
        if parent is not None:
            messages[(current, parent)] = send_message(clique_graph, clique_potentials, current, parent)

    return messages

# Step 3: Downward Pass (Distribute Evidence)
def distribute_evidence(junction_tree, clique_potentials, messages, root_node):
    """
    Distribute evidence by passing messages downward from the root to the leaves.
    """
    beliefs = {v: clique_potentials[v].copy() for v in junction_tree.vertices()}  # Initialize beliefs
    stack = [(root_node, None)]  # Stack for DFS traversal (reverse direction)

    # Downward message passing
    while stack:
        current, parent = stack.pop()
        for neighbor in current.all_neighbors():
            if neighbor != parent:
                stack.append((neighbor, current))
                # Update belief at the neighbor node using the incoming message from current
                belief_update = send_message(junction_tree, clique_potentials, current, neighbor)
                for state, log_prob in belief_update.items():
                    beliefs[neighbor][state] += log_prob

    return beliefs

# Step 4: Exact Inference
def exact_inference_junction_tree(clique_graph):
    """
    Perform exact inference using belief propagation on a junction tree.

    Parameters:
        junction_tree (Graph): The constructed junction tree with weight properties on vertices.

    Returns:
        beliefs (dict): Dictionary of posterior beliefs for each node in the junction tree.
    """
    # Step 1: Initialize clique potentials (log-space)
    clique_potentials = {}
    for v in clique_graph.vertices():
        # Potential function for this clique
        clique_weights = clique_graph.vp["v_weights"][v]['weight']
        clique_potentials[v] = {state: np.log(prob) for state, prob in clique_weights.items()}  # Convert to log


    # for v in clique_graph.vertices():
    #     print(v, clique_graph.vp['v_weights'][v]['weight'])
    
    sorted_labels = sort_nodes([clique_graph.vertex_properties["v_label"][v] for v in clique_graph.vertices()])
    root_node = sorted_labels[0]
    root_node_id = clique_vertex_reversed[root_node]

    # Step 2: Upward pass (collect evidence)
    messages_up = collect_evidence(clique_graph, clique_potentials, root_node_id)

    # Step 3: Downward pass (distribute evidence)
    beliefs = distribute_evidence(clique_graph, clique_potentials, messages_up, root_node)

    # Step 4: Normalize beliefs
    for v in beliefs:
        log_probs = np.array(list(beliefs[v].values()))
        log_norm_const = logsumexp(log_probs)
        beliefs[v] = {state: np.exp(log_prob - log_norm_const) for state, log_prob in beliefs[v].items()}  # Normalize

    return beliefs




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

    frag_graph = fragment_model.graph

    chordal_frag = chordal_contraction_fragment_graph(frag_graph)
    print('Chordal Fragment Graph constructed.')

    ch_e_labels = chordal_frag.edge_properties["e_label"]
    ch_v_labels = chordal_frag.vertex_properties["v_label"]
    gt.graph_draw(chordal_frag, output_size=(500, 500), vertex_text=ch_v_labels, edge_text=ch_e_labels, vertex_font_size=14,  
    edge_font_size=12)


    # # create quotient graph
    # quotient_g = QuotientGraph(fragment_model)
    # quotient_g.construct(input_handler, config)

    # e_labels_q = quotient_g.graph.edge_properties["e_label"]
    # v_labels_q = quotient_g.graph.vertex_properties["v_label"]
    # gt.graph_draw(quotient_g.graph, output_size=(500, 500), vertex_text=v_labels_q, edge_text=e_labels_q, vertex_font_size=14,  
    # edge_font_size=12)

    # quotient_g_v_label_reversed = quotient_g.v_label_reversed

    # edges_map_quotient = {}
    # for k in quotient_g.e_label_reversed.keys():
    #     edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]


    # chordal_graph = chordal_contraction_graph_tool(quotient_g.graph, config, fragment_model)

    # v_labels_c = chordal_graph.vertex_properties["v_label"]
    # e_labels_c = chordal_graph.edge_properties["e_label"]
    # gt.graph_draw(chordal_graph, output_size=(500, 500), vertex_text=v_labels_c, edge_text=e_labels_c, vertex_font_size=14,  
    # edge_font_size=12)

    true_haplotypes = pd.read_csv(genotype_path).T
    all_poss_genotype = input_handler.genotype[1:]
    # Step 1: Find maximal cliques using graph_tool
    cliques = []
    for clique in gt.max_cliques(chordal_frag):
        cliques.append(set(clique))


    # Step 2: Create the clique graph
    clique_graph = gt.Graph(directed=False)
    clique_map = {i: clique for i, clique in enumerate(cliques)}
    clique_indices = list(clique_map.keys())
    clique_vertices = {}
    
    vertex_labels = clique_graph.new_vertex_property("string")
    clique_vertex_reversed = {}
    for idx in clique_indices:
        indices_in_quotient = list(clique_map[idx])
        v_names_quo = [frag_graph.vertex_properties["v_label"][ii] for ii in indices_in_quotient]
        v_names_quo2 = [ii.split('-') for ii in v_names_quo]
        flattened_list = sorted(set([item for sublist in v_names_quo2 for item in sublist]))
        v_name = '-'.join([str(nnn) for nnn in flattened_list])
        clique_vertices[idx] = clique_graph.add_vertex()
        vertex_labels[clique_vertices[idx]] = v_name
        clique_vertex_reversed[v_name] = clique_vertices[idx]

    # Add edges weighted by intersection size
    edge_weights = clique_graph.new_edge_property("double")
    for i in clique_indices:
        for j in clique_indices:
            if i < j:
                weight = len(clique_map[i] & clique_map[j])
                if weight > 0:  # Only add edges for cliques with shared variables
                    e = clique_graph.add_edge(clique_vertices[i], clique_vertices[j])
                    edge_weights[e] = weight

    clique_graph.ep['intersection'] = edge_weights
    clique_graph.vp['v_label'] = vertex_labels

    # clique_e_labels = clique_graph.edge_properties["e_label"]
    clique_v_labels = clique_graph.vertex_properties["v_label"]
    gt.graph_draw(clique_graph, output_size=(500, 500), vertex_text=clique_v_labels, vertex_font_size=14,  
    edge_font_size=12)


    clique_graph.clear_filters()
    inverse_weight = clique_graph.new_edge_property("double")
    v_weights = clique_graph.new_vertex_property("object")

    fragment_list_positions = fragment_model.fragment_list[::2] # [1::2] # [::2]

    for e in clique_graph.edges():
        source, target = e.source(), e.target()
        source_label = clique_graph.vp['v_label'][source]
        target_label = clique_graph.vp['v_label'][target]
        poss = sorted(set([int(nn) for nn in source_label.split('-')] + [int(nn) for nn in target_label.split('-')]))
        shared = set([int(nn) for nn in source_label.split('-')]).intersection(set([int(nn) for nn in target_label.split('-')]))
        relev_reads = [r for r in fragment_list_positions if len(set(poss) & set(r)) > len(shared)]
        inverse_weight[e] = 1/(1 + len(relev_reads))

    clique_graph.ep['inverse_weight'] = inverse_weight
    # Step 3: Compute the Maximum Spanning inverse_weight (MST)
    tree = gt.min_spanning_tree(clique_graph, weights=inverse_weight)

    # Create a set of MST edges (for efficient lookup)
    mst_edges_set = set((int(e.source()), int(e.target())) for e in clique_graph.edges() if tree[e])

    # Remove edges that are not part of the MST
    edges_to_remove = [e for e in clique_graph.edges() if (int(e.source()), int(e.target())) not in mst_edges_set]

    for e in edges_to_remove:
        clique_graph.remove_edge(e)


    for v in clique_graph.vertices():
        
        label = clique_graph.vp['v_label'][v]
        poss = sorted(set([int(nn) for nn in label.split('-')]))
        poss_id = [pos -1 for pos in poss]
        this_gen = ''.join([str(nn) for nn in list(true_haplotypes[poss_id].to_numpy().sum(axis=0))])
        possible_phasings = generate_phasings_ploidy_long(ploidy, this_gen)
        
        # relevant_reads = [r for r in fragment_list_positions if len(set(poss) & set(r)) > len(shared)]
        
        # 1. Extract positions and haplotypes lists
        fragment_list_positions = fragment_model.fragment_list[::2]  # Odd indices (positions)
        fragment_list_haplotypes = fragment_model.fragment_list[1::2]  # Even indices (haplotypes)

        # 2. Extract relevant reads, corresponding haplotypes, and indices
        relevant_reads = []
        relevant_haplotypes = []
        read_shared_indices = []  # To store indices of shared positions in the read
        poss_shared_indices = []  # To store indices of shared positions in `poss`

        for positions, haps in zip(fragment_list_positions, fragment_list_haplotypes):
            # Shared positions with their indices
            shared_positions = list(set(poss) & set(positions))
            if len(shared_positions) >= 2:  # Criteria for relevance
                relevant_reads.append(positions)
                relevant_haplotypes.append(haps)
                
                # Indices in the read and poss for shared positions
                read_indices = [positions.index(p) for p in shared_positions]  # Indices in the read
                poss_indices = [poss.index(p) for p in shared_positions]  # Indices in poss
                
                read_shared_indices.append(read_indices)
                poss_shared_indices.append(poss_indices)
        
        if v not in v_weights:
            v_weights[v] = {}
        
        if 'weight' not in v_weights[v]:
            v_weights[v]['weight'] = {}

        for phas in possible_phasings:
            wei = 0
            for elem_id in range(len(relevant_reads)):
                obs = relevant_haplotypes[elem_id]
                obs_pos = read_shared_indices[elem_id]
                phas_pos = poss_shared_indices[elem_id]
                this_phas_read_weight = compute_likelihood_generalized_plus(np.array(obs), phas, obs_pos, phas_pos, config.error_rate)
                wei += this_phas_read_weight
            v_weights[v]['weight'][phas_2_str(phas)] = wei
    clique_graph.vp['v_weights'] = v_weights

    for v in clique_graph.vertices():
        print(v, clique_graph.vp['v_weights'][v]['weight'])

    normalize_potentials(clique_graph)













from graph_tool import draw

# Set node fill color to red, border color to black, and edge color to black
node_fill_color = fragment_model.graph.new_vertex_property("vector<double>")
node_border_color = fragment_model.graph.new_vertex_property("vector<double>")
edge_color = fragment_model.graph.new_edge_property("vector<double>")

# Set red for node fill color (RGB: (1, 0, 0))
for v in fragment_model.graph.vertices():
    node_fill_color[v] = (1, 0, 0, 1)  # RGBA, where A=1 means fully opaque
    node_border_color[v] = (0, 0, 0, 1)  # Black border

# Set edge color to black
for e in fragment_model.graph.edges():
    edge_color[e] = (0, 0, 0, 1)  # Black

# Draw the graph with the updated styles
gt.graph_draw(
    fragment_model.graph,
    output_size=(500, 500),
    vertex_text=v_labels,
    edge_text=e_labels,
    vertex_font_size=18,
    edge_font_size=16,
    vertex_fill_color=node_fill_color,
    vertex_color=node_border_color,  # Border color
    edge_color=edge_color,
    vertex_pen_width=2  # Thickness of node border
)