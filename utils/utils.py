import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import itertools
import os
import subprocess
import shutil
import pandas as pd
# import networkit as nk
import pickle
import graph_tool.all as gt
from scipy.stats import entropy


# @profile
def get_matching_reads_for_positions(pos, fragment_list):
    # pos = [1, 2, 3]
    
    obs_positions = fragment_list[::2]
    obs_reads = fragment_list[1::2]
    
    matching_obs = [(obs_positions[idx], obs_reads[idx]) for idx in range(len(obs_positions)) if
                    len(set(pos).intersection(set(obs_positions[idx]))) > 1]
    
    # for idx in range(len(obs_positions)):
    #     if len(set(pos).intersection(set(obs_positions[idx]))) > 1:
    #         print(idx, obs_positions[idx], obs_reads[idx])


    matches = []
    for mat in matching_obs:
        # print(mat)
        shared_idx = [i for i, element in enumerate(mat[0]) if element in pos]
        # print([shared_idx, mat[0], mat[1]])
        matches.append([shared_idx, mat[0], mat[1]])
        # print('-------------------------------------')
        # [mat[0][i] for i in shared_idx]
        # [mat[1][i] for i in shared_idx]
    
    return matches


def str_2_phas(phasings, ploidy):
    return [np.array([int(p) for p in [*phas]]).reshape(ploidy, -1) for phas in phasings]


def str_2_phas_1(phasing, ploidy):
    return np.array([int(p) for p in [*phasing]]).reshape(ploidy, -1)


def phas_2_str(phas):
    return ''.join([str(ph) for ph in list(np.ravel(phas))])


def plot_graph(graph):
    pos = nx.spring_layout(graph, weight=None, seed=9)
    nx.draw(graph, with_labels=True, node_color='tab:pink', node_size=500, pos=pos)
    plt.show()


def is_chordal(G):
    """ Check if the graph is chordal """
    return nx.is_chordal(G)


def contract_edge(G, u, v):
    """ Contract edge (u, v) in graph G """
    G = nx.contracted_edge(G, (u, v), self_loops=False)
    return G


def find_cycle_without_chord(G):
    """ Find a cycle without a chord in the graph """
    for cycle in nx.cycle_basis(G):
        if len(cycle) >= 4:
            # Check if the cycle has a chord
            has_chord = any((u, v) in G.edges() for u, v in itertools.combinations(cycle, 2) if u not in cycle or v not in cycle)
            if not has_chord:
                return cycle
    return None


def convert_to_chordal(G):
    """ Convert graph G to a chordal graph using edge contraction """
    while not is_chordal(G):
        cycle = find_cycle_without_chord(G)
        if cycle:
            # first edge or sample
            G = contract_edge(G, cycle[0], cycle[1])
        else:
            break
    return G


def longest_common_substring(s1, s2):
    m, n = len(s1), len(s2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    length = 0
    end_pos = 0

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
                if dp[i][j] > length:
                    length = dp[i][j]
                    end_pos = i
            else:
                dp[i][j] = 0

    return s1[end_pos - length: end_pos]


def merge_phasings_dict(u_phase, v_phase, v, shared_dict):
    shared_idx = shared_dict[f'{v}']['shared_idx'][0][1:]
    non_redundants = shared_dict[f'{v}']['non_redundant_sites']
    u_sorted = u_phase[u_phase[:, shared_idx[0]].argsort()]
    v_sorted = v_phase[v_phase[:, shared_idx[1]].argsort()]
    matx = np.hstack([u_sorted, v_sorted])
    nr_matx = matx[:, non_redundants]
    return nr_matx


def merge_one_row2(samples, columns, row, shared_dict2):
    u = columns[0]
    u_phase = str_2_phas([samples.loc[row, u]], 3)[0]
    for co in range(1, len(columns)):
        v = columns[co]
        v_phase = str_2_phas([samples.loc[row, v]], 3)[0]
        u_phase = merge_phasings_dict(u_phase, v_phase, v, shared_dict2)
    return phas_2_str(u_phase)


def merge_one_row(columns, row, shared_dict):
    u = columns[0]
    u_phase = str_2_phas([row[u]], 3)[0]
    for co in range(1, len(columns)):
        v = columns[co]
        v_phase = str_2_phas([row[v]], 3)[0]
        u_phase = merge_phasings_dict(u_phase, v_phase, v, shared_dict)
    return phas_2_str(u_phase)

    
    
def make_shared_dict_general(samples_df):
    columns = list(samples_df.columns.values)
    shared_dict = {}
    for co in range(len(columns) - 1):
        col1, col2 = columns[co], columns[co + 1]
        shared_dict[f'{col1}:{col2}'] = {}
        corresponding_indices = [(item, i, col2.split('-').index(item)) for i, item in enumerate(col1.split('-')) if
                                 item in col2.split('-')]
        shared_dict[f'{col1}:{col2}']['shared_idx'] = corresponding_indices

        u_sites = [(i, elem) for i, elem in enumerate(col1.split('-'))]
        v_sites = [(i + len(u_sites), elem) for i, elem in enumerate(col2.split('-'))]

        redundant_sites = [elem[0] for elem in v_sites if elem[1] in col1.split('-')]

        non_redundant_sites = [col for col in range(len(u_sites) + len(v_sites)) if col not in redundant_sites]
        # shared_dict[f'{col1}:{col2}']['redundant_sites'] = redundant_sites
        shared_dict[f'{col1}:{col2}']['non_redundant_sites'] = non_redundant_sites
    return shared_dict


def make_shared_dict(samples_df):
    columns = list(samples_df.columns.values)
    col1 = columns[0]
    shared_dict = {}
    for co in range(1, len(columns)):
        col2 = columns[co]
        shared_dict[f'{col2}'] = {}
        corresponding_indices = [(item, i, col2.split('-').index(item)) for i, item in enumerate(col1.split('-')) if
                                 item in col2.split('-')]
        
        shared_dict[f'{col2}']['shared_idx'] = corresponding_indices
        
        u_sites = [(i, elem) for i, elem in enumerate(col1.split('-'))]
        v_sites = [(i + len(u_sites), elem) for i, elem in enumerate(col2.split('-'))]
        
        redundant_sites = [elem[0] for elem in v_sites if elem[1] in col1.split('-')]
        
        non_redundant_sites = [col for col in range(len(u_sites) + len(v_sites)) if col not in redundant_sites]
        # shared_dict[f'{col1}:{col2}']['redundant_sites'] = redundant_sites
        shared_dict[f'{col2}']['non_redundant_sites'] = non_redundant_sites
        col1 = '-'.join(sorted(list(set(col1.split('-') + col2.split('-')))))
        
    return shared_dict

def wsl_available() -> bool:
    """
    heuristic to detect if Windows Subsystem for Linux is available.

    Uses presence of /etc/os-release in the WSL image to say Linux is there.
    This is a de facto file standard across Linux distros.
    """
    print('before3')
    print(os.name)
    if os.name == "nt":
        print('before4')
        wsl = shutil.which("wsl")
        if not wsl:
            return False
        # can't read this file or test with
        # pathlib.Path('//wsl$/Ubuntu/etc/os-release').
        # A Python limitation?
        ret = subprocess.run(["wsl", "test", "-f", "/etc/os-release"])
        return ret.returncode == 0
    print('before5')
    return False


def give_marginals(factor_graph, qg, beliefs):
    variable_nodes = [node for node in factor_graph.nodes() if isinstance(node, str)]
    
    marginals = {}
    max_phasings = {}
    for variable in variable_nodes:
        marginal = beliefs.query(variables=[variable], show_progress=False)
        marginals[variable] = marginal
        
        max_index = list(marginal.values).index(max(list(marginal.values)))
        max_phase = list(qg.nodes[variable]['weight'])[max_index]
        max_phasings[variable] = max_phase
        # print(variable, ':', max_phase)
    return marginals, max_phasings


def creat_vcf(max_phase, positions, config):
    max_ph_np = str_2_phas([max_phase], config.ploidy)[0]
    h_df = pd.DataFrame(columns=['Position', 'phasing'], index=range(len(positions)))
    for id, pos in enumerate(positions):
        h_df.loc[id, :] = pos, '|'.join([str(i) for i in list(max_ph_np[:,id])])
    return h_df


# def dfs(graph, start, visited=None, path=None):
#     if visited is None:
#         visited = set()
#     if path is None:
#         path = [start]
#     visited.add(start)
#     for neighbor in graph.iterNeighbors(start):
#         if neighbor in visited:
#             if neighbor == path[-2]:
#                 # Skip the immediate parent node in path
#                 continue
#             # Cycle found
#             cycle = path[path.index(neighbor):] + [neighbor]
#             yield cycle
#         else:
#             yield from dfs(graph, neighbor, visited, path + [neighbor])
#     visited.remove(start)


# def networkit_find_cycles(graph):
#     visited = set()
#     for node in graph.iterNodes():
#         if node not in visited:
#             yield from dfs(graph, node, visited)
#
#
# def cycle_is_chordless(graph, cycle):
#     for i in range(len(cycle)):
#         for j in range(i+2, len(cycle) - (1 if i == 0 else 0)):
#             if graph.hasEdge(cycle[i], cycle[j]):
#                 return False
#     return True


def networkit_is_chordal(graph):
    edges = list(graph.iterEdges())
    # nodes = list(graph.iterNodes())
    tempnx = nx.Graph()
    tempnx.add_edges_from(edges)
    chordless_cycles = [cyc for cyc in list(nx.chordless_cycles(tempnx)) if len(cyc) > 3]
    if len(chordless_cycles) > 0:
        return False
    else:
        return True


def nx2nk(graphnx):
    nk_node_map = {node: i for i, node in enumerate(graphnx.nodes())}
    reverse_map = {i: node for node, i in nk_node_map.items()}
    graphnk = nk.Graph(graphnx.number_of_nodes(), weighted=False, directed=graphnx.is_directed())
    for u, v in graphnx.edges():
        graphnk.addEdge(nk_node_map[u], nk_node_map[v])
    return graphnk, reverse_map


def nk2nx(graphnk, reverse_map):
    graphnx = nx.Graph()
    nodes = list(reverse_map.values())
    graphnx.add_nodes_from(nodes)
    for u, v in graphnk.iterEdges():
        graphnx.add_edge(reverse_map[u], reverse_map[v])
    return graphnx


def networkit_find_cliques(graphnk):
    cliqueFinder = nk.clique.MaximalCliques(graphnk)
    cliqueFinder.run()
    cliques = cliqueFinder.getCliques()
    return cliques

# def cycle_is_chordless(graph, cycle):
#     cycle_length = len(cycle)
#     for i in range(cycle_length):
#         for j in range(i + 2, cycle_length + (i > 0)):
#             node1 = cycle[i]
#             node2 = cycle[j % cycle_length]
#             if graph.hasEdge(node1, node2):
#                 return False
#     return True

# def is_hole(graph, cycle):
#     cycle_length = len(cycle)
#     if cycle_length < 4:
#         return False  # A hole must have at least 4 nodes
#     for i in range(cycle_length):
#         for j in range(i + 2, cycle_length + (i > 0)):
#             node1 = cycle[i]
#             node2 = cycle[j % cycle_length]
#             # Check if there's an edge between non-consecutive nodes
#             if graph.hasEdge(node1, node2):
#                 return False  # Not a hole if there's a chord
#     return True
#
# def is_hole(graph, cycle):
#     cycle_length = len(cycle)
#     for i in range(cycle_length):
#         for j in range(i + 2, cycle_length + (i > 0) - 2):
#             # print(i, j)
#             node1 = cycle[i]
#             node2 = cycle[j % cycle_length]
#             if not graph.hasEdge(node1, node2):
#                 return True
#     return False
#
# def simple_cycles(graph, start, visited=None, path=None, cycles=None):
#     if visited is None:
#         visited = set()
#     if path is None:
#         path = []
#     if cycles is None:
#         cycles = []
#     visited.add(start)
#     path.append(start)
#     for neighbor in graph.iterNeighbors(start):
#         if neighbor not in visited:
#             simple_cycles(graph, neighbor, visited, path, cycles)
#         elif len(path) > 2 and neighbor == path[-3]:
#             # Found a cycle
#             cycle = path[path.index(neighbor):]
#             cycles.append(cycle)
#     path.pop()
#     if len(path) == 0:  # Reset visited set when unwinding back to start
#         visited.clear()
#     return cycles
#
#
# def dfs(graph, node, start, visited, stack, cycles):
#     visited[node] = True
#     stack.append(node)
#     for neighbor in graph.iterNeighbors(node):
#         if neighbor == start and len(stack) > 2:
#             # Found a cycle, adding to cycles list
#             cycles.append(stack.copy())
#         elif not visited[neighbor]:
#             dfs(graph, neighbor, start, visited, stack, cycles)
#     stack.pop()  # Remove current node from stack on backtrack
#     visited[node] = False  # Mark node as unvisited on backtrack
#
# def find_simple_cycles(graph):
#     cycles = []
#     visited = [False] * graph.numberOfNodes()
#     stack = []
#     for node in graph.iterNodes():
#         dfs(graph, node, node, visited, stack, cycles)
#         visited[node] = True  # Ensure each cycle is only found once
#     return cycles
#
def find_simple_cycles2(graphnk):
    edges = list(graphnk.iterEdges())
    nodes = list(graphnk.iterNodes())
    tempnx = nx.Graph()
    tempnx.add_edges_from(edges)
    nx.chordless_cycles(tempnx)

def nk2nx_simple(graphnk):
    edges = list(graphnk.iterEdges())
    tempnx = nx.Graph()
    tempnx.add_edges_from(edges)
    return tempnx


def convert_to_int_list(s):
    """Function to convert a string like '5-12' to a list of integers [5, 12]"""
    return [int(x) for x in s.split('-')]


def read_fragment_graph(main_path, contig, coverage, frag_file):
    fragment_graph_path = os.path.join(main_path, 'fragment_graphs', contig, coverage)
    reversed_map_path = os.path.join(main_path, 'reverse_maps', contig, coverage)

    file_name = frag_file.split('.')[0]

    frag_graph_path = os.path.join(fragment_graph_path, file_name + '.gt.gz')

    fragment_v_label_revered_path = os.path.join(reversed_map_path, 'fg_v_label_' + file_name + '.pkl')
    fragment_e_label_revered_path = os.path.join(reversed_map_path, 'fg_e_label_' + file_name + '.pkl')

    with open(fragment_v_label_revered_path, "rb") as f:
        v_label_reversed = pickle.load(f)

    with open(fragment_e_label_revered_path, "rb") as f:
        e_label_reversed = pickle.load(f)
    
    g_loaded = gt.load_graph(frag_graph_path)
    return g_loaded, v_label_reversed, e_label_reversed

# @profile
def read_quotient_graph(main_path, contig, coverage, frag_file):
    quotient_graph_path = os.path.join(main_path, 'quotient_graphs', contig, coverage)
    reversed_map_path = os.path.join(main_path, 'reverse_maps', contig, coverage)

    file_name = frag_file.split('.')[0]

    quotient_v_label_revered_path = os.path.join(reversed_map_path, 'qg_v_label_' + file_name + '.pkl')
    quotient_e_label_revered_path = os.path.join(reversed_map_path, 'qg_e_label_' + file_name + '.pkl')

    quot_graph_path = os.path.join(quotient_graph_path, file_name + '.gt.gz')

    with open(quotient_v_label_revered_path, "rb") as f:
        v_label_reversed = pickle.load(f)

    with open(quotient_e_label_revered_path, "rb") as f:
        e_label_reversed = pickle.load(f)

    g_loaded = gt.load_graph(quot_graph_path)
    return g_loaded, v_label_reversed, e_label_reversed



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


def find_common_element_and_index_generalized(node1, node2):
    # Split the node names into components (e.g., '1-2' -> ['1', '2'])
    node1_parts = node1.split('-')
    node2_parts = node2.split('-')
    
    # Find the common elements and their indices in both nodes
    common_indices_node1 = []
    common_indices_node2 = []
    
    for i, part in enumerate(node1_parts):
        if part in node2_parts:
            # Collect all occurrences of this common element in both lists
            common_indices_node1.append(i)
            common_indices_node2.append(node2_parts.index(part))
    
    if not common_indices_node1:
        raise ValueError(f"No common elements found between {node1} and {node2}")
    
    return common_indices_node1, common_indices_node2


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


def find_matchings_generalized(part1, part2):

    # Sort both parts and remember the original indices.
    sorted_part1 = sorted(enumerate([''.join([str(s) for s in list(ss)]) for ss in part1]), key=lambda x: x[1])
    sorted_part2 = sorted(enumerate([''.join([str(s) for s in list(ss)]) for ss in part2]), key=lambda x: x[1])
    
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


# @profile
def shortest_path_excluding_node(graph, source, target, exclude_node):
    # Create a vertex filter property map (all True initially)
    v_filter = graph.new_vertex_property("bool", val=True)
    
    # Exclude the specified node by setting its filter to False
    v_filter[graph.vertex(exclude_node)] = False
    
    # Create a filtered view of the graph
    subgraph_2 = gt.GraphView(graph, vfilt=v_filter)
    
    # Find the shortest path in the filtered graph
    path_vertices, path_edges = gt.shortest_path(subgraph_2, source=source, target=target)

    # Convert vertices and edges to lists for readability
    # path_v_list = [int(v) for v in path_vertices]
    # path_e_list = [(int(e.source()), int(e.target())) for e in path_edges]
    # path_e_list = path_edges + [source, target]
    # Clear the filter after use (restore original graph state)
    graph.set_vertex_filter(None)

    return path_edges


def plot_subgraph_graph_tool(new_graph, included_vertices):
    v_filter = new_graph.new_vertex_property("bool")

    # Set the filter to True for vertices you want to keep
    for v in new_graph.vertices():
        if int(v) in included_vertices:
            v_filter[v] = True

    # Create a graph view using the vertex filter
    subgraph = gt.GraphView(new_graph, vfilt=v_filter)

    # Optionally: Convert the graph view into a new standalone graph
    # subgraph_copy = gt.Graph(subgraph, prune=True)

    e_labels = subgraph.edge_properties["e_label"]
    v_labels = subgraph.vertex_properties["v_label"]
    # gt.graph_draw(subgraph, output_size=(500, 500), vertex_text=subgraph.vertex_index, edge_text=e_labels, vertex_font_size=16,  
    # edge_font_size=12)
    gt.graph_draw(subgraph, output_size=(500, 500), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=14,  
    edge_font_size=12)

# @profile
def compute_phasings_and_weights(e_poss, input_handler, config, fragment_model):
    poss_genotype = input_handler.get_genotype_positions([int(p) for p in e_poss])
    all_phasings = generate_phasings_ploidy_long(config.ploidy, poss_genotype, allel_set=[0, 1])
    matches = get_matching_reads_for_positions([int(i) for i in e_poss], fragment_model.fragment_list)
    weights = {phas_2_str(phas): 0 for phas in all_phasings}
    for phas in all_phasings:
        for indc, this_po, obs in matches:
            # print(indc, this_po, obs)
            # for obs in all_obs:
            #     obs_np = np.array([int(po) for po in obs])
            #     weights[phas_2_str(phas)] += compute_likelihood(obs_np, phas, error_rate)
            weights[phas_2_str(phas)] += compute_likelihood_generalized_plus(np.array(obs), phas, indc,
                                                                                list(range(len(indc))),
                                                                                config.error_rate)
    entr = entropy(list(weights.values()), base=10)
    weight =  {"weight": weights, "entropy": entr}
    return weight


def generate_phasings_ploidy_local():
        source = edge.source()
        target = edge.target()
        source_weights = quotient_graph.vertex_properties["v_weights"][source]['weight']
        target_weights = quotient_graph.vertex_properties["v_weights"][target]['weight']
        source_label = quotient_graph.vertex_properties["v_label"][source]
        target_label = quotient_graph.vertex_properties["v_label"][target]
        # e_label = quotient_graph.edge_properties["e_label"][edge]
        e_label = '--'.join(sorted([source_label, target_label]))
        e_weights = quotient_graph.edge_properties["e_weights"][edge]['weight']

        
        common_ff, common_sf = find_common_element_and_index(source_label, target_label)
        source_phasings = list(source_weights.keys())
        target_phasings = list(target_weights.keys())
        # transitions_dict = {'source': source_phasings, 'target': target_phasings}
        transitions_mtx = np.zeros((len(source_phasings), len(target_phasings)))
        for i, ffstr in enumerate(source_phasings):
            for j, sfstr in enumerate(target_phasings):
                
                matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, 3), str_2_phas_1(sfstr, 3), common_ff, common_sf, source_label, target_label)
                sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)
                matched_phasings_str = [phas_2_str(pm) for pm in sorted_phasings] 
                this_weight = np.sum([e_weights[pm] for pm in matched_phasings_str if pm in e_weights.keys()])
                transitions_mtx[i, j] = this_weight

# @profile
def compute_edge_weight(new_vertex_name, v_label, source_phasings, target_phasings, fragment_model, config):
    possitions = sorted(set([int(nn) for nn in new_vertex_name.split('-')] + [int(nn) for nn in v_label.split('-')]))
    common_ff, common_sf = find_common_element_and_index(new_vertex_name, v_label)
    all_phasings =[]
    for ffstr in source_phasings:
        for sfstr in target_phasings:
            
            matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, 3), str_2_phas_1(sfstr, 3), common_ff, common_sf, new_vertex_name, v_label)
            
            sorted_phasings = []
            for mtx in matched_phasings:
                sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                sorted_phasings.append(sorted_matrix)
            matched_phasings_str = [phas_2_str(pm) for pm in sorted_phasings]
            all_phasings += matched_phasings_str
    matches = get_matching_reads_for_positions(possitions, fragment_model.fragment_list)
    weights = {phas_2_str(phas): 0 for phas in all_phasings}
    for phas in all_phasings:
        for indc, this_po, obs in matches:
            
            weights[phas_2_str(phas)] += compute_likelihood_generalized_plus(np.array(obs), str_2_phas_1(phas, config.ploidy), indc,
                                                                                list(range(len(indc))),
                                                                                config.error_rate)
    entr = entropy(list(weights.values()), base=10)
    final_weight = {"weight": weights, "entropy": entr}
    return final_weight

# @profile
def sort_strings(strings):
    # Sort the list of strings using the custom comparison logic
    return sorted(strings, key=lambda x: list(map(int, x.split('-'))))


def get_minimum_spanning_tree(quotient_graph):
    quotient_graph.clear_filters()
    e_entropy = quotient_graph.new_edge_property("double")
    for e in quotient_graph.edges():
        e_entropy[e] = quotient_graph.ep['e_weights'][e]['entropy']
    quotient_graph.ep['e_entropy'] = e_entropy
    tree = gt.min_spanning_tree(quotient_graph, weights=e_entropy)
    
    mst_graph = gt.GraphView(quotient_graph, efilt=tree)
    non_mst_graph = gt.GraphView(quotient_graph, efilt=lambda e: not tree[e])

    return mst_graph, non_mst_graph, tree


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


def compute_likelihood_generalized(observed, phasing, used_pos, error_rate):
    """This likelihood computation can accept shorter observed and longer or equal length phasing, then the used pos
    are the positions on the phasing that should be used"""
    assert len(used_pos) == len(observed)
    new_phasing = phasing[:, used_pos]
    y = np.tile(observed, (new_phasing.shape[0], 1))
    diff = y - new_phasing
    diff[diff != 0] = 1
    comp_diff = 1 - diff
    term1 = diff * error_rate
    term2 = comp_diff * (1 - error_rate)
    terms = term1 + term2
    probs = np.prod(terms, axis=1)
    likelihood = np.mean(probs)
    return likelihood

# @profile
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


def plot_subgraph_graph_tool_simple(new_graph, included_vertices):
    v_filter = new_graph.new_vertex_property("bool")

    if not included_vertices is None:
        for v in new_graph.vertices():
            if int(v) in included_vertices:
                v_filter[v] = True
        subgraph = gt.GraphView(new_graph, vfilt=v_filter)

        gt.graph_draw(subgraph, output_size=(500, 500), vertex_font_size=14,edge_font_size=12, vetex_size=10)
    else:
        gt.graph_draw(new_graph, output_size=(500, 500), vertex_font_size=14, edge_font_size=12, vetex_size=10)