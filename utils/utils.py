import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import itertools
import os
import subprocess
import shutil
import pandas as pd
# import networkit as nk


def get_matching_reads_for_positions(pos, fragment_list):
    # pos = [1, 2, 3]
    
    obs_positions = fragment_list[::2]
    obs_reads = fragment_list[1::2]
    
    matching_obs = [(obs_positions[idx], obs_reads[idx]) for idx in range(len(obs_positions)) if
                    len(set(pos).intersection(set(obs_positions[idx]))) > 1]
    
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
def find_simple_cycles2(graph):
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

