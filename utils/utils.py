import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import itertools


def get_matching_reads_for_positions(pos, fragment_list):

    obs_positions = fragment_list[::2]
    obs_reads = fragment_list[1::2]

    matching_obs = [obs_reads[idx] for idx in range(len(obs_positions)) if len(set(pos).intersection(set(obs_positions[idx])))==len(pos)]
    return matching_obs

def str_2_phas(phasings, ploidy):
    return [np.array([int(p) for p in [*phas]]).reshape(ploidy, -1) for phas in phasings]


def phas_2_str(phas):
    return ''.join([str(ph) for ph in list(np.ravel(phas))])


def plot_graph(graph):
    pos = nx.spring_layout(graph, seed=9)
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
