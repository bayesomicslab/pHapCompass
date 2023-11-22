import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def phas_2_str(phas):
    return ''.join([str(ph) for ph in list(np.ravel(phas))])


def plot_graph(graph):
    pos = nx.spring_layout(graph, seed=9)
    nx.draw(graph, with_labels=True, node_color='tab:pink', node_size=500, pos=pos)
    plt.show()
