import networkx as nx
import itertools


class QuotientGraph:
    def __init__(self, frag_graph):
        self.fragment_graph = frag_graph

    def construct(self):
        output_graph = nx.Graph()
        for u, v in self.fragment_graph.edges():
            # print(u, v)
            output_graph.add_node(f"{u}-{v}")
    
        for node in self.fragment_graph:
            node_list = []
            for nbr, datadict in self.fragment_graph.adj[node].items():
                # print('-'.join(sorted([str(node), str(nbr)])))
                node_list.append('-'.join(sorted([str(node), str(nbr)])))
            for ed in list(itertools.combinations(node_list, 2)):
                if ed not in output_graph.edges():
                    output_graph.add_edge(ed[0], ed[1])
    
        return output_graph
    