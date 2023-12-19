import networkx as nx
import itertools
import numpy as np
from utils.utils import phas_2_str, get_matching_reads_for_positions
from algorithm.haplotype_assembly_helper import generate_phasings_ploidy_long, compute_likelihood
from scipy.stats import entropy


class QuotientGraph:
    def __init__(self, frag_graph):
        self.fragment_graph = frag_graph
        
    
    def construct2(self):
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
    
    def construct(self, fragment_list, ploidy, error_rate):
        output_graph = nx.Graph()
        for u, v, attrs in self.fragment_graph.edges(data=True):
            # print(u, v, attrs)
            output_graph.add_node(f"{u}-{v}", weight=attrs['weight'])
        # plot_graph(output_graph)
        for node in self.fragment_graph:
            node_list = []
            for nbr, datadict in self.fragment_graph.adj[node].items():
                # print('-'.join(sorted([str(node), str(nbr)])))
                node_list.append('-'.join(sorted([str(node), str(nbr)])))
            for ed in list(itertools.combinations(node_list, 2)):
                if ed not in output_graph.edges():
                    poss = sorted(list(set(ed[0].split('-') + ed[1].split('-'))))
                    poss_genotype = ''.join(['1' for i in poss])
                    all_phasings = generate_phasings_ploidy_long(ploidy, poss_genotype, allel_set=[0, 1])
                    all_obs = get_matching_reads_for_positions([int(i) for i in poss], fragment_list)
                    weights = {phas_2_str(phas): 0 for phas in all_phasings}
                    for phas in all_phasings:
                        for obs in all_obs:
                            obs_np = np.array([int(po) for po in obs])
                            weights[phas_2_str(phas)] += compute_likelihood(obs_np, phas, error_rate)
                    entr = entropy(list(weights.values()), base=2)
                    output_graph.add_edge(ed[0], ed[1], weight=weights, entropy=entr)
        return output_graph
