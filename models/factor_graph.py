import itertools
import numpy as np
from pgmpy.models import FactorGraph
from pgmpy.factors.discrete import DiscreteFactor
from algorithm.haplotype_assembly_helper import generate_phasings_ploidy, compute_likelihood
from utils.utils import phas_2_str


class Factorgraph:
    def __init__(self, ploidy, error_rate, epsilon):
        self.ploidy = ploidy
        self.error_rate = error_rate
        self.epsilon = epsilon
        
    def get_weights(self, fragment_list, quotientg, conf, inp_hand):
        # ploidy = 3
        # error_rate = 0.01
        nodes = {}
        edges = {}
        for i in range(int(len(fragment_list) / 2)):
            # print(i)
            idx = i * 2
            positions = fragment_list[idx]
            pos_alleles = fragment_list[idx + 1]
            comb2 = zip(list(itertools.combinations(positions, 2)), list(itertools.combinations(pos_alleles, 2)))
            for pp in comb2:
                # print(pp)
                node_name = '-'.join([str(p) for p in pp[0]])
                obs = ''.join([str(p) for p in pp[1]])
                node_genotype = '11'
                node_phasings = list(conf.global_likelihoods[2][node_genotype].keys())
                if node_name not in nodes.keys():
                    nodes[node_name] = {ph: conf.global_likelihoods[2][node_genotype][ph][obs] for ph in node_phasings}
                else:
                    for ph in node_phasings:
                        nodes[node_name][ph] += conf.global_likelihoods[2][node_genotype][ph][obs]
        
            if len(positions) >= 3:
                comb2 = zip(list(itertools.combinations(positions, 2)),
                            list(itertools.combinations(pos_alleles, 2)))
                for pp in list(itertools.combinations(comb2, 2)):
                
                    if len(set(pp[0][0]).intersection(set(pp[1][0]))) > 0:
                        pp_nodes = list([pp[0][0], pp[1][0]])
                        pp_allels = list([pp[0][1], pp[1][1]])
                        print('found', pp_nodes)
                        poss_factor_name = '-'.join(['-'.join([str(q) for q in po]) for po in pp_nodes])
                    
                        positions_genotypes = [''.join([str(q) for q in inp_hand.get_genotype_positions(po)]) for po in
                                               pp_nodes]
                        phasings = [generate_phasings_ploidy(self.ploidy, pg) for pg in positions_genotypes]
                    
                        minp = int(min(min(pp[0][0]), min(pp[1][0])))
                        observed = np.zeros(3, dtype=int)
                        phasing_np3 = np.zeros([self.ploidy, 3], dtype=int)
                    
                        for ii in range(2):
                            for jj in range(2):
                                observed[pp_nodes[ii][jj] - minp] = pp_allels[ii][jj]
                        if poss_factor_name not in edges.keys():
                            edges[poss_factor_name] = {}
                            for ii in range(len(phasings[0])):
                                for jj in range(len(phasings[1])):
                                    ph_name = '-'.join([phas_2_str(phasings[0][ii]), phas_2_str(phasings[1][jj])])
                                    phasing_np3[:, np.array(pp_nodes[0]) - minp] = phasings[0][ii]
                                    phasing_np3[:, np.array(pp_nodes[1]) - minp] = phasings[1][jj]
                                    thisl = compute_likelihood(observed, phasing_np3, self.error_rate)
                                    edges[poss_factor_name][ph_name] = thisl
                        else:
                            for ii in range(len(phasings[0])):
                                for jj in range(len(phasings[1])):
                                    ph_name = '-'.join([phas_2_str(phasings[0][ii]), phas_2_str(phasings[1][jj])])
                                    phasing_np3[:, np.array(pp_nodes[0]) - minp] = phasings[0][ii]
                                    phasing_np3[:, np.array(pp_nodes[1]) - minp] = phasings[1][jj]
                                    thisl = compute_likelihood(observed, phasing_np3, self.error_rate)
                                    edges[poss_factor_name][ph_name] += thisl
        
            for ed in list(quotientg.edges()):
                edge_name = '-'.join(ed)
                if edge_name not in edges.keys():
                    print('not found', edge_name)
                    edges[edge_name] = {}
                    pp_nodes_str = [elem.split('-') for elem in ed]
                    pp_nodes = [(int(elem[0]), int(elem[1])) for elem in pp_nodes_str]
                    positions_genotypes = [''.join([str(q) for q in inp_hand.get_positions_genotypes(po)]) for po in
                                           pp_nodes]
                    phasings = [generate_phasings_ploidy(self.ploidy, pg) for pg in positions_genotypes]
                    for ii in range(len(phasings[0])):
                        for jj in range(len(phasings[1])):
                            ph_name = '-'.join([phas_2_str(phasings[0][ii]), phas_2_str(phasings[1][jj])])
                            edges[edge_name][ph_name] = self.epsilon
    
        return nodes, edges

    def construct_factor_graph(self, g_nodes, g_edges):
        
        factor_graph = FactorGraph()
        for node in g_nodes.keys():
            # print(node)
            # phasings = global_phasings[2]['11']
            # add variable
            factor_graph.add_node(str(node))
        
            #  potentials for node
            node_likelihood = [val if val != 0 else self.epsilon for val in g_nodes[node].values()]
            disf = DiscreteFactor(variables=[str(node)], cardinality=[len(node_likelihood)], values=node_likelihood)
        
            # Add node potentials to graph
            factor_graph.add_factors(disf)
            factor_graph.add_edges_from([(str(node), disf)])
    
        for edge in g_edges.keys():
            # print(edge)
            node1 = '-'.join(edge.split('-')[0:2])
            node2 = '-'.join(edge.split('-')[2:])
            # print(node1, node2)
            # factor_node = (u, v)
            edge_likelihood = [val if val != 0 else self.epsilon for val in g_edges[edge].values()]
            edge_potential = DiscreteFactor(variables=[node1, node2],
                                            cardinality=[len(g_nodes[node1]), len(g_nodes[node2])],
                                            values=edge_likelihood)
        
            # Add edge potentials to graph
            factor_graph.add_factors(edge_potential)
        
            # Add factor nodes and connect them to variable nodes
            factor_graph.add_edges_from([(node1, edge_potential), (node2, edge_potential)])
    
        for factor in factor_graph.factors:
            # print(factor)
            # print('=================================')
            for variable in factor.scope():
                assert (variable, factor) in factor_graph.edges(), f"Missing edge for {factor} and {variable}"
    
        print('No missing edges between factors and variables.')
    
        for edge in factor_graph.edges():
            if isinstance(edge[1], DiscreteFactor):
                assert edge[1] in factor_graph.factors, f"Missing factor for edge {edge}"
        print('No missing factor for edges.')
    
        return factor_graph



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

# Example usage
seq1 = "ABCBDAB"
seq2 = "BDCAB"
print(longest_common_substring(seq1, seq2))