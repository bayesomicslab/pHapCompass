import numpy as np
from pgmpy.models import FactorGraph
from pgmpy.factors.discrete import DiscreteFactor
from algorithm.haplotype_assembly_helper import compute_likelihood_generalized_plus
from utils.utils import get_matching_reads_for_positions, str_2_phas


class Factorgraph:
    def __init__(self, ploidy, error_rate, epsilon):
        self.ploidy = ploidy
        self.error_rate = error_rate
        self.epsilon = epsilon
       
    def construct(self, chordal_g, fragment_list):
        chordal_nodes = chordal_g.nodes(data=True)
        chordal_edges = chordal_g.edges(data=True)
        
        factor_graph = FactorGraph()
        for node, attr in chordal_nodes:
            print(node, attr)
            factor_graph.add_node(str(node))
            node_likelihoods = [val if val != 0 else self.epsilon for val in attr['weight'].values()]
            disf = DiscreteFactor(variables=[str(node)], cardinality=[len(node_likelihoods)], values=node_likelihoods)
            factor_graph.add_factors(disf)
            factor_graph.add_edges_from([(str(node), disf)])
        
        for u, v, attr in chordal_edges:
            # shared_sites = sorted([i for i in u.split('-') if i in v.split('-')])
            # shares = [(u.split('-').index(sha), v.split('-').index(sha), sha) for sha in shared_sites]
            all_sites = sorted(list(set(u.split('-') + v.split('-'))))
            matches = get_matching_reads_for_positions([int(p) for p in all_sites], fragment_list)
            
            u_sites = [(i, elem) for i, elem in enumerate(u.split('-'))]
            v_sites = [(i + len(u_sites), elem) for i, elem in enumerate(v.split('-'))]
            redundant_sites = [elem[0] for elem in v_sites if elem[1] in u.split('-')]
            non_redundant_sites = [col for col in range(len(u_sites) + len(v_sites)) if col not in redundant_sites]
            
            # edge_weights = list(chordal_nodes[u]['weight'].values()) + list(chordal_nodes[v]['weight'].values())
            
            edge_weights = []
            for phas1 in chordal_nodes[u]['weight'].keys():
                for phas2 in chordal_nodes[v]['weight'].keys():
                    matx1 = str_2_phas([phas1], self.ploidy)[0]
                    matx2 = str_2_phas([phas2], self.ploidy)[0]
                    matx = np.hstack([matx1, matx2])
                    nr_matx = matx[:, non_redundant_sites]  # this_phase
                    
                    this_marg = 0
                    for obs_indc, this_po, obs in matches:
                        phase_indc = [i for i, elem in enumerate(all_sites) if int(elem) in this_po]
                        this_marg += compute_likelihood_generalized_plus(np.array(obs), nr_matx, obs_indc, phase_indc,
                                                                         self.error_rate)
                    edge_weights.append(this_marg)
            
            edge_likelihoods = [val if val != 0 else self.epsilon for val in edge_weights]
            edge_potential = DiscreteFactor(variables=[u, v],
                                            cardinality=[len(chordal_nodes[u]['weight']), len(chordal_nodes[v]['weight'])],
                                            values=edge_likelihoods)
            # Add edge potentials to graph
            factor_graph.add_factors(edge_potential)
            
            # Add factor nodes and connect them to variable nodes
            factor_graph.add_edges_from([(u, edge_potential), (v, edge_potential)])
        
        for factor in factor_graph.factors:
            # print(factor)
            # print('=================================')
            for variable in factor.scope():
                assert (variable, factor) in factor_graph.edges(), f"Missing edge for {factor} and {variable}"
        
        print('All edges are included.')
        
        for edge in factor_graph.edges():
            if isinstance(edge[1], DiscreteFactor):
                assert edge[1] in factor_graph.factors, f"Missing factor for edge {edge}"
        print('All factors are included.')
        
        return factor_graph
