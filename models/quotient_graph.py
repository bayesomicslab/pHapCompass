import networkx as nx
import itertools
import numpy as np
from utils.utils import phas_2_str, get_matching_reads_for_positions, convert_to_int_list, sort_strings
from algorithm.haplotype_assembly_helper import generate_phasings_ploidy_long, compute_likelihood_generalized_plus
from scipy.stats import entropy
import graph_tool.all as gt



class QuotientGraph:
    def __init__(self, frag_graph):
        self.fragment_graph = frag_graph
        self.graph = gt.Graph(directed=False)
        self.v_label = self.graph.new_vertex_property("string")
        self.v_weights = self.graph.new_vertex_property("object")
        self.e_label = self.graph.new_edge_property("string")
        self.e_weights = self.graph.new_edge_property("object")
        self.v_label_reversed = {}
        self.e_label_reversed = {}

    def construct(self, fragment_list, inpt_handler, config):
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
                    # poss_genotype = ''.join(['1' for i in poss])
                    poss_genotype = inpt_handler.get_genotype_positions([int(p) for p in poss])
                    all_phasings = generate_phasings_ploidy_long(config.ploidy, poss_genotype, allel_set=[0, 1])
                    matches = get_matching_reads_for_positions([int(i) for i in poss], fragment_list)
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
                    output_graph.add_edge(ed[0], ed[1], weight=weights, entropy=entr)
        return output_graph


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


    def construct3(self, inpt_handler, config):
        # output_graph = nx.Graph()
        # for u, v, attrs in self.fragment_graph.edges(data=True):
        # self.fragment_graph.e_label_reversed
        for edge_name in self.fragment_graph.e_label_reversed.keys():
            # print(u, v, attrs)
            # node_name = f"{u}-{v}", weight=attrs['weight'])
            edge_id = self.fragment_graph.e_label_reversed[edge_name]
            wei = self.fragment_graph.e_weights[edge_id]['weight']
            v1 = self.graph.add_vertex()
            self.v_label[v1] = edge_name
            self.v_label_reversed[edge_name] = int(v1)
            # self.current_vertex_index = self.current_vertex_index + 1
            self.current_vertex_index = int(v1) + 1
            self.v_weights[v1] = {'weight': wei}
            

        # plot_graph(output_graph)
        for node_id in self.fragment_graph.graph.vertices():
            
            node_label = self.fragment_graph.v_label[node_id]
            # print(node_id, node_label)
            
            node_list = []
            for nbr in self.fragment_graph.graph.iter_all_neighbors(node_id):
                # print(quotient_g.fragment_graph.v_label[nbr])
                nbr_name = self.fragment_graph.v_label[nbr]
                sorted_node_name = '-'.join([str(nnn) for nnn in sorted([int(node_label), int(nbr_name)])])
                # print(sorted_node_name)
                node_list.append(sorted_node_name)


            # for nbr, datadict in self.fragment_graph.adj[node].items():
            #     # print('-'.join(sorted([str(node), str(nbr)])))
            #     node_list.append('-'.join(sorted([str(node), str(nbr)])))
            edges_to_add = list(itertools.combinations(node_list, 2))
            for ed in edges_to_add:
                # print(node_id, ed)
                poss = sorted(list(set(ed[0].split('-') + ed[1].split('-'))))
                # edge_label = '-'. join([str(elem) for elem in sorted([int(p) for p in poss])])
                sorted_labels = sort_strings([ed[0], ed[1]])
                edge_label = '--'.join(sorted_labels)


                if ed not in self.e_label_reversed.keys():
                # if ed not in output_graph.edges():
                    # poss = sorted(list(set(ed[0].split('-') + ed[1].split('-'))))
                    # poss_genotype = ''.join(['1' for i in poss])
                    poss_genotype = inpt_handler.get_genotype_positions([int(p) for p in poss])
                    all_phasings = generate_phasings_ploidy_long(config.ploidy, poss_genotype, allel_set=[0, 1])
                    matches = get_matching_reads_for_positions([int(i) for i in poss], self.fragment_graph.fragment_list)
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
                    # output_graph.add_edge(ed[0], ed[1], weight=weights, entropy=entr)
                    sorted_ed = sorted(ed, key=convert_to_int_list)

                    v1 = self.v_label_reversed[sorted_ed[0]]
                    v2 = self.v_label_reversed[sorted_ed[1]]
                    edge = self.graph.add_edge(v1, v2)

                    # Assign a label to the edge
                    self.e_label[edge] = edge_label
                    self.e_label_reversed[edge_label] = edge
                    # Assign weights to the edge (as a dictionary)
                    self.e_weights[edge] = {"weight": weights, "entropy": entr}

        self.graph.ep['e_weights'] = self.e_weights
        self.graph.ep['e_label'] = self.e_label
        self.graph.vp['v_weights'] = self.v_weights
        self.graph.vp['v_label'] = self.v_label
