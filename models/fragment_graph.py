import networkx as nx
import numpy as np
import graph_tool.all as gt

class FragmentGraph:
    def __init__(self, data_path, genotype_path, ploidy, alleles):
        self.ploidy = ploidy
        self.alleles = alleles
        self.data_path = data_path
        self.genotype_path = genotype_path
        self.current_vertex_index = 0
        self.current_edge_index = 0
        self.graph = gt.Graph(directed=False)
        self.v_label = self.graph.new_vertex_property("string")
        self.v_weights = self.graph.new_vertex_property("object")
        self.e_label = self.graph.new_edge_property("string")
        self.e_weights = self.graph.new_edge_property("object")
        self.v_label_reversed = {}
        # self.v_weights_reversed = {}
        self.e_label_reversed = {}
        # self.e_weights_reversed = {}
        self.fragment_list = []

    def create_graph(self):
        # Logic to create the graph
        pass

    # def construct_graph(self, inp_hand, conf):
    #     '''
    #     Loads in a fragment file created by extract-poly
    #     :param fragment_filename: The file path to the fragment file
    #     :return: A graph representation of the fragments
    #     '''
    #     # G = nk.Graph(vcf_df.shape[0], weighted=True, directed=False)
    #     graph = nx.Graph()
    #     fragment_list = []
    #     # for fragment in open(self.data_path, 'r'):
    #     for fragment in open(self.data_path, 'r'):
    #         # print(fragment)
    #         # read index when there are pairs, readname, start pos (in variants), allele sequence, alleles
    #         parts = fragment.split()
    #         positions = []
    #         alleles = []
    #         for iii in range(int(parts[0])):
    #             # process i+2,i+3.... i+4,i+5...
    #             start_idx_of_read = iii * 2 + 3
    #             seq_len = len(parts[start_idx_of_read])
    #             positions.extend(
    #                 list(range(int(parts[start_idx_of_read - 1]), int(parts[start_idx_of_read - 1]) + seq_len)))
    #             [alleles.append(int(a)) for a in parts[start_idx_of_read]]
    #             self.update_fragment_graph(graph, positions, alleles, inp_hand, conf)
    #             fragment_list.append(positions)
    #             fragment_list.append(alleles)
    
    #     return graph, fragment_list
        
    # @profile
    def construct(self, inp_hand, conf):
        '''
        Loads in a fragment file created by extract-poly
        :param fragment_filename: The file path to the fragment file
        :return: A graph representation of the fragments
        '''
        # G = nk.Graph(vcf_df.shape[0], weighted=True, directed=False)
        # graph = gt.Graph(directed=False)
        # v_label = graph.new_vertex_property("string")
        # v_weights = graph.new_vertex_property("object")
        # e_label = graph.new_edge_property("string")
        # e_weights = graph.new_edge_property("object")
        # fragment_list = []
        # for fragment in open(self.data_path, 'r'):
        for fragment in open(self.data_path, 'r'):
            print(fragment)
            # read index when there are pairs, readname, start pos (in variants), allele sequence, alleles
            parts = fragment.split()
            positions = []
            alleles = []
            for iii in range(int(parts[0])):
                # process i+2,i+3.... i+4,i+5...
                start_idx_of_read = iii * 2 + 3
                seq_len = len(parts[start_idx_of_read])
                positions.extend(
                    list(range(int(parts[start_idx_of_read - 1]), int(parts[start_idx_of_read - 1]) + seq_len)))
                [alleles.append(int(a)) for a in parts[start_idx_of_read]]
                self.update_fragment_graph(positions, alleles, inp_hand, conf)
                self.fragment_list.append(positions)
                self.fragment_list.append(alleles)
                
        self.graph.ep['e_weights'] = self.e_weights
        self.graph.ep['e_label'] = self.e_label
        self.graph.vp['v_weights'] = self.v_weights
        self.graph.vp['v_label'] = self.v_label
        # print_graph(graph, v_label, v_weights, e_label, e_weights)
        # return graph, fragment_list


    # def update_fragment_graph(self, graph, pos, alle, inp_hand, conf):
    #     n_allels = 2
    #     n_positions = 2
        
    #     genotypes = inp_hand.get_genotype_positions(pos)
    #     # conf.global_phasings
    #     # create nodes if they don't exist
    #     for i in range(len(pos)):
    #         # genotype is an array, genotypes (count of minor allele) by evidence
    #         if not graph.has_node(pos[i]):
    #             graph.add_node(pos[i], genotype=genotypes[i], weight=np.zeros(n_allels))
    #             graph.nodes[pos[i]]['weight'][alle[i]] += 1
    #         else:
    #             graph.nodes[pos[i]]['weight'][alle[i]] += 1
    
    #     # change node and edge potentials
    #     for i in range(len(pos)):
    #         # add evidence for edge potentials
    #         for j in range(i):
    #             str_genotypes = ''.join([str(gn) for gn in [genotypes[j], genotypes[i]]])
    #             # print(conf.global_likelihoods[n_positions])
    #             out_comes = conf.global_likelihoods[n_positions][str_genotypes]
    #             if not graph.has_edge(pos[j], pos[i]):
    #                 # phasings,weights = generate_phasings(ploidy)
    #                 # all_phasings = generate_phasings(ploidy, genotype=str_genotypes)
    #                 # str_phasings = global_phasings[2][str_genotypes]
    #                 edge_weights = {outk: 0 for outk in out_comes.keys()}
    #                 graph.add_edge(pos[j], pos[i], weight=edge_weights)
                    
    #                 # update likelihood
    #                 observation = np.array([alle[j], alle[i]])
    #                 str_observation = ''.join([str(ob) for ob in observation])
    #                 # updating_weights = [k for k in out_comes.keys() if k[0:2] == str_observation]
    #                 # updating_weights = [k for k in out_comes.keys() if k[0:2] == str_observation]
    #                 # print('aaaaaaaaaaaa', len(updating_weights), updating_weights)
    #                 for wei in out_comes.keys():
    #                     # graph.get_edge_data(pos[j], pos[i])['weight'][wei] += \
    #                     # conf.global_likelihoods[n_positions][str_genotypes][wei]
    #                     graph.get_edge_data(pos[j], pos[i])['weight'][wei] += conf.global_likelihoods[n_positions][str_genotypes][wei][str_observation]
                        
    #             else:
    #                 observation = np.array([alle[j], alle[i]])
    #                 str_observation = ''.join([str(ob) for ob in observation])
    #                 # updating_weights = [k for k in out_comes.keys() if k[0:2] == str_observation]
    #                 for wei in out_comes.keys():
    #                     graph.get_edge_data(pos[j], pos[i])['weight'][wei] += \
    #                     conf.global_likelihoods[n_positions][str_genotypes][wei][str_observation]
                        
    #     # return G


    def update_fragment_graph(self, pos, alle, inp_hand, conf):
        n_allels = 2
        n_positions = 2
        
        genotypes = inp_hand.get_genotype_positions(pos)
        # conf.global_phasings
        # create nodes if they don't exist
        for i in range(len(pos)):
            # genotype is an array, genotypes (count of minor allele) by evidence
            # if not graph.has_node(pos[i]):
            # has_n, tv = graph_has_node(graph, str(pos[i]), v_label)
            node_name = str(pos[i])
            if node_name not in self.v_label_reversed.keys():
            # if not has_n:
                # graph.add_node(pos[i], genotype=genotypes[i], weight=np.zeros(n_allels))
                v1 = self.graph.add_vertex()
                self.v_label[v1] = node_name
                self.v_label_reversed[node_name] = int(v1)
                # self.current_vertex_index = self.current_vertex_index + 1
                self.current_vertex_index = int(v1) + 1
                self.v_weights[v1] = {}
                self.v_weights[v1]['weight'] = {al: 0 for al in conf.alleles}
                self.v_weights[v1]['weight'][alle[i]] += 1
                # graph.nodes[pos[i]]['weight'][alle[i]] += 1
            else:
                self.v_weights[self.v_label_reversed[node_name]]['weight'][alle[i]] += 1
                # graph.nodes[pos[i]]['weight'][alle[i]] += 1
        
        # change node and edge potentials
        for i in range(len(pos)):
            # add evidence for edge potentials
            for j in range(i):
                str_genotypes = ''.join([str(gn) for gn in [genotypes[j], genotypes[i]]])
                # print(conf.global_likelihoods[n_positions])
                out_comes = conf.global_likelihoods[n_positions][str_genotypes]
                edge_label = '-'.join([str(pos[j]), str(pos[i])])
                if edge_label not in self.e_label_reversed.keys():
                # has_e, te = graph_has_edge(graph, edge_label, e_label)
                # if not graph.has_edge(pos[j], pos[i]):
                # if not has_e:
                    # phasings,weights = generate_phasings(ploidy)
                    # all_phasings = generate_phasings(ploidy, genotype=str_genotypes)
                    # str_phasings = global_phasings[2][str_genotypes]
                    # v1 = find_vertex_by_label(graph, str(pos[j]), v_label)
                    # v2 = find_vertex_by_label(graph, str(pos[i]), v_label)
                    v1 = self.graph.vertex(self.v_label_reversed[str(pos[j])])
                    v2 = self.graph.vertex(self.v_label_reversed[str(pos[i])])
                    e1 = self.graph.add_edge(v1, v2)
                    self.e_label[e1] = edge_label
                    self.e_label_reversed[edge_label] = e1
                    self.current_edge_index = self.current_edge_index + 1
                    # graph.add_edge(pos[j], pos[i])
                    edge_weights = {outk: 0 for outk in out_comes.keys()}
                    # graph.add_edge(pos[j], pos[i], weight=edge_weights)
                    self.e_weights[e1] = {'weight': edge_weights}

                    # update likelihood
                    observation = np.array([alle[j], alle[i]])
                    str_observation = ''.join([str(ob) for ob in observation])
                    # updating_weights = [k for k in out_comes.keys() if k[0:2] == str_observation]
                    # updating_weights = [k for k in out_comes.keys() if k[0:2] == str_observation]
                    # print('aaaaaaaaaaaa', len(updating_weights), updating_weights)
                    for wei in out_comes.keys():
                        # graph.get_edge_data(pos[j], pos[i])['weight'][wei] += \
                        # conf.global_likelihoods[n_positions][str_genotypes][wei]
                        # graph.get_edge_data(pos[j], pos[i])['weight'][wei] += conf.global_likelihoods[n_positions][str_genotypes][wei][str_observation]
                        self.e_weights[e1]['weight'][wei] += conf.global_likelihoods[n_positions][str_genotypes][wei][str_observation]
                else:
                    te = self.e_label_reversed[edge_label]
                    observation = np.array([alle[j], alle[i]])
                    str_observation = ''.join([str(ob) for ob in observation])
                    # updating_weights = [k for k in out_comes.keys() if k[0:2] == str_observation]
                    for wei in out_comes.keys():
                        self.e_weights[te]['weight'][wei] += conf.global_likelihoods[n_positions][str_genotypes][wei][str_observation]

                        # graph.get_edge_data(pos[j], pos[i])['weight'][wei] += \
                        # conf.global_likelihoods[n_positions][str_genotypes][wei][str_observation]
                        
        # return G


    def add_node(self, node):
        # Logic to add a node
        pass

    def add_edge(self, edge):
        # Logic to add an edge
        pass

def graph_has_node(g, label, v_label):
    has_node = False
    target_v = None
    for v in g.vertices():
        if v_label[v] == label:
            has_node = True
            target_v = v
            break
    return has_node, target_v

def graph_has_edge(g, label, e_label):
    has_edge = False
    target_e = None
    for e in g.edges():
        if e_label[e] == label:
            has_edge = True
            target_e = e
            break
    return has_edge, target_e


def find_vertex_by_label(g, label, v_label):
    for v in g.vertices():
        if v_label[v] == label:
            return v
    return None

def print_graph(graph, v_label, v_weights, e_label, e_weights):
    print("Vertices:")
    for v in graph.vertices():
        print(f"Vertex {int(v)}: Label = {v_label[v]}, Weights = {v_weights[v]}")

    # Print all edges with their labels and weights
    print("\nEdges:")
    for e in graph.edges():
        print(f"Edge from vertex {int(e.source())} to vertex {int(e.target())}: Label = {e_label[e]}, Weights = {e_weights[e]}")