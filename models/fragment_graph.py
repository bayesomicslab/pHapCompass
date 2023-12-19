import networkx as nx
import numpy as np

class FragmentGraph:
    def __init__(self, data_path, genotype_path, ploidy, alleles):
        self.ploidy = ploidy
        self.alleles = alleles
        self.data_path = data_path
        self.genotype_path = genotype_path
        

    def create_graph(self):
        # Logic to create the graph
        pass

    def construct_graph(self, inp_hand, conf):
        '''
        Loads in a fragment file created by extract-poly
        :param fragment_filename: The file path to the fragment file
        :return: A graph representation of the fragments
        '''
        # G = nk.Graph(vcf_df.shape[0], weighted=True, directed=False)
        graph = nx.Graph()
        fragment_list = []
        for fragment in open(self.data_path, 'r'):
            # print(fragment)
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
                self.update_fragment_graph(graph, positions, alleles, inp_hand, conf)
                fragment_list.append(positions)
                fragment_list.append(alleles)
    
        return graph, fragment_list
    
    def update_fragment_graph(self, graph, pos, alle, inp_hand, conf):
        n_allels = 2
        n_positions = 2
        
        genotypes = inp_hand.get_genotype_positions(pos)
        # conf.global_phasings
        # create nodes if they don't exist
        for i in range(len(pos)):
            # genotype is an array, genotypes (count of minor allele) by evidence
            if not graph.has_node(pos[i]):
                graph.add_node(pos[i], genotype=genotypes[i], weight=np.zeros(n_allels))
                graph.nodes[pos[i]]['weight'][alle[i]] += 1
            else:
                graph.nodes[pos[i]]['weight'][alle[i]] += 1
    
        # change node and edge potentials
        for i in range(len(pos)):
            # add evidence for edge potentials
            for j in range(i):
                str_genotypes = ''.join([str(gn) for gn in [genotypes[j], genotypes[i]]])
                out_comes = conf.global_likelihoods[n_positions][str_genotypes]
                if not graph.has_edge(pos[j], pos[i]):
                    # phasings,weights = generate_phasings(ploidy)
                    # all_phasings = generate_phasings(ploidy, genotype=str_genotypes)
                    # str_phasings = global_phasings[2][str_genotypes]
                    edge_weights = {outk: 0 for outk in out_comes.keys()}
                    graph.add_edge(pos[j], pos[i], weight=edge_weights)
                    
                    # update likelihood
                    observation = np.array([alle[j], alle[i]])
                    str_observation = ''.join([str(ob) for ob in observation])
                    # updating_weights = [k for k in out_comes.keys() if k[0:2] == str_observation]
                    # updating_weights = [k for k in out_comes.keys() if k[0:2] == str_observation]
                    # print('aaaaaaaaaaaa', len(updating_weights), updating_weights)
                    for wei in out_comes.keys():
                        # graph.get_edge_data(pos[j], pos[i])['weight'][wei] += \
                        # conf.global_likelihoods[n_positions][str_genotypes][wei]
                        graph.get_edge_data(pos[j], pos[i])['weight'][wei] += conf.global_likelihoods[n_positions][str_genotypes][wei][str_observation]
                        
                else:
                    observation = np.array([alle[j], alle[i]])
                    str_observation = ''.join([str(ob) for ob in observation])
                    # updating_weights = [k for k in out_comes.keys() if k[0:2] == str_observation]
                    for wei in out_comes.keys():
                        graph.get_edge_data(pos[j], pos[i])['weight'][wei] += \
                        conf.global_likelihoods[n_positions][str_genotypes][wei][str_observation]
                        
        # return G

    def add_node(self, node):
        # Logic to add a node
        pass

    def add_edge(self, edge):
        # Logic to add an edge
        pass

