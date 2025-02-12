
import os
import argparse
import sys
import networkx as nx
from data.input_handler import InputHandler
from data.configuration import Configuration
from algorithm.haplotype_assembly import HaplotypeAssembly
from models.fragment_graph import FragmentGraph
from utils.utils import *
from algorithm.inference import *
from algorithm.chordal_contraction import *
from multiprocessing import Pool
from algorithm.haplotype_assembly_helper import *

# @profile
def get_top_k_weights(mydict, k):
    # Extract the 'weight' dictionary
    weights = mydict['weight']
    
    # Sort the weights dictionary by value in descending order and take the top k items
    top_k_weights = dict(sorted(weights.items(), key=lambda item: item[1], reverse=True)[:k])
    
    # Return a new dictionary with the top k weights and the same entropy
    return {'weight': top_k_weights, 'entropy': mydict['entropy']}


# @profile
def chordal_contraction_graph_tool_approx(inp):
    save_path, subg_id, subg, config, fragment_model, k = inp
    this_path = os.path.join(save_path, 'chordal_sub_' + str(subg_id) + '.gt.gz')
    print('Working on', this_path)
    new_graph = subg.copy()
    e_weights = new_graph.edge_properties["e_weights"]
    # new_graph.clear_filters()
    e_entropy = new_graph.new_edge_property("double")

    # Loop over edges and assign entropy from the e_weights property
    for e in new_graph.edges():
        e_entropy[e] = e_weights[e]['entropy']

    new_graph.ep['e_entropy'] = e_entropy

    chordless_cycles = get_chordless_cycles(new_graph)

    to_be_removed_nodes = []

    for cyc_id, cyc in enumerate(chordless_cycles):
        
        # print(cyc_id)
        
        edges = [new_graph.edge(cyc[-1], cyc[0])]
        for i in range(len(cyc) - 1):
            edges += [new_graph.edge(cyc[i], cyc[i+1])]
        edges = [x for x in edges if x is not None]
        while len(edges) > 3:
            min_edge = min(edges, key=lambda e: new_graph.ep['e_entropy'][e])
            source_label = new_graph.vp['v_label'][min_edge.source()]
            target_label = new_graph.vp['v_label'][min_edge.target()]
            # new node positions
            poss = sorted(set([int(nn) for nn in source_label.split('-')] + [int(nn) for nn in target_label.split('-')]))
            # new vertex properties:
            new_vertex_name = '-'.join([str(nnn) for nnn in poss])
            vertex_weights = new_graph.ep['e_weights'][min_edge]
            vertex_weights_appr = get_top_k_weights(vertex_weights, k)

            new_graph.vertex_properties["v_weights"][min_edge.source()] = vertex_weights_appr
            new_graph.vertex_properties["v_label"][min_edge.source()] = new_vertex_name

            source_nbrs = [n for n in min_edge.source().all_neighbors() if n != min_edge.target()]
            target_nbrs = [n for n in min_edge.target().all_neighbors() if n != min_edge.source()]
            common_nbrs = set(source_nbrs).intersection(set(target_nbrs))

            for n in common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]
                # e_poss = sorted(set([int(nn) for nn in v_label.split('-')] + poss))
                # print(len(e_poss))
                # new_edge_name = '-'.join([str(nnn) for nnn in e_poss])
                sorted_labels = sort_nodes([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)
                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)
                final_weights_appr = get_top_k_weights(final_weight, k)

                e1 = new_graph.edge(min_edge.source(), n)
                e2 = new_graph.edge(min_edge.target(), n)
                
                new_graph.edge_properties["e_weights"][e1] = final_weights_appr
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weights_appr['entropy']
                new_graph.remove_edge(e2)

            for n in set(source_nbrs)-common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]

                sorted_labels = sort_nodes([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)
                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)
                final_weights_appr = get_top_k_weights(final_weight, k)


                e1 = new_graph.edge(min_edge.source(), n)
                # e2 = new_graph.edge(min_edge.target(), n)
                new_graph.edge_properties["e_weights"][e1] = final_weights_appr
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weights_appr['entropy']
                # new_graph.edge_properties["e_weights"][e2]

            
            for n in set(target_nbrs)-common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]
                sorted_labels = sort_nodes([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)

                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)
                final_weights_appr = get_top_k_weights(final_weight, k)

                e2 = new_graph.edge(min_edge.target(), n)
                new_graph.remove_edge(e2)
                e1 = new_graph.add_edge(min_edge.source(), n)
                new_graph.edge_properties["e_weights"][e1] = final_weights_appr
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weights_appr['entropy']
            
            # to_be_removed_nodes += [min_edge.target()]
            to_be_removed_nodes.append(min_edge.target())
            new_graph.remove_edge(min_edge)
            edges.remove(min_edge)

    new_graph.remove_vertex(to_be_removed_nodes)
    new_graph.save(this_path)
    print('[Done]', this_path)
    

class Args:
    def __init__(self):
        self.vcf_path = 'example/62_ID0.vcf'
        # self.data_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/c2/ART_0.frag.txt'
        self.data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data/Contig1_k3/c2/ART_0.frag.txt'
        self.bam_path = 'example/example.bam'
        # self.genotype_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/real_haps_contig1_k3.txt'
        self.genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data/Contig1_k3/real_haps_contig1_k3.txt'
        self.ploidy = 3
        self.error_rate = 0.001
        self.epsilon = 0.0001
        self.output_path = 'output'
        self.root_dir = 'D:/UCONN/HaplOrbit'
        self.alleles = [0, 1]

# Create the mock args object
args = Args()

# Initialize classes with parsed arguments
input_handler = InputHandler(args)

config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)
fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)

# frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)
# frag_graph, fragment_list = fragment_model.construct2(input_handler, config)
fragment_model.construct(input_handler, config)

# main_path = '/home/mok23003/BML/HaplOrbit/old_simulated_data_graphs'
main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/old_simulated_data_graphs'

contig = 'Contig1_k3'
coverage = 'c2'
frag_file = 'ART_0.frag.txt'
m = 100
k = 5
quotient_graph, v_label_reversed, e_label_reversed = read_quotient_graph(main_path, contig, coverage, frag_file)


save_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_graphs/test_chordal_graphs2'
if not os.path.exists(save_path):
    os.makedirs(save_path)

subgraphs = divide_graph_by_labels(quotient_graph, m)
subgraph_copies = [gt.Graph(subgraph, prune=True) for subgraph in subgraphs]

inputs = [[save_path, subg_id, subg, config, fragment_model, k] for subg_id, subg in enumerate(subgraph_copies)]
print('Number of subgraphs:', len(inputs))
pool = Pool(2)
pool.map(chordal_contraction_graph_tool_approx, inputs)
# chordal_contraction_graph_tool_approx(inputs[0])
# for inp in inputs:
#     chordal_contraction_graph_tool_approx(inp)
