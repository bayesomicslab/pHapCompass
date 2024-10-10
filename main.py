import os
import argparse
import sys
import networkx as nx
from data.input_handler import InputHandler
from data.configuration import Configuration
from algorithm.haplotype_assembly import HaplotypeAssembly
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from models.factor_graph import Factorgraph
from algorithm.chordal_contraction import chordal_contraction_cycle_base, chordal_contraction
from utils.utils import *
from algorithm.inference import *
from test.FFBS import generate_hmm_with_weights_and_emissions




def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Haplotype Assembly and Phasing")
    parser.add_argument("-d", "--data_path", type=str, required=False, help="Path to the input data", default=None)
    parser.add_argument("-p", "--ploidy", type=int, required=True, help="Ploidy of the organism", default=3)
    parser.add_argument("-g", "--genotype_path", type=str, required=True, help="Path to the genotype data",
                        default='example/genotype.txt')
    parser.add_argument("-a", "--alleles", required=False, nargs='*', help="List of alleles (optional)", default={0, 1})
    parser.add_argument("--error_rate", type=float, required=True, help="Error rate", default=0.001)
    parser.add_argument("--epsilon", type=float, required=True, help="epsilon", default=0.0001)
    parser.add_argument("-v", "--vcf_path", required=True, help="VCF file for called variants (string)")
    parser.add_argument("-b", "--bam_path", help="sam or bam file input (string)")
    parser.add_argument("-o", "--output_path", required=True, help="output path")
    parser.add_argument("-r", "--root_dir", required=True, help="root directory")
    
    # parser.add_argument("--epsilon", help="epsilon in computing prob.")

    args = parser.parse_args()


    class Args:
        def __init__(self):
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = '/home/mok23003/BML/HaplOrbit/example/Contig1_k3/c2/ART_0.frag.txt'
            self.bam_path = 'example/example.bam'
            self.genotype_path = '/home/mok23003/BML/HaplOrbit/example/Contig1_k3/real_haps_contig1_k3.txt'
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
    fragment_model.construct2(input_handler, config)
    
    # fragment_model = FragmentGraph(args.data_path, args.genotype_path, args.ploidy, input_handler.alleles)
    # frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)

    # plot_graph(frag_graph)
    print('Fragment Graph constructed.')

    # quotient_g = QuotientGraph(fragment_model.graph).construct(fragment_model.fragment_list, input_handler, config)
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct3(input_handler, config)
    # plot_graph(quotient_g)
    print('Quotient Graph constructed.')


    qg = chordal_contraction_cycle_base(quotient_g, fragment_list, input_handler, config)

    # qg = chordal_contraction(quotient_g, fragment_list, input_handler, config)
    plot_graph(qg)    
    print('Chordal Graph constructed.')


    error_rate = 0.001
    # state_names, transition_matrix, emission_prob_matrix, emission_index_map = generate_hmm_with_weights_and_emissions2(qg, error_rate)


    # for nn in qg.nodes(data=True):
    #     print(nn)
    # Save the graph in graph-tool binary format (preserves all properties)
#     quotient_g.graph.save("/home/mok23003/BML/HaplOrbit/example/Contig1_k3/graphs/c2/ART_0_quotient.gt")

#     # Save the graph in GraphML format (to preserve properties for external tools)
#     quotient_g.graph.save("/home/mok23003/BML/HaplOrbit/example/Contig1_k3/graphs/c2/ART_0_quotient.graphml", fmt="xml")
    
#     v_label_loaded = g_loaded.vertex_properties["v_label"]  # Corrected name for vertex label
#     e_weights_loaded = g_loaded.edge_properties["e_weights"]  # Corrected name for edge weights
#     e_label_loaded = g_loaded.edge_properties["e_label"]  # Corrected name for edge label

#     # Print the vertex labels
#     print("Vertex Labels:")
#     for v in g_loaded.vertices():
#         print(f"Vertex {int(v)}: Label = {v_label_loaded[v]}")

#     # Print the edge labels and weights
#     print("\nEdge Labels and Weights:")
#     for e in g_loaded.edges():
#         print(f"Edge {int(e.source())} -> {int(e.target())}: Label = {e_label_loaded[e]}, Weight = {e_weights_loaded[e]}")



def generate_quotient_graph(fragment_model, input_handler, config):
    fragment_files_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3'
    fragment_graphs_path = os.path.join(fragment_files_path, 'fragment_graphs')
    quotient_graph_path = os.path.join(fragment_files_path, 'quotient_graphs')


    import graph_tool.all as gt
    quotient_g.graph.save("/home/mok23003/BML/HaplOrbit/graphs/c2/ART_0.gt.gz")
    g_loaded = gt.load_graph("/home/mok23003/BML/HaplOrbit/graphs/c2/ART_0.gt.gz")
    
    v_label_loaded = g_loaded.vertex_properties["v_label"]  
    v_weights_loaded = g_loaded.vertex_properties["v_weights"]  

    e_label_loaded = g_loaded.edge_properties["e_label"]
    e_weight_loaded = g_loaded.edge_properties["e_weight"]



    # Access the correct vertex and edge properties
    v_label_loaded = g_loaded.vertex_properties["v_label"]   # Corrected from "e_label" to "v_label"
    v_weights_loaded = g_loaded.vertex_properties["v_weights"]  # Corrected from "e_weights" to "v_weights"

    e_label_loaded = g_loaded.edge_properties["e_label"]  # Edge property, correctly loaded from edge_properties
    e_weight_loaded = g_loaded.edge_properties["e_weights"]  # Edge property, correctly loaded from edge_properties

    # Print the vertex labels and weights
    print("Vertex Labels:")
    for v in g_loaded.vertices():
        print(f"Vertex {int(v)}: Label = {v_label_loaded[v]} | Weight = {v_weights_loaded[v]}")

    # Print the edge labels and weights
    print("\nEdge Labels and Weights:")
    for e in g_loaded.edges():
        print(f"Edge {int(e.source())} -> {int(e.target())}: Name = {e_label_loaded[e]} | Weight = {e_weight_loaded[e]}")



with open("vertex_dict.pkl", "wb") as f:
    pickle.dump(vertex_dict, f)
