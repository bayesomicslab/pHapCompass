import argparse
import networkx as nx
from data.input_handler import InputHandler
from data.configuration import Configuration
from algorithm.haplotype_assembly import HaplotypeAssembly
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from models.factor_graph import Factorgraph
from algorithm.chordal_contraction import chordal_contraction_networkit, chordal_contraction
from utils.utils import *
from algorithm.inference import *


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

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)
    
    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
    frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)
    
    # fragment_model = FragmentGraph(args.data_path, args.genotype_path, args.ploidy, input_handler.alleles)
    # frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)

    # plot_graph(frag_graph)
    print('Fragment Graph constructed.')

    quotient_g = QuotientGraph(frag_graph).construct(fragment_list, input_handler, config)
    # plot_graph(quotient_g)

    # quotient_g.nodes(data=True)
    # quotient_g.edges(data=True)

    #
    #
    # interpreter = '/home/FCAM/mhosseini/anaconda3/envs/t2t/bin/python3'
    # import networkx as nx
    # import networkit as nk
    # nodes = ['1-2', '1-3', '2-3', '3-4', '4-5', '4-7', '5-6', '6-8', '7-8']
    # edges = [('1-2', '1-3'), ('1-2', '2-3'), ('1-3', '2-3'), ('1-3', '3-4'), ('2-3', '3-4'), ('3-4', '4-5'), ('3-4', '4-7'), ('4-5', '4-7'), ('4-5', '5-6'), ('4-7', '7-8'), ('5-6', '6-8'), ('6-8', '7-8')]
    # graphnx = nx.Graph()
    # graphnx.add_nodes_from(nodes)
    # graphnx.add_edges_from(edges)
    #
    # graphnk, reverse_map = nx2nk(graphnx)
    # cliques = networkit_find_cliques(graphnk)
    # cycles_g4 = [cyc for cyc in list(simple_cycles(graphnk, 0)) if len(cyc) > 4]
    # cycles_g4_unique = [list(x) for x in set(tuple(x) for x in cycles_g4)]
    #
    # chordless_cycles = list(nx.chordless_cycles(tempnx))
    #
    #
    #
    #
    # nodes_dict = dict(quotient_g.nodes(data=True))
    # edges_list = list(quotient_g.edges(data=True))
    # edges_dict = {str(item[0:2]): item[2] for item in edges_list}

    # nodes = list(quotient_g.nodes())
    # edges = list(quotient_g.edges())
    # new_g = nx.Graph()
    # new_g.add_nodes_from(nodes)
    # new_g.add_edges_from(edges)
    # graphnk, reverse_map = nx2nk(new_g)
    #

    # edges_w_data = list(quotient_g.edges(data=True))
    # dict(quotient_g.edges(data=True))

    print('Quotient Graph constructed.')
    
    qg = chordal_contraction(quotient_g, fragment_list, input_handler, config)
    print('Chordal Graph constructed.')
    # plot_graph(qg)
    # qg = quotient_g.copy()
    factor_graph = Factorgraph(config.ploidy, config.error_rate, config.epsilon).construct(qg, fragment_list)

    beliefs = factor_graph_inference(factor_graph)

    marginals, max_phasings = give_marginals(factor_graph, qg, beliefs)

    max_phase, positions = query_paths_gibbs_max(fragment_list, qg, beliefs, n_samples=1000)
    h_df = creat_vcf(max_phase, positions, config)
    print(h_df)
    
    # va_inference = VariableElimination(factor_graph)
    # result = va_inference.query(variables=['1-2'])
    

if __name__ == "__main__":
    main()
