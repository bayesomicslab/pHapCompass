import argparse

import pandas as pd

from data.input_handler import InputHandler
from data.configuration import Configuration
from algorithm.haplotype_assembly import HaplotypeAssembly
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from models.factor_graph import Factorgraph
from algorithm.chordal_contraction import chordal_contraction
from utils.utils import *
from algorithm.inference import *

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Haplotype Assembly and Phasing")
    parser.add_argument("-d", "--data_path", type=str, required=True, help="Path to the input data",
                        default='example/test.txt')
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
    
    fragment_model = FragmentGraph(args.data_path, args.genotype_path, args.ploidy, input_handler.alleles)
    frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)
    
    # fragment_model = FragmentGraph(args.data_path, args.genotype_path, args.ploidy, input_handler.alleles)
    # frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)

    plot_graph(frag_graph)
    print('Fragment Graph constructed.')

    quotient_g = QuotientGraph(frag_graph).construct(fragment_list, input_handler, config)
    plot_graph(quotient_g)
    print('Quotient Graph constructed.')
    
    qg = chordal_contraction(quotient_g, fragment_list, input_handler, config)
    print('Chordal Graph constructed.')
    plot_graph(qg)

    factor_graph = Factorgraph(args.ploidy, args.error_rate, args.epsilon).construct(qg, fragment_list)

    beliefs = factor_graph_inference(factor_graph)

    marginals, max_phasings = give_marginals(factor_graph, qg, beliefs)

    max_phase, positions = query_paths_gibbs_max(fragment_list, qg, beliefs, n_samples=1000)
    h_df = creat_vcf(max_phase, positions, config)

    va_inference = VariableElimination(factor_graph)
    result = va_inference.query(variables=['1-2'])
    
    


if __name__ == "__main__":
    main()
