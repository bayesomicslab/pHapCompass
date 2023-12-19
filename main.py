import argparse
from data.input_handler import InputHandler
from data.configuration import Configuration
from algorithm.haplotype_assembly import HaplotypeAssembly
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from models.factor_graph import Factorgraph
from algorithm.chordal_contraction import chordal_contraction

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Haplotype Assembly and Phasing")
    parser.add_argument("--data_path", type=str, required=True, help="Path to the input data")
    parser.add_argument("--ploidy", type=int, required=True, help="Ploidy of the organism")
    parser.add_argument("--genotype_path", type=str, required=True, help="Path to the genotype data")
    parser.add_argument("--error_rate", type=float, required=True, help="Error rate")
    parser.add_argument("--epsilon", type=float, required=True, help="epsilon")
    parser.add_argument("--alleles", nargs='*', help="List of alleles (optional)")
    # parser.add_argument("--epsilon", help="epsilon in computing prob.")

    args = parser.parse_args()

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args.data_path, args.genotype_path, args.ploidy, args.alleles)
    config = Configuration(args.ploidy, args.error_rate, args.epilon, input_handler.alleles)
    fragment_model = FragmentGraph(args.data_path, args.genotype_path, args.ploidy, input_handler.alleles)
    frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)
    
    quotient_g = QuotientGraph(frag_graph).construct(fragment_list, args.ploidy, args.error_rate)
    qg = chordal_contraction(quotient_g, args.ploidy, fragment_list, args.error_rate)
    
    # todo: after this should be fixed
    fg = Factorgraph(args.ploidy, args.error_rate, epsilon=0.0001)
    fg_nodes, fg.edges = fg.get_weights(fragment_list, quotient_g, config, input_handler)
    factor_graph = fg.construct_factor_graph(fg_nodes, fg.edges)
    print(factor_graph.nodes())
    
    # data_manager = DataManager()
    haplotype_assembly = HaplotypeAssembly(ploidy=args.ploidy)

    # Perform initial computations or setup
    haplotype_assembly.initialize_assembly()

    # Further processing...


if __name__ == "__main__":
    main()
