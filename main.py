import argparse
from data.input_handler import InputHandler
from data.configuration import Configuration
from algorithm.haplotype_assembly import HaplotypeAssembly
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from models.factor_graph import Factorgraph
from algorithm.chordal_contraction import chordal_contraction
from utils.utils import plot_graph
from algorithm.inference import *

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Haplotype Assembly and Phasing")
    parser.add_argument("-d", "--data_path", type=str, required=True, help="Path to the input data",
                        default='example/test.txt')
    parser.add_argument("-p", "--ploidy", type=int, required=True, help="Ploidy of the organism", default=3)
    parser.add_argument("-g", "--genotype_path", type=str, required=True, help="Path to the genotype data",
                        default='example/genotype.txt')
    parser.add_argument("--error_rate", type=float, required=True, help="Error rate", default=0.001)
    parser.add_argument("--epsilon", type=float, required=True, help="epsilon", default=0.0001)
    parser.add_argument("-a", "--alleles", required=False, nargs='*', help="List of alleles (optional)", default={0, 1})
    # parser.add_argument("--epsilon", help="epsilon in computing prob.")

    args = parser.parse_args()

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args.data_path, args.genotype_path, args.ploidy, args.alleles)
    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)
    fragment_model = FragmentGraph(args.data_path, args.genotype_path, args.ploidy, input_handler.alleles)
    frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)

    plot_graph(frag_graph)
    
    quotient_g = QuotientGraph(frag_graph).construct(fragment_list, input_handler, config)
    plot_graph(quotient_g)
    
    qg = chordal_contraction(quotient_g, fragment_list, input_handler, config)
    plot_graph(qg)

    factor_graph = Factorgraph(args.ploidy, args.error_rate, args.epsilon).construct(qg, fragment_list)

    beliefs = factor_graph_inference(factor_graph)
    source, target = '1-2', '4-7-8'
    samples = sample_gibbs(source, target, beliefs, qg)

    from utils.utils import str_2_phas
    row = 0
    columns = list(samples.columns.values)

    # shared_points =


    for i in range(len(columns) - 1):
        u, v = columns[i], columns[i + 1]
        u_phase = str_2_phas([samples.loc[row, u]], 3)[0]
        v_phase = str_2_phas([samples.loc[row, v]], 3)[0]
        # u_sites = [(i, elem) for i, elem in enumerate(u.split('-'))]
        # v_sites = [(i, elem) for i, elem in enumerate(v.split('-'))]
        u_sites = [(i, elem) for i, elem in enumerate(u.split('-'))]
        v_sites = [(i + len(u_sites), elem) for i, elem in enumerate(v.split('-'))]
    
        redundant_sites = [elem[0] for elem in v_sites if elem[1] in u.split('-')]
    
        non_redundant_sites = [col for col in range(len(u_sites) + len(v_sites)) if col not in redundant_sites]
        u_shared_site = [elem[0] for elem in u_sites if int(elem[1]) in redundant_sites]
        v_shared_site = [elem[0] - len(u_sites) for elem in v_sites if int(elem[1]) in redundant_sites]
    
        for u_idx in u_shared_site:
            u_phase = u_phase[u_phase[:, u_idx].argsort()]
    
        for v_idx in v_shared_site:
            v_phase = v_phase[v_phase[:, v_idx].argsort()]
    
        matx = np.hstack([u_phase, v_phase])
    
        nr_matx = matx[:, non_redundant_sites]


if __name__ == "__main__":
    main()
