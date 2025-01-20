import os
from data.input_handler import InputHandler
from data.configuration import Configuration
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from models.factor_graph import Factorgraph
from utils.utils import *
from algorithm.inference import *
from algorithm.chordal_contraction import *
import graph_tool.all as gt



def main():
    # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/test.frag'
    frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_100/ploidy_8/cov_50/frag/93.frag'
    ploidy= 8
    # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/haplotypes.csv'
    genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_100/ploidy_8/haplotypes.csv'

    class Args:
        def __init__(self):
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = frag_path
            self.bam_path = 'example/example.bam'
            self.genotype_path = genotype_path
            self.ploidy = ploidy
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

    fragment_model.construct(input_handler, config)
    print('Fragment Graph constructed.')

    fragment_list = fragment_model.fragment_list

    obs_positions = fragment_list[::2]
    obs_reads = fragment_list[1::2]