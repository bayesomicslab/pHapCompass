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
import pickle
from multiprocessing import Pool


def generate_quotient_graph_make_input():
    inputs = []
    contigs = ['Contig1_k3']
    simulated_data_path = '/home/mok23003/BML/HaplOrbit/simulated_data'
    graph_path = '/home/mok23003/BML/HaplOrbit/simulated_data_graphs/'
    for cont in contigs:
        fragment_files_path = os.path.join(simulated_data_path, cont)
        coverges = sorted([d for d in os.listdir(fragment_files_path) if os.path.isdir(os.path.join(fragment_files_path, d))])
        for coverage in coverges:
            this_frag_path = os.path.join(fragment_files_path, coverage)
            this_fragment_coverage_path = os.path.join(graph_path, 'fragment_graphs', cont, coverage)
            this_quotient_coverage_path = os.path.join(graph_path, 'quotient_graphs', cont, coverage)
            this_reverse_maps_path = os.path.join(graph_path, 'reverse_maps', cont, coverage)

            if not os.path.exists(this_quotient_coverage_path):
                os.makedirs(this_quotient_coverage_path)
            if not os.path.exists(this_fragment_coverage_path):
                os.makedirs(this_fragment_coverage_path)
            if not os.path.exists(this_reverse_maps_path):
                os.makedirs(this_reverse_maps_path)

            frag_files = [f for f in os.listdir(this_frag_path) if '.frag.txt' in f]
        
            for frag_file in frag_files:
                # frag_file_path = os.path.join(this_frag_path, frag_file)
                inputs.append([this_frag_path, this_fragment_coverage_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file])

    return inputs


def generate_quotient_graph(inp):
    this_frag_path, this_fragment_coverage_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file = inp
    file_name = frag_file.split('.')[0]
    print('Working on', os.path.join(this_frag_path, frag_file))
    fragment_v_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_v_label_' + file_name + '.pkl')
    fragment_e_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_e_label_' + file_name + '.pkl')
    quotient_v_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + file_name + '.pkl')
    quotient_e_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + file_name + '.pkl')


    class Args:
        def __init__(self):
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join(this_frag_path, frag_file)
            self.bam_path = 'example/example.bam'
            self.genotype_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/real_haps_contig1_k3.txt'
            self.ploidy = 3
            self.error_rate = 0.001
            self.epsilon = 0.0001
            self.output_path = 'output'
            self.root_dir = 'D:/UCONN/HaplOrbit'
            self.alleles = [0, 1]


    args = Args()
    input_handler = InputHandler(args)
    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)
    
    # create fragment graph
    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
    fragment_model.construct2(input_handler, config)

    # save fragment graph
    frag_graph_path = os.path.join(this_fragment_coverage_path, file_name + '.gt.gz')
    fragment_model.graph.save(frag_graph_path)

    with open(fragment_v_label_revered_path, "wb") as f:
        pickle.dump(fragment_model.v_label_reversed, f)

    edges_map_fragment = {}
    for k in fragment_model.e_label_reversed.keys():
        edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]

    with open(fragment_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_fragment, f)


    # create quotient graph
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct3(input_handler, config)

    # save quotient graph
    quot_graph_path = os.path.join(this_quotient_coverage_path, file_name + '.gt.gz')
    quotient_g.graph.save(quot_graph_path)

    with open(quotient_v_label_revered_path, "wb") as f:
        pickle.dump(quotient_g.v_label_reversed, f)

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    with open(quotient_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_quotient, f)

    print('[Done]', os.path.join(this_frag_path, frag_file))


def generate_graphs():
    inputs = generate_quotient_graph_make_input()
    print('Number of inputs:', len(inputs))
    pool = Pool(10)
    pool.map(generate_quotient_graph, inputs)


if __name__ == '__main__':
    generate_graphs()
