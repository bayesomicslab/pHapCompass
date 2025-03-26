import os
import sys
import pickle
from algorithm.FFBS import *
from algorithm.chordal_contraction import *
from data.input_handler import InputHandler
from data.configuration import Configuration
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from evaluation.evaluation import *
from utils.utils import *
import mrftools
from algorithm.SMBP import create_markov_net, update_entropies, update_transitions, update_vertices
from algorithm.viterbi import *


def make_inputs_for_generate_qoutient_graph(simulator):
    inputs = []
    simulator.contig_lens = [100]
    simulator.ploidies = [3, 4, 6, 8]
    for contig_len in simulator.contig_lens:
        for ploidy in simulator.ploidies:
            # stop
            genotype_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
            # genotype_df = pd.read_csv(genotype_path)
            for coverage in simulator.coverages:
                this_cov_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                frag_path = os.path.join(this_cov_path, 'frag')
                frag_graph_path = os.path.join(this_cov_path, 'fgraph')
                quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
                qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')

                if not os.path.exists(frag_graph_path):
                    os.makedirs(frag_graph_path)
                if not os.path.exists(quotient_graph_path):
                    os.makedirs(quotient_graph_path)
                if not os.path.exists(qgraph_reverse_maps_path):
                    os.makedirs(qgraph_reverse_maps_path)

                existing_files_qg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_e_label' in ff]
                existing_files_qg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_v_label' in ff]
                existing_files_fg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_e_label' in ff]
                existing_files_fg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_v_label' in ff]
                existing_fg = [ff for ff in os.listdir(frag_graph_path) if '.gt.gz' in ff]
                existing_qg = [ff for ff in os.listdir(quotient_graph_path) if '.gt.gz' in ff]

                for rd in range(simulator.n_samples):
                    if 'qg_e_label_' + str(rd).zfill(2) + '.pkl' not in existing_files_qg_e or 'qg_v_label_' + str(rd).zfill(2) + '.pkl' not in existing_files_qg_v or 'fg_e_label_' + str(rd).zfill(2) + '.pkl' not in existing_files_fg_e or 'fg_v_label_' + str(rd).zfill(2) + '.pkl' not in existing_files_fg_v or str(rd).zfill(2) + '.gt.gz' not in existing_fg or str(rd).zfill(2) + '.gt.gz' not in existing_qg:
                        inp = [frag_path, frag_graph_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path]
                        inputs.append(inp)
    return inputs


def make_inputs_for_run_count(simulator):
    simulator.contig_lens = [100]
    inputs = []
    for contig_len in simulator.contig_lens:
        for ploidy in simulator.ploidies:
            # stop
            genotype_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
            # genotype_df = pd.read_csv(genotype_path)
            for coverage in simulator.coverages:
                this_cov_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                frag_path = os.path.join(this_cov_path, 'frag')
                frag_graph_path = os.path.join(this_cov_path, 'fgraph')
                quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
                qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')
                results_path = os.path.join(this_cov_path, 'results_counts')

                if not os.path.exists(frag_graph_path):
                    os.makedirs(frag_graph_path)
                if not os.path.exists(quotient_graph_path):
                    os.makedirs(quotient_graph_path)
                if not os.path.exists(qgraph_reverse_maps_path):
                    os.makedirs(qgraph_reverse_maps_path)
                if not os.path.exists(results_path):
                    os.makedirs(results_path)

                # existing_files_qg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_e_label' in ff]
                # existing_files_qg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_v_label' in ff]
                # existing_files_fg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_e_label' in ff]
                # existing_files_fg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_v_label' in ff]
                # existing_fg = [ff for ff in os.listdir(frag_graph_path) if '.gt.gz' in ff]
                # existing_qg = [ff for ff in os.listdir(quotient_graph_path) if '.gt.gz' in ff]
                existing_results = [ff for ff in os.listdir(results_path) if 'FFBS' in ff]
                # existing_results = []
                for rd in range(simulator.n_samples):
                    if 'FFBS_{}.pkl'.format(str(rd).zfill(2)) not in existing_results and \
                        '{}.gt.gz'.format(str(rd).zfill(2)) in os.listdir(quotient_graph_path) and \
                        'qg_e_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path) and \
                        'qg_v_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path):
                        inp = [frag_path, frag_graph_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path, results_path]
                        # inp = [frag_path, frag_graph_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path]

                        inputs.append(inp)
        inputs = sorted(inputs, key=lambda x: x[0])
    return inputs


# def make_inputs_for_run_likelihood(simulator):
#     simulator.contig_lens = [100]
#     simulator.ploidies = [3, 4, 6]
#     inputs = []
#     for contig_len in simulator.contig_lens:
#         for ploidy in simulator.ploidies:
#             # stop
#             genotype_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
#             # genotype_df = pd.read_csv(genotype_path)
#             for coverage in simulator.coverages:
#                 this_cov_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
#                 frag_path = os.path.join(this_cov_path, 'frag')
#                 frag_graph_path = os.path.join(this_cov_path, 'fgraph')
#                 quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
#                 qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')
#                 results_path = os.path.join(this_cov_path, 'results_likelihood')

#                 if not os.path.exists(frag_graph_path):
#                     os.makedirs(frag_graph_path)
#                 if not os.path.exists(quotient_graph_path):
#                     os.makedirs(quotient_graph_path)
#                 if not os.path.exists(qgraph_reverse_maps_path):
#                     os.makedirs(qgraph_reverse_maps_path)
#                 if not os.path.exists(results_path):
#                     os.makedirs(results_path)

#                 # existing_files_qg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_e_label' in ff]
#                 # existing_files_qg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_v_label' in ff]
#                 # existing_files_fg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_e_label' in ff]
#                 # existing_files_fg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_v_label' in ff]
#                 # existing_fg = [ff for ff in os.listdir(frag_graph_path) if '.gt.gz' in ff]
#                 # existing_qg = [ff for ff in os.listdir(quotient_graph_path) if '.gt.gz' in ff]
#                 # existing_results = [ff for ff in os.listdir(results_path) if 'FFBS' in ff]
#                 existing_results = []
#                 # for rd in range(90, 100):
#                 for rd in range(simulator.n_samples):
#                     if 'FFBS_{}.pkl'.format(str(rd).zfill(2)) not in existing_results and \
#                         '{}.gt.gz'.format(str(rd).zfill(2)) in os.listdir(quotient_graph_path) and \
#                         'qg_e_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path) and \
#                         'qg_v_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path):
#                         print(frag_path)
#                         inp = [frag_path, frag_graph_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path, results_path]
#                         # inp = [frag_path, frag_graph_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path]

#                         inputs.append(inp)
#         inputs = sorted(inputs, key=lambda x: x[0], reverse=True)
#     return inputs


def make_inputs_for_running_FFBS(simulator):
    simulator.contig_lens = [10]
    inputs = []
    for contig_len in simulator.contig_lens:
        for ploidy in simulator.ploidies:
            # stop
            genotype_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
            # genotype_df = pd.read_csv(genotype_path)
            for coverage in simulator.coverages:
                this_cov_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                frag_path = os.path.join(this_cov_path, 'frag')
                frag_graph_path = os.path.join(this_cov_path, 'fgraph')
                quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
                qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')
                results_path = os.path.join(this_cov_path, 'results_algorithm_v2')

                if not os.path.exists(frag_graph_path):
                    os.makedirs(frag_graph_path)
                if not os.path.exists(quotient_graph_path):
                    os.makedirs(quotient_graph_path)
                if not os.path.exists(qgraph_reverse_maps_path):
                    os.makedirs(qgraph_reverse_maps_path)
                if not os.path.exists(results_path):
                    os.makedirs(results_path)

                # existing_files_qg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_e_label' in ff]
                # existing_files_qg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_v_label' in ff]
                # existing_files_fg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_e_label' in ff]
                # existing_files_fg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_v_label' in ff]
                # existing_fg = [ff for ff in os.listdir(frag_graph_path) if '.gt.gz' in ff]
                # existing_qg = [ff for ff in os.listdir(quotient_graph_path) if '.gt.gz' in ff]
                # existing_results = [ff for ff in os.listdir(results_path) if 'FFBS' in ff]
                existing_results = []
                for rd in range(simulator.n_samples):
                    if 'FFBS_{}.pkl'.format(str(rd).zfill(2)) not in existing_results and \
                        '{}.gt.gz'.format(str(rd).zfill(2)) in os.listdir(quotient_graph_path) and \
                        'qg_e_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path) and \
                        'qg_v_label_' + str(rd).zfill(2) + '.pkl' in os.listdir(qgraph_reverse_maps_path):
                        inp = [frag_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path, results_path]
                        inputs.append(inp)
        inputs = sorted(inputs, key=lambda x: x[0])
    return inputs


def make_inputs_for_chordal_contraction(simulator):
    simulator.contig_lens = [100]
    k = 10
    inputs = []
    for contig_len in simulator.contig_lens:
        for ploidy in simulator.ploidies:
            # stop
            genotype_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
            # genotype_df = pd.read_csv(genotype_path)
            for coverage in simulator.coverages:
                this_cov_path = os.path.join(simulator.main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                frag_path = os.path.join(this_cov_path, 'frag')
                frag_graph_path = os.path.join(this_cov_path, 'fgraph')
                quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
                qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')
                chordal_graph_path = os.path.join(this_cov_path, 'chgraph')

                if not os.path.exists(frag_graph_path):
                    os.makedirs(frag_graph_path)
                if not os.path.exists(quotient_graph_path):
                    os.makedirs(quotient_graph_path)
                if not os.path.exists(qgraph_reverse_maps_path):
                    os.makedirs(qgraph_reverse_maps_path)
                if not os.path.exists(chordal_graph_path):
                    os.makedirs(chordal_graph_path)

                # existing_files_qg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_e_label' in ff]
                # existing_files_qg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_v_label' in ff]
                # existing_files_fg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_e_label' in ff]
                # existing_files_fg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_v_label' in ff]
                # existing_fg = [ff for ff in os.listdir(frag_graph_path) if '.gt.gz' in ff]
                # existing_qg = [ff for ff in os.listdir(quotient_graph_path) if '.gt.gz' in ff]
                existing_results = [ff for ff in os.listdir(chordal_graph_path) if '.gt.gz' in ff]

                for rd in range(simulator.n_samples):
                    if '{}.gt.gz'.format(str(rd).zfill(2)) not in existing_results and \
                        '{}.gt.gz'.format(str(rd).zfill(2)) in os.listdir(quotient_graph_path):
                        inp = [frag_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), chordal_graph_path, genotype_path, ploidy, k]
                        
                        inputs.append(inp)
        inputs = sorted(inputs, key=lambda x: x[0])
    return inputs


def make_inputs_for_run_likelihood():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_dir = '/mnt/research/aguiarlab/proj/HaplOrbit/inputs'
    contig_lens = [100]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    inputs = []
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # stop
            genotype_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
            # genotype_df = pd.read_csv(genotype_path)
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                frag_path = os.path.join(this_cov_path, 'frag')
                frag_graph_path = os.path.join(this_cov_path, 'fgraph')
                quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
                qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')
                # results_path = os.path.join(this_cov_path, 'results_FFBS_Multiple')
                results_path = os.path.join(this_cov_path, 'results_FFBS_single_mv_samples')

                if not os.path.exists(frag_graph_path):
                    os.makedirs(frag_graph_path)
                if not os.path.exists(quotient_graph_path):
                    os.makedirs(quotient_graph_path)
                if not os.path.exists(qgraph_reverse_maps_path):
                    os.makedirs(qgraph_reverse_maps_path)
                if not os.path.exists(results_path):
                    os.makedirs(results_path)

                # existing_files_qg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_e_label' in ff]
                # existing_files_qg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'qg_v_label' in ff]
                # existing_files_fg_e = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_e_label' in ff]
                # existing_files_fg_v = [ff for ff in os.listdir(os.path.join(qgraph_reverse_maps_path)) if 'fg_v_label' in ff]
                # existing_fg = [ff for ff in os.listdir(frag_graph_path) if '.gt.gz' in ff]
                # existing_qg = [ff for ff in os.listdir(quotient_graph_path) if '.gt.gz' in ff]
                existing_results = [ff for ff in os.listdir(results_path) if 'FFBS' in ff]
                # existing_results = []
                # for rd in range(90, 100):
                for rd in range(n_samples):
                    
                    if 'FFBS_{}.pkl'.format(str(rd).zfill(2)) not in existing_results:
                        print(frag_path)
                        inp = [frag_path, frag_graph_path, quotient_graph_path, qgraph_reverse_maps_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path, results_path]
                        inputs.append(inp)
    # inputs = sorted(inputs, key=lambda x: x[0], reverse=True)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, inp in enumerate(inputs):
        input_file = os.path.join(output_dir, f"input_{i}.pkl")
        with open(input_file, "wb") as f:
            pickle.dump(inp, f)
    print(f"Saved {len(inputs)} inputs to {output_dir}")


def make_inputs_for_checking_pairwise_samples():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_dir = '/mnt/research/aguiarlab/proj/HaplOrbit/inputs_pairwise'
    contig_lens = [100]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    inputs = []
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # stop
            genotype_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
            # genotype_df = pd.read_csv(genotype_path)
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                frag_path = os.path.join(this_cov_path, 'frag')
                frag_graph_path = os.path.join(this_cov_path, 'fgraph')
                quotient_graph_path = os.path.join(this_cov_path, 'qgraph')
                qgraph_reverse_maps_path = os.path.join(this_cov_path, 'reverse_maps')
                results_path = os.path.join(this_cov_path, 'pairwise_results')

                if not os.path.exists(frag_graph_path):
                    os.makedirs(frag_graph_path)
                if not os.path.exists(quotient_graph_path):
                    os.makedirs(quotient_graph_path)
                if not os.path.exists(qgraph_reverse_maps_path):
                    os.makedirs(qgraph_reverse_maps_path)
                if not os.path.exists(results_path):
                    os.makedirs(results_path)

                existing_results = []
                # for rd in range(90, 100):
                for rd in range(n_samples):
                    
                    if 'FFBS_{}.pkl'.format(str(rd).zfill(2)) not in existing_results:
                        print(frag_path)
                        inp = [frag_path, '{}.frag'.format(str(rd).zfill(2)), ploidy, genotype_path, results_path]
                        inputs.append(inp)
    inputs = sorted(inputs, key=lambda x: x[0], reverse=True)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, inp in enumerate(inputs):
        input_file = os.path.join(output_dir, f"input_{i}.pkl")
        with open(input_file, "wb") as f:
            pickle.dump(inp, f)
    print(f"Saved {len(inputs)} inputs to {output_dir}")


def save_inputs(inputs, output_dir):
    """
    Save each input as a separate pickle file in the specified output directory.
    """
    output_dir = '/mnt/research/aguiarlab/proj/HaplOrbit/inputs100_2'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, inp in enumerate(inputs):
        input_file = os.path.join(output_dir, f"input_{i}.pkl")
        with open(input_file, "wb") as f:
            pickle.dump(inp, f)
    print(f"Saved {len(inputs)} inputs to {output_dir}")


def run_chordal_contraction_graph_tool_top_k(inp):
    this_frag_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, chordal_graph_path, genotype_path, ploidy, k = inp

    # save_path, subg_id, subg, config, fragment_model, k = inp
    # this_path = os.path.join(save_path, 'chordal_sub_' + str(subg_id) + '.gt.gz')

    print('Working on:', os.path.join(this_frag_path, frag_file))
    chordal_v_label_revered_path = os.path.join(this_reverse_maps_path, 'ch_v_label_' + frag_file.split('.')[0] + '.pkl')
    chordal_e_label_revered_path = os.path.join(this_reverse_maps_path, 'ch_e_label_' + frag_file.split('.')[0] + '.pkl')


    class Args:
        def __init__(self):
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join(this_frag_path, frag_file)
            # self.data_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/c2/ART_90.frag.txt'
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

    quotient_g_path = os.path.join(this_quotient_coverage_path, frag_file.split('.')[0] + '.gt.gz')
    quotient_g = gt.load_graph(quotient_g_path)

    new_graph = quotient_g.copy()
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

    e_labels_ch = new_graph.edge_properties["e_label"]
    v_labels_ch = new_graph.vertex_properties["v_label"]
        
    v_label_reversed = {}
    for v in new_graph.vertices():
        v_label = v_labels_ch[v]
        v_label_reversed[v_label] = int(v)
    
    e_label_reversed = {}
    for e in new_graph.edges():
        e_label = e_labels_ch[e]
        v1_label, v2_label = e_label.split('--')
        e_label_reversed[e_label] = [v_label_reversed[v1_label], v_label_reversed[v2_label]]
        

    this_path = os.path.join(chordal_graph_path, frag_file.split('.')[0] + '.gt.gz')
    new_graph.save(this_path)

    with open(chordal_v_label_revered_path, 'wb') as f:
        pickle.dump(v_label_reversed, f)

    with open(chordal_e_label_revered_path, 'wb') as f:
        pickle.dump(e_label_reversed, f)

    print('[Done]', this_path)


def run_FFBS_quotient(inp):
    this_frag_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
    print('Working on:', os.path.join(this_frag_path, frag_file))
    # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
    # frag_path = '/labs/Aguiar/pHapCompass/FFBS/test2.frag'
    # ploidy= 3
    # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
    # genotype_path = '/labs/Aguiar/pHapCompass/FFBS/haplotypes.csv'

    class Args:
        def __init__(self):
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join(this_frag_path, frag_file)
            # self.data_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/c2/ART_90.frag.txt'
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

    # frag_graph_path = os.path.join(this_frag_graph_path, frag_file.split('.')[0] + '.gt.gz')
    # frag_graph = gt.load_graph(frag_graph_path)

    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
    fragment_model.construct(input_handler, config)

    quotient_g_path = os.path.join(this_quotient_coverage_path, frag_file.split('.')[0] + '.gt.gz')
    quotient_g = gt.load_graph(quotient_g_path)

    edge_map_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + frag_file.split('.')[0] + '.pkl')
    with open(edge_map_path, 'rb') as f:
        edges_map_quotient = pickle.load(f)
        
    quotient_g_v_label_reversed_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + frag_file.split('.')[0] + '.pkl')
    with open(quotient_g_v_label_reversed_path, 'rb') as f:
        quotient_g_v_label_reversed = pickle.load(f)

    start_time = time.time()

    transitions_dict, transitions_dict_extra = transition_matrices_v2(quotient_g, edges_map_quotient, ploidy, config, fragment_model)
    emission_dict = emissions_v2(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, _ =  assign_slices_and_interfaces(nodes, edges)

    assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, args.data_path)

    forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, args.data_path)

    # backward_messages = compute_backward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, args.data_path)

    samples = sample_states_book(slices, edges, forward_messages, transitions_dict)
    # samples = sample_states(slices, edges, forward_messages, transitions_dict)
    # for k in samples.keys():
    #     kedges = samples[k].keys()
    #     for e in kedges:
    #         print(e, samples[k][e])
    
    predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config)
    end_time = time.time()
    elapsed_time = round(end_time - start_time, 2)


    # print('Predicted Haplotypes:\n', predicted_haplotypes)
    # print('\nTrue Haplotypes:\n', pd.read_csv(genotype_path).T)
    sampled_positions = [c for c in predicted_haplotypes.columns.values if np.nan not in list(predicted_haplotypes[c].values)]

    predicted_haplotypes_np = predicted_haplotypes[sampled_positions].to_numpy()
    true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()[:, sampled_positions]

    vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_np, true_haplotypes)
    accuracy, _ = calculate_accuracy(predicted_haplotypes_np, true_haplotypes)
    mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_np, true_haplotypes)
    mec_ = mec(predicted_haplotypes_np, fragment_model.fragment_list)
    results_name = 'FFBS_{}.pkl'.format(frag_file.split('.')[0])
    results = {}
    results['evaluation'] = {'vector_error_rate': vector_error_rate, 'vector_error': vector_error, 'backtracking_steps': backtracking_steps, 
                             'dp_table': dp_table, 'accuracy': accuracy, 'mismatch_error': mismatch_error, 'mec': mec_}
    results['predicted_haplotypes'] = predicted_haplotypes_np
    results['true_haplotypes'] = true_haplotypes
    results['forward_messages'] = forward_messages
    # results['backward_messages'] = backward_messages
    results['transitions_dict'] = transitions_dict
    results['transitions_dict_extra'] = transitions_dict_extra
    results['emission_dict'] = emission_dict
    results['assignment_dict'] = assignment_dict
    results['samples'] = samples
    results['slices'] = slices
    results['best_permutation'] = best_permutation
    results['fragment_list'] = fragment_model.fragment_list
    results['time'] = elapsed_time

    with open(os.path.join(results_path, results_name), 'wb') as f:
        pickle.dump(results, f)

    print('Saved results in {}.'.format(os.path.join(results_path, results_name)), 'vector_error_rate', vector_error_rate, 'accuracy', accuracy, 'mismatch_error', mismatch_error, 'mec', mec_)
    

def run_FFBS_quotient_count(inp):
    this_frag_path, this_fragment_coverage_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
    file_name = frag_file.split('.')[0]

    fragment_v_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_v_label_' + file_name + '.pkl')
    fragment_e_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_e_label_' + file_name + '.pkl')
    quotient_v_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + file_name + '.pkl')
    quotient_e_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + file_name + '.pkl')
    print('Working on:', os.path.join(this_frag_path, frag_file))

    class Args:
        def __init__(self):
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join(this_frag_path, frag_file)
            # self.data_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/c2/ART_90.frag.txt'
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

    start_time = time.time()

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
    fragment_model.construct(input_handler, config)

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
    quotient_g.construct(input_handler, config)

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

    quotient_g_v_label_reversed = quotient_g.v_label_reversed

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model, config)
    emission_dict = emissions(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, interfaces =  assign_slices_and_interfaces(nodes, edges)

    assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, input_handler.data_path)

    forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)

    samples = sample_states_book(slices, edges, forward_messages, transitions_dict)

    predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="counts")

    end_time = time.time()

    elapsed_time = round(end_time - start_time, 2)

    sampled_positions = [c for c in predicted_haplotypes.columns.values if np.nan not in list(predicted_haplotypes[c].values)]

    predicted_haplotypes_np = predicted_haplotypes[sampled_positions].to_numpy()
    true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()[:, sampled_positions]

    vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_np, true_haplotypes)
    accuracy, _ = calculate_accuracy(predicted_haplotypes_np, true_haplotypes)
    mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_np, true_haplotypes)
    mec_ = mec(predicted_haplotypes_np, fragment_model.fragment_list)
    results_name = 'FFBS_{}.pkl'.format(frag_file.split('.')[0])
    results = {}
    results['evaluation'] = {'vector_error_rate': vector_error_rate, 'vector_error': vector_error, 'backtracking_steps': backtracking_steps, 
                             'dp_table': dp_table, 'accuracy': accuracy, 'mismatch_error': mismatch_error, 'mec': mec_}
    results['predicted_haplotypes'] = predicted_haplotypes_np
    results['true_haplotypes'] = true_haplotypes
    results['forward_messages'] = forward_messages
    results['transitions_dict'] = transitions_dict
    results['transitions_dict_extra'] = transitions_dict_extra
    results['emission_dict'] = emission_dict
    results['assignment_dict'] = assignment_dict
    results['samples'] = samples
    results['slices'] = slices
    results['best_permutation'] = best_permutation
    results['fragment_list'] = fragment_model.fragment_list
    results['time'] = elapsed_time
    # print('Results:', results['evaluation'])

    with open(os.path.join(results_path, results_name), 'wb') as f:
        pickle.dump(results, f)

    print('Saved results in {}.'.format(os.path.join(results_path, results_name)), 'vector_error_rate', vector_error_rate, 'accuracy', accuracy, 'mismatch_error', mismatch_error, 'mec', mec_)
    

def run_FFBS_quotient_likelihood(inp):
    this_frag_path, this_fragment_coverage_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
    file_name = frag_file.split('.')[0]

    fragment_v_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_v_label_' + file_name + '.pkl')
    fragment_e_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_e_label_' + file_name + '.pkl')
    quotient_v_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + file_name + '.pkl')
    quotient_e_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + file_name + '.pkl')
    print('Working on:', os.path.join(this_frag_path, frag_file))

    class Args:
        def __init__(self):
            # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
            # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join(this_frag_path, frag_file)
            # self.data_path = frag_path
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

    start_time = time.time()

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
    fragment_model.construct(input_handler, config)

    frag_graph_path = os.path.join(this_fragment_coverage_path, file_name + '.gt.gz')
    fragment_model.graph.save(frag_graph_path)

    frag_graph_plot_path = os.path.join(this_fragment_coverage_path, file_name + '.png')

    e_labels = fragment_model.graph.edge_properties["e_label"]
    v_labels = fragment_model.graph.vertex_properties["v_label"]
    gt.graph_draw(fragment_model.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=frag_graph_plot_path)

    with open(fragment_v_label_revered_path, "wb") as f:
        pickle.dump(fragment_model.v_label_reversed, f)

    edges_map_fragment = {}
    for k in fragment_model.e_label_reversed.keys():
        edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]

    with open(fragment_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_fragment, f)

    # create quotient graph
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct(input_handler, config)

    # save quotient graph
    quot_graph_path = os.path.join(this_quotient_coverage_path, file_name + '.gt.gz')
    quotient_g.graph.save(quot_graph_path)


    quotient_graph_plot_path = os.path.join(this_quotient_coverage_path, file_name + '.png')
    e_labels = quotient_g.graph.edge_properties["e_label"]
    v_labels = quotient_g.graph.vertex_properties["v_label"]
    gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=quotient_graph_plot_path)

    # gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    # edge_font_size=10)
    gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, vertex_font_size=24)

    with open(quotient_v_label_revered_path, "wb") as f:
        pickle.dump(quotient_g.v_label_reversed, f)

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    with open(quotient_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_quotient, f)

    quotient_g_v_label_reversed = quotient_g.v_label_reversed

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model, config)
    emission_dict = emissions(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, interfaces =  assign_slices_and_interfaces(nodes, edges)

    assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, input_handler.data_path)

    forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)

    # backward_messages = compute_backward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)   


    samples = sample_states_book(slices, edges, forward_messages, transitions_dict)
    # samples = sample_states_book_multiple_times(slices, edges, forward_messages, transitions_dict, n=100)
    # samples = sample_states_ground_truth(slices, nodes, genotype_path)
    # fragment_list = fragment_model.fragment_list
    # reads_dict = calculate_pair_counts(fragment_list)

    ffbs_acc = evaulate_ffbs_acc_sample(genotype_path, samples, ploidy)
    # print('FFBS Accuracy:', ffbs_acc)
    predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")

    # for _ in range(10):
    #     samples = sample_states_book(slices, edges, forward_messages, transitions_dict)
    #     predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")
    #     ffbs_acc = evaulate_ffbs_acc_sample(genotype_path, samples, ploidy)
    #     this_phasing_likelihood = compute_global_phasing_likelihood(predicted_haplotypes, fragment_model, config)
    #     print(ffbs_acc, this_phasing_likelihood)
    #     compute_vector_error_rate(predicted_haplotypes.to_numpy(), true_haplotypes.to_numpy())

    end_time = time.time()

    elapsed_time = round(end_time - start_time, 2)

    true_haplotypes = pd.read_csv(genotype_path).T

    block_info, components = get_block_info(quotient_g, predicted_haplotypes, true_haplotypes, fragment_model)

    sampled_positions = [c for c in predicted_haplotypes.columns.values if np.nan not in list(predicted_haplotypes[c].values)]

    predicted_haplotypes_np = predicted_haplotypes[sampled_positions].to_numpy()
    # true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()[:, sampled_positions]
    
    true_haplotypes_np = true_haplotypes.to_numpy()[:, sampled_positions]

    vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_np, true_haplotypes_np)
    accuracy, _ = calculate_accuracy(predicted_haplotypes_np, true_haplotypes_np)
    mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_np, true_haplotypes_np)
    mec_ = mec(predicted_haplotypes_np, fragment_model.fragment_list)
    results_name = 'FFBS_{}.pkl'.format(frag_file.split('.')[0])
    results = {}
    results['block_evaluation'] = block_info
    results['components'] = components
    results['n_blocks'] = len(components.keys())
    results['average_block_size'] = block_info['average_block_size']
    results['length_phased'] = len(sampled_positions)
    results['evaluation'] = {'vector_error_rate': vector_error_rate, 'vector_error': vector_error, 'backtracking_steps': backtracking_steps, 
                             'dp_table': dp_table, 'accuracy': accuracy, 'mismatch_error': mismatch_error, 'mec': mec_, 'ffbs_acc': ffbs_acc}
    results['predicted_haplotypes'] = predicted_haplotypes
    results['true_haplotypes'] = pd.read_csv(genotype_path).T
    results['forward_messages'] = forward_messages
    results['transitions_dict'] = transitions_dict
    results['transitions_dict_extra'] = transitions_dict_extra
    results['emission_dict'] = emission_dict
    results['assignment_dict'] = assignment_dict
    results['samples'] = samples
    results['slices'] = slices
    results['best_permutation'] = best_permutation
    results['fragment_list'] = fragment_model.fragment_list
    results['time'] = elapsed_time
    # print('Results:', results['evaluation'])

    with open(os.path.join(results_path, results_name), 'wb') as f:
        pickle.dump(results, f)

    print('Saved results in {}.'.format(os.path.join(results_path, results_name)), 'vector_error_rate', vector_error_rate, 'vector_error', vector_error, 'mismatch_error', mismatch_error, 'mec', mec_, 'ffbs_acc', ffbs_acc)
    

def run_FFBS_Single(inp):
    this_frag_path, this_fragment_coverage_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
    file_name = frag_file.split('.')[0]

    fragment_v_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_v_label_' + file_name + '.pkl')
    fragment_e_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_e_label_' + file_name + '.pkl')
    quotient_v_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + file_name + '.pkl')
    quotient_e_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + file_name + '.pkl')
    print('Working on:', os.path.join(this_frag_path, frag_file))

    class Args:
        def __init__(self):
            # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
            # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join(this_frag_path, frag_file)
            # self.data_path = frag_path
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

    start_time = time.time()

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
    fragment_model.construct(input_handler, config)

    edges_map_fragment = {}
    for k in fragment_model.e_label_reversed.keys():
        edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]

    # save fragment graph data:
    frag_graph_path = os.path.join(this_fragment_coverage_path, file_name + '.gt.gz')
    fragment_model.graph.save(frag_graph_path)

    frag_graph_plot_path = os.path.join(this_fragment_coverage_path, file_name + '.png')

    e_labels = fragment_model.graph.edge_properties["e_label"]
    v_labels = fragment_model.graph.vertex_properties["v_label"]
    gt.graph_draw(fragment_model.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=frag_graph_plot_path)

    with open(fragment_v_label_revered_path, "wb") as f:
        pickle.dump(fragment_model.v_label_reversed, f)

    with open(fragment_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_fragment, f)

    # create quotient graph
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct(input_handler, config)

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    quotient_g_v_label_reversed = quotient_g.v_label_reversed

    # save quotient graph
    quot_graph_path = os.path.join(this_quotient_coverage_path, file_name + '.gt.gz')
    quotient_g.graph.save(quot_graph_path)

    quotient_graph_plot_path = os.path.join(this_quotient_coverage_path, file_name + '.png')
    e_labels = quotient_g.graph.edge_properties["e_label"]
    v_labels = quotient_g.graph.vertex_properties["v_label"]
    gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=quotient_graph_plot_path)

    # gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    # edge_font_size=10)
    # gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, vertex_font_size=24)

    with open(quotient_v_label_revered_path, "wb") as f:
        pickle.dump(quotient_g.v_label_reversed, f)

    with open(quotient_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_quotient, f)


    transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model, config)
    emission_dict = emissions(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, _ =  assign_slices_and_interfaces(nodes, edges)

    assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, input_handler.data_path)

    forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)

    samples = sample_states_book(slices, edges, forward_messages, transitions_dict)

    ffbs_acc = evaulate_ffbs_acc_sample(genotype_path, samples, ploidy)
    # ffbs_acc, stats_dict = evaluate_ffbs_acc_sample_with_uncertainty(genotype_path, consensus_samples, ploidy)
    predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")

    end_time = time.time()

    elapsed_time = round(end_time - start_time, 2)

    true_haplotypes = pd.read_csv(genotype_path).T

    sampled_positions = [c for c in predicted_haplotypes.columns.values if np.nan not in list(predicted_haplotypes[c].values)]


    predicted_haplotypes_np = predicted_haplotypes.to_numpy()

    true_haplotypes_np = true_haplotypes.to_numpy()


    vector_error_rate, _, _ = compute_vector_error_rate_with_missing_positions(true_haplotypes_np, predicted_haplotypes_np)


    results_name = 'FFBS_{}.pkl'.format(frag_file.split('.')[0])
    results = {'evaluation': {}}
    results['evaluation']['length_phased'] = len(sampled_positions)
    results['evaluation']['vector_error_rate'] = vector_error_rate
    results['ffbs_acc'] = ffbs_acc
    results['true_haplotypes'] = pd.read_csv(genotype_path).T
    results['forward_messages'] = forward_messages
    results['transitions_dict'] = transitions_dict
    results['transitions_dict_extra'] = transitions_dict_extra
    results['emission_dict'] = emission_dict
    results['assignment_dict'] = assignment_dict
    results['samples'] = samples
    results['fragment_list'] = fragment_model.fragment_list
    results['time'] = elapsed_time
    # print('Results:', results['evaluation'])

    with open(os.path.join(results_path, results_name), 'wb') as f:
        pickle.dump(results, f)

    print('Saved results in {}.'.format(os.path.join(results_path, results_name))) #, 'vector_error_rate', vector_error_rate, 'ffbs_acc', ffbs_acc)


def load_test_data():
    # Load test data
    frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
    genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
    ploidy = 3
    class Args:
        def __init__(self):
            # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
            # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = frag_path
            # self.data_path = frag_path
            self.bam_path = 'example/example.bam'
            self.genotype_path = genotype_path
            self.ploidy = ploidy
            self.error_rate = 0.001
            self.epsilon = 0.0001
            self.output_path = 'output'
            self.root_dir = 'D:/UCONN/HaplOrbit'
            self.alleles = [0, 1]
    args = Args()
    return args


def run_viterbi_Single(inp):
    this_frag_path, this_fragment_coverage_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
    file_name = frag_file.split('.')[0]

    fragment_v_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_v_label_' + file_name + '.pkl')
    fragment_e_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_e_label_' + file_name + '.pkl')
    quotient_v_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + file_name + '.pkl')
    quotient_e_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + file_name + '.pkl')
    print('Working on:', os.path.join(this_frag_path, frag_file))

    class Args:
        def __init__(self):
            # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
            # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join(this_frag_path, frag_file)
            # self.data_path = frag_path
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

    # args = load_test_data()

    start_time = time.time()

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
    fragment_model.construct(input_handler, config)

    edges_map_fragment = {}
    for k in fragment_model.e_label_reversed.keys():
        edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]

    # save fragment graph data:
    frag_graph_path = os.path.join(this_fragment_coverage_path, file_name + '.gt.gz')
    fragment_model.graph.save(frag_graph_path)

    frag_graph_plot_path = os.path.join(this_fragment_coverage_path, file_name + '.png')

    e_labels = fragment_model.graph.edge_properties["e_label"]
    v_labels = fragment_model.graph.vertex_properties["v_label"]
    gt.graph_draw(fragment_model.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=frag_graph_plot_path)

    with open(fragment_v_label_revered_path, "wb") as f:
        pickle.dump(fragment_model.v_label_reversed, f)

    with open(fragment_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_fragment, f)

    # create quotient graph
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct(input_handler, config)

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    quotient_g_v_label_reversed = quotient_g.v_label_reversed

    # save quotient graph
    quot_graph_path = os.path.join(this_quotient_coverage_path, file_name + '.gt.gz')
    quotient_g.graph.save(quot_graph_path)

    quotient_graph_plot_path = os.path.join(this_quotient_coverage_path, file_name + '.png')
    e_labels = quotient_g.graph.edge_properties["e_label"]
    v_labels = quotient_g.graph.vertex_properties["v_label"]
    gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=quotient_graph_plot_path)

    # gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    # edge_font_size=10)
    # gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, vertex_font_size=24)

    with open(quotient_v_label_revered_path, "wb") as f:
        pickle.dump(quotient_g.v_label_reversed, f)

    with open(quotient_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_quotient, f)


    transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, args.ploidy, fragment_model, config)
    emission_dict = emissions(args.ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, _ =  assign_slices_and_interfaces(nodes, edges)

    assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, input_handler.data_path)

    samples = run_viterbi(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)


    ffbs_acc = evaulate_ffbs_acc_sample(genotype_path, samples, args.ploidy)
    # ffbs_acc, stats_dict = evaluate_ffbs_acc_sample_with_uncertainty(genotype_path, consensus_samples, args.ploidy)
    predicted_haplotypes = predict_haplotypes(nodes, edges, samples, args.ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")

    end_time = time.time()

    elapsed_time = round(end_time - start_time, 2)

    true_haplotypes = pd.read_csv(genotype_path).T

    sampled_positions = [c for c in predicted_haplotypes.columns.values if np.nan not in list(predicted_haplotypes[c].values)]


    predicted_haplotypes_np = predicted_haplotypes.to_numpy()

    true_haplotypes_np = true_haplotypes.to_numpy()


    vector_error_rate, _, _ = compute_vector_error_rate_with_missing_positions(true_haplotypes_np, predicted_haplotypes_np)


    results_name = 'FFBS_{}.pkl'.format(frag_file.split('.')[0])
    results = {'evaluation': {}}
    results['evaluation']['length_phased'] = len(sampled_positions)
    results['evaluation']['vector_error_rate'] = vector_error_rate
    results['ffbs_acc'] = ffbs_acc
    results['true_haplotypes'] = pd.read_csv(genotype_path).T
    results['forward_messages'] = forward_messages
    results['transitions_dict'] = transitions_dict
    results['transitions_dict_extra'] = transitions_dict_extra
    results['emission_dict'] = emission_dict
    results['assignment_dict'] = assignment_dict
    results['samples'] = samples
    results['fragment_list'] = fragment_model.fragment_list
    results['time'] = elapsed_time
    # print('Results:', results['evaluation'])

    with open(os.path.join(results_path, results_name), 'wb') as f:
        pickle.dump(results, f)

    print('Saved results in {}.'.format(os.path.join(results_path, results_name))) #, 'vector_error_rate', vector_error_rate, 'ffbs_acc', ffbs_acc)


def run_FFBS_single_multi_samles(inp):
    this_frag_path, this_fragment_coverage_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
    file_name = frag_file.split('.')[0]

    fragment_v_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_v_label_' + file_name + '.pkl')
    fragment_e_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_e_label_' + file_name + '.pkl')
    quotient_v_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + file_name + '.pkl')
    quotient_e_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + file_name + '.pkl')
    print('Working on:', os.path.join(this_frag_path, frag_file))

    class Args:
        def __init__(self):
            # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
            # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join(this_frag_path, frag_file)
            # self.data_path = frag_path
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

    start_time = time.time()

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
    fragment_model.construct(input_handler, config)

    edges_map_fragment = {}
    for k in fragment_model.e_label_reversed.keys():
        edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]

    # save fragment graph data:
    frag_graph_path = os.path.join(this_fragment_coverage_path, file_name + '.gt.gz')
    fragment_model.graph.save(frag_graph_path)

    frag_graph_plot_path = os.path.join(this_fragment_coverage_path, file_name + '.png')

    e_labels = fragment_model.graph.edge_properties["e_label"]
    v_labels = fragment_model.graph.vertex_properties["v_label"]
    gt.graph_draw(fragment_model.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=frag_graph_plot_path)

    with open(fragment_v_label_revered_path, "wb") as f:
        pickle.dump(fragment_model.v_label_reversed, f)

    with open(fragment_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_fragment, f)

    # create quotient graph
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct(input_handler, config)

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    quotient_g_v_label_reversed = quotient_g.v_label_reversed

    # save quotient graph
    quot_graph_path = os.path.join(this_quotient_coverage_path, file_name + '.gt.gz')
    quotient_g.graph.save(quot_graph_path)

    quotient_graph_plot_path = os.path.join(this_quotient_coverage_path, file_name + '.png')
    e_labels = quotient_g.graph.edge_properties["e_label"]
    v_labels = quotient_g.graph.vertex_properties["v_label"]
    gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=quotient_graph_plot_path)

    with open(quotient_v_label_revered_path, "wb") as f:
        pickle.dump(quotient_g.v_label_reversed, f)

    with open(quotient_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_quotient, f)


    transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model, config)
    emission_dict = emissions(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, _ =  assign_slices_and_interfaces(nodes, edges)

    assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, input_handler.data_path)

    predicted_haplotypes_single_list = []
    predicted_haplotypes_mv_list = []
    ffbs_acc_list = []
    likelihood_single_list = []
    likelihood_mv_list = []

    for i in range(10):
        forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)
        samples = sample_states_book(slices, edges, forward_messages, transitions_dict)
        ffbs_acc = evaulate_ffbs_acc_sample(genotype_path, samples, ploidy)

        predicted_haplotypes_single = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")
        predicted_haplotypes_np_single = predicted_haplotypes_single.to_numpy()
        likelihood_single = compute_solution_likelihood(predicted_haplotypes_np_single, fragment_model, config)
        likelihood_single_list.append(likelihood_single)

        predicted_haplotypes_mv = predict_haplotypes_multiple_variants(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")
        predicted_haplotypes_np_mv = predicted_haplotypes_mv.to_numpy()
        likelihood_mv = compute_solution_likelihood(predicted_haplotypes_np_mv, fragment_model, config)
        likelihood_mv_list.append(likelihood_mv)

        predicted_haplotypes_single_list.append(predicted_haplotypes_single)
        predicted_haplotypes_mv_list.append(predicted_haplotypes_mv)
        ffbs_acc_list.append(ffbs_acc)

    
    # Find the index of the highest likelihood in single likelihood list
    best_single_idx = np.argmax(likelihood_single_list)
    best_predicted_haplotypes_single = predicted_haplotypes_single_list[best_single_idx]
    best_ffbs_acc_single = ffbs_acc_list[best_single_idx]

    # Find the index of the highest likelihood in mv likelihood list
    best_mv_idx = np.argmax(likelihood_mv_list)
    best_predicted_haplotypes_mv = predicted_haplotypes_mv_list[best_mv_idx]
    best_ffbs_acc_mv = ffbs_acc_list[best_mv_idx]

    # Compute entropy of likelihood lists
    entropy_single = entropy(likelihood_single_list)
    entropy_mv = entropy(likelihood_mv_list)

    # Compute confidence of the best sample (portion of elements that are the same as max)
    max_single = max(likelihood_single_list)
    max_mv = max(likelihood_mv_list)

    confidence_single = likelihood_single_list.count(max_single) / len(likelihood_single_list)
    confidence_mv = likelihood_mv_list.count(max_mv) / len(likelihood_mv_list)

    end_time = time.time()

    elapsed_time = round(end_time - start_time, 2)

    true_haplotypes = pd.read_csv(genotype_path).T

    sampled_positions_single = [c for c in best_predicted_haplotypes_single.columns.values if np.nan not in list(best_predicted_haplotypes_single[c].values)]
    sampled_positions_mv = [c for c in best_predicted_haplotypes_mv.columns.values if np.nan not in list(best_predicted_haplotypes_mv[c].values)]

    best_predicted_haplotypes_single_np = best_predicted_haplotypes_single.to_numpy()
    best_predicted_haplotypes_mv_np = best_predicted_haplotypes_mv.to_numpy()

    true_haplotypes_np = true_haplotypes.to_numpy()

    vector_error_rate_single, _, _ = compute_vector_error_rate_with_missing_positions(true_haplotypes_np, best_predicted_haplotypes_single_np)
    vector_error_rate_mv, _, _ = compute_vector_error_rate_with_missing_positions(true_haplotypes_np, best_predicted_haplotypes_mv_np)
    
    results_name = 'FFBS_{}.pkl'.format(frag_file.split('.')[0])
    results = {'Single Variant': {}, 'Multi Variant': {}}
    
    results['Single Variant']['length_phased'] = len(sampled_positions_single)
    results['Multi Variant']['length_phased'] = len(sampled_positions_mv)

    results['Single Variant']['vector_error_rate'] = vector_error_rate_single
    results['Multi Variant']['vector_error_rate'] = vector_error_rate_mv

    results['Single Variant']['ffbs_acc'] = best_ffbs_acc_single
    results['Multi Variant']['ffbs_acc'] = best_ffbs_acc_mv

    results['Single Variant']['solution_confidence'] = confidence_single
    results['Multi Variant']['solution_confidence'] = confidence_mv

    results['Single Variant']['solution_entropy'] = entropy_single
    results['Multi Variant']['solution_entropy'] = entropy_mv

    results['true_haplotypes'] = pd.read_csv(genotype_path).T
    results['fragment_list'] = fragment_model.fragment_list
    results['time'] = elapsed_time
    # print('Results:', results['evaluation'])

    with open(os.path.join(results_path, results_name), 'wb') as f:
        pickle.dump(results, f)

    print('Saved results in {}.'.format(os.path.join(results_path, results_name))) #, 'vector_error_rate', vector_error_rate, 'ffbs_acc', ffbs_acc)



def run_LBP_FFBS_Single(inp):
    this_frag_path, this_fragment_coverage_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
    file_name = frag_file.split('.')[0]

    fragment_v_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_v_label_' + file_name + '.pkl')
    fragment_e_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_e_label_' + file_name + '.pkl')
    quotient_v_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + file_name + '.pkl')
    quotient_e_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + file_name + '.pkl')
    print('Working on:', os.path.join(this_frag_path, frag_file))

    class Args:
        def __init__(self):
            # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
            # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join(this_frag_path, frag_file)
            # self.data_path = frag_path
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

    start_time = time.time()

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
    fragment_model.construct(input_handler, config)

    edges_map_fragment = {}
    for k in fragment_model.e_label_reversed.keys():
        edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]

    # save fragment graph data:
    frag_graph_path = os.path.join(this_fragment_coverage_path, file_name + '.gt.gz')
    fragment_model.graph.save(frag_graph_path)

    frag_graph_plot_path = os.path.join(this_fragment_coverage_path, file_name + '.png')

    e_labels = fragment_model.graph.edge_properties["e_label"]
    v_labels = fragment_model.graph.vertex_properties["v_label"]
    gt.graph_draw(fragment_model.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=frag_graph_plot_path)

    with open(fragment_v_label_revered_path, "wb") as f:
        pickle.dump(fragment_model.v_label_reversed, f)

    with open(fragment_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_fragment, f)

    # create quotient graph
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct(input_handler, config)

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    quotient_g_v_label_reversed = quotient_g.v_label_reversed

    # save quotient graph
    quot_graph_path = os.path.join(this_quotient_coverage_path, file_name + '.gt.gz')
    quotient_g.graph.save(quot_graph_path)

    quotient_graph_plot_path = os.path.join(this_quotient_coverage_path, file_name + '.png')
    e_labels = quotient_g.graph.edge_properties["e_label"]
    v_labels = quotient_g.graph.vertex_properties["v_label"]
    gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=quotient_graph_plot_path)

    # gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    # edge_font_size=10)
    # gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, vertex_font_size=24)

    with open(quotient_v_label_revered_path, "wb") as f:
        pickle.dump(quotient_g.v_label_reversed, f)

    with open(quotient_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_quotient, f)


    transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model, config)
    emission_dict = emissions(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, _ =  assign_slices_and_interfaces(nodes, edges)

    assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, input_handler.data_path)

    # mn = create_markov_net(quotient_g.graph, transitions=transitions_dict) # Creates mrftools TorchMarkovNetwork

    # torch_bp = mrftools.TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=False, var_on=False)
    # torch_bp.infer(display="full")
    # torch_bp.load_beliefs()
    # update_entropies(quotient_g.graph, torch_bp)

    # transitions_dict = update_transitions(quotient_g.graph, transitions_dict, torch_bp)
    # # print("Transitions Updated")

    # vertices_dict = update_vertices(quotient_g_v_label_reversed, emission_dict, slices, torch_bp)

    forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)

    consensus_samples, samples = sample_states_book_multiple_times(slices, edges, forward_messages, transitions_dict, n=20)

    ffbs_acc = evaulate_ffbs_acc_sample(genotype_path, consensus_samples, ploidy)
    # ffbs_acc, stats_dict = evaluate_ffbs_acc_sample_with_uncertainty(genotype_path, consensus_samples, ploidy)
    predicted_haplotypes_single_count = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")
    predicted_haplotypes_single_count = predict_haplotypes_entropy_based(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")

    # predicted_haplotypes_single_evidence = predict_haplotypes_single_evidence_count(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, assignment_dict, config, priority="probabilities")
    # predicted_haplotypes_mv = predict_haplotypes_multiple_variants(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")





    end_time = time.time()

    elapsed_time = round(end_time - start_time, 2)

    true_haplotypes = pd.read_csv(genotype_path).T

    # block_info, components = get_block_info(quotient_g, predicted_haplotypes, true_haplotypes, fragment_model)

    sampled_positions_count = [c for c in predicted_haplotypes_single_count.columns.values if np.nan not in list(predicted_haplotypes_single_count[c].values)]
    sampled_positions_evidence = [c for c in predicted_haplotypes_single_evidence.columns.values if np.nan not in list(predicted_haplotypes_single_evidence[c].values)]
    sampled_positions_mv = [c for c in predicted_haplotypes_mv.columns.values if np.nan not in list(predicted_haplotypes_mv[c].values)]

    predicted_haplotypes_np_count = predicted_haplotypes_single_count.to_numpy()
    predicted_haplotypes_np_evidence = predicted_haplotypes_single_evidence.to_numpy()
    predicted_haplotypes_np_mv = predicted_haplotypes_mv.to_numpy()
    # predicted_haplotypes_np = predicted_haplotypes[sampled_positions].to_numpy()
    # true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()[:, sampled_positions]
    
    # true_haplotypes_np = true_haplotypes.to_numpy()[:, sampled_positions]
    true_haplotypes_np = true_haplotypes.to_numpy()

    # vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_np, true_haplotypes_np)
    # accuracy, _ = calculate_accuracy(predicted_haplotypes_np, true_haplotypes_np)
    # mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_np, true_haplotypes_np)
    # mec_ = mec(predicted_haplotypes_np, fragment_model.fragment_list)
    vector_error_rate_count, _, _ = compute_vector_error_rate_with_missing_positions(true_haplotypes_np, predicted_haplotypes_np_count)
    vector_error_rate_evidence, _, _ = compute_vector_error_rate_with_missing_positions(true_haplotypes_np, predicted_haplotypes_np_evidence)
    vector_error_rate_mv, _, _ = compute_vector_error_rate_with_missing_positions(true_haplotypes_np, predicted_haplotypes_np_mv)
    
    results_name = 'FFBS_{}.pkl'.format(frag_file.split('.')[0])
    results = {'Edge count': {}, 'Evidence count': {}, 'M.V.': {}}
    # results['block_evaluation'] = block_info
    # results['components'] = components
    # results['n_blocks'] = len(components.keys())
    # results['average_block_size'] = block_info['average_block_size']
    
    results['Edge count']['length_phased'] = len(sampled_positions_count)
    results['Evidence count']['length_phased'] = len(sampled_positions_evidence)
    results['M.V.']['length_phased'] = len(sampled_positions_mv)
    results['Edge count']['vector_error_rate'] = vector_error_rate_count
    results['Evidence count']['vector_error_rate'] = vector_error_rate_evidence
    results['M.V.']['vector_error_rate'] = vector_error_rate_mv
    results['ffbs_acc'] = ffbs_acc
    results['ffbs_stats_dict'] = stats_dict
    # results['evaluation'] = {'vector_error_rate': vector_error_rate, 'ffbs_acc': ffbs_acc, 'ffbs_stats_dict': stats_dict}
    # results['predicted_haplotypes'] = predicted_haplotypes
    results['true_haplotypes'] = pd.read_csv(genotype_path).T
    results['forward_messages'] = forward_messages
    results['transitions_dict'] = transitions_dict
    results['transitions_dict_extra'] = transitions_dict_extra
    results['emission_dict'] = emission_dict
    results['assignment_dict'] = assignment_dict
    results['samples'] = samples
    # results['slices'] = slices
    # results['best_permutation'] = best_permutation
    results['fragment_list'] = fragment_model.fragment_list
    results['time'] = elapsed_time
    # print('Results:', results['evaluation'])

    with open(os.path.join(results_path, results_name), 'wb') as f:
        pickle.dump(results, f)

    print('Saved results in {}.'.format(os.path.join(results_path, results_name))) #, 'vector_error_rate', vector_error_rate, 'ffbs_acc', ffbs_acc)


def run_FFBS_multiple(inp):
    this_frag_path, this_fragment_coverage_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
    file_name = frag_file.split('.')[0]

    fragment_v_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_v_label_' + file_name + '.pkl')
    fragment_e_label_revered_path = os.path.join(this_reverse_maps_path, 'fg_e_label_' + file_name + '.pkl')
    quotient_v_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + file_name + '.pkl')
    quotient_e_label_revered_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + file_name + '.pkl')
    print('Working on:', os.path.join(this_frag_path, frag_file))

    class Args:
        def __init__(self):
            # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
            # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join(this_frag_path, frag_file)
            # self.data_path = frag_path
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

    start_time = time.time()

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
    fragment_model.construct(input_handler, config)

    edges_map_fragment = {}
    for k in fragment_model.e_label_reversed.keys():
        edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]

    # save fragment graph data:
    frag_graph_path = os.path.join(this_fragment_coverage_path, file_name + '.gt.gz')
    fragment_model.graph.save(frag_graph_path)

    frag_graph_plot_path = os.path.join(this_fragment_coverage_path, file_name + '.png')

    e_labels = fragment_model.graph.edge_properties["e_label"]
    v_labels = fragment_model.graph.vertex_properties["v_label"]
    gt.graph_draw(fragment_model.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=frag_graph_plot_path)

    with open(fragment_v_label_revered_path, "wb") as f:
        pickle.dump(fragment_model.v_label_reversed, f)

    with open(fragment_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_fragment, f)

    # create quotient graph
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct(input_handler, config)

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    quotient_g_v_label_reversed = quotient_g.v_label_reversed

    # save quotient graph
    quot_graph_path = os.path.join(this_quotient_coverage_path, file_name + '.gt.gz')
    quotient_g.graph.save(quot_graph_path)

    quotient_graph_plot_path = os.path.join(this_quotient_coverage_path, file_name + '.png')
    e_labels = quotient_g.graph.edge_properties["e_label"]
    v_labels = quotient_g.graph.vertex_properties["v_label"]
    gt.graph_draw(quotient_g.graph, output_size=(1000, 1000), vertex_text=v_labels, edge_text=e_labels, vertex_font_size=16,  
    edge_font_size=10, output=quotient_graph_plot_path)

    with open(quotient_v_label_revered_path, "wb") as f:
        pickle.dump(quotient_g.v_label_reversed, f)

    with open(quotient_e_label_revered_path, "wb") as f:
        pickle.dump(edges_map_quotient, f)

    transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model, config)
    emission_dict = emissions(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, _ =  assign_slices_and_interfaces(nodes, edges)

    assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, input_handler.data_path)

    forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)
    samples = sample_states_book_multiple_times(slices, edges, forward_messages, transitions_dict, n=20)
    # samples = sample_states_book(slices, edges, forward_messages, transitions_dict)

    ffbs_acc = evaulate_ffbs_acc_sample(genotype_path, samples, ploidy)
    # predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")
    predicted_haplotypes = predict_haplotypes_multiple_variants(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")

    end_time = time.time()

    elapsed_time = round(end_time - start_time, 2)

    true_haplotypes = pd.read_csv(genotype_path).T

    # block_info, components = get_block_info(quotient_g, predicted_haplotypes, true_haplotypes, fragment_model)

    sampled_positions = [c for c in predicted_haplotypes.columns.values if np.nan not in list(predicted_haplotypes[c].values)]

    predicted_haplotypes_np = predicted_haplotypes.to_numpy()
    # predicted_haplotypes_np = predicted_haplotypes[sampled_positions].to_numpy()
    # true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()[:, sampled_positions]
    
    # true_haplotypes_np = true_haplotypes.to_numpy()[:, sampled_positions]
    true_haplotypes_np = true_haplotypes.to_numpy()

    # vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_np, true_haplotypes_np)
    # accuracy, _ = calculate_accuracy(predicted_haplotypes_np, true_haplotypes_np)
    # mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_np, true_haplotypes_np)
    # mec_ = mec(predicted_haplotypes_np, fragment_model.fragment_list)
    vector_error_rate, _, _ = compute_vector_error_rate_with_missing_positions(true_haplotypes_np, predicted_haplotypes_np)
    results_name = 'FFBS_{}.pkl'.format(frag_file.split('.')[0])
    results = {}
    # results['block_evaluation'] = block_info
    # results['components'] = components
    # results['n_blocks'] = len(components.keys())
    # results['average_block_size'] = block_info['average_block_size']
    results['length_phased'] = len(sampled_positions)
    results['evaluation'] = {'vector_error_rate': vector_error_rate, 'ffbs_acc': ffbs_acc}
    results['predicted_haplotypes'] = predicted_haplotypes
    results['true_haplotypes'] = pd.read_csv(genotype_path).T
    results['forward_messages'] = forward_messages
    results['transitions_dict'] = transitions_dict
    results['transitions_dict_extra'] = transitions_dict_extra
    results['emission_dict'] = emission_dict
    results['assignment_dict'] = assignment_dict
    results['samples'] = samples
    # results['slices'] = slices
    # results['best_permutation'] = best_permutation
    results['fragment_list'] = fragment_model.fragment_list
    results['time'] = elapsed_time
    # print('Results:', results['evaluation'])

    with open(os.path.join(results_path, results_name), 'wb') as f:
        pickle.dump(results, f)

    print('Saved results in {}.'.format(os.path.join(results_path, results_name)), 'vector_error_rate', vector_error_rate, 'ffbs_acc', ffbs_acc)
    

def compute_pairwise_samples_accuracy(inp):
    this_frag_path, frag_file, ploidy, genotype_path, results_path = inp
    # file_name = frag_file.split('.')[0]
    print('Working on:', os.path.join(this_frag_path, frag_file))

    class Args:
        def __init__(self):
            # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
            # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join(this_frag_path, frag_file)
            # self.data_path = frag_path
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

    edges_map_fragment = {}
    for k in fragment_model.e_label_reversed.keys():
        edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]

    # create quotient graph
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct(input_handler, config)

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    quotient_g_v_label_reversed = quotient_g.v_label_reversed

    transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model, config)
    emission_dict = emissions(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, _ =  assign_slices_and_interfaces(nodes, edges)

    assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, input_handler.data_path)

    forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)
    samples_ffbs = sample_states_book(slices, edges, forward_messages, transitions_dict)
    samples_likelihood = sample_states_max_likelihood(slices, nodes, quotient_g, quotient_g_v_label_reversed)
    
    ffbs_acc1 = evaulate_ffbs_acc_sample(genotype_path, samples_ffbs, ploidy)
    likelihood_acc2 = evaulate_ffbs_acc_sample(genotype_path, samples_likelihood, ploidy)

    mn = create_markov_net(quotient_g.graph, transitions=transitions_dict) # Creates mrftools TorchMarkovNetwork
    torch_bp = mrftools.TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=False, var_on=False)
    torch_bp.infer(display="full")
    torch_bp.load_beliefs()
    update_entropies(quotient_g.graph, torch_bp)
    transitions_dict = update_transitions(quotient_g.graph, transitions_dict, torch_bp)
    # print("Transitions Updated")
    vertices_dict = update_vertices(quotient_g_v_label_reversed, emission_dict, slices, torch_bp)
    highest_lbp_samples = sample_states_max_lbp(slices, vertices_dict)
    highest_lbp_acc3 = evaulate_ffbs_acc_sample(genotype_path, highest_lbp_samples, ploidy)

    lbp_samples_no_forwarding = sample_states_book(slices, edges, vertices_dict, transitions_dict)
    lbp_no_forwarding_ffbs_acc4 = evaulate_ffbs_acc_sample(genotype_path, lbp_samples_no_forwarding, ploidy)

    forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, input_handler.data_path)

    lbp_ffbs_samples = sample_states_book(slices, edges, forward_messages, transitions_dict)
    lbp_ffbs_acc5 = evaulate_ffbs_acc_sample(genotype_path, lbp_ffbs_samples, ploidy)

    results_name = 'samples_accuracy_{}.pkl'.format(frag_file.split('.')[0])
    results = {'FFBS': ffbs_acc1, 'Likelihood': likelihood_acc2, 'Beliefs': highest_lbp_acc3, 'LBP Samples': lbp_no_forwarding_ffbs_acc4, 'LBP + FFBS': lbp_ffbs_acc5}


    with open(os.path.join(results_path, results_name), 'wb') as f:
        pickle.dump(results, f)

    print('Saved results in {}.'.format(os.path.join(results_path, results_name)))
    

def run_pHapcompass_from_input(input_file):
    with open(input_file, "rb") as f:
        inp = pickle.load(f)

    # run_FFBS_quotient_likelihood(inp)
    # run_LBP_FFBS_Single(inp)
    # compute_pairwise_samples_accuracy(inp)
    # run_FFBS_multiple(inp)
    run_FFBS_single_multi_samles(inp)

if __name__ == '__main__':

    # simulate_na12878()
    # simulate_awri()

    if len(sys.argv) != 2:
        print("Usage: python3 simulator_paper.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    run_pHapcompass_from_input(input_file)
