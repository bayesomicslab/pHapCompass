import os
import pandas as pd
import random
import pysam
import sys
import numpy as np
# sys.path.append('/home/mok23003/BML/HaplOrbit')
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import pickle
# from ..generate_simulated_graphs import generate_quotient_graph
from data.input_handler import InputHandler
from data.configuration import Configuration
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from multiprocessing import Pool
import matplotlib.pyplot as plt
import graph_tool.all as gt
from FFBS.FFBS_quotient_graph import *
import time


def generate_scripts():
    sh_path = '/labs/Aguiar/pHapCompass/scripts/run_method'
    chromosomes = ['CP142453.1', 'CP142454.1', 'CP142455.1', 'CP142456.1', 'CP142457.1', 'CP142458.1', 'CP142459.1', 'CP142460.1', 'CP142461.1', 'CP142462.1'] # 'CP142451.1', 'CP142452.1', 
    for chrom in chromosomes:
        to_print = '#!/bin/bash\n'\
                    '#BATCH --job-name=like{}\n'\
                    '#SBATCH -N 1\n'\
                    '#SBATCH -n 1\n'\
                    '#SBATCH -c 1\n'\
                    '#SBATCH -x xanadu-[46,49,25,10]\n'\
                    '#SBATCH --partition=general\n'\
                    '#SBATCH --qos=general\n'\
                    '#SBATCH --mail-type=END\n'\
                    '#SBATCH --mem=50G\n'\
                    '#SBATCH --mail-user=marjan.hosseini@uconn.edu\n'\
                    '#SBATCH -o like{}.out\n'\
                    '#SBATCH -e like{}.err\n\n'.format(chrom[6:8], chrom[6:8], chrom[6:8])
        to_print += 'echo `hostname`\n\n'
        command = '/home/FCAM/mhosseini/miniforge3/envs/hap/bin/python3 /labs/Aguiar/pHapCompass/run_sassafras.py  --chrom "{}"'.format(chrom)
        to_print += command
        with open(os.path.join(sh_path, 'sass{}.sh'.format(chrom[6:8])), 'w') as f:
            f.write(to_print)


def run_pHapcompass_sassafras_chr(chrom):
    # chrom = 'CP142451.1'
    frag_path = '/labs/Aguiar/pHapCompass/datasets/sassafras/frag/{}_RagTag.frag'.format(chrom)
    ploidy= 4
    genotype_path = '/labs/Aguiar/pHapCompass/datasets/sassafras/vcf/haplotype_filtered_{}_RagTag.csv'.format(chrom)
    main_results_path = '/labs/Aguiar/pHapCompass/datasets/sassafras/results'
    this_results_path = os.path.join(main_results_path, chrom)
    if not os.path.exists(this_results_path):
        os.makedirs(this_results_path)

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

    start_time = time.time()
    # Initialize classes with parsed arguments
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)

    fragment_model.construct(input_handler, config)
    fragment_model.graph.save(os.path.join(this_results_path, 'fragment_graph_{}.gt'.format(chrom)))

    with open(os.path.join(this_results_path, 'fragment_graph_map_v_label_{}.pkl').format(chrom), "wb") as f:
        pickle.dump(fragment_model.v_label_reversed, f)
    print('Fragment Graph constructed.')

    edges_map_fragment = {}
    for k in fragment_model.e_label_reversed.keys():
        edges_map_fragment[k] = [int(fragment_model.e_label_reversed[k].source()), int(fragment_model.e_label_reversed[k].target())]

    with open(os.path.join(this_results_path, 'fragment_graph_map_e_label_{}.pkl').format(chrom), "wb") as f:
        pickle.dump(edges_map_fragment, f)

    # create quotient graph
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct(input_handler, config)
    quotient_g.graph.save(os.path.join(this_results_path, 'quotient_graph_{}.gt'.format(chrom)))

    with open(os.path.join(this_results_path, 'quotient_graph_map_v_label_{}.pkl').format(chrom), "wb") as f:
        pickle.dump(quotient_g.v_label_reversed, f)


    quotient_g_v_label_reversed = quotient_g.v_label_reversed

    edges_map_quotient = {}
    for k in quotient_g.e_label_reversed.keys():
        edges_map_quotient[k] = [int(quotient_g.e_label_reversed[k].source()), int(quotient_g.e_label_reversed[k].target())]

    with open(os.path.join(this_results_path, 'quotient_graph_map_e_label_{}.pkl').format(chrom), "wb") as f:
        pickle.dump(edges_map_quotient, f)


    transitions_dict, transitions_dict_extra = transition_matrices(quotient_g, edges_map_quotient, ploidy, fragment_model, config)
    emission_dict = emissions(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, _ =  assign_slices_and_interfaces(nodes, edges)

    assignment_dict = assign_evidence_to_states_and_transitions(nodes, edges, frag_path)

    forward_messages = compute_forward_messages(slices, edges, assignment_dict, emission_dict, transitions_dict, frag_path)

    samples = sample_states_book(slices, edges, forward_messages, transitions_dict)

    predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config, priority="probabilities")
    end_time = time.time()
    elapsed_time = round(end_time - start_time, 2)

    sampled_positions = [c for c in predicted_haplotypes.columns.values if np.nan not in list(predicted_haplotypes[c].values)]

    predicted_haplotypes_np = predicted_haplotypes[sampled_positions].to_numpy()
    true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()[:, sampled_positions]


    vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_np, true_haplotypes)
    accuracy, best_permutation = calculate_accuracy(predicted_haplotypes_np, true_haplotypes)
    mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_np, true_haplotypes)
    mec_ = mec(predicted_haplotypes_np, fragment_model.fragment_list)
    print('Vector Error Rate:', vector_error_rate, '\nVector Error:', vector_error, '\nAccuracy:', accuracy, '\nMismatch Error:', mismatch_error, '\nMEC:', mec_)
    
    results_name = 'likelihood_FFBS_{}.pkl'.format(chrom)
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

    with open(os.path.join(this_results_path, results_name), 'wb') as f:
        pickle.dump(results, f)

    print('Saved results in {}.'.format(os.path.join(this_results_path, results_name)), 'vector_error_rate', vector_error_rate, 'vector_error', vector_error, 'mismatch_error', mismatch_error, 'mec', mec_)
    
if __name__ == '__main__':
    # chrom = sys.argv[1]
    
    parser = argparse.ArgumentParser(description="Haplotype Assembly for sassafras")
    parser.add_argument("-c", "--chrom", type=str, required=False, help="chromosome name", default='None')

    mainargs = parser.parse_args()
    chrom = mainargs.chrom
    run_pHapcompass_sassafras_chr(chrom)