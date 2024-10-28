import os
import pandas as pd
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import graph_tool.all as gt
import numpy as np
import pickle
from algorithm.chordal_contraction import *
from utils.utils import *
from data.configuration import *
from data.input_handler import InputHandler
from models.fragment_graph import FragmentGraph
from evaluation.metrics import *


def generate_ffbs_input():
    inputs = []
    contigs = ['Contig1_k3']
    # contigs = ['Contig1_k4', 'Contig1_k5', 'Contig1_k6']
    # contigs = ['Contig2_k3']
    # simulated_data_path = '/home/mok23003/BML/HaplOrbit/simulated_data'
    # main_path = '/home/mok23003/BML/HaplOrbit/simulated_data_graphs/'
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_graphs/'
    # '/home/mok23003/BML/HaplOrbit/simulated_data_graphs/quotient_graphs/Contig1_k3/c2'
    for cont in contigs:
        cont_path = os.path.join(main_path, 'quotient_graphs', cont)
        coverges = sorted([d for d in os.listdir(cont_path) if os.path.isdir(os.path.join(cont_path, d))])
        coverges = ['c6']
        for coverage in coverges:
            ploidy = int(cont.split('_k')[1])
            this_cont_coverage_path = os.path.join(cont_path, coverage)
            this_save_path = os.path.join(main_path, 'ffbs', cont, coverage)

            if not os.path.exists(this_save_path):
                os.makedirs(this_save_path)
            existing_files = [ff for ff in os.listdir(this_save_path) if '.pkl.gz' in ff]

            graph_files = [f for f in os.listdir(this_cont_coverage_path) if '.gt.gz' in f]

            frag_files = [f.split('.gt.gz')[0] + '.frag.txt' for f in graph_files if not f.split('.gt.gz')[0] + '.pkl.gz' in existing_files]
            for frag_file in frag_files:
                this_input = [main_path, cont, coverage, frag_file, this_save_path, ploidy]
                inputs.append(this_input)
    return inputs


def retrive_graphs(inp):
    main_path, contig, coverage, frag_file, this_save_path, ploidy = inp
    genotype_path = os.path.join('/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data/', contig, 'real_haps_' + contig.lower() +  '.txt')
    # ffbs_path = os.path.join(this_save_path, frag_file.split('.frag.txt')[0] + '.pkl.gz')
    # print('Working on', ffbs_path)
    
    class Args:
        def __init__(self):
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = os.path.join('/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data/', contig, coverage, frag_file)
            # self.data_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/c2/ART_90.frag.txt'
            self.bam_path = 'example/example.bam'
            self.genotype_path = genotype_path
            self.ploidy = ploidy
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
    
    
    quotient_graph, v_label_reversed, e_label_reversed = read_quotient_graph(main_path, contig, coverage, frag_file)
    # fragment_graph, v_label_reversed_fragment, e_label_reversed_fragment = read_fragment_graph(main_path, contig, coverage, frag_file)
    return quotient_graph, v_label_reversed, e_label_reversed, fragment_model, input_handler, config, ploidy
    

def run_ffbs():
    inputs = generate_ffbs_input()
    inp = inputs[0]
    quotient_graph, v_label_reversed, e_label_reversed, fragment_model, input_handler, config, ploidy = retrive_graphs(inp)


    
    snp_list = []
    for vvv in quotient_graph.vertex_properties["v_label"]:
        snp_list += [int(v) for v in vvv.split('-')]
    
    min_snp = min(snp_list)
    max_snp = max(snp_list)
    
    all_vertices = [vvv for vvv in quotient_graph.vertex_properties["v_label"]]


    min_vertices = []
    max_vertices = []
    for vvv in quotient_graph.vertex_properties["v_label"]:
        if min_snp in [int(v) for v in vvv.split('-')]:
            min_vertices.append(vvv)
        if max_snp in [int(v) for v in vvv.split('-')]:
            max_vertices.append(vvv)

    source_node = min_vertices[0]
    target_node = max_vertices[0]
    target_node = '149-182'
    source = v_label_reversed[source_node]
    target = v_label_reversed[target_node]


    
    mst_graph, non_mst_graph, tree = get_minimum_spanning_tree(quotient_graph)
    mst_graph = gt.Graph(mst_graph, prune=True)

    path_vertices, path_edges = gt.shortest_path(mst_graph, source, target)
    alpha = forward_sum(quotient_graph, path_vertices, path_edges, ploidy, fragment_model, config)
    seq_len = len(path_vertices)

    sampled_states = ffbs(alpha, seq_len, path_vertices, path_edges, quotient_graph, ploidy, fragment_model, config)


    # for v_id in range(len(path_vertices)):
    #     vertex = path_vertices[v_id]
    #     print(quotient_graph.vertex_properties["v_label"][vertex], sampled_states[v_id])


    phasings = []

    vertex = path_vertices[0]
    source_label = quotient_graph.vertex_properties["v_label"][vertex]
    source_sampled_phasing = str_2_phas([sampled_states[0]], ploidy)
    # cites = [int(v) for v in source_label.split('-')]

    for v_id in range(1, len(path_vertices)):
        # if v_id == 11:
        #     stop
        # vertex = path_vertices[v_id]
        next_vertex = path_vertices[v_id]
        # print(quotient_graph.vertex_properties["v_label"][vertex], sampled_states[v_id])
        # source_label = quotient_graph.vertex_properties["v_label"][vertex]
        target_label = quotient_graph.vertex_properties["v_label"][next_vertex]
        # common_ff, common_sf = find_common_element_and_index(source_label, target_label)
        common_ff, common_sf = find_common_element_and_index_generalized(source_label, target_label)
        # source_sampled_phasing = str_2_phas_1(sampled_states[v_id], ploidy)
        target_sampled_phasing = str_2_phas_1(sampled_states[v_id], ploidy)
        phasings = []
        
        for sp in source_sampled_phasing:
            
        # for sp in source_sampled_phasing:
            # print(sp)
            matched_phasings = find_phasings_matches_generalized(sp, target_sampled_phasing, 
                                                                    common_ff, common_sf, 
                                                                    source_label, target_label)
            target_positions = [int(tl) for tl in target_label.split('-')]
            # all([tl in [int(tl) for tl in source_label.split('-')] for tl in target_positions])
            if len(matched_phasings) == 0 and all([tl in [int(tl) for tl in source_label.split('-')] for tl in target_positions]):
                # print(target_positions)
                continue
            else:
                # phasings = []
                sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)
                matched_phasings_str = list(set([phas_2_str(pm) for pm in sorted_phasings]))
                # sorted_phasings = []
                # for mtx in matched_phasings:
                #     sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                #     sorted_phasings.append(sorted_matrix)
                # matched_phasings_str = list(set([phas_2_str(pm) for pm in sorted_phasings]))
                
                phasings += matched_phasings_str

        phasings = list(set(phasings))
        
        # print(v_id, len(phasings), len(phasings[0])/ploidy)
        current_vertex = next_vertex
        source_label = '-'.join([str(s) for s in sorted(set([int(i) for i in source_label.split('-')] + [int(i) for i in target_label.split('-')]))])
        if len(phasings) != 0:
            source_sampled_phasing = str_2_phas(phasings, ploidy)
        # print(v_id, len(source_sampled_phasing), source_label, target_label)
        print(v_id, len(source_sampled_phasing), source_sampled_phasing[0].shape[1], source_label, target_label)
        # else:
        #     source_sampled_phasing = source_sampled_phasing
        # cites = [int(v) for v in source_label.split('-')]
        # input_handler.get_genotype_positions([96, 101])
        # np.sum(source_sampled_phasing[0], axis=1)

    dfs, gen_df = prepare_data_for_eval(source_label, source_sampled_phasing)
    true_hap_np = gen_df.to_numpy()
    true_hap_sorted = true_hap_np[np.argsort([''.join(map(str, row)) for row in true_hap_np])]
    constructed_haps = [df.to_numpy() for df in dfs]

    mec_scores = []
    vector_error_rates = []
    correct_phasing_rates = []
    perfect_solution_rates = []
    
    for haplotypes in constructed_haps:
        mec_score = calculate_mec_score(true_hap_sorted, haplotypes)
        vector_error_rate = calculate_vector_error_rate(haplotypes, true_hap_sorted)
        correct_phasing_rate = calculate_correct_phasing_rate(haplotypes, true_hap_sorted)
        perfect_solution_rate = calculate_perfect_solution_rate(haplotypes, true_hap_sorted)
        mec_scores.append(mec_score)
        vector_error_rates.append(vector_error_rate)
        correct_phasing_rates.append(correct_phasing_rate)
        perfect_solution_rates.append(perfect_solution_rate)


    
    categories = ['vector error rate', 'correct phasing rate', 'perfect solution rate']
    results_df = pd.DataFrame({
        'Value': vector_error_rates + correct_phasing_rates + perfect_solution_rates,
        'Category': ([categories[0]] * len(vector_error_rates)) +
                    ([categories[1]] * len(correct_phasing_rates)) +
                    ([categories[2]] * len(perfect_solution_rates))
    })

    # Create the box plot
    plt.figure(figsize=(8, 6))
    sns.boxplot(x='Category', y='Value', data=results_df)

    # Customize plot
    plt.title('Results on Contig1_k3')
    plt.show()






def prepare_data_for_eval(source_label, source_sampled_phasing, input_handler):
    columns = source_label.split('-')
    positions = [int(c) for c in columns]
    sorted_phasings = []
    for mtx in source_sampled_phasing:
        sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
        sorted_phasings.append(sorted_matrix)

    dfs = [pd.DataFrame(ssp, columns=columns) for ssp in sorted_phasings]
    gen_df = input_handler.get_haplotype()
    gen_df = gen_df.loc[positions, :]
    gen_df = gen_df.T

    return dfs, gen_df


def forward_sum(quotient_graph, path_vertices, path_edges, ploidy, fragment_model, config):
    alpha = {vi: np.zeros(len(quotient_graph.vertex_properties["v_weights"][v]['weight'].keys())) for vi, v in enumerate(path_vertices)}
    alpha[0] = np.array(list(quotient_graph.vertex_properties["v_weights"][path_vertices[0]]['weight'].values()))
    alpha[0] = torch.tensor(alpha[0] / np.sum(alpha[0]))

    seq_len = len(path_vertices)
    for sl in range(1, seq_len):
        # print(sl)
        source = path_vertices[sl-1]
        target = path_vertices[sl]
        source_weights = quotient_graph.vertex_properties["v_weights"][source]['weight']
        target_weights = quotient_graph.vertex_properties["v_weights"][target]['weight']
        source_label = quotient_graph.vertex_properties["v_label"][source]
        target_label = quotient_graph.vertex_properties["v_label"][target]
        e_label = quotient_graph.edge_properties["e_label"][path_edges[sl-1]]
        e_weights = quotient_graph.edge_properties["e_weights"][path_edges[sl-1]]['weight']

        common_ff, common_sf = find_common_element_and_index(source_label, target_label)
        source_phasings = list(source_weights.keys())
        target_phasings = list(target_weights.keys())
        transitions_dict = {'source': source_phasings, 'target': target_phasings}
        transitions_mtx = np.zeros((len(source_phasings), len(target_phasings)))
        for i, ffstr in enumerate(source_phasings):
            for j, sfstr in enumerate(target_phasings):
                
                matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, ploidy), str_2_phas_1(sfstr, ploidy), common_ff, common_sf, source_label, target_label)
                sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)
                
                matched_phasings_str = list(set([phas_2_str(pm) for pm in sorted_phasings]))
                poss = sorted(list(set([int(ss) for ss in source_label.split('-')] + [int(tt) for tt in target_label.split('-')])))
                match_reads = get_matching_reads_for_positions([int(i) for i in poss], fragment_model.fragment_list)
                wei = 0
                for phas in matched_phasings_str:
                    for indc, this_po, obs in match_reads:
                        # print(indc, this_po, obs)
                        # for obs in all_obs:
                        #     obs_np = np.array([int(po) for po in obs])
                        #     weights[phas_2_str(phas)] += compute_likelihood(obs_np, phas, error_rate)
                        wei += compute_likelihood_generalized_plus(np.array(obs), str_2_phas_1(phas, ploidy), indc, list(range(len(indc))), 
                                                                   config.error_rate)


                # this_weight = np.sum([e_weights[pm] for pm in matched_phasings_str if pm in e_weights.keys()])
                transitions_mtx[i, j] = wei

        transitions_mtx = transitions_mtx / transitions_mtx.sum(axis=1, keepdims=True)

        dp_eq = alpha[sl-1][:, np.newaxis] * transitions_mtx
        alpha[sl] = torch.tensor(dp_eq.sum(axis=0))
    return alpha


def backward_sum(quotient_graph, spath_vertices, spath_edges):
    beta = {vi: np.zeros(len(quotient_graph.vertex_properties["v_weights"][v]['weight'].keys())) for vi, v in enumerate(spath_vertices)}
    # last_vertex = spath_vertices[-1]
    beta[seq_len - 1] = torch.ones(len(quotient_graph.vertex_properties["v_weights"][last_vertex]['weight'].keys()))

    for sl in range(0, seq_len-1)[::-1]:
        source = spath_vertices[sl]
        target = spath_vertices[sl + 1]
        source_weights = quotient_graph.vertex_properties["v_weights"][source]['weight']
        target_weights = quotient_graph.vertex_properties["v_weights"][target]['weight']
        source_label = quotient_graph.vertex_properties["v_label"][source]
        target_label = quotient_graph.vertex_properties["v_label"][target]
        e_label = quotient_graph.edge_properties["e_label"][spath_edges[sl]]
        e_weights = quotient_graph.edge_properties["e_weights"][spath_edges[sl]]['weight']
        for i, ffstr in enumerate(source_phasings):
            for j, sfstr in enumerate(target_phasings):
                matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, 3), str_2_phas_1(sfstr, 3), common_ff, common_sf, source_label, target_label)
                sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)

                matched_phasings_str = [phas_2_str(pm) for pm in sorted_phasings] 
                this_weight = np.sum([e_weights[pm] for pm in matched_phasings_str if pm in e_weights.keys()])
                transitions_mtx[i, j] = this_weight

        transitions_mtx = transitions_mtx / transitions_mtx.sum(axis=1, keepdims=True)

        dp_eq = torch.tensor(transitions_mtx) * beta[sl + 1][np.newaxis, :]
        beta[sl] = dp_eq.sum(axis=1)


def ffbs(alpha, seq_len, path_vertices, path_edges, quotient_graph, ploidy, fragment_model, config):
    # Initialize a list to store the sampled states
    sampled_states = [None] * seq_len

    # Step 1: Sample the last state x_T
    last_vertex = path_vertices[seq_len - 1]
    target_weights = quotient_graph.vertex_properties["v_weights"][last_vertex]['weight']
    target_phasings = list(target_weights.keys())

    # Convert alpha_T to probabilities
    probs = alpha[seq_len - 1].numpy()
    probs = probs / probs.sum()

    # Sample the last state based on alpha_T
    x_T_index = np.random.choice(len(target_phasings), p=probs)
    x_T = target_phasings[x_T_index]
    sampled_states[seq_len - 1] = x_T

    # Step 2: Backward Sampling for t from T-1 down to 0
    for t in range(seq_len - 2, -1, -1):
        source = path_vertices[t]
        target = path_vertices[t + 1]
        source_weights = quotient_graph.vertex_properties["v_weights"][source]['weight']
        target_weights = quotient_graph.vertex_properties["v_weights"][target]['weight']
        source_label = quotient_graph.vertex_properties["v_label"][source]
        target_label = quotient_graph.vertex_properties["v_label"][target]
        
        e_label = quotient_graph.edge_properties["e_label"][path_edges[t]]
        e_weights = quotient_graph.edge_properties["e_weights"][path_edges[t]]['weight']
        
        # Find common elements and indices
        common_ff, common_sf = find_common_element_and_index(source_label, target_label)
        source_phasings = list(source_weights.keys())
        target_phasings = list(target_weights.keys())
        transitions_mtx = np.zeros((len(source_phasings), len(target_phasings)))
        print(t, source_label, source_phasings, target_label, target_label)
        # Recompute the transition matrix at time t
        for i, ffstr in enumerate(source_phasings):
            for j, sfstr in enumerate(target_phasings):
                matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, ploidy), str_2_phas_1(sfstr, ploidy), common_ff, common_sf, source_label, target_label)
                sorted_phasings = []

                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)
                
                matched_phasings_str = list(set([phas_2_str(pm) for pm in sorted_phasings]))
                poss = sorted(list(set([int(ss) for ss in source_label.split('-')] + [int(tt) for tt in target_label.split('-')])))
                match_reads = get_matching_reads_for_positions([int(i) for i in poss], fragment_model.fragment_list)
                wei = 0
                for phas in matched_phasings_str:
                    for indc, this_po, obs in match_reads:
                        # print(indc, this_po, obs)
                        # for obs in all_obs:
                        #     obs_np = np.array([int(po) for po in obs])
                        #     weights[phas_2_str(phas)] += compute_likelihood(obs_np, phas, error_rate)
                        wei += compute_likelihood_generalized_plus(np.array(obs), str_2_phas_1(phas, ploidy), indc, list(range(len(indc))), 
                                                                   config.error_rate)


                # this_weight = np.sum([e_weights[pm] for pm in matched_phasings_str if pm in e_weights.keys()])
                transitions_mtx[i, j] = wei

        
        # Normalize the transition matrix
        transitions_mtx = transitions_mtx / transitions_mtx.sum(axis=1, keepdims=True)
        
        # Get the index of the next sampled state
        x_t1 = sampled_states[t + 1]
        x_t1_index = target_phasings.index(x_t1)
        
        # Compute the probabilities for the current state
        alpha_t = alpha[t].numpy()
        probs = alpha_t * transitions_mtx[:, x_t1_index]
        probs = probs / probs.sum()
        
        # Sample the current state
        x_t_index = np.random.choice(len(source_phasings), p=probs)
        x_t = source_phasings[x_t_index]
        sampled_states[t] = x_t
    return sampled_states


