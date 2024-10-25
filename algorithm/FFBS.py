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
    
    min_vertices = []
    max_vertices = []
    for vvv in quotient_graph.vertex_properties["v_label"]:
        if min_snp in [int(v) for v in vvv.split('-')]:
            min_vertices.append(vvv)
        if max_snp in [int(v) for v in vvv.split('-')]:
            max_vertices.append(vvv)

    source_node = min_vertices[0]
    target_node = max_vertices[0]
    source = v_label_reversed[source_node]
    target = v_label_reversed[target_node]

    mst_graph, non_mst_graph, tree = get_minimum_spanning_tree(quotient_graph)
    mst_graph = gt.Graph(mst_graph, prune=True)

    path_vertices, path_edges = gt.shortest_path(mst_graph, source, target)
    alpha = forward_sum(quotient_graph, path_vertices, path_edges, ploidy, fragment_model)
    seq_len = len(path_vertices)

    sampled_states = ffbs(alpha, seq_len, path_vertices, path_edges, quotient_graph, ploidy, fragment_model, config)
    revsered_sampled_states = sampled_states[::-1]

    connections = []
    for v_id in range(len(path_vertices) - 1):
        vertex = path_vertices[v_id]
        next_vertex = path_vertices[v_id + 1]
        print(quotient_graph.vertex_properties["v_label"][vertex], revsered_sampled_states[v_id])
        source_label = quotient_graph.vertex_properties["v_label"][vertex]
        target_label = quotient_graph.vertex_properties["v_label"][next_vertex]
        common_ff, common_sf = find_common_element_and_index(source_label, target_label)
        source_sampled_phasing = str_2_phas_1(revsered_sampled_states[v_id], ploidy)
        target_sampled_phasing = str_2_phas_1(revsered_sampled_states[v_id + 1], ploidy)
        matched_phasings = find_phasings_matches(source_sampled_phasing, target_sampled_phasing, common_ff, common_sf, source_label, target_label)

        input_handler.get_genotype_positions([3, 4, 6])
        fragment_model.v_label_reversed['3']
        fragment_model.v_label_reversed['4']
        fragment_model.v_label_reversed['6']

        fragment_model.graph.vertex_properties['v_weights'][fragment_model.graph.vertex(fragment_model.v_label_reversed['3'])]
        fragment_model.graph.vertex_properties['v_weights'][fragment_model.graph.vertex(fragment_model.v_label_reversed['4'])]
        fragment_model.graph.vertex_properties['v_weights'][fragment_model.graph.vertex(fragment_model.v_label_reversed['6'])]


        edge_34 = fragment_model.graph.edge(fragment_model.v_label_reversed['3'], fragment_model.v_label_reversed['4'])
        edge_36 = fragment_model.graph.edge(fragment_model.v_label_reversed['3'], fragment_model.v_label_reversed['6'])
        edge_46 = fragment_model.graph.edge(fragment_model.v_label_reversed['4'], fragment_model.v_label_reversed['6'])


        frag_weights34 = fragment_model.graph.edge_properties["e_weights"][edge_34]['weight'].keys()
        frag_weights36 = fragment_model.graph.edge_properties["e_weights"][edge_36]['weight'].keys()
        frag_weights46 = fragment_model.graph.edge_properties["e_weights"][edge_46]['weight'].keys()


        print(quotient_graph.vertex_properties["v_weights"][quotient_graph.vertex(v_label_reversed['3-4'])]['weight'].keys())
        print(quotient_graph.vertex_properties["v_weights"][quotient_graph.vertex(v_label_reversed['3-6'])]['weight'].keys())
        print(quotient_graph.vertex_properties["v_weights"][quotient_graph.vertex(v_label_reversed['4-6'])]['weight'].keys())



        q_edge_3436 = quotient_graph.edge(quotient_graph.vertex(v_label_reversed['3-4']), quotient_graph.vertex(v_label_reversed['3-6']))
        q_edge_3446 = quotient_graph.edge(quotient_graph.vertex(v_label_reversed['3-4']), quotient_graph.vertex(v_label_reversed['4-6']))
        q_edge_3646 = quotient_graph.edge(quotient_graph.vertex(v_label_reversed['3-6']), quotient_graph.vertex(v_label_reversed['4-6']))

        quotient_graph.vertex_properties["v_weights"][quotient_graph.vertex(v_label_reversed['3-4'])]['weight']
        quotient_graph.vertex_properties["v_weights"][quotient_graph.vertex(v_label_reversed['3-6'])]['weight']
        quotient_graph.vertex_properties["v_weights"][quotient_graph.vertex(v_label_reversed['4-6'])]['weight']


        print(quotient_graph.edge_properties["e_weights"][q_edge_3436]['weight'].keys())
        print(quotient_graph.edge_properties["e_weights"][q_edge_3446]['weight'].keys())
        print(quotient_graph.edge_properties["e_weights"][q_edge_3646]['weight'].keys())




        sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)
                
                matched_phasings_str = list(set([phas_2_str(pm) for pm in sorted_phasings]))




        phase = revsered_sampled_states[v_id]
        start, end = map(int, pos.split('-'))
        connections.append((start, end, phase))
    



def forward_sum(quotient_graph, path_vertices, path_edges, ploidy, fragment_model):
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


