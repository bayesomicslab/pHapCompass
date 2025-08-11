import graph_tool.all as gt
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from data.input_handler import InputHandler
from data.configuration import Configuration
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
import matplotlib.pyplot as plt
import torch
import numpy as np
import time
import gzip
import pickle
import scipy
from algorithm.FFBS import *
from evaluation.evaluation import *
sys.path.append(os.path.abspath("/mnt/research/aguiarlab/proj/HaplOrbit/algorithm/wSMBP"))
from belief_propagator import TorchMatrixBeliefPropagator
from factor_graph import FactorGraph

# Create an mrftools Markov Network
def create_markov_net(q_graph, transitions):
    # initialize markov net object
    markov_net = FactorGraph(is_cuda=False)
    for vertex in q_graph.iter_vertices():
        # label = q_graph.vp['v_label'][vertex]
        weight_dict = q_graph.vp['v_weights'][vertex]['weight']
        weights = list(weight_dict.values())
        markov_net.set_unary_factor(vertex, torch.tensor(weights))

    # print("Vertices encoded")

    added_edges = []
    for s, t in q_graph.iter_edges():
        edge = (s, t)
        # s_size = len(q_graph.vp['e_weights'][s]['weight'].values())
        s_label = q_graph.vp['v_label'][s]
        # t_size = len(q_graph.vp['e_weights'][t]['weight'].values())
        t_label = q_graph.vp['v_label'][t]
        if (sorted((s_label, t_label))) in added_edges:
            continue
        e_label = '--'.join([s_label, t_label])
        weight_matrix = transitions[e_label]
        markov_net.set_edge_factor((s, t), torch.from_numpy(weight_matrix))
        added_edges.append(sorted((s_label, t_label)))

    # print("Edge factors encoded")

    return markov_net

# Update edge entropies based on the computed BP beliefs
def update_entropies(q_graph, bp):
    for s, t in q_graph.iter_edges():
        transition_matrix = bp.pair_beliefs[(s, t)]
        entropy_array = scipy.stats.entropy(transition_matrix)
        entropy = np.mean(entropy_array)
        q_graph.ep['e_weights'][(s, t)]['entropy'] = entropy
    return


# update the transition matrix with computed BP beliefs
def update_transitions(q_graph, transition_dict, bp):
    added_edges = []
    for s, t in q_graph.iter_edges():
        s_label = q_graph.vp['v_label'][s]
        # t_size = len(q_graph.vp['e_weights'][t]['weight'].values())
        t_label = q_graph.vp['v_label'][t]
        if (sorted((s_label, t_label))) in added_edges:
            continue
        e_label = '--'.join([s_label, t_label])
        probs = np.exp(bp.pair_beliefs[(s, t)].numpy())
        probs = np.clip(probs, 0, 1)  # Ensure probabilities are in [0, 1]
        if np.isnan(probs).any():
            # print(f"Warning: Probabilities for edge {e_label} contain NaN values. Setting to zero.")
            probs = np.nan_to_num(probs)
        # if any probs are negative or greater than 1, set to zero
        probs[probs < 0] = 0
        probs[probs > 1] = 0
        # if all probs are zero, make edge uniform
        if np.all(probs == 0):
            probs = np.ones_like(probs) / len(probs)
        # normalize the probabilities
        probs /= np.sum(probs)
        transition_dict[e_label] = probs
        added_edges.append(sorted((s_label, t_label)))
    return transition_dict

# update the vertex marginals according to slices with computed BP beliefs
def update_vertices(rev_label, emission_dict, slices, bp):
    vertices_dict = {}
    for v_slice in slices:
        vertices_dict[v_slice] = {}
        vertices = slices[v_slice]
        for vertex in vertices:
            vertices_dict[v_slice][vertex] = {}
            state_space = list(emission_dict[vertex].keys())
            i = 0
            for phase in state_space:
                prob = np.exp(bp.var_beliefs[rev_label[vertex]][i].numpy())
                if prob < 0 or prob > 1 or np.isnan(prob):
                    # print(f"Warning: Probability {prob} for vertex {vertex}, phase {phase} is out of bounds [0, 1]. Setting to 0.")
                    prob = 0.0
                vertices_dict[v_slice][vertex][phase] = prob
                i +=1
            probs = np.array(list(vertices_dict[v_slice][vertex].values()))
            if np.sum(probs) < 1e-10 or np.sum(probs) > 1 or np.isnan(probs).any():
                # print(f"Warning: Probabilities for vertex {vertex} sum to zero. Setting to uniform distribution.")
                probs = np.ones_like(probs) / len(probs)
                for phase in state_space:
                    vertices_dict[v_slice][vertex][phase] = probs[0]

    return vertices_dict


# Each of These Take in an MN and output MN w/ computed marginals ---------------------------------------------------

def torch_belief_propagation(markov_network, is_cuda=True):
    return

def pgmpy_belief_propagation(markov_network):
    return

def independent_belief_propagation(markov_network, is_cuda=True):
    return

# ------------------------------------------------------------------------------------------------------------------

# Take in Quotient Graph, Transitions Matrix, and maybe a Method, and output --- factor graph with computed marginals post BP?
def inference(q_graph, transitions, method="smbp"):
    cuda_avail = torch.cuda.is_available()
    if method == "smbp":
        mn = create_markov_net(q_graph, transitions=transitions)
        bp = torch_belief_propagation(mn, is_cuda=cuda_avail)
        # sample this graph with updated marginals
    elif method == "indep":
        mn = create_markov_net(q_graph, transitions=transitions)
        bp = independent_belief_propagation(mn, is_cuda=cuda_avail)
    elif method == "pgmpy":
        # some markov network creation function
        bp = pgmpy_belief_propagation(mn)
    return

def pipeline(ploidy, contig, data, coverage, sample, weighting=None):
    print(f"Beginning pipeline for ploidy {ploidy}, coverage {coverage}, sample {sample}, and weighting {str(weighting)}.")
    t0 = time.time()
    init_path = f"/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_{data}/contig_{contig}/ploidy_{ploidy}/"

    genotype_path = init_path + "haplotypes.csv"
    fgraph_path = init_path + f"cov_{coverage}/fgraph/{sample}.gt.gz"
    qgraph_path = init_path + f"cov_{coverage}/qgraph/{sample}.gt.gz"
    edge_map_path = init_path + f"cov_{coverage}/reverse_maps/qg_e_label_{sample}.pkl"
    vertex_reversed_path = init_path + f"cov_{coverage}/reverse_maps/qg_v_label_{sample}.pkl"
    results_path = init_path + f"cov_{coverage}/results_wSMBP/weighting_{str(weighting)}/"
    if not os.path.exists(results_path):
        os.makedirs(results_path, exist_ok=True)
    results_file = results_path + f"wSMBP_{sample}.pkl"

    with open(edge_map_path, 'rb') as f:
        edges_map_quotient = pickle.load(f)

    class Args:
        def __init__(self):
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = init_path + f'cov_{coverage}/frag/{sample}.frag'
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

    fragment_model = FragmentGraph(input_handler.data_path, genotype_path, ploidy=ploidy, alleles=[0, 1])
    fragment_model.construct(input_handler, config)

    quotient_g = gt.load_graph(qgraph_path)

    with open(vertex_reversed_path, 'rb') as f:
        quotient_g_v_label_reversed = pickle.load(f)

    transitions_dict, transitions_dict_extra = transition_matrices_v2(quotient_g, edges_map_quotient, ploidy=ploidy, config=config, fragment_model=fragment_model)
    emission_dict = emissions_v2(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, _ =  assign_slices_and_interfaces(nodes, edges)

    # Initialize Markov Network
    mn = create_markov_net(quotient_g, transitions=transitions_dict) # Creates mrftools TorchMarkovNetwork

    # print("Markov Network Created")
    # print("Unary matrix: ")
    # print(mn.unary_mat)
    # print("Edge factors: ")
    # print(mn.edge_pot_tensor)

    # Belief Propagation Start
    torch_bp = TorchMatrixBeliefPropagator(
        markov_net=mn, is_cuda=False, weighting=weighting
    )
    torch_bp.infer(display="iter")
    torch_bp.load_beliefs()

    transitions_dict = update_transitions(quotient_g, transitions_dict, torch_bp)

    vertices_dict = update_vertices(quotient_g_v_label_reversed, emission_dict, slices, torch_bp)

    samples = sample_states_book(slices, edges, vertices_dict, transitions_dict)
    ffbs_acc = evaulate_ffbs_acc_sample(genotype_path, samples, ploidy)


    predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config)
    # print(predicted_haplotypes)
    t1 = time.time()
    elapsed_time = t1 - t0

    true_haplotypes = pd.read_csv(genotype_path).T

    sampled_positions = [c for c in predicted_haplotypes.columns.values if np.nan not in list(predicted_haplotypes[c].values)]

    predicted_haplotypes_filtered = predicted_haplotypes[sampled_positions].to_numpy()
    true_haplotypes_filtered = true_haplotypes[sampled_positions].to_numpy()

    # block_info, components = get_block_info(quotient_g, predicted_haplotypes, true_haplotypes, fragment_model)

    predicted_haplotypes_filtered = np.array(predicted_haplotypes_filtered.astype(int))
    true_haplotypes_filtered = np.array(true_haplotypes_filtered.astype(int))

    predicted_haplotypes_numeric = predicted_haplotypes.apply(pd.to_numeric, errors='coerce')
    true_haplotypes_numeric = true_haplotypes.apply(pd.to_numeric, errors='coerce')

    predicted_haplotypes_np = predicted_haplotypes_numeric.to_numpy()
    true_haplotypes_np = true_haplotypes_numeric.to_numpy()

    # predicted_haplotypes_np = predicted_haplotypes_np.astype(int)
    # true_haplotypes_np = true_haplotypes_np.astype(int)

    # print(f"Predicted haplotypes np: {predicted_haplotypes_np}")
    # print(f"True haplotypes np: {true_haplotypes_np}")

    # try:
    #     frag_list = fragment_model.fragment_list
    #     original_indices = set()
    #     for i in range(0, len(frag_list), 2):
    #         for idx in frag_list[i]:
    #             original_indices.add(idx)
        
    #     sorted_indices = sorted(list(original_indices))
    #     index_mapping = {old_idx: new_idx for new_idx, old_idx in enumerate(sorted_indices)}

    #     for i in range(0, len(frag_list), 2):
    #         frag_list[i] = [index_mapping[idx] for idx in frag_list[i]]

    # except Exception as e:
    #     print(f"Error processing fragment list: {e}")
    #     frag_list = None

    vector_error_rate, vector_error, backtracking_steps, dp_table = compute_vector_error_rate(predicted_haplotypes_filtered, true_haplotypes_filtered)
    accuracy, _ = calculate_accuracy(predicted_haplotypes_filtered, true_haplotypes_filtered)
    mismatch_error, best_permutation = calculate_mismatch_error(predicted_haplotypes_filtered, true_haplotypes_filtered)

    try:
        mec_ = mec_full(predicted_haplotypes_np, true_haplotypes_np, fragment_model.fragment_list, probabalistic=True)
    except Exception as e:
        print(f"Error calculating MEC: {e}")
        mec_ = None
    results = {}
    # results['block_evaluation'] = block_info
    # results['components'] = components
    # results['n_blocks'] = len(components.keys())
    # results['average_block_size'] = block_info['average_block_size']
    results['length_phased'] = len(sampled_positions)
    results['evaluation'] = {'vector_error_rate': vector_error_rate, 'vector_error': vector_error, 'backtracking_steps': backtracking_steps, 
                             'dp_table': dp_table, 'accuracy': accuracy, 'mismatch_error': mismatch_error, 'mec': mec_, 'ffbs_acc': ffbs_acc}
    results['predicted_haplotypes'] = predicted_haplotypes
    results['true_haplotypes'] = pd.read_csv(genotype_path).T
    results['transitions_dict'] = transitions_dict
    results['transitions_dict_extra'] = transitions_dict_extra
    results['emission_dict'] = emission_dict
    results['samples'] = samples
    results['slices'] = slices
    results['best_permutation'] = best_permutation
    results['fragment_list'] = fragment_model.fragment_list
    results['time'] = elapsed_time

    with open(results_file, 'wb') as f:
        pickle.dump(results, f)
    try:
        print(f"Ploidy: {ploidy}, Contig: {contig}, Coverage: {coverage}, Sample: {sample}, Weighting: {str(weighting)} - MEC: {mec_:.2f}, Vector Error Rate: {vector_error_rate:.2f}, Pairwise Phasing Accuracy: {ffbs_acc:.2f}, Time: {elapsed_time:.2f} seconds")
    except Exception as e:
        print(f"Error printing results: {e}")
    return


if __name__ == "__main__":

    print("Starting")
    t0 = time.time()

    ploidy = 3
    contig = 100
    data = 'NA12878'
    coverage = 10
    sample = '00'

    init_path = f"/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_{data}/contig_{contig}/ploidy_{ploidy}/"

    frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
    genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
    ploidy=3
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
    args = Args()
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)

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

    transitions_dict, transitions_dict_extra = transition_matrices_v2(quotient_g, edges_map_quotient, ploidy=ploidy, config=config, fragment_model=fragment_model)
    emission_dict = emissions_v2(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    print("Transition Matrices Constructed")

    nodes = list(emission_dict.keys())
    edges = [(e.split('--')[0], e.split('--')[1]) for e in list(transitions_dict.keys())]

    slices, _ =  assign_slices_and_interfaces(nodes, edges)

    # Initialize Markov Network
    mn = create_markov_net(quotient_g, transitions=transitions_dict) # Creates mrftools TorchMarkovNetwork
    t1 = time.time()
    print("Markov Network Created")
    print(f"Time elapsed: {t1-t0}")

    # Belief Propagation Start
    t0 = time.time()
    print("Starting BP")
    torch_bp = mrftools.TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=False, var_on=False)
    print("BP initialized")
    torch_bp.infer(display="full")
    print("BP inference done")
    torch_bp.load_beliefs()

    t1 = time.time()
    print(f"Time elapsed: {t1-t0}")

    update_entropies(quotient_g, torch_bp)

    transitions_dict = update_transitions(quotient_g, transitions_dict, torch_bp)
    print("Transitions Updated")

    vertices_dict = update_vertices(quotient_g_v_label_reversed, emission_dict, slices, torch_bp)
    print("Vertices Updated")

    samples = sample_states_book(slices, edges, vertices_dict, transitions_dict)

    predicted_haplotypes = predict_haplotypes(nodes, edges, samples, ploidy, genotype_path, fragment_model, transitions_dict_extra, config)

    print(predicted_haplotypes)
    print("done")