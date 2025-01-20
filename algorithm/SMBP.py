import mrftools
import IndependentTorchMatrixBeliefPropagator
import graph_tool.all as gt
import pgmpy
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
from Simulator.simulator_paper import transition_matrices_v2, emissions_v2

# Create an mrftools Markov Network
def create_markov_net(q_graph, transitions):
    # initialize markov net object
    markov_net = mrftools.TorchMarkovNet(is_cuda=False, var_on=False)
    for vertex in q_graph.iter_vertices():
        # label = q_graph.vp['v_label'][vertex]
        weight_dict = q_graph.vp['e_weights'][vertex]['weight']
        weights = list(weight_dict.values())
        markov_net.set_unary_factor(vertex, torch.tensor(weights))  

    print("Vertices encoded")

    for s, t in q_graph.iter_edges():
        edge = (s, t)
        # s_size = len(q_graph.vp['e_weights'][s]['weight'].values())
        s_label = q_graph.vp['v_label'][s]
        # t_size = len(q_graph.vp['e_weights'][t]['weight'].values())
        t_label = q_graph.vp['v_label'][t]
        e_label = '--'.join(sorted([s_label, t_label]))        
        weight_matrix = transitions[e_label]['transitions']
        markov_net.set_edge_factor((s, t), torch.from_numpy(weight_matrix))

    print("Edge factors encoded")

    return markov_net

# Update edge entropies based on the computed BP beliefs
def update_entropies(q_graph, bp):
    for s, t in q_graph.iter_edges():
        transition_matrix = bp.pair_beliefs[(s, t)]
        entropy_array = scipy.stats.entropy(transition_matrix)
        entropy = np.mean(entropy_array)
        q_graph.ep['e_weights'][(s, t)]['entropy'] = entropy
    return

# Produce graph-tool edge index from given source-target nodes
def to_edge(source, target):
    return

# Given the assignments of each node, build possible (1 phasing if no ambiguity, else some branching phasing)
def build_phasing(q_graph, node_assignment_dict):
    v_label = q_graph.vp['v_label'][0]
    phasing_list = q_graph.vp['v_weights'][0]['weight'].keys()
    local_phasing = phasing_list[state]
    positions = v_label.split('-')
    ploidy = len(positions)
    phasings = {}
    for node, state in node_assignment_dict.items():
        v_label = q_graph.vp['v_label'][node]
        phasing_list = q_graph.vp['v_weights'][node]['weight'].keys()
        local_phasing = phasing_list[state]
        positions = v_label.split('-')
        local_haplotypes = [local_phasing[i:i+len(positions)] for i in range(0, len(positions)*ploidy, positions)]
        for i, position in enumerate(positions):
            if phasings[position]:
                # do some matching
                pass
            else:
                phasing = ""
                for j in range(ploidy):
                    phasing += (local_haplotypes[j][i])
                
    return

# Iterator through all nodes in a spanning tree graph
def node_search(tree_map, prev_node, curr_node, node_assignments, bp):
    prev_node_state = node_assignments[prev_node]
    marginal = bp.pair_beliefs[(prev_node, curr_node)]
    curr_node_state = np.argmax(marginal[prev_node_state])
    node_assignments[curr_node]=curr_node_state
    for neighbor in curr_node.all_neighbors():
        if neighbor == prev_node:
            pass
        if tree_map[to_edge(curr_node, neighbor)]:
            node_search(tree_map=tree_map, prev_node=curr_node, curr_node=neighbor, node_assignments=node_assignments, bp=bp)
    return

# return a single proportionately-sampled phasing from a quotient graph and computed BP
def single_sample_phasings(q_graph, bp):
    # build spanning tree
    tree_map = gt.random_spanning_tree(q_graph)

    # sample edge indices according to spanning tree

    # build phasings
    return

# Return max probability samples phasings over a spanning tree given Quotient Graph and Computed BP graph
def max_phasings(q_graph, bp):
    tree_map = gt.random_spanning_tree(q_graph)
    node_assignments = {}
    starting_node = 0
    node_belief = np.argmax(bp.var_beliefs[starting_node])
    node_assignments[starting_node]=node_belief
    for neighbor in starting_node.all_neighbors():
        if tree_map[to_edge(starting_node, neighbor)]:
            node_search(tree_map=tree_map, prev_node=starting_node, curr_node=neighbor, node_assignments=node_assignments, bp=bp)
    
    build_phasing(node_assignments)

    return


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


if __name__ == "__main__":

    t0 = time.time()
    this_frag_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_3/cov_10/frag"
    this_quotient_coverage_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_3/cov_10/qgraph"
    this_reverse_maps_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_3/cov_10/reverse_maps"
    frag_file = "00.frag"
    ploidy = 3
    genotype_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_3/haplotypes.csv"
    results_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_3/cov_10/results"

    # this_frag_path, this_quotient_coverage_path, this_reverse_maps_path, frag_file, ploidy, genotype_path, results_path = inp
    print('Working on:', os.path.join(this_frag_path, frag_file))
    # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/test.frag'
    # frag_path = '/labs/Aguiar/pHapCompass/test/test2.frag'
    # ploidy= 3
    # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/test/haplotypes.csv'
    # genotype_path = '/labs/Aguiar/pHapCompass/test/haplotypes.csv'

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

    edge_map_path = os.path.join(this_reverse_maps_path, 'qg_e_label_' + frag_file.split('.')[0] + '.pkl')
    with open(edge_map_path, 'rb') as f:
        edges_map_quotient = pickle.load(f)
        
    quotient_g_v_label_reversed_path = os.path.join(this_reverse_maps_path, 'qg_v_label_' + frag_file.split('.')[0] + '.pkl')
    with open(quotient_g_v_label_reversed_path, 'rb') as f:
        quotient_g_v_label_reversed = pickle.load(f)

    transitions_dict, transitions_dict_extra = transition_matrices_v2(quotient_g, edges_map_quotient, ploidy, config, fragment_model)
    emission_dict = emissions_v2(ploidy, quotient_g, quotient_g_v_label_reversed, config.error_rate)

    # frag_graph_path = os.path.join(this_frag_graph_path, frag_file.split('.')[0] + '.gt.gz')
    # frag_graph = gt.load_graph(frag_graph_path)

    # Initialize Markov Network
    mn = create_markov_net(quotient_g, transitions=transitions_dict) # Creates mrftools TorchMarkovNetwork
    t1 = time.time()
    print("Markov Network Created")
    print(f"Time elapsed: {t1-t0}")

    # Belief Propagation Start
    t0 = time.time()
    print("Starting BP")
    torch_bp = mrftools.TorchMatrixBeliefPropagator(
        markov_net=mn, is_cuda=False, var_on=False
    )
    print("BP initialized")
    torch_bp.infer(display="full")
    print("BP inference done")
    torch_bp.load_beliefs()
    t1 = time.time()
    print(f"Time elapsed: {t1-t0}")

    # print("Independent Inference")
    # t0 = time.time()
    # color_dict = gt.sequential_vertex_coloring(quotient_graph)
    # try:
    #     itmbp = IndependentTorchMatrixBeliefPropagator.IndependentTorchMatrixBeliefPropagator(markov_net=mn, is_cuda=False, var_on=False, color_dict=color_dict)
    #     itmbp.independent_inference(display="full")
    #     itmbp.load_beliefs()
    # except Exception as e:
    #     print(e)
    # t1 = time.time()
    # print(f"Time elapsed: {t1-t0}")

    # translate mrftools beliefs into graph tool graph ?

    update_entropies(quotient_g, torch_bp)

    # get path to sample
    print("done")