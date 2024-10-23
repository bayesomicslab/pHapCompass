import mrftools
import IndependentTorchMatrixBeliefPropagator
import graph_tool.all as gt
import pgmpy
import torch
import numpy as np
import time
import gzip
import pickle as cPickle
import scipy

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

def update_entropies(q_graph, bp):
    for s, t in q_graph.iter_edges():
        transition_matrix = bp.pair_beliefs[(s, t)]
        entropy_array = scipy.stats.entropy(transition_matrix)
        entropy = np.mean(entropy_array)
        q_graph.ep['e_weights'][(s, t)]['entropy'] = entropy
    return

def to_edge(source, target):
    return

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

def single_sample_phasings(q_graph, bp):
    # build spanning tree
    tree_map = gt.random_spanning_tree(q_graph)

    # sample edge indices according to spanning tree

    # build phasings
    return

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

def torch_belief_propagation(markov_network, is_cuda=True):
    return

def pgmpy_belief_propagation(markov_network):
    return

def independent_belief_propagation(markov_network, is_cuda=True):
    return

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
    quotient_graph = gt.Graph()
    print("Q Graph initialized")
    # quotient_graph.load("/home/tgbergendahl/research/hap_assembly/HaplOrbit/example/ART_0_quotient_graph.gt.gz")
    quotient_graph.load("/Users/tgbergendahl/Research/HaplOrbit/example/ART_0_quotient_graph.gt.gz")
    print("Q Graph loaded")
    transitions_path = '/Users/tgbergendahl/Research/HaplOrbit/example/ART_0.pkl.gz'
    with gzip.open(transitions_path, 'rb') as f:
        transitions = cPickle.load(f)
    print("Transitions Loaded")
    mn = create_markov_net(quotient_graph, transitions=transitions)
    t1 = time.time()
    print("Markov Network Created")
    print(f"Time elapsed: {t1-t0}")
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
    print("done")