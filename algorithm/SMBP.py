import mrftools
import IndependentTorchMatrixBeliefPropagator
import graph_tool.all as gt
import pgmpy
import torch
import numpy as np
import time
import gzip
import pickle as cPickle

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
        s_size = len(q_graph.vp['e_weights'][s]['weight'].values())
        # s_label = q_graph.vp['v_label'][s]
        t_size = len(q_graph.vp['e_weights'][t]['weight'].values())
        # t_label = q_graph.vp['v_label'][t]
        e_label = q_graph.ep['e_label'][edge]
        weight_matrix = transitions[e_label]['transitions']
        markov_net.set_edge_factor((s, t), torch.from_numpy(weight_matrix))

    print("Edge factors encoded")

    return markov_net

def sample_phasings(q_graph, markov_net):
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