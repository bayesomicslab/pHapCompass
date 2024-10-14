import mrftools
import IndependentTorchMatrixBeliefPropagator
import graph_tool.all as gt
import pgmpy
import torch
import numpy as np
import time

def create_markov_net(q_graph):
    # initialize markov net object
    markov_net = mrftools.TorchMarkovNet(is_cuda=True, var_on=False)

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
        weight_dict = q_graph.ep['e_weights'][edge]['weight']
        weights = list(weight_dict.values())
        weight_array = np.array(weights)
        markov_net.set_edge_factor((s, t), torch.from_numpy(np.reshape(weight_array, (s_size, t_size))))

    return markov_net

if __name__ == "__main__":
    t0 = time.time()
    quotient_graph = gt.Graph()
    print("Q Graph initialized")
    quotient_graph.load("/home/tgbergendahl/research/hap_assembly/HaplOrbit/example/ART_0_quotient_graph.gt.gz")
    print("Q Graph loaded")
    mn = create_markov_net(quotient_graph)
    t1 = time.time()
    print(f"Time elapsed: {t1-t0}")
    print("done")