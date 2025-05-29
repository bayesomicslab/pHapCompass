import mrftoolstest
import numpy as np
import pandas as pd
import torch

if __name__ == "__main__":
    a_weight = 100000
    d_weight = 1
    epsilon = 0.0001

    nodes = ["A", "B", "C", "D"]
    node_dict = {
        "A": torch.tensor([a_weight, epsilon]),
        "B": torch.tensor([epsilon, epsilon]),
        "C": torch.tensor([epsilon, epsilon]),
        "D": torch.tensor([epsilon, d_weight])
    }
    edges = [("A", "B"), ("A", "C"), ("B", "C"), ("B", "D"), ("C", "D")]

    edge_transition = [[0.99, 0.01], [0.01, 0.99]]
    edge_dict = {
        ("A", "B"): torch.tensor(edge_transition),
        ("C", "D"): torch.tensor(edge_transition),
        ("B", "C"): torch.tensor(edge_transition),
        ("A", "C"): torch.tensor(edge_transition),
        ("B", "D"): torch.tensor(edge_transition)
    }

    weight_map = ()

    map_indices = torch.tensor([
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        [1, 0, 2, 0, 2, 1, 3, 1, 3, 2]
    ])

    map_weights = torch.DoubleTensor(
        [a_weight, 1, a_weight, 1, 1, 1, 1, d_weight, 1, d_weight]
    )

    map_size = torch.Size([2*len(edges), 4])

    weight_map = torch.sparse.DoubleTensor(
        map_indices,
        map_weights,
        map_size
    )

    markov_net = mrftoolstest.TorchMarkovNet(is_cuda=False, var_on=False, weight_map=weight_map)
    for node in nodes:
        markov_net.set_unary_factor(node, node_dict[node])
    for edge in edges:
        markov_net.set_edge_factor(edge, edge_dict[edge])
    
    torch_bp = mrftoolstest.TorchMatrixBeliefPropagator(
        markov_net, is_cuda=False, var_on=False
    )

    torch_bp.infer(display="full")

    torch_bp.load_beliefs()

    for key, value in torch_bp.var_beliefs.items():
        print(key)
        print(np.exp(value))

    for key, value in torch_bp.pair_beliefs.items():
        print(key)
        print(np.exp(value))

    print("Done")