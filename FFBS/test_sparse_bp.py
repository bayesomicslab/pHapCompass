from mrftools import *
import unittest
import time
import torch as t
import numpy as np
from pgmpy.models import MarkovNetwork
from pgmpy.factors.discrete import DiscreteFactor
from pgmpy.inference.ExactInference import BeliefPropagation

# EXAMPLE 1
#
#   A ----- B ------- C -------- D
#   |       |                    |
#   E ----- F ------- G -------- H
#

# EXAMPLE 2
#
#  A ---- B ---- C ------ D ----- E
#

np.random.seed(1)


def example_1():
    nodes = ["A", "B", "C", "D", "E", "F", "G", "H"]
    edges = [
        ("A", "B"),
        ("B", "C"),
        ("C", "D"),
        ("A", "E"),
        ("E", "F"),
        ("B", "F"),
        ("F", "G"),
        ("G", "H"),
        ("D", "H"),
    ]
    pgmpy_model = MarkovNetwork(edges)
    torch_model = TorchMarkovNet(is_cuda=False, var_on=False)

    for node in nodes:
        factor_values = np.random.rand(2)
        factor = DiscreteFactor([node], cardinality=[2], values=factor_values)
        pgmpy_model.add_factors(factor)
        torch_model.set_unary_factor(node, t.from_numpy(factor_values))

    for edge in edges:
        factor_values = np.random.rand(4)
        factor = DiscreteFactor(edge, cardinality=[2, 2], values=factor_values)
        pgmpy_model.add_factors(factor)
        torch_model.set_edge_factor(
            edge, t.from_numpy(np.reshape(factor_values, (2, 2)))
        )

    print("Torch Edge Potentials ---------------")
    for key, value in torch_model.edge_potentials.items():
        print("Node {}".format(key))
        print(value)
    print("\n")

    print("PGMpy Edge Potentials ---------------")
    for factor in pgmpy_model.get_factors():
        print(factor)
    print("\n")

    bp = BeliefPropagation(pgmpy_model)
    bp.calibrate()

    torch_model.create_matrices()

    torch_bp = TorchMatrixBeliefPropagator(
        markov_net=torch_model, is_cuda=False, var_on=False
    )
    torch_bp.infer(display="full")
    torch_bp.load_beliefs()

    for i in torch_model.variables:
        print(
            "Belief prop unary marginal of %s: %s"
            % (i, repr(t.exp(torch_bp.var_beliefs[i])))
        )

    # print("PGMpy Implementation Dictionary ------------------------------")
    # print(bp.query(variables=nodes))
    print("PGMpy Unary Factors --------------------------------")
    query_dict = bp.query(variables=nodes, joint=False)
    for key, value in query_dict.items():
        print("Factor {}".format(key))
        print(value)
    print("PGMpy Implementation MAP ----------------------------------")
    print(bp.map_query(variables=nodes))

def example_2():
    nodes = ["A", "B", "C"]
    edges = [
        ("A", "B"),
        ("B", "C"),
        ("A", "C")
    ]

    pgmpy_model = MarkovNetwork(edges)
    torch_model = TorchMarkovNet(is_cuda=False, var_on=False)
    for node in nodes:
        if node=="A":
            factor_values = np.array([1, 1000])
        else:
            factor_values = np.array([9, 1])
        factor = DiscreteFactor([node], cardinality=[2], values=factor_values)
        pgmpy_model.add_factors(factor)
        torch_model.set_unary_factor(node, t.from_numpy(factor_values))

    for edge in edges:
        factor_values = np.array([1, 9, 9, 1])
        factor = DiscreteFactor(edge, cardinality=[2, 2], values=factor_values)
        pgmpy_model.add_factors(factor)
        torch_model.set_edge_factor(
            edge, t.from_numpy(np.reshape(factor_values, (2, 2)))
        )

    print("Torch Edge Potentials ---------------")
    for key, value in torch_model.edge_potentials.items():
        print("Node {}".format(key))
        print(value)
    print("\n")

    print("PGMpy Edge Potentials ---------------")
    for factor in pgmpy_model.get_factors():
        print(factor)
    print("\n")

    bp = BeliefPropagation(pgmpy_model)
    bp.calibrate()

    torch_model.create_matrices()

    torch_bp = TorchMatrixBeliefPropagator(
        markov_net=torch_model, is_cuda=False, var_on=False
    )
    torch_bp.infer(display="full")
    torch_bp.load_beliefs()

    for i in torch_model.variables:
        print(
            "Belief prop unary marginal of %s: %s"
            % (i, repr(t.exp(torch_bp.var_beliefs[i])))
        )

    # print("PGMpy Implementation Dictionary ------------------------------")
    # print(bp.query(variables=nodes))
    print("PGMpy Unary Factors --------------------------------")
    query_dict = bp.query(variables=nodes, joint=False)
    for key, value in query_dict.items():
        print("Factor {}".format(key))
        print(value)
    print("PGMpy Implementation MAP ----------------------------------")
    print(bp.map_query(variables=nodes))

    return

import networkx as nx
import matplotlib.pyplot as plt

def graph_coloring():
    G = nx.complete_graph(12)
    color_dict = nx.greedy_color(G)
    nx.draw(G, node_color=[x for x in color_dict.values()])
    plt.draw()
    plt.show()
    return color_dict

def torch_test():
    available = t.cuda.is_available()
    count = t.cuda.device_count()
    curr = t.cuda.current_device()
    name = t.cuda.get_device_name(curr)
    print(f"Cuda is available: {available}")
    print(f"Device count: {count}")
    print(f"Current device: {name}")
    return

if __name__ == "__main__":
    example_2()
    print("done")
