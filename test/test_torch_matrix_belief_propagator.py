"SOURCED FROM https://github.com/VTCSML/mrftools/tree/pytorch-test"

"""Tests for Matrix belief propagation"""

from mrftools import *
import unittest
import time
import torch as t
import numpy as np

try:
    import opengm
except:
    print("Could not load opengm package. At least one test will fail.")


class TestTorchMatrixBeliefPropagator(unittest.TestCase):
    """Test class for TorchMatrixBeliefPropagator"""

    def create_chain_model(self, is_cuda):
        """Create chain MRF with variable of different cardinalities."""
        mn = TorchMarkovNet(is_cuda=is_cuda, var_on=False)

        np.random.seed(1)

        k = [4, 3, 6, 2, 5]

        if is_cuda:
            mn.set_unary_factor(0, t.from_numpy(np.random.randn(k[0])).cuda())
            mn.set_unary_factor(1, t.from_numpy(np.random.randn(k[1])).cuda())
            mn.set_unary_factor(2, t.from_numpy(np.random.randn(k[2])).cuda())
            mn.set_unary_factor(3, t.from_numpy(np.random.randn(k[3])).cuda())
        else:
            mn.set_unary_factor(0, t.from_numpy(np.random.randn(k[0])))
            mn.set_unary_factor(1, t.from_numpy(np.random.randn(k[1])))
            mn.set_unary_factor(2, t.from_numpy(np.random.randn(k[2])))
            mn.set_unary_factor(3, t.from_numpy(np.random.randn(k[3])))

        if is_cuda:
            factor4 = t.from_numpy(np.random.randn(k[4])).cuda()
        else:
            factor4 = t.from_numpy(np.random.randn(k[4]))
        factor4[2] = -float("inf")

        mn.set_unary_factor(4, factor4)

        if is_cuda:
            mn.set_edge_factor((0, 1), t.from_numpy(np.random.randn(k[0], k[1])).cuda())
            mn.set_edge_factor((1, 2), t.from_numpy(np.random.randn(k[1], k[2])).cuda())
            mn.set_edge_factor((2, 3), t.from_numpy(np.random.randn(k[2], k[3])).cuda())
            mn.set_edge_factor((3, 4), t.from_numpy(np.random.randn(k[3], k[4])).cuda())
        else:
            mn.set_edge_factor((0, 1), t.from_numpy(np.random.randn(k[0], k[1])))
            mn.set_edge_factor((1, 2), t.from_numpy(np.random.randn(k[1], k[2])))
            mn.set_edge_factor((2, 3), t.from_numpy(np.random.randn(k[2], k[3])))
            mn.set_edge_factor((3, 4), t.from_numpy(np.random.randn(k[3], k[4])))

        mn.create_matrices()

        return mn

    def create_chain_model_old(self):
        """Create chain MRF with variable of different cardinalities."""
        mn = MarkovNet()

        np.random.seed(1)

        k = [4, 3, 6, 2, 5]

        mn.set_unary_factor(0, np.random.randn(k[0]))
        mn.set_unary_factor(1, np.random.randn(k[1]))
        mn.set_unary_factor(2, np.random.randn(k[2]))
        mn.set_unary_factor(3, np.random.randn(k[3]))

        factor4 = np.random.randn(k[4])
        factor4[2] = -float("inf")

        mn.set_unary_factor(4, factor4)

        mn.set_edge_factor((0, 1), np.random.randn(k[0], k[1]))
        mn.set_edge_factor((1, 2), np.random.randn(k[1], k[2]))
        mn.set_edge_factor((2, 3), np.random.randn(k[2], k[3]))
        mn.set_edge_factor((3, 4), np.random.randn(k[3], k[4]))
        mn.create_matrices()

        return mn

    def create_loop_model(self, is_cuda):
        """Create a loop-structured MRF"""
        mn = self.create_chain_model(is_cuda=is_cuda)

        k = [4, 3, 6, 2, 5]

        if is_cuda:
            mn.set_edge_factor((3, 0), t.from_numpy(np.random.randn(k[3], k[0])).cuda())
        else:
            mn.set_edge_factor((3, 0), t.from_numpy(np.random.randn(k[3], k[0])))
        mn.create_matrices()
        return mn

    def create_grid_model(self, is_cuda, my_l, my_k):
        """Create a grid-structured MRF"""
        mn = TorchMarkovNet(is_cuda=is_cuda, var_on=False)
        np.random.seed(1)

        length = my_l

        k = my_k

        for x in range(length):
            for y in range(length):
                if is_cuda:
                    mn.set_unary_factor(
                        (x, y), t.from_numpy(np.random.random(k)).cuda()
                    )
                else:
                    mn.set_unary_factor((x, y), t.from_numpy(np.random.random(k)))

        for x in range(length - 1):
            for y in range(length):
                if is_cuda:
                    mn.set_edge_factor(
                        ((x, y), (x + 1, y)),
                        t.from_numpy(np.random.random((k, k))).cuda(),
                    )
                    mn.set_edge_factor(
                        ((y, x), (y, x + 1)),
                        t.from_numpy(np.random.random((k, k))).cuda(),
                    )
                else:
                    mn.set_edge_factor(
                        ((x, y), (x + 1, y)), t.from_numpy(np.random.random((k, k)))
                    )
                    mn.set_edge_factor(
                        ((y, x), (y, x + 1)), t.from_numpy(np.random.random((k, k)))
                    )

        return mn

    def create_grid_model_old(self, my_l, my_k):
        """Create a grid-structured MRF"""
        mn = MarkovNet()
        np.random.seed(1)

        length = my_l

        k = my_k

        for x in range(length):
            for y in range(length):
                mn.set_unary_factor((x, y), np.random.random(k))

        for x in range(length - 1):
            for y in range(length):
                mn.set_edge_factor(((x, y), (x + 1, y)), np.random.random((k, k)))
                mn.set_edge_factor(((y, x), (y, x + 1)), np.random.random((k, k)))

        return mn

    def create_grid_model_simple_edges(self, is_cuda):
        """Create a grid-structured MRFs with edge potentials that are attractive."""
        mn = TorchMarkovNet(is_cuda=is_cuda, var_on=False)
        np.random.seed(1)

        length = 2

        k = 8

        for x in range(length):
            for y in range(length):
                if is_cuda:
                    mn.set_unary_factor(
                        (x, y), t.from_numpy(np.random.random(k)).cuda()
                    )
                else:
                    mn.set_unary_factor((x, y), t.from_numpy(np.random.random(k)))

        for x in range(length - 1):
            for y in range(length):
                if is_cuda:
                    mn.set_edge_factor(((x, y), (x + 1, y)), t.eye(k).cuda())
                    mn.set_edge_factor(((y, x), (y, x + 1)), t.eye(k).cuda())
                else:
                    mn.set_edge_factor(((x, y), (x + 1, y)), t.eye(k))
                    mn.set_edge_factor(((y, x), (y, x + 1)), t.eye(k))

        return mn

    def create_grid_model_simple_edges_old(self):
        """Create a grid-structured MRFs with edge potentials that are attractive."""
        mn = MarkovNet()
        np.random.seed(1)

        length = 2

        k = 8

        for x in range(length):
            for y in range(length):
                mn.set_unary_factor((x, y), np.random.random(k))

        for x in range(length - 1):
            for y in range(length):
                mn.set_edge_factor(((x, y), (x + 1, y)), np.eye(k))
                mn.set_edge_factor(((y, x), (y, x + 1)), np.eye(k))

        return mn

    def test_exactness(self):
        """Test that Matrix BP produces the true marginals in a chain model."""
        is_cuda = False
        mn = self.create_chain_model(is_cuda=is_cuda)
        bp = TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=is_cuda, var_on=False)
        bp.infer(display="full")
        bp.load_beliefs()

        mn_old = self.create_chain_model_old()
        bf = BruteForce(mn_old)

        for i in mn.variables:
            print(
                "Brute force unary marginal of %d: %s" % (i, repr(bf.unary_marginal(i)))
            )
            print(
                "Belief prop unary marginal of %d: %s"
                % (i, repr(t.exp(bp.var_beliefs[i])))
            )
            if is_cuda:
                assert np.allclose(
                    bf.unary_marginal(i), t.exp(bp.var_beliefs[i]).cpu().numpy()
                ), "beliefs aren't exact on chain model"
            else:
                assert np.allclose(
                    bf.unary_marginal(i), t.exp(bp.var_beliefs[i]).numpy()
                ), "beliefs aren't exact on chain model"

        print("Brute force pairwise marginal: " + repr(bf.pairwise_marginal(0, 1)))
        print("Belief prop pairwise marginal: " + repr(t.exp(bp.pair_beliefs[(0, 1)])))

        print("Bethe energy functional: %f" % bp.compute_energy_functional())

        print("Brute force log partition function: %f" % np.log(bf.compute_z()))

        assert np.allclose(
            np.log(bf.compute_z()), bp.compute_energy_functional()
        ), "log partition function is not exact on chain model"

    def test_consistency(self):
        """Test that loopy matrix BP infers marginals that are locally consistent."""
        is_cuda = False
        mn = self.create_loop_model(is_cuda=is_cuda)

        bp = TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=is_cuda, var_on=False)
        bp.infer(display="full")

        bp.load_beliefs()

        for var in mn.variables:
            unary_belief = t.exp(bp.var_beliefs[var])
            for neighbor in mn.get_neighbors(var):
                pair_belief = t.sum(t.exp(bp.pair_beliefs[(var, neighbor)]), 1)
                print(pair_belief, unary_belief)
                if is_cuda:
                    assert np.allclose(
                        pair_belief.cpu().numpy(), unary_belief.cpu().numpy()
                    ), "unary and pairwise beliefs are inconsistent"
                else:
                    assert np.allclose(
                        pair_belief.numpy(), unary_belief.numpy()
                    ), "unary and pairwise beliefs are inconsistent"

    def test_normalization(self):
        """Test that the unary and pairwise beliefs properly sum to 1.0"""
        is_cuda = False
        mn = self.create_loop_model(is_cuda=is_cuda)

        bp = TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=is_cuda, var_on=False)
        bp.infer(display="full")

        bp.load_beliefs()

        for var in mn.variables:
            unary_belief = t.exp(bp.var_beliefs[var])
            assert np.allclose(
                t.sum(unary_belief), 1.0
            ), "unary belief is not normalized"
            for neighbor in mn.get_neighbors(var):
                pair_belief = t.exp(bp.pair_beliefs[(var, neighbor)])
                assert np.allclose(
                    t.sum(pair_belief), 1.0
                ), "pairwise belief is not normalize"

    def test_conditioning(self):
        """Test that conditioning on variable properly sets variables to conditioned state"""
        is_cuda = False
        mn = self.create_loop_model(is_cuda=is_cuda)

        bp = TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=is_cuda, var_on=False)

        bp.condition(2, 0)

        bp.infer()
        bp.load_beliefs()

        assert np.allclose(
            bp.var_beliefs[2][0], 0
        ), "Conditioned variable was not set to correct state"

        beliefs0 = bp.var_beliefs[0]

        bp.condition(2, 1)
        bp.infer()
        bp.load_beliefs()
        beliefs1 = bp.var_beliefs[0]

        assert not np.allclose(
            beliefs0.numpy(), beliefs1.numpy()
        ), "Conditioning var 2 did not change beliefs of var 0"

    def test_overflow(self):
        """Test that MatrixBP does not fail when given very large, poorly scaled factors"""
        is_cuda = False
        mn = self.create_chain_model(is_cuda=is_cuda)

        # set a really large factor
        mn.set_unary_factor(0, t.FloatTensor([1000, 2000, 3000, 4000]))

        mn.create_matrices()

        bp = TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=is_cuda, var_on=False)

        with np.errstate(all="raise"):
            bp.infer()
            bp.load_beliefs()

    def test_grid_consistency(self):
        """Test that matrix BP infers consistent marginals on a grid MRF"""
        is_cuda = False
        mn = self.create_grid_model(is_cuda=is_cuda, my_l=64, my_k=8)
        bp = TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=is_cuda, var_on=False)
        bp.infer(display="full")

        bp.load_beliefs()

        for var in mn.variables:
            unary_belief = np.exp(bp.var_beliefs[var])
            for neighbor in mn.get_neighbors(var):
                pair_belief = t.sum(t.exp(bp.pair_beliefs[(var, neighbor)]), 1)
                # print pair_belief, unary_belief
                assert np.allclose(
                    pair_belief.numpy(), unary_belief.numpy()
                ), "unary and pairwise beliefs are inconsistent"

    def test_speedup(self):
        """Test that matrix BP is faster than loop-based BP"""
        is_cuda = False
        mn = self.create_grid_model(is_cuda=is_cuda, my_l=64, my_k=8)
        old_mn = self.create_grid_model_old(my_l=64, my_k=8)

        bp = TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=is_cuda, var_on=False)
        old_bp = MatrixBeliefPropagator(old_mn)

        bp.set_max_iter(1000)
        old_bp.set_max_iter(1000)

        t0 = time.time()
        bp.infer(display="final")
        t1 = time.time()

        bp_time = t1 - t0

        t0 = time.time()
        old_bp.infer(display="final")
        t1 = time.time()

        old_bp_time = t1 - t0

        print(
            "Torch Matrix BP took %f, Matrix BP took %f. Speedup was %f"
            % (bp_time, old_bp_time, old_bp_time / bp_time)
        )
        assert bp_time > old_bp_time, "Torch Matrix form was slower than Matrix BP"

        # check marginals
        bp.load_beliefs()
        old_bp.load_beliefs()

        for var in mn.variables:
            assert np.allclose(
                bp.var_beliefs[var].numpy(), old_bp.var_beliefs[var]
            ), "unary beliefs don't agree"
            for neighbor in mn.get_neighbors(var):
                edge = (var, neighbor)
                assert np.allclose(
                    bp.pair_beliefs[edge].numpy(), old_bp.pair_beliefs[edge]
                ), (
                    "pairwise beliefs don't agree"
                    + "\n"
                    + repr(bp.pair_beliefs[edge])
                    + "\n"
                    + repr(old_bp.pair_beliefs[edge])
                )

    def test_belief_propagator_messages(self):
        """Test that matrix BP and loop-based BP calculate the same messages and beliefs each iteration of inference"""
        is_cuda = False
        model_old = self.create_grid_model_simple_edges_old()
        model = self.create_grid_model_simple_edges(is_cuda=is_cuda)
        bp = BeliefPropagator(model_old)
        bp.load_beliefs()

        mat_bp = TorchMatrixBeliefPropagator(
            markov_net=model, is_cuda=is_cuda, var_on=False
        )
        mat_bp.load_beliefs()

        for i in range(4):
            for var in sorted(bp.mn.variables):
                for neighbor in sorted(bp.mn.get_neighbors(var)):
                    edge = (var, neighbor)
                    bp_message = bp.messages[edge]

                    if edge in mat_bp.mn.message_index:
                        edge_index = mat_bp.mn.message_index[edge]
                    else:
                        edge_index = (
                            mat_bp.mn.message_index[(edge[1], edge[0])]
                            + mat_bp.mn.num_edges
                        )

                    # ugly transition to and from numpy because of ravel
                    mat_bp_message = t.from_numpy(
                        mat_bp.message_mat[:, edge_index].numpy().ravel()
                    )

                    assert np.allclose(bp_message, mat_bp_message.numpy()), (
                        "BP and matBP did not agree on message for edge %s in iter %d"
                        % (repr(edge), i)
                        + "\nBP: "
                        + repr(bp_message)
                        + "\nmatBP: "
                        + repr(mat_bp_message)
                    )

                    # print "Message %s is OK" % repr(edge)

                    assert np.allclose(
                        bp.pair_beliefs[edge], mat_bp.pair_beliefs[edge].numpy()
                    ), (
                        "BP and matBP did not agree on pair beliefs after %d message updates"
                        % i
                    )

                assert np.allclose(
                    bp.var_beliefs[var], mat_bp.var_beliefs[var].numpy()
                ), (
                    "BP and matBP did not agree on unary beliefs after %d message updates"
                    % i
                )

            bp.update_messages()
            bp.load_beliefs()
            mat_bp.update_messages()
            mat_bp.load_beliefs()

    def test_torch_logsumexp(self):
        ins_try_mat2 = t.DoubleTensor(
            [
                [1.62434536, -2.50015019, -3.89192232, -0.98001738, -1.7067365],
                [-1.37871011, -3.49762318, -3.663615, 1.13376944, -0.17242821],
                [-1.37513475, 1.47182113, -2.63786005, -float("inf"), -float("inf")],
                [-1.72311984, -float("inf"), 1.29355548, -float("inf"), -2.9843092],
                [-float("inf"), -float("inf"), -4.25291378, -float("inf"), -1.37131876],
                [
                    -float("inf"),
                    -float("inf"),
                    -1.56686887,
                    -float("inf"),
                    -float("inf"),
                ],
            ]
        )

        res_try_mat2 = t.DoubleTensor(
            [1.75064451, 1.49727762, 1.38283992, 1.24779407, 0.28323888]
        )

        ins_try_mat3 = t.DoubleTensor(
            [
                [
                    [
                        -1.35561165,
                        -1.92600932,
                        -2.51375888,
                        -1.20307934,
                        -2.3394558,
                        -3.5011019,
                        -1.41583884,
                        -1.97543269,
                    ],
                    [
                        -1.98859991,
                        -2.74598077,
                        -2.8762014,
                        1.32163145,
                        -3.26105924,
                        -2.57888885,
                        -0.55043521,
                        -1.85899228,
                    ],
                    [
                        -2.53039981,
                        -0.30862473,
                        -2.51098691,
                        -float("inf"),
                        0.54215823,
                        -2.97333881,
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -2.42195099,
                        -float("inf"),
                        0.24836905,
                        -float("inf"),
                        -float("inf"),
                        -0.60327896,
                        -float("inf"),
                        -1.42182921,
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -3.82727569,
                        -float("inf"),
                        -float("inf"),
                        -5.23017336,
                        -float("inf"),
                        0.32597336,
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -3.12953745,
                        -float("inf"),
                        -float("inf"),
                        -2.14072381,
                        -float("inf"),
                        -float("inf"),
                    ],
                ],
                [
                    [
                        0.88973124,
                        -2.08404226,
                        -3.16617905,
                        -2.014102,
                        -0.73634228,
                        -1.15412702,
                        -2.85852736,
                        -0.96854569,
                    ],
                    [
                        -1.5902383,
                        -3.66373879,
                        -2.08302593,
                        1.98739004,
                        -3.504927,
                        -0.99163904,
                        -0.54752809,
                        0.62467596,
                    ],
                    [
                        -3.34327902,
                        -0.99642814,
                        -1.71546949,
                        -float("inf"),
                        -1.04316035,
                        -2.49392823,
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -3.64396721,
                        -float("inf"),
                        0.34048772,
                        -float("inf"),
                        -float("inf"),
                        -1.47824333,
                        -float("inf"),
                        -2.7903713,
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -3.22772185,
                        -float("inf"),
                        -float("inf"),
                        -4.8600213,
                        -float("inf"),
                        -1.31012189,
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -1.12227135,
                        -float("inf"),
                        -float("inf"),
                        -0.6826849,
                        -float("inf"),
                        -float("inf"),
                    ],
                ],
                [
                    [
                        0.64659825,
                        -1.91008275,
                        -float("inf"),
                        -1.22749101,
                        -1.36172685,
                        -2.76312145,
                        -1.92490339,
                        -float("inf"),
                    ],
                    [
                        -3.17482211,
                        -4.5976185,
                        -float("inf"),
                        -1.5946817,
                        -5.34155238,
                        -2.37067885,
                        0.38843783,
                        -float("inf"),
                    ],
                    [
                        -2.67539767,
                        -0.168517,
                        -float("inf"),
                        -float("inf"),
                        -0.62732057,
                        -2.11117719,
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -3.34905998,
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        1.50964077,
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -3.99270406,
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -1.7579364,
                        -float("inf"),
                        -float("inf"),
                    ],
                ],
                [
                    [
                        -float("inf"),
                        -1.25150122,
                        -float("inf"),
                        -1.79158087,
                        -0.70848116,
                        -float("inf"),
                        -0.87702575,
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -5.29341191,
                        -float("inf"),
                        -1.64229917,
                        -5.0974437,
                        -float("inf"),
                        0.73291672,
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        1.74082265,
                        -float("inf"),
                        -float("inf"),
                        -0.75618602,
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                ],
                [
                    [
                        -float("inf"),
                        -2.35614697,
                        -float("inf"),
                        -0.58437977,
                        -float("inf"),
                        -float("inf"),
                        -1.43042184,
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -5.15294124,
                        -float("inf"),
                        -0.70265123,
                        -float("inf"),
                        -float("inf"),
                        0.6869558,
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -0.23927354,
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                ],
                [
                    [
                        -float("inf"),
                        -1.00442093,
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -2.4704071,
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -2.71332835,
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        1.05468279,
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        0.25777062,
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                    [
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                    ],
                ],
            ]
        )

        res_try_mat3 = t.DoubleTensor(
            [
                [
                    [
                        -0.5743833,
                        -0.05723698,
                        0.44848604,
                        1.39866793,
                        0.61758397,
                        -0.1898393,
                        -0.19915766,
                        0.65294042,
                    ]
                ],
                [
                    [
                        0.99318571,
                        -0.65535341,
                        0.75019678,
                        2.00551315,
                        -0.14936972,
                        0.41206502,
                        -0.45297991,
                        0.94719271,
                    ]
                ],
                [
                    [
                        0.72019508,
                        0.00305616,
                        -float("inf"),
                        -0.70117939,
                        -0.22937903,
                        1.60814951,
                        0.48277494,
                        -float("inf"),
                    ]
                ],
                [
                    [
                        -float("inf"),
                        1.79061412,
                        -float("inf"),
                        -1.0210098,
                        -0.03256725,
                        -float("inf"),
                        0.9151542,
                        -float("inf"),
                    ]
                ],
                [
                    [
                        -float("inf"),
                        -0.11904638,
                        -float("inf"),
                        0.05137918,
                        -float("inf"),
                        -float("inf"),
                        0.80059409,
                        -float("inf"),
                    ]
                ],
                [
                    [
                        -float("inf"),
                        0.54616246,
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        -float("inf"),
                        1.08370665,
                        -float("inf"),
                    ]
                ],
            ]
        )

        ins_except = t.DoubleTensor(
            [
                [
                    1.00000000e03,
                    8.65407629e-01,
                    -7.61206901e-01,
                    -3.84054355e-01,
                    -1.09989127e00,
                ],
                [
                    2.00000000e03,
                    -2.30153870e00,
                    3.19039096e-01,
                    1.13376944e00,
                    -1.72428208e-01,
                ],
                [
                    3.00000000e03,
                    1.74481176e00,
                    -2.49370375e-01,
                    -float("inf"),
                    -float("inf"),
                ],
                [
                    4.00000000e03,
                    -float("inf"),
                    1.46210794e00,
                    -float("inf"),
                    4.22137467e-02,
                ],
                [
                    -float("inf"),
                    -float("inf"),
                    -2.06014071e00,
                    -float("inf"),
                    5.82815214e-01,
                ],
                [
                    -float("inf"),
                    -float("inf"),
                    -3.22417204e-01,
                    -float("inf"),
                    -float("inf"),
                ],
            ]
        )

        res_except = t.DoubleTensor(
            [4.00000000e03, 2.10424425e00, 2.05272230e00, 1.33195481e00, 1.38847124e00]
        )

        ins_try_dim2 = t.DoubleTensor(
            [
                [
                    [-3.67236818, -6.33946165, -5.65306655, -2.65163243],
                    [-4.21779995, -5.61582591, -4.49341063, -2.53519202],
                    [1.0693653, -5.07447092, -float("inf"), -float("inf")],
                    [-float("inf"), -2.29221174, -float("inf"), -2.09802895],
                    [-float("inf"), -7.01803503, -float("inf"), -0.35022638],
                    [-float("inf"), -3.62502753, -float("inf"), -float("inf")],
                ],
                [
                    [-4.30535645, -5.59186408, -6.69107919, 0.02220336],
                    [-6.69776948, -5.62795342, -4.08582763, 1.61542501],
                    [-2.75205506, -6.19443765, -float("inf"), -float("inf")],
                    [-float("inf"), -4.76655342, -float("inf"), -1.79962225],
                    [-float("inf"), -8.24726029, -float("inf"), -0.31937284],
                    [-float("inf"), -3.76636593, -float("inf"), -float("inf")],
                ],
                [
                    [-4.84715635, -3.42749867, -5.60475731, -float("inf")],
                    [-8.45081021, -3.2336334, -2.9971638, -float("inf")],
                    [-2.25263062, -2.03832678, -float("inf"), -float("inf")],
                    [-float("inf"), 1.99469051, -float("inf"), -float("inf")],
                    [-float("inf"), -3.60658321, -float("inf"), -float("inf")],
                    [-float("inf"), -1.0682576, -float("inf"), -float("inf")],
                ],
                [
                    [-4.73870753, -float("inf"), -0.98300863, -float("inf")],
                    [-8.75149839, -float("inf"), 0.92118612, -float("inf")],
                    [-2.92629294, -float("inf"), -float("inf"), -float("inf")],
                    [-float("inf"), -float("inf"), -float("inf"), -float("inf")],
                    [-float("inf"), -float("inf"), -float("inf"), -float("inf")],
                    [-float("inf"), -float("inf"), -float("inf"), -float("inf")],
                ],
                [
                    [-float("inf"), -float("inf"), -7.01178459, -float("inf")],
                    [-float("inf"), -float("inf"), -4.60015467, -float("inf")],
                    [-float("inf"), -float("inf"), -float("inf"), -float("inf")],
                    [-float("inf"), -float("inf"), -float("inf"), -float("inf")],
                    [-float("inf"), -float("inf"), -float("inf"), -float("inf")],
                    [-float("inf"), -float("inf"), -float("inf"), -float("inf")],
                ],
                [
                    [-float("inf"), -float("inf"), -5.72345593, -float("inf")],
                    [-float("inf"), -float("inf"), -1.90411375, -float("inf")],
                    [-float("inf"), -float("inf"), -float("inf"), -float("inf")],
                    [-float("inf"), -float("inf"), -float("inf"), -float("inf")],
                    [-float("inf"), -float("inf"), -float("inf"), -float("inf")],
                    [-float("inf"), -float("inf"), -float("inf"), -float("inf")],
                ],
            ]
        )

        res_try_dim2 = t.DoubleTensor([1.16561449, 2.09239801, 1.14294035, 2.06957452])

        np_exp = np.exp(ins_try_mat3.numpy())
        torch_exp = t.exp(ins_try_mat3)

        np_sum = np.sum(np_exp, 1, keepdims=True)
        torch_sum = t.sum(torch_exp, dim=1, keepdim=True)

        np_log = np.log(np_sum)
        torch_log = t.log(torch_sum)

        np_max = np.nan_to_num(np.max(ins_try_mat3.numpy(), axis=1, keepdims=True))
        torch_max = torch_nan_to_num(t.max(ins_try_mat3, dim=1, keepdim=True)[0])

        np_max_tuple = np.nan_to_num(
            np.max(ins_try_mat3.numpy(), axis=(0, 1), keepdims=True)
        )
        torch_max_tuple = ins_try_mat3.max(dim=0, keepdim=True)[0].max(
            dim=1, keepdim=True
        )[0]

        np_sum_tuple = (
            np.sum(np.exp(ins_try_mat3.numpy() - np_max_tuple), (0, 1), keepdims=True)
            + np_max_tuple
        )
        torch_sum_tuple = (
            t.exp(ins_try_mat3 - torch_max_tuple)
            .sum(dim=0, keepdim=True)
            .sum(dim=1, keepdim=True)
            + torch_max_tuple
        )

        assert np.allclose(
            np_exp, torch_exp.numpy()
        ), "exponential on tensors is different between numpy and torch"
        assert np.allclose(
            np_sum, torch_sum.numpy()
        ), "sum on tensors is different between numpy and torch"
        assert np.allclose(
            np_log, torch_log.numpy()
        ), "logarithms on tensors is different between numpy and torch"
        assert np.allclose(
            np_max, torch_max.numpy()
        ), "maxes on tensors is different between numpy and torch"
        assert np.allclose(
            res_try_mat2.numpy(), torch_logsumexp(ins_try_mat2, 0).numpy()
        ), "torch_logsumexp failed for a 2 dimensional tensor"
        assert np.allclose(
            res_try_mat3.numpy(), torch_logsumexp(ins_try_mat3, 1).numpy()
        ), "torch_logsumexp failed for a 3 dimensional tensor"
        assert np.allclose(
            res_except.numpy(), torch_logsumexp(ins_except, 0).numpy()
        ), "torch_logsumexp failed when max_val was needed"
        assert np.allclose(
            np_max_tuple, torch_max_tuple.numpy()
        ), "maxes on tensors is different between numpy and torch when using tuples"
        assert np.allclose(
            np_sum_tuple, torch_sum_tuple.numpy()
        ), "sums on tensors is different between numpy and torch when using tuples"
        assert np.allclose(
            res_try_dim2.numpy(), torch_logsumexp(ins_try_dim2, (0, 1)).numpy()
        ), "torch_logsumexp failed when it was given a tuple"

    def test_gpu_speedup(self):
        if not t.cuda.is_available():
            print("CUDA is not available.")
            return

        self.my_l = 8
        self.my_k = 8
        slow = True

        # CUDA must be "started up" in order to achieve max speeds
        # Approximate loss of X seconds on first run otherwise
        init_mn = self.create_grid_model(is_cuda=True, my_l=2, my_k=2)
        init_bp = TorchMatrixBeliefPropagator(
            markov_net=init_mn, is_cuda=True, var_on=False
        )
        init_bp.set_max_iter(1000)
        init_bp.infer(display="off")

        t_prime0 = time.time()

        cuda_mn = self.create_grid_model(is_cuda=True, my_l=self.my_l, my_k=self.my_k)
        cuda_bp = TorchMatrixBeliefPropagator(
            markov_net=cuda_mn, is_cuda=True, var_on=False
        )
        cuda_bp.set_max_iter(1000)
        t0 = time.time()
        cuda_bp.infer(display="off")
        t1 = time.time()
        cuda_bp_time = t1 - t0

        mn = self.create_grid_model(is_cuda=False, my_l=self.my_l, my_k=self.my_k)
        bp = TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=False, var_on=False)
        bp.set_max_iter(1000)
        t0 = time.time()
        bp.infer(display="off")
        t1 = time.time()
        bp_time = t1 - t0

        old_mn = self.create_grid_model_old(my_l=self.my_l, my_k=self.my_k)
        old_bp = MatrixBeliefPropagator(old_mn)
        old_bp.set_max_iter(1000)
        t0 = time.time()
        old_bp.infer(display="off")
        t1 = time.time()
        old_bp_time = t1 - t0

        if slow:
            slow_bp = BeliefPropagator(old_mn)
            slow_bp.set_max_iter(1000)
            t0 = time.time()
            slow_bp.infer(display="off")
            t1 = time.time()
            slow_bp_time = t1 - t0

        t_prime1 = time.time()

        start_time = (t_prime1 - t_prime0) - cuda_bp_time - bp_time - old_bp_time
        if slow:
            start_time = start_time - slow_bp_time

        print("Build time took %f" % (start_time))
        print("CUDA Torch Matrix BP took %f" % (cuda_bp_time))
        print("Torch Matrix BP took %f" % (bp_time))
        print("Sparse Matrix BP took %f" % (old_bp_time))
        if slow:
            print("Loop Matrix BP took %f" % (slow_bp_time))
        # assert bp_time > old_bp_time, "Torch Matrix BP was faster than Matrix BP"
        # assert cuda_bp_time < bp_time, "CUDA Torch Matrix BP was slower than Torch Matrix BP"

        # check marginals
        cuda_bp.load_beliefs()
        old_bp.load_beliefs()

        for var in mn.variables:
            assert np.allclose(
                cuda_bp.var_beliefs[var].cpu().numpy(), old_bp.var_beliefs[var]
            ), "unary beliefs don't agree"
            for neighbor in mn.get_neighbors(var):
                edge = (var, neighbor)
                assert np.allclose(
                    cuda_bp.pair_beliefs[edge].cpu().numpy(), old_bp.pair_beliefs[edge]
                ), (
                    "pairwise beliefs don't agree"
                    + "\n"
                    + repr(cuda_bp.pair_beliefs[edge])
                    + "\n"
                    + repr(old_bp.pair_beliefs[edge])
                )

    def validate_params(self, unary_potentials, pairwise_params, edges):
        n_states = unary_potentials.shape[-1]
        if pairwise_params.shape == (n_states, n_states):
            # only one matrix given
            pairwise_potentials = np.repeat(
                pairwise_params[np.newaxis, :, :], edges.shape[0], axis=0
            )
        else:
            if pairwise_params.shape != (edges.shape[0], n_states, n_states):
                raise ValueError(
                    "Expected pairwise_params either to "
                    "be of shape n_states x n_states "
                    "or n_edges x n_states x n_states, but"
                    " got shape %s" % repr(pairwise_params.shape)
                )
            pairwise_potentials = pairwise_params
        return n_states, pairwise_potentials

    def init_ogm(
        self,
        unary_potentials,
        pairwise_potentials,
        edges,
        reserveNumFactorsPerVariable=4,
        **kwargs
    ):
        """Inference with OpenGM backend.

        Parameters
        ----------
        unary_potentials : nd-array, shape (n_nodes, n_states)
            Unary potentials of energy function.

        pairwise_potentials : nd-array, shape (n_states, n_states) or (n_states, n_states, n_edges).
            Pairwise potentials of energy function.
            If the first case, edge potentials are assumed to be the same for all edges.
            In the second case, the sequence needs to correspond to the edges.

        edges : nd-array, shape (n_edges, 2)
            Graph edges for pairwise potentials, given as pair of node indices. As
            pairwise potentials are not assumed to be symmetric, the direction of
            the edge matters.

        init : nd-array
            Initial solution for starting inference (ignored by some algorithms).

        reserveNumFactorsPerVariable :
            reserve a certain number of factors for each variable can speed up
            the building of a graphical model.
            ( For a 2d grid with second order factors one should set this to 5
             4 2-factors and 1 unary factor for most pixels )

        """
        n_states, pairwise_potentials = self.validate_params(
            unary_potentials, pairwise_potentials, edges
        )
        n_nodes = len(unary_potentials)

        gm = opengm.gm(np.ones(n_nodes, dtype=opengm.label_type) * n_states)

        nFactors = int(n_nodes + edges.shape[0])
        gm.reserveFactors(nFactors)
        gm.reserveFunctions(nFactors, "explicit")

        # all unaries as one numpy array
        # (opengm's value_type == float64 but all types are accepted)
        unaries = np.require(unary_potentials, dtype=opengm.value_type) * -1.0
        # add all unart functions at once
        fidUnaries = gm.addFunctions(unaries)
        visUnaries = np.arange(n_nodes, dtype=opengm.label_type)
        # add all unary factors at once
        gm.addFactors(fidUnaries, visUnaries)

        # add all pariwise functions at once
        # - first axis of secondOrderFunctions iterates over the function)

        secondOrderFunctions = -np.require(pairwise_potentials, dtype=opengm.value_type)
        fidSecondOrder = gm.addFunctions(secondOrderFunctions)
        gm.addFactors(fidSecondOrder, edges.astype(np.uint64))

        return gm

    def create_grid_model_opengm(self, my_l, my_k):
        np.random.seed(1)
        length = my_l
        k = my_k

        unary_potentials = np.zeros((length * length, k))
        for x in range(length * length):
            unary_potentials[x] = np.random.random(k)

        n_edges = 0
        for x in range(length - 1):
            for y in range(length):
                n_edges += 2
        edges = np.zeros((n_edges, 2))
        count = 0
        for x in range((length - 1) * length):
            edges[count] = [x, x + 1]
            count += 1
            edges[count] = [x + 1, x]
            count += 1

        pairwise_potentials = np.zeros((n_edges, k, k))
        for x in range(n_edges):
            pairwise_potentials[x] = np.random.random((k, k))

        gm = self.init_ogm(
            unary_potentials=unary_potentials,
            pairwise_potentials=pairwise_potentials,
            edges=np.sort(edges),
        )

        return gm

    def test_speedup_multi(self):
        self.my_l = 8
        self.my_k = 8
        slow = True

        # CUDA must be "started up" in order to achieve max speeds
        # Approximate loss of X seconds on first run otherwise
        if t.cuda.is_available():
            init_mn = self.create_grid_model(is_cuda=True, my_l=2, my_k=2)
            init_bp = TorchMatrixBeliefPropagator(
                markov_net=init_mn, is_cuda=True, var_on=False
            )
            init_bp.set_max_iter(1000)
            init_bp.infer(display="off")

        t_total0 = time.time()
        print("k = %d" % self.my_k)
        if slow:
            print(
                "length\tTorch-CUDA\tTorch-Py\tSparse-Py\tLoop-Py \tOpenGM\tSeq-Loop-BP\tBuild Time"
            )
        else:
            print(
                "length\tTorch-CUDA\tTorch-Py\tSparse-Py\tOpenGM\tSeq-Loop-BP\tBuild Time"
            )

        while self.my_k <= 64:
            t_prime0 = time.time()

            if t.cuda.is_available():
                cuda_mn = self.create_grid_model(
                    is_cuda=True, my_l=self.my_l, my_k=self.my_k
                )
                cuda_bp = TorchMatrixBeliefPropagator(
                    markov_net=cuda_mn, is_cuda=True, var_on=False
                )
                cuda_bp.set_max_iter(1000)

                t0 = time.time()
                cuda_bp.infer(display="off")
                t1 = time.time()
                cuda_bp_time = t1 - t0
            else:
                cuda_bp_time = np.nan  # if cuda is not available

            mn = self.create_grid_model(is_cuda=False, my_l=self.my_l, my_k=self.my_k)
            bp = TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=False, var_on=False)
            bp.set_max_iter(1000)
            t0 = time.time()
            bp.infer(display="off")
            t1 = time.time()
            bp_time = t1 - t0

            old_mn = self.create_grid_model_old(my_l=self.my_l, my_k=self.my_k)
            old_bp = MatrixBeliefPropagator(old_mn)
            old_bp.set_max_iter(1000)
            t0 = time.time()
            old_bp.infer(display="off")
            t1 = time.time()
            old_bp_time = t1 - t0

            t0 = time.time()
            # infer with asynchronous BP
            seq_bp = BeliefPropagator(old_mn)
            seq_bp.set_max_iter(1000)
            t0 = time.time()
            seq_bp.sequential_infer(display="off")
            t1 = time.time()
            sequential_bp_time = t1 - t0

            opengm_mn = self.create_grid_model_opengm(my_l=self.my_l, my_k=self.my_k)
            t0 = time.time()
            # infer with parallel opengm
            inf = opengm.inference.BeliefPropagation(
                opengm_mn,
                parameter=opengm.InfParam(
                    steps=1000, damping=0, convergenceBound=1e-8, isAcyclic=False
                ),
            )
            inf.infer()
            t1 = time.time()
            opengm_time = t1 - t0

            if slow:
                slow_bp = BeliefPropagator(old_mn)
                slow_bp.set_max_iter(1000)
                t0 = time.time()
                slow_bp.infer(display="off")
                t1 = time.time()
                slow_bp_time = t1 - t0

            t_prime1 = time.time()

            start_time = (
                (t_prime1 - t_prime0)
                - cuda_bp_time
                - bp_time
                - old_bp_time
                - opengm_time
                - sequential_bp_time
            )
            if slow:
                start_time = start_time - slow_bp_time

            if slow:
                print(
                    "%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f"
                    % (
                        self.my_l,
                        cuda_bp_time,
                        bp_time,
                        old_bp_time,
                        slow_bp_time,
                        opengm_time,
                        sequential_bp_time,
                        start_time,
                    )
                )
            else:
                print(
                    "%d\t%f\t%f\t%f\t%f\t%f\t%f"
                    % (
                        self.my_l,
                        cuda_bp_time,
                        bp_time,
                        old_bp_time,
                        opengm_time,
                        sequential_bp_time,
                        start_time,
                    )
                )
            self.my_l *= 2
            """
            try:
                self.create_grid_model(is_cuda=True, my_l=self.my_l, my_k=self.my_k)
            except Exception:
                self.my_l = 8
                self.my_k *= 2
                print("k = %d" % self.my_k)
            """
            if (self.my_l * self.my_k) >= 8192:
                self.my_l = 8
                self.my_k *= 2
                print("k = %d" % self.my_k)
        t_total1 = time.time()
        total_time = t_total1 - t_total0
        print("Total runtime took %f" % (total_time))

    def test_cuda_time_loss(self):
        if not t.cuda.is_available():
            print("CUDA is not available.")
            return

        self.my_l = 8
        self.my_k = 8
        init_mn = self.create_grid_model(is_cuda=True, my_l=self.my_l, my_k=self.my_k)
        init_bp = TorchMatrixBeliefPropagator(
            markov_net=init_mn, is_cuda=True, var_on=False
        )
        init_bp.set_max_iter(1000)
        init_time0 = time.time()
        init_bp.infer(display="off")
        init_time1 = time.time()

        init_time = init_time1 - init_time0

        cuda_mn = self.create_grid_model(is_cuda=True, my_l=self.my_l, my_k=self.my_k)
        cuda_bp = TorchMatrixBeliefPropagator(
            markov_net=cuda_mn, is_cuda=True, var_on=False
        )
        cuda_bp.set_max_iter(1000)
        t0 = time.time()
        cuda_bp.infer(display="off")
        t1 = time.time()
        cuda_bp_time = t1 - t0

        time_diff = init_time - cuda_bp_time
        print("Total time lost to CUDA initialization %f" % (time_diff))
        print(
            "In this case, this is a factor of %f times slower"
            % (time_diff / cuda_bp_time)
        )

    def test_autograd(self):
        if not t.cuda.is_available():
            print("CUDA is not available.")
            return

        is_cuda = True
        mn = self.create_grid_model(is_cuda=is_cuda, my_l=64, my_k=8)
        bp = TorchMatrixBeliefPropagator(markov_net=mn, is_cuda=is_cuda, var_on=False)
        bp.infer(display="full")

        bp.load_beliefs()
        print(bp.autograd_unary())


if __name__ == "__main__":
    unittest.main()
