import os
import mrftools
import torch as t
from torch.autograd import Variable
import numpy as np
from scipy.sparse import coo_matrix

class IndependentMarkovNet():
    def __init__(self, markov_net):
        self.mn = markov_net
    
    def create_matrices(self, rev_color_dict):
        """
        Create matrix representations of the MRF structure and potentials to allow inference to be done via 
        matrix operations
        :return: None
        """
        self.mn.matrix_mode = True

        self.mn.max_states = max([len(x) for x in self.mn.unary_potentials.values()])
        self.mn.unary_mat = -np.inf * np.ones((self.mn.max_states, len(self.mn.variables)), dtype=self.mn.dtype)
        self.mn.degrees = np.zeros(len(self.mn.variables), dtype=self.mn.dtype)

        # var_index allows looking up the numerical index of a variable by its hashable name
        self.mn.var_index = dict()
        self.mn.var_list = []

        message_num = 0
        for var in self.mn.variables:
            potential = self.mn.unary_potentials[var]
            self.mn.unary_mat[0:len(potential), message_num] = potential
            self.mn.var_index[var] = message_num
            self.mn.var_list.append(var)
            self.mn.degrees[message_num] = len(self.mn.neighbors[var])
            message_num += 1

        # set up pairwise tensor
        self.mn.num_edges = 0
        for var in self.mn.variables:
            for neighbor in self.mn.neighbors[var]:
                if var < neighbor:
                    self.mn.num_edges += 1

        self.mn.edge_pot_tensor = -np.inf * np.ones((self.mn.max_states, self.mn.max_states, 2 * self.mn.num_edges), dtype=self.mn.dtype)
        self.mn.message_index = {}

        # set up sparse matrix representation of adjacency
        from_rows = []
        from_cols = []
        to_rows = []
        to_cols = []

        message_num = 0
        for var in self.mn.variables:
            for neighbor in self.mn.neighbors[var]:
                if var < neighbor:
                    # for each unique edge

                    potential = self.mn.get_potential((var, neighbor))
                    dims = potential.shape

                    # store copies of the potential for each direction messages can travel on the edge
                    # forward
                    self.edge_pot_tensor[0:dims[1], 0:dims[0], message_num] = potential.T
                    # and backward
                    self.mn.edge_pot_tensor[0:dims[0], 0:dims[1], message_num + self.mn.num_edges] = potential

                    # get numerical index of var and neighbor
                    var_i = self.mn.var_index[var]
                    neighbor_i = self.mn.var_index[neighbor]

                    # store that the forward slice represents a message from var
                    from_rows.append(message_num)
                    from_cols.append(var_i)

                    # store that the backward slice represents a message from neighbor
                    from_rows.append(message_num + self.mn.num_edges)
                    from_cols.append(neighbor_i)

                    # store that the forward slice represents a message to neighbor
                    to_rows.append(message_num)
                    to_cols.append(neighbor_i)

                    # store that the backward slice represents a message to var
                    to_rows.append(message_num + self.mn.num_edges)
                    to_cols.append(var_i)

                    self.mn.message_index[(var, neighbor)] = message_num

                    message_num += 1

        # generate a sparse matrix representation of the message indices to variables that receive messages
        self.mn.message_to_map = coo_matrix((np.ones(len(to_rows), dtype=self.mn.dtype), (to_rows, to_cols)),
                                         (2 * self.mn.num_edges, len(self.mn.variables)))

        self.mapping_matrices = {}
        self.message_to_matrices = {}
        self.message_from_matrices = {}

        for color in rev_color_dict.keys():
            for var in rev_color_dict[color]:
                for neighbor in self.mn.neighbors[var]:
                    if var < neighbor:
                        # for each unique edge

                        potential = self.mn.get_potential((var, neighbor))
                        dims = potential.shape

                        # store copies of the potential for each direction messages can travel on the edge
                        # forward
                        self.mn.edge_pot_tensor[0:dims[1], 0:dims[0], message_num] = potential.T
                        # and backward
                        self.mn.edge_pot_tensor[0:dims[0], 0:dims[1], message_num + self.mn.num_edges] = potential

                        # get numerical index of var and neighbor
                        var_i = self.mn.var_index[var]
                        neighbor_i = self.mn.var_index[neighbor]

                        # store that the forward slice represents a message from var
                        from_rows.append(message_num)
                        from_cols.append(var_i)

                        # store that the backward slice represents a message from neighbor
                        from_rows.append(message_num + self.mn.num_edges)
                        from_cols.append(neighbor_i)

                        # store that the forward slice represents a message to neighbor
                        to_rows.append(message_num)
                        to_cols.append(neighbor_i)

                        # store that the backward slice represents a message to var
                        to_rows.append(message_num + self.mn.num_edges)
                        to_cols.append(var_i)

                        self.mn.message_index[(var, neighbor)] = message_num

                        message_num += 1
            self.mapping_matrices[color] = coo_matrix((np.ones(len(to_rows), dtype=self.mn.dtype), (to_rows, to_cols)),
                                         (2 * self.mn.num_edges, len(self.mn.variables)))
            self.message_to_matrices[color] = np.zeros(2 * self.mn.num_edges, dtype=np.intp)
            self.message_to_matrices[color][to_rows] = to_cols
            self.message_from_matrices[color] = np.zeros(2 * self.mn.num_edges, dtype=np.intp)
            self.message_from_matrices[color][from_rows] = from_cols


        # store an array that lists which variable each message is sent to
        self.mn.message_to = np.zeros(2 * self.mn.num_edges, dtype=np.intp)
        self.mn.message_to[to_rows] = to_cols

        # store an array that lists which variable each message is received from
        self.mn.message_from = np.zeros(2 * self.mn.num_edges, dtype=np.intp)
        self.mn.message_from[from_rows] = from_cols


class IndependentTorchMatrixBeliefPropagator(mrftools.TorchMatrixBeliefPropagator):
    def __init__(self, markov_net, is_cuda, var_on, color_dict):
        super().__init__(markov_net, is_cuda, var_on)
        self.color_dict = color_dict
        # self.rev_colors = {v: k for k, v in color_dict.items()}
        self.rev_colors = {}
        for v in self.mn.variables:
            if color_dict[v] in self.rev_colors:
                self.rev_colors[color_dict[v]] = self.rev_colors[color_dict[v]].append(v)
            else:
                self.rev_colors[color_dict[v]] = [v]
        self.ind_mn = IndependentMarkovNet(self.mn)
        
        self.ind_mn.create_matrices(self.rev_colors)

    def compute_independent_beliefs(self, color):
        if not self.fully_conditioned:
            self.belief_mat = self.mn.unary_mat + self.augmented_mat
            self.belief_mat += mrftools.TorchMatrixBeliefPropagator.sparse_dot(self.message_mat, self.mn.mapping_matrices[color])

            self.belief_mat -= mrftools.TorchMatrixBeliefPropagator.torch_logsumexp(self.belief_mat, 0)
        return
        
    def update_independent_messages(self, color):
        self.compute_independent_beliefs(color)

        # Using the beliefs as the sum of all incoming log messages, subtract the outgoing messages and add the edge
        # potential.
        adjusted_message_prod = self.mn.edge_pot_tensor + self.belief_mat[:, self.mn.message_from_matrices[color]] \
                                - t.cat((self.message_mat[:, self.mn.num_edges:],
                                             self.message_mat[:, :self.mn.num_edges]), 1)

        # Have to convert NaNs to negative infinity at this point
        adjusted_message_prod = mrftools.TorchMatrixBeliefPropagator.torch_nan_to_neginf(adjusted_message_prod)

        messages = t.squeeze(mrftools.TorchMatrixBeliefPropagator.torch_logsumexp(adjusted_message_prod, 1))
        messages = mrftools.TorchMatrixBeliefPropagator.torch_nan_to_zero(messages - t.max(messages, 0)[0])

        change = t.sum(mrftools.TorchMatrixBeliefPropagator.torch_nan_to_zero(t.abs(messages - self.message_mat)))

        self.message_mat = messages

        return change
    
    def independent_inference(self, tolerance=1e-8, display='iter'):
        change = float('inf')
        iteration = 0
#        while change > tolerance and iteration < self.max_iter:
        while (iteration < self.max_iter) and (change > tolerance) :
            for color in self.node_colors.values():
                change = self.update_independent_messages(color)
            if display == "full":
                energy_func = self.compute_energy_functional()
                disagreement = self.compute_inconsistency()
                dual_obj = self.compute_dual_objective()
                print("Iteration %d, change in messages %f. Calibration disagreement: %f, "
                      "energy functional: %f, dual obj: %f" % (iteration, change, disagreement, energy_func, dual_obj))
            elif display == "iter":
                print("Iteration %d, change in messages %f." % (iteration, change))
            iteration += 1
        if display == 'final' or display == 'full' or display == 'iter':
            print("Belief propagation finished in %d iterations." % iteration)
        return

# Graph Coloring on Chordal Graph
#
# Greedy Coloring given by Networkx -> color_dict = greedy_color(G)
#
# Manual Coding
# def color_graph(G):
#     degree = G.degree
#     colors = range(degree)
#     for node in list(G.nodes):
#         pass

#     return G

# Inference

if __name__ =="__main__":
    print("Done")