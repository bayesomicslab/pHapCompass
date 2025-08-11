import torch
import warnings

class FactorGraph:
    def __init__(self, is_cuda=False):
        self.edge_potentials = dict()
        self.unary_potentials = dict()
        self.neighbors = dict()
        self.variables = set()
        self.num_states = dict()

        self.tree_probabilities = dict()

        self.max_states = None
        self.message_to_map = None
        self.message_to = None
        self.message_from = None
        self.var_index = None
        self.var_list = None
        self.unary_mat = None
        self.edge_pot_tensor = None
        self.num_edges = None
        self.message_index = None
        self.degrees = None
        self.is_cuda = is_cuda

        self.matrix_mode = False

    def set_unary_factor(self, variable, potential):
        """
        Set the potential function for the unary factor. Implicitly declare variable.
        Must be called before setting edge factors.

        :param variable: name of the variable (can be any hashable object)
        :param potential: length-k vector of log potential values for the respective k states
        :return: None
        """
        self.unary_potentials[variable] = potential
        if variable not in self.variables:
            self.declare_variable(variable, potential.size()[0])
    
    def declare_variable(self, variable, num_states):
        """
        Indicate the existence of a variable

        :param variable: name of the variable (can be any hashable object)
        :param num_states: integer number of states the variable can take
        :return: None
        """
        if variable not in self.variables:
            self.variables.add(variable)
            self.neighbors[variable] = set()
            self.num_states[variable] = num_states
        else:
            print("Warning: declaring a variable %s that was previously declared." % repr(variable))

    def set_edge_factor(self, edge, potential):
        """
        Set a factor by inputting the involved variables then the potential function.
        The potential function should be a np matrix.

        :param edge: 2-tuple of the variables in the edge. Can be in any order.
        :param potential: k1 by k2 matrix of potential values for the joint state of the two variables
        :return: None
        """
        assert potential.size() == (len(self.unary_potentials[edge[0]]), len(self.unary_potentials[edge[1]])), \
            "potential size %d, %d incompatible with unary sizes %d, %d" % \
            (potential.size()[0], potential.size()[1], len(self.unary_potentials[edge[0]]),
             len(self.unary_potentials[edge[1]]))

        if edge[0] < edge[1]:
            self.edge_potentials[edge] = potential
        else:
            self.edge_potentials[(edge[1], edge[0])] = potential.t()

        self.neighbors[edge[0]].add(edge[1])
        self.neighbors[edge[1]].add(edge[0])

    def get_potential(self, key):
        """
        Return the potential between pair[0] and pair[1]. If (pair[1], pair[0]) is in our dictionary instead,
        return the transposed potential.

        :param key: name of the key whose potential to get. Can either be a variable name or a pair of variables (edge)
        :return potential table for the key
        """
        if key in self.edge_potentials:
            return self.edge_potentials[key]
        else:
            return self.edge_potentials[(key[1], key[0])].T

    def get_neighbors(self, variable):
        """
        Return the neighbors of variable.

        :param variable: name of variable
        :return: set of neighboring variables connected in MRF
        """
        return self.neighbors[variable]

    def evaluate_state(self, states):
        """
        Evaluate the energy of a state. states should be a dictionary of variable: state (int) pairs.

        :param states: dictionary of variable states with a key-value pair for each variable
        :return: MRF energy value for the state as a float
        """
        energy = 0.0
        for var in self.variables:
            energy += self.unary_potentials[var][states[var]]

            for neighbor in self.neighbors[var]:
                if var < neighbor:
                    energy += self.get_potential((var, neighbor))[states[var], states[neighbor]]

        return energy

    def set_unary_mat(self, unary_mat):
        """
        Set the matrix representation of the unary potentials

        :param unary_mat: (num vars) by (max states) of unary potentials
        :return: None
        """
        assert torch.eq(self.unary_mat.size(), unary_mat.size())
        self.unary_mat[:, :] = unary_mat
        if self.is_cuda:
            self.unary_mat = self.unary_mat.cuda()

    def set_edge_tensor(self, edge_tensor):
        """
        Set the tensor representation of the edge potentials
        :param edge_tensor: (max states) by (max states) by (num edges) tensor of the edge potentials
        :return: None
        """
        if torch.eq(self.edge_pot_tensor.size(), edge_tensor.size()):
            self.edge_pot_tensor[:, :, :] = edge_tensor
        else:
            mirrored_edge_tensor = torch.cat((edge_tensor, edge_tensor.permute(1, 0, 2)), 2)
            assert torch.eq(self.edge_pot_tensor.size(), mirrored_edge_tensor.size())

            self.edge_pot_tensor[:, :, :] = mirrored_edge_tensor
    
    def create_matrices(self):
        """
        Create matrix representations of the MRF structure and potentials to allow inference to be done via
        matrix operations
        :return: None
        """
        self.matrix_mode = True

        self.max_states = max([len(x) for x in self.unary_potentials.values()])
        self.unary_mat = -float('inf') * torch.ones(self.max_states, len(self.variables)).double()

        if self.is_cuda:
            self.unary_mat = self.unary_mat.cuda()

        self.degrees = torch.DoubleTensor(len(self.variables)).zero_()

        if self.is_cuda:
            self.degrees = self.degrees.cuda()

        # var_index allows looking up the numerical index of a variable by its hashable name
        self.var_index = dict()
        self.var_list = []

        message_num = 0
        for var in self.variables:
            potential = self.unary_potentials[var]
            self.unary_mat[0:len(potential), message_num] = potential
            self.var_index[var] = message_num
            self.var_list.append(var)
            self.degrees[message_num] = len(self.neighbors[var])
            message_num += 1

        # set up pairwise tensor
        self.num_edges = 0
        for var in self.variables:
            for neighbor in self.neighbors[var]:
                if var < neighbor:
                    self.num_edges += 1

        self.edge_pot_tensor = -float('inf') * torch.ones(self.max_states, self.max_states, 2 * self.num_edges).double()

        if self.is_cuda:
            self.edge_pot_tensor = self.edge_pot_tensor.cuda()

        self.message_index = {}

        # set up sparse matrix representation of adjacency
        self.from_rows = []
        self.from_cols = []
        self.to_rows = []
        self.to_cols = []

        message_num = 0
        for var in self.variables:
            for neighbor in self.neighbors[var]:
                if var < neighbor:
                    # for each unique edge

                    potential = self.get_potential((var, neighbor))
                    dims = potential.shape

                    # store copies of the potential for each direction messages can travel on the edge
                    # forward
                    self.edge_pot_tensor[0:dims[1], 0:dims[0], message_num] = potential.t()
                    # and backward
                    self.edge_pot_tensor[0:dims[0], 0:dims[1], message_num + self.num_edges] = potential

                    # get numerical index of var and neighbor
                    var_i = self.var_index[var]
                    neighbor_i = self.var_index[neighbor]

                    # store that the forward slice represents a message from var
                    self.from_rows.append(message_num)
                    self.from_cols.append(var_i)

                    # store that the backward slice represents a message from neighbor
                    self.from_rows.append(message_num + self.num_edges)
                    self.from_cols.append(neighbor_i)

                    # store that the forward slice represents a message to neighbor
                    self.to_rows.append(message_num)
                    self.to_cols.append(neighbor_i)

                    # store that the backward slice represents a message to var
                    self.to_rows.append(message_num + self.num_edges)
                    self.to_cols.append(var_i)

                    self.message_index[(var, neighbor)] = message_num

                    message_num += 1
        # generate a sparse matrix representation of the message indices to variables that receive messages
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.message_to_map = torch.sparse.DoubleTensor(
                torch.LongTensor([self.to_rows, self.to_cols]),
                torch.ones(len(self.to_rows)).double(),
                torch.Size([2 * self.num_edges, len(self.variables)])
            )

        if self.is_cuda:
            self.message_to_map = self.message_to_map.cuda()
        
        # store an array that lists which variable each message is sent to
        self.message_to = torch.LongTensor(2 * self.num_edges).zero_()
        self.message_to[torch.LongTensor(self.to_rows)] = torch.LongTensor(self.to_cols)

        if self.is_cuda:
            self.message_to = self.message_to.cuda()
            self.message_to[torch.LongTensor(self.to_rows).cuda()] = torch.LongTensor(self.to_cols).cuda()

        # store an array that lists which variable each message is received from
        self.message_from = torch.LongTensor(2 * self.num_edges).zero_()
        self.message_from[torch.LongTensor(self.from_rows)] = torch.LongTensor(self.from_cols)

        if self.is_cuda:
            self.message_from = self.message_from.cuda()
            self.message_from[torch.LongTensor(self.from_rows).cuda()] = torch.LongTensor(self.from_cols).cuda()

if __name__ == "__main__":
    # Example usage
    fg = FactorGraph()
    fg.declare_variable('A', 3)
    fg.declare_variable('B', 2)
    fg.set_unary_factor('A', torch.tensor([0.1, 0.2, 0.3]))
    fg.set_unary_factor('B', torch.tensor([0.4, 0.5]))
    fg.set_edge_factor(('A', 'B'), torch.tensor([[0.6, 0.7], [0.8, 0.9], [1.0, 1.1]]))
    
    print(fg.evaluate_state({'A': 1, 'B': 0}))  # Example state evaluation
    print(fg.get_neighbors('A'))  # Should return neighbors of variable 'A'
    print(fg.get_potential(('A', 'B')))  # Should return the potential
    print(fg.get_potential(('B', 'A')))  # Should return

    fg.create_matrices()  # Create matrix representations
    print(fg.unary_mat)  # Should print the unary matrix
