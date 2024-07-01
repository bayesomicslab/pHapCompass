import argparse
import networkx as nx
from data.input_handler import InputHandler
from data.configuration import Configuration
from algorithm.haplotype_assembly import HaplotypeAssembly
from models.fragment_graph import FragmentGraph
from models.quotient_graph import QuotientGraph
from models.factor_graph import Factorgraph
from algorithm.chordal_contraction import chordal_contraction_cycle_base, chordal_contraction
from utils.utils import *
from algorithm.inference import *


def test():
    inp = 1


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Haplotype Assembly and Phasing")
    parser.add_argument("-d", "--data_path", type=str, required=False, help="Path to the input data", default=None)
    parser.add_argument("-p", "--ploidy", type=int, required=True, help="Ploidy of the organism", default=3)
    parser.add_argument("-g", "--genotype_path", type=str, required=True, help="Path to the genotype data",
                        default='example/genotype.txt')
    parser.add_argument("-a", "--alleles", required=False, nargs='*', help="List of alleles (optional)", default={0, 1})
    parser.add_argument("--error_rate", type=float, required=True, help="Error rate", default=0.001)
    parser.add_argument("--epsilon", type=float, required=True, help="epsilon", default=0.0001)
    parser.add_argument("-v", "--vcf_path", required=True, help="VCF file for called variants (string)")
    parser.add_argument("-b", "--bam_path", help="sam or bam file input (string)")
    parser.add_argument("-o", "--output_path", required=True, help="output path")
    parser.add_argument("-r", "--root_dir", required=True, help="root directory")
    
    # parser.add_argument("--epsilon", help="epsilon in computing prob.")

    args = parser.parse_args()

    # Initialize classes with parsed arguments
    input_handler = InputHandler(args)

    config = Configuration(args.ploidy, args.error_rate, args.epsilon, input_handler.alleles)
    
    fragment_model = FragmentGraph(input_handler.data_path, input_handler.genotype_path, input_handler.ploidy, input_handler.alleles)
    frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)
    
    # fragment_model = FragmentGraph(args.data_path, args.genotype_path, args.ploidy, input_handler.alleles)
    # frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)

    plot_graph(frag_graph)
    print('Fragment Graph constructed.')

    quotient_g = QuotientGraph(frag_graph).construct(fragment_list, input_handler, config)
    plot_graph(quotient_g)
    print('Quotient Graph constructed.')

    # qg = chordal_contraction_cycle_base(quotient_g, fragment_list, input_handler, config)
    qg = chordal_contraction(quotient_g, fragment_list, input_handler, config)
    plot_graph(qg)
    
    print('Chordal Graph constructed.')
    
    plot_graph(qg)
    
    # for edge in quotient_g.edges(data=True):
    #     print(edge)
    # quotient_g.nodes(data=True)
    #
    # G = nx.Graph()
    #
    # entropies = np.random.uniform(0, 1, size=14)
    # G.add_weighted_edges_from([('1', '2', {'original_order': ('1', '2'), 'endtropy': entropies[0]}),
    #                            ('1', '5', {'original_order': ('1', '5'), 'endtropy': entropies[1]}),
    #                            ('2', '3', {'original_order': ('2', '3'), 'endtropy': entropies[2]}),
    #                            ('2', '6', {'original_order': ('2', '6'), 'endtropy': entropies[3]}),
    #                            ('3', '7', {'original_order': ('3', '7'), 'endtropy': entropies[4]}),
    #                            ('3', '4', {'original_order': ('3', '4'), 'endtropy': entropies[5]}),
    #                            ('4', '5', {'original_order': ('4', '5'), 'endtropy': entropies[6]}),
    #                            ('4', '8', {'original_order': ('4', '8'), 'endtropy': entropies[7]}),
    #                            ('5', '9', {'original_order': ('5', '9'), 'endtropy': entropies[8]}),
    #                            ('6', '7', {'original_order': ('6', '7'), 'endtropy': entropies[9]}),
    #                            ('7', '8', {'original_order': ('7', '8'), 'endtropy': entropies[10]}),
    #                            ('8', '9', {'original_order': ('8', '9'), 'endtropy': entropies[11]}),
    #                            ('5', '10', {'original_order': ('5', '10'), 'endtropy': entropies[12]}),
    #                            ('9', '10', {'original_order': ('9', '10'), 'endtropy': entropies[13]})])
    # plot_graph(G)
    # G.edges()

    # interpreter = '/home/FCAM/mhosseini/anaconda3/envs/t2t/bin/python3'
    # import networkx as nx
    # import networkit as nk
    # nodes = ['1-2', '1-3', '2-3', '3-4', '4-5', '4-7', '5-6', '6-8', '7-8']
    # edges = [('1-2', '1-3'), ('1-2', '2-3'), ('1-3', '2-3'), ('1-3', '3-4'), ('2-3', '3-4'), ('3-4', '4-5'),
    # ('3-4', '4-7'), ('4-5', '4-7'), ('4-5', '5-6'), ('4-7', '7-8'), ('5-6', '6-8'), ('6-8', '7-8')]
    # graphnx = nx.Graph()
    # graphnx.add_nodes_from(nodes)
    # graphnx.add_edges_from(edges)
    #
    # graphnk, reverse_map = nx2nk(graphnx)
    # cliques = networkit_find_cliques(graphnk)
    # cycles_g4 = [cyc for cyc in list(simple_cycles(graphnk, 0)) if len(cyc) > 4]
    # cycles_g4_unique = [list(x) for x in set(tuple(x) for x in cycles_g4)]
    #
    # chordless_cycles = list(nx.chordless_cycles(tempnx))
    #
    #
    # nodes_dict = dict(quotient_g.nodes(data=True))
    # edges_list = list(quotient_g.edges(data=True))
    # edges_dict = {str(item[0:2]): item[2] for item in edges_list}

    # nodes = list(quotient_g.nodes())
    # edges = list(quotient_g.edges())
    # new_g = nx.Graph()
    # new_g.add_nodes_from(nodes)
    # new_g.add_edges_from(edges)
    # graphnk, reverse_map = nx2nk(new_g)
    
    # edges_w_data = list(quotient_g.edges(data=True))
    # dict(quotient_g.edges(data=True))
    
    # plot_graph(qg)
    # qg = quotient_g.copy()
    
    factor_graph = Factorgraph(config.ploidy, config.error_rate, config.epsilon).construct(qg, fragment_list)

    

    beliefs = factor_graph_inference(factor_graph)

    variables = [node for node in factor_graph.nodes() if not isinstance(node, DiscreteFactor)]

    # Step 3: Perform inference
    for var in variables:
        print(f"Belief for {var}:")
        result = beliefs.query(variables=list(var))
        for state_index, prob in enumerate(result.values):
            print(f"Value {result.variables} = {state_index}: {prob}")


    for var in variables:
        print(f"Belief for {var}:")
        print(result[var])

    print(result)
    for var in variables:
        print(f"Belief for {var}:")
        print(f"Values: {result[var].values}")
        print(f"Variables: {result[var].variables}")
    # beliefs = BeliefPropagation(factor_graph)

    result12 = beliefs.query(variables=[list(variables)[0]])
    print(result12)
    # marginals, max_phasings = give_marginals(factor_graph, qg, beliefs)
    #
    # max_phase, positions = query_paths_gibbs_max(fragment_list, qg, beliefs, n_samples=1000)
    # h_df = creat_vcf(max_phase, positions, config)
    # print(h_df)
    
    # va_inference = VariableElimination(factor_graph)
    # result = va_inference.query(variables=['1-2'])
    
    from hmmlearn import hmm
    
    # Define states and observations
    states = ["Rainy", "Sunny"]
    observations = ["Walk", "Shop", "Clean"]
    n_states = len(states)
    n_observations = len(observations)
    
    # Example emission probabilities
    # The rows correspond to states, and the columns correspond to observations
    emission_probabilities = np.array([
        [0.1, 0.4, 0.5],  # Emission probabilities for "Rainy"
        [0.6, 0.3, 0.1]   # Emission probabilities for "Sunny"
    ])
    
    # Generate a random sequence of observations (for example purposes)
    # Normally, you would use your actual observation sequence data
    np.random.seed(42)
    sequence_length = 100
    X = np.random.choice(n_observations, sequence_length).reshape(-1, 1)
    
    # Initialize the HMM
    model = hmm.MultinomialHMM(n_components=n_states, n_iter=100, random_state=42)
    
    # Set the emission probabilities
    model.emissionprob_ = emission_probabilities
    
    # Fit the model to the observation sequence
    model.fit(X)
    
    # Extract transition probabilities
    transition_probabilities = model.transmat_

    import torch
    from torch import nn
    import torch_treecrf

    # Define the number of states and possible observations
    num_states = 3  # Number of hidden states
    num_observations = 3  # Number of observation symbols

    # Example transition and emission probabilities (replace with actual data)
    transition_probs = torch.tensor([[0.5, 0.2, 0.3],
                                     [0.3, 0.5, 0.2],
                                     [0.2, 0.3, 0.5]], dtype=torch.float)

    emission_probs = torch.tensor([[0.6, 0.3, 0.1],
                                   [0.2, 0.5, 0.3],
                                   [0.1, 0.3, 0.6]], dtype=torch.float)

    # Define the CRF model
    class HMMCRF(nn.Module):
        def __init__(self, num_states, num_observations):
            super(HMMCRF, self).__init__()
            self.num_states = num_states
            self.num_observations = num_observations
            self.crf = CRF(num_tags=num_states, batch_first=True)
        
            # Initialize the emission probabilities as a learnable parameter
            self.emission_probs = nn.Parameter(emission_probs)
    
        def forward(self, observations):
            # Convert observations to one-hot encoding
            one_hot_observations = nn.functional.one_hot(observations, num_classes=self.num_observations).float()
        
            # Calculate emission scores
            emission_scores = torch.matmul(one_hot_observations, self.emission_probs.t())
        
            # Use the CRF to get the most likely state sequence
            best_path = self.crf.decode(emission_scores)
            return best_path
    
        def log_likelihood(self, observations, states):
            one_hot_observations = nn.functional.one_hot(observations, num_classes=self.num_observations).float()
            emission_scores = torch.matmul(one_hot_observations, self.emission_probs.t())
            log_likelihood = self.crf(emission_scores, states)
            return log_likelihood

    # Create example observation sequences and state sequences
    observations = torch.tensor([[0, 1, 2], [2, 0, 1]], dtype=torch.long)  # Example observations
    states = torch.tensor([[0, 1, 2], [2, 0, 1]], dtype=torch.long)  # Example true state sequences

    # Initialize the HMM-CRF model
    model = HMMCRF(num_states=num_states, num_observations=num_observations)

    # Compute the most likely state sequences
    with torch.no_grad():
        best_paths = model(observations)
        print("Best paths:", best_paths)

    # Compute the log likelihood of the true state sequences
    log_likelihood = model.log_likelihood(observations, states)
    print("Log likelihood:", log_likelihood)

    # Example training loop
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)
    num_epochs = 100

    for epoch in range(num_epochs):
        optimizer.zero_grad()
        loss = -model.log_likelihood(observations, states).mean()  # Negative log likelihood
        loss.backward()
        optimizer.step()
    
        if epoch % 10 == 0:
            print(f"Epoch {epoch}, Loss: {loss.item()}")
    

