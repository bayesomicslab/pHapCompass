import os
import argparse
import sys
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
from test.FFBS import generate_hmm_with_weights_and_emissions




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


    class Args:
        def __init__(self):
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = 'example/test.txt'
            self.bam_path = 'example/example.bam'
            self.genotype_path = 'example/genotype.txt'
            self.ploidy = 3
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




    error_rate = 0.001
    # state_names, transition_matrix, emission_prob_matrix, emission_index_map = generate_hmm_with_weights_and_emissions2(qg, error_rate)


    # for nn in qg.nodes(data=True):
    #     print(nn)



#     import numpy as np
#     import itertools


# def find_matchings(nodes_part1, nodes_part2):
#     # Sort both parts and remember the original indices.
#     sorted_part1 = sorted(enumerate(nodes_part1), key=lambda x: x[1])
#     sorted_part2 = sorted(enumerate(nodes_part2), key=lambda x: x[1])
    
#     # Split nodes by type and collect their original indices.
#     def split_by_type(sorted_nodes):
#         grouped = {}
#         for idx, t in sorted_nodes:
#             if t not in grouped:
#                 grouped[t] = []
#             grouped[t].append(idx)
#         return grouped
    
#     grouped_part1 = split_by_type(sorted_part1)
#     grouped_part2 = split_by_type(sorted_part2)
#     if grouped_part1.keys() != grouped_part2.keys():
#         return []
#     if any([len((grouped_part1[i])) != len((grouped_part2[i])) for i in grouped_part1.keys()]):
#         return []
#     # Start with a single empty matching.
#     matchings = [[]]
#     for node_type, indices1 in grouped_part1.items():
#         indices2 = grouped_part2[node_type]
        
#         # For each current matching, extend it with all possible permutations for the current type.
#         new_matchings = []
#         for perm in itertools.permutations(indices2, len(indices2)):
#             for current_matching in matchings:
#                 # Add new matching to the results only if it doesn't conflict with the current matching.
#                 if all((i1, i2) not in current_matching for i1, i2 in zip(indices1, perm)):
#                     new_matchings.append(current_matching + list(zip(indices1, perm)))
#         matchings = new_matchings
    
#     return matchings


# def str_2_phas_1(phasing, ploidy):
#     return np.array([int(p) for p in [*phasing]]).reshape(ploidy, -1)


# def phas_2_str(phas):
#     return ''.join([str(ph) for ph in list(np.ravel(phas))])


# def find_phasings_matches(ff, sf, common_ff, common_sf):
#     templates = []
#     all_local = find_matchings(list(ff[:, -1]), list(sf[:, 0]))
#     for al in all_local:
#         ff_ordering = [ii[0] for ii in al]
#         sf_ordering = [ii[1] for ii in al]
#         assert any(ff[ff_ordering, common_ff] == sf[sf_ordering, common_sf])
#         temp = np.hstack([ff[ff_ordering, :], sf[sf_ordering, 1:]])
#         byte_set = {a.tobytes() for a in templates}
#         if temp.tobytes() not in byte_set:
#             templates.append(temp)
#     return templates

# # Function to find the common element and its index between two nodes
# def find_common_element_and_index(node1, node2):
#     # Split the node names into components (e.g., '1-2' -> ['1', '2'])
#     node1_parts = node1.split('-')
#     node2_parts = node2.split('-')
    
#     # Find the common element between node1 and node2
#     common_element = None
#     for part in node1_parts:
#         if part in node2_parts:
#             common_element = part
#             break
    
#     if common_element is None:
#         raise ValueError(f"No common element found between {node1} and {node2}")
    
#     # Find the index of the common element in both nodes
#     common_ff = node1_parts.index(common_element)
#     common_sf = node2_parts.index(common_element)
    
#     return common_ff, common_sf


# def permute_rows(A):
#     # Generate all permutations of the rows
#     perms = list(itertools.permutations(A))
    
#     # Convert each permutation into a NumPy array and return as a list of matrices
#     permuted_matrices = [np.array(p) for p in perms]
    
#     return permuted_matrices



# # Function to generate state names and transition matrix with weighted transitions
# def generate_hmm_with_weights(qg, error_rate=0.001):
#     state_names = []  # To store the names of states in hmmg
#     state_index_map = {}  # Map state names to their index in the transition matrix
    
#     # Step 1: Generate the state names
#     for node, node_data in qg.nodes(data=True):
#         for state_key in node_data['weight'].keys():
#             state_name = f"{node}-{state_key}"
#             state_names.append(state_name)
    
#     # Create an index for each state
#     for i, state_name in enumerate(state_names):
#         state_index_map[state_name] = i

#     # Step 2: Initialize transition matrix
#     num_states = len(state_names)
#     transition_matrix = np.zeros((num_states, num_states))  # Start with a zero matrix

#     # Step 3: Define transitions between neighboring nodes with weights
#     for node1, node2 in qg.edges():
#         # Get the states of node1
#         for state_key1 in qg.nodes[node1]['weight'].keys():
#             state1 = f"{node1}-{state_key1}"
#             state1_idx = state_index_map[state1]
#             # Get the states of node2
#             for state_key2 in qg.nodes[node2]['weight'].keys():
#                 state2 = f"{node2}-{state_key2}"
#                 state2_idx = state_index_map[state2]
#                 # print(state1, state2)
#                 # Step 4: Compute the transition weight based on phasings
#                 # 1. Convert the state keys to phasings
#                 ff = str_2_phas_1(state_key1, 3)
#                 sf = str_2_phas_1(state_key2, 3)
                
#                 # 2. Find the common index (e.g., common_ff and common_sf)
#                 common_ff, common_sf = find_common_element_and_index(node1, node2)
#                 # print('common indices:', common_ff, common_sf)

#                 # 3. Get the matching phasings
#                 phasings = find_phasings_matches(ff, sf, common_ff, common_sf)
#                 all_phasings_str = []
#                 for phas in phasings:
#                     all_phasings = permute_rows(phas)
#                     phasings_str = [phas_2_str(phasss) for phasss in all_phasings]
#                     all_phasings_str += phasings_str
#                 # print('matched phasings', phasings_str)

#                 # 4. Calculate the transition weight by summing the values of matching keys
#                 transition_weight = 0
#                 for phasing_key in phasings_str:
#                     if phasing_key in qg[node1][node2]['weight']:
#                         transition_weight += qg[node1][node2]['weight'][phasing_key]
#                     # else:

#                 # print('edge weights:', qg[node1][node2]['weight'])
#                 # print('transition prob.', transition_weight)
#                 # 5. Set the transition weight in the matrix
#                 transition_matrix[state1_idx][state2_idx] = transition_weight
    
#     # transition_matrix = transition_matrix / transition_matrix.sum(axis=1, keepdims=True)
#     # Step 5: Normalize the transition matrix row-wise, skipping zero-sum rows
#     for i in range(num_states):
#         row_sum = transition_matrix[i].sum()
#         if row_sum > 0:
#             transition_matrix[i] /= row_sum
#         else:
#             # Optionally, set a uniform distribution for zero-sum rows
#             transition_matrix[i] = np.zeros(num_states)  # or set to uniform probabilities
#             # Uncomment the following line if you prefer uniform probabilities instead of zeros
#             # transition_matrix[i] = np.full(num_states, 1/num_states)


#     return state_names, transition_matrix

#     state_names, transition_matrix = generate_hmm_with_weights(qg)





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


    for nn in factor_graph.nodes(data=True):
        print(nn)
    
    for edge in factor_graph.edges(data=True):
        print(edge)



    beliefs = factor_graph_inference(factor_graph)
    # print("")

    # ve_map = var_elim(factor_graph)
    # print("Variable Elimination----------------------------------------------")
    # print("Max Prob Phasing: {} | Marginal: {}".format(ve_map.map_query(), ve_map.max_marginal()))

    # print("")

    # variables = [node for node in factor_graph.nodes() if not isinstance(node, DiscreteFactor)]
    #
    # # Step 3: Perform inference
    # for var in variables:
    #     print(f"Belief for {var}:")
    #     result = beliefs.query(variables=list(var))
    #     for state_index, prob in enumerate(result.values):
    #         print(f"Value {result.variables} = {state_index}: {prob}")
    #
    #
    # for var in variables:
    #     print(f"Belief for {var}:")
    #     print(result[var])
    #
    # print(result)
    # for var in variables:
    #     print(f"Belief for {var}:")
    #     print(f"Values: {result[var].values}")
    #     print(f"Variables: {result[var].variables}")
    # # beliefs = BeliefPropagation(factor_graph)
    #
    # result12 = beliefs.query(variables=[list(variables)[0]])
    # print(result12)
    # # marginals, max_phasings = give_marginals(factor_graph, qg, beliefs)
    #
    # #
    # # max_phase, positions = query_paths_gibbs_max(fragment_list, qg, beliefs, n_samples=1000)
    # # h_df = creat_vcf(max_phase, positions, config)
    # # print(h_df)
    #
    # # va_inference = VariableElimination(factor_graph)
    # # result = va_inference.query(variables=['1-2'])
    #
    # from hmmlearn import hmm
    #
    # # Define states and observations
    # states = ["Rainy", "Sunny"]
    # observations = ["Walk", "Shop", "Clean"]
    # n_states = len(states)
    # n_observations = len(observations)
    #
    # # Example emission probabilities
    # # The rows correspond to states, and the columns correspond to observations
    # emission_probabilities = np.array([
    #     [0.1, 0.4, 0.5],  # Emission probabilities for "Rainy"
    #     [0.6, 0.3, 0.1]   # Emission probabilities for "Sunny"
    # ])
    #
    # # Generate a random sequence of observations (for example purposes)
    # # Normally, you would use your actual observation sequence data
    # np.random.seed(42)
    # sequence_length = 100
    # X = np.random.choice(n_observations, sequence_length).reshape(-1, 1)
    #
    # # Initialize the HMM
    # model = hmm.MultinomialHMM(n_components=n_states, n_iter=100, random_state=42)
    #
    # # Set the emission probabilities
    # model.emissionprob_ = emission_probabilities
    #
    # # Fit the model to the observation sequence
    # model.fit(X)
    #
    # Extract transition probabilities
    # transition_probabilities = model.transmat_

    

