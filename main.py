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
from algorithm.chordal_contraction import *
# from test.FFBS import generate_hmm_with_weights_and_emissions




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
            self.data_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/c2/ART_0.frag.txt'
            self.bam_path = 'example/example.bam'
            self.genotype_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/real_haps_contig1_k3.txt'
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

    # frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)
    # frag_graph, fragment_list = fragment_model.construct2(input_handler, config)
    fragment_model.construct2(input_handler, config)
    
    # fragment_model = FragmentGraph(args.data_path, args.genotype_path, args.ploidy, input_handler.alleles)
    # frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)

    # plot_graph(frag_graph)
    print('Fragment Graph constructed.')

    # quotient_g = QuotientGraph(fragment_model.graph).construct(fragment_model.fragment_list, input_handler, config)
    quotient_g = QuotientGraph(fragment_model)
    quotient_g.construct3(input_handler, config)
    # plot_graph(quotient_g)
    print('Quotient Graph constructed.')


    qg = chordal_contraction_graph_tool(quotient_graph, input_handler, config, fragment_model)
    print('Chordal Graph constructed.')





class Args:
    def __init__(self):
        self.vcf_path = 'example/62_ID0.vcf'
        self.data_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/c2/ART_0.frag.txt'
        self.bam_path = 'example/example.bam'
        self.genotype_path = '/home/mok23003/BML/HaplOrbit/simulated_data/Contig1_k3/real_haps_contig1_k3.txt'
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

# frag_graph, fragment_list = fragment_model.construct_graph(input_handler, config)
# frag_graph, fragment_list = fragment_model.construct2(input_handler, config)
fragment_model.construct2(input_handler, config)


main_path = '/home/mok23003/BML/HaplOrbit/old_simulated_data_graphs'

contig = 'Contig1_k3'
coverage = 'c2'
frag_file = 'ART_0.frag.txt'


quotient_graph, v_label_reversed, e_label_reversed = read_quotient_graph(main_path, contig, coverage, frag_file)


def chordal_contraction_graph_tool2(quotient_graph, input_handler, config, fragment_model):
    new_graph = quotient_graph.copy()
    e_weights = new_graph.edge_properties["e_weights"]
    new_graph.clear_filters()
    e_entropy = new_graph.new_edge_property("double")

    # Loop over edges and assign entropy from the e_weights property
    for e in new_graph.edges():
        e_entropy[e] = e_weights[e]['entropy']

    new_graph.ep['e_entropy'] = e_entropy

    chordless_cycles = get_chordless_cycles(new_graph)

    to_be_removed_nodes = []

    for cyc_id, cyc in enumerate(chordless_cycles):
        
        print(cyc_id)
        
        edges = [new_graph.edge(cyc[-1], cyc[0])]
        for i in range(len(cyc) - 1):
            edges += [new_graph.edge(cyc[i], cyc[i+1])]
        edges = [x for x in edges if x is not None]
        while len(edges) > 3:
            min_edge = min(edges, key=lambda e: new_graph.ep['e_entropy'][e])
            source_label = new_graph.vp['v_label'][min_edge.source()]
            target_label = new_graph.vp['v_label'][min_edge.target()]
            # new node positions
            poss = sorted(set([int(nn) for nn in source_label.split('-')] + [int(nn) for nn in target_label.split('-')]))
            # new vertex properties:
            new_vertex_name = '-'.join([str(nnn) for nnn in poss])
            vertex_weights = new_graph.ep['e_weights'][min_edge]
            
            new_graph.vertex_properties["v_weights"][min_edge.source()] = vertex_weights
            new_graph.vertex_properties["v_label"][min_edge.source()] = new_vertex_name

            source_nbrs = [n for n in min_edge.source().all_neighbors() if n != min_edge.target()]
            target_nbrs = [n for n in min_edge.target().all_neighbors() if n != min_edge.source()]
            common_nbrs = set(source_nbrs).intersection(set(target_nbrs))

            for n in common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]
                # e_poss = sorted(set([int(nn) for nn in v_label.split('-')] + poss))
                # print(len(e_poss))
                # new_edge_name = '-'.join([str(nnn) for nnn in e_poss])
                sorted_labels = sort_strings([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)
                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)

                e1 = new_graph.edge(min_edge.source(), n)
                e2 = new_graph.edge(min_edge.target(), n)
                
                new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
                new_graph.remove_edge(e2)

            for n in set(source_nbrs)-common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]

                sorted_labels = sort_strings([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)
                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)


                e1 = new_graph.edge(min_edge.source(), n)
                # e2 = new_graph.edge(min_edge.target(), n)
                new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
                # new_graph.edge_properties["e_weights"][e2]

            
            for n in set(target_nbrs)-common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]
                sorted_labels = sort_strings([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)

                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)

                e2 = new_graph.edge(min_edge.target(), n)
                new_graph.remove_edge(e2)
                e1 = new_graph.add_edge(min_edge.source(), n)
                new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
            
            # to_be_removed_nodes += [min_edge.target()]
            to_be_removed_nodes.append(min_edge.target())
            new_graph.remove_edge(min_edge)
            edges.remove(min_edge)

    new_graph.remove_vertex(to_be_removed_nodes)
    return new_graph
    


def chordal_contraction_graph_tool_approx(quotient_graph, input_handler, config, fragment_model):
    new_graph = quotient_graph.copy()
    e_weights = new_graph.edge_properties["e_weights"]
    new_graph.clear_filters()
    e_entropy = new_graph.new_edge_property("double")

    # Loop over edges and assign entropy from the e_weights property
    for e in new_graph.edges():
        e_entropy[e] = e_weights[e]['entropy']

    new_graph.ep['e_entropy'] = e_entropy

    chordless_cycles = get_chordless_cycles(new_graph)

    to_be_removed_nodes = []

    for cyc_id, cyc in enumerate(chordless_cycles):
        
        print(cyc_id)
        
        edges = [new_graph.edge(cyc[-1], cyc[0])]
        for i in range(len(cyc) - 1):
            edges += [new_graph.edge(cyc[i], cyc[i+1])]
        edges = [x for x in edges if x is not None]
        while len(edges) > 3:
            min_edge = min(edges, key=lambda e: new_graph.ep['e_entropy'][e])
            source_label = new_graph.vp['v_label'][min_edge.source()]
            target_label = new_graph.vp['v_label'][min_edge.target()]
            # new node positions
            poss = sorted(set([int(nn) for nn in source_label.split('-')] + [int(nn) for nn in target_label.split('-')]))
            # new vertex properties:
            new_vertex_name = '-'.join([str(nnn) for nnn in poss])
            vertex_weights = new_graph.ep['e_weights'][min_edge]
            
            new_graph.vertex_properties["v_weights"][min_edge.source()] = vertex_weights
            new_graph.vertex_properties["v_label"][min_edge.source()] = new_vertex_name

            source_nbrs = [n for n in min_edge.source().all_neighbors() if n != min_edge.target()]
            target_nbrs = [n for n in min_edge.target().all_neighbors() if n != min_edge.source()]
            common_nbrs = set(source_nbrs).intersection(set(target_nbrs))

            for n in common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]
                # e_poss = sorted(set([int(nn) for nn in v_label.split('-')] + poss))
                # print(len(e_poss))
                # new_edge_name = '-'.join([str(nnn) for nnn in e_poss])
                sorted_labels = sort_strings([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)
                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)

                e1 = new_graph.edge(min_edge.source(), n)
                e2 = new_graph.edge(min_edge.target(), n)
                
                new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
                new_graph.remove_edge(e2)

            for n in set(source_nbrs)-common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]

                sorted_labels = sort_strings([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)
                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)


                e1 = new_graph.edge(min_edge.source(), n)
                # e2 = new_graph.edge(min_edge.target(), n)
                new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
                # new_graph.edge_properties["e_weights"][e2]

            
            for n in set(target_nbrs)-common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]
                sorted_labels = sort_strings([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)

                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)

                e2 = new_graph.edge(min_edge.target(), n)
                new_graph.remove_edge(e2)
                e1 = new_graph.add_edge(min_edge.source(), n)
                new_graph.edge_properties["e_weights"][e1] = final_weight
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weight['entropy']
            
            # to_be_removed_nodes += [min_edge.target()]
            to_be_removed_nodes.append(min_edge.target())
            new_graph.remove_edge(min_edge)
            edges.remove(min_edge)

    new_graph.remove_vertex(to_be_removed_nodes)
    return new_graph
    







for vvv in quotient_graph.vertex_properties["v_label"]:
    print(vvv)
vertex1 = '125-126'
vertex2 = '324-342'


spath_edges = gt.shortest_path(quotient_graph, source=quotient_graph.vertex(v_label_reversed[vertex1]), target=quotient_graph.vertex(v_label_reversed[vertex2]))[1]
spath_vertices = gt.shortest_path(quotient_graph, source=quotient_graph.vertex(v_label_reversed[vertex1]), target=quotient_graph.vertex(v_label_reversed[vertex2]))[0]


def forward_sum(quotient_graph, spath_vertices, spath_edges):
    alpha = {vi: np.zeros(len(quotient_graph.vertex_properties["v_weights"][v]['weight'].keys())) for vi, v in enumerate(spath_vertices)}
    alpha[0] = np.array(list(quotient_graph.vertex_properties["v_weights"][spath_vertices[0]]['weight'].values()))
    alpha[0] = torch.tensor(alpha[0] / np.sum(alpha[0]))

    seq_len = len(spath_vertices)
    for sl in range(1, seq_len):
        source = spath_vertices[sl-1]
        target = spath_vertices[sl]
        source_weights = quotient_graph.vertex_properties["v_weights"][source]['weight']
        target_weights = quotient_graph.vertex_properties["v_weights"][target]['weight']
        source_label = quotient_graph.vertex_properties["v_label"][source]
        target_label = quotient_graph.vertex_properties["v_label"][target]
        e_label = quotient_graph.edge_properties["e_label"][spath_edges[sl-1]]
        e_weights = quotient_graph.edge_properties["e_weights"][spath_edges[sl-1]]['weight']

        common_ff, common_sf = find_common_element_and_index(source_label, target_label)
        source_phasings = list(source_weights.keys())
        target_phasings = list(target_weights.keys())
        transitions_dict = {'source': source_phasings, 'target': target_phasings}
        transitions_mtx = np.zeros((len(source_phasings), len(target_phasings)))
        for i, ffstr in enumerate(source_phasings):
            for j, sfstr in enumerate(target_phasings):
                matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, 3), str_2_phas_1(sfstr, 3), common_ff, common_sf)
                sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)

                matched_phasings_str = [phas_2_str(pm) for pm in sorted_phasings] 
                this_weight = np.sum([e_weights[pm] for pm in matched_phasings_str if pm in e_weights.keys()])
                transitions_mtx[i, j] = this_weight

        transitions_mtx = transitions_mtx / transitions_mtx.sum(axis=1, keepdims=True)

        dp_eq = alpha[sl-1][:, np.newaxis] * transitions_mtx
        alpha[sl] = torch.tensor(dp_eq.sum(axis=0))
    return alpha


def backward_sum(quotient_graph, spath_vertices, spath_edges):
    beta = {vi: np.zeros(len(quotient_graph.vertex_properties["v_weights"][v]['weight'].keys())) for vi, v in enumerate(spath_vertices)}
    # last_vertex = spath_vertices[-1]
    beta[seq_len - 1] = torch.ones(len(quotient_graph.vertex_properties["v_weights"][last_vertex]['weight'].keys()))

    for sl in range(0, seq_len-1)[::-1]:
        source = spath_vertices[sl]
        target = spath_vertices[sl + 1]
        source_weights = quotient_graph.vertex_properties["v_weights"][source]['weight']
        target_weights = quotient_graph.vertex_properties["v_weights"][target]['weight']
        source_label = quotient_graph.vertex_properties["v_label"][source]
        target_label = quotient_graph.vertex_properties["v_label"][target]
        e_label = quotient_graph.edge_properties["e_label"][spath_edges[sl]]
        e_weights = quotient_graph.edge_properties["e_weights"][spath_edges[sl]]['weight']
        for i, ffstr in enumerate(source_phasings):
            for j, sfstr in enumerate(target_phasings):
                matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, 3), str_2_phas_1(sfstr, 3), common_ff, common_sf)
                sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)

                matched_phasings_str = [phas_2_str(pm) for pm in sorted_phasings] 
                this_weight = np.sum([e_weights[pm] for pm in matched_phasings_str if pm in e_weights.keys()])
                transitions_mtx[i, j] = this_weight

        transitions_mtx = transitions_mtx / transitions_mtx.sum(axis=1, keepdims=True)

        dp_eq = torch.tensor(transitions_mtx) * beta[sl + 1][np.newaxis, :]
        beta[sl] = dp_eq.sum(axis=1)


def ffbs(alpha, seq_len, spath_vertices, spath_edges, quotient_graph):

    
    # Initialize a list to store the sampled states
    sampled_states = [None] * seq_len

    # Step 1: Sample the last state x_T
    last_vertex = spath_vertices[seq_len - 1]
    target_weights = quotient_graph.vertex_properties["v_weights"][last_vertex]['weight']
    target_phasings = list(target_weights.keys())

    # Convert alpha_T to probabilities
    probs = alpha[seq_len - 1].numpy()
    probs = probs / probs.sum()

    # Sample the last state based on alpha_T
    x_T_index = np.random.choice(len(target_phasings), p=probs)
    x_T = target_phasings[x_T_index]
    sampled_states[seq_len - 1] = x_T

    # Step 2: Backward Sampling for t from T-1 down to 0
    for t in range(seq_len - 2, -1, -1):
        source = spath_vertices[t]
        target = spath_vertices[t + 1]
        source_weights = quotient_graph.vertex_properties["v_weights"][source]['weight']
        target_weights = quotient_graph.vertex_properties["v_weights"][target]['weight']
        source_label = quotient_graph.vertex_properties["v_label"][source]
        target_label = quotient_graph.vertex_properties["v_label"][target]
        e_label = quotient_graph.edge_properties["e_label"][spath_edges[t]]
        e_weights = quotient_graph.edge_properties["e_weights"][spath_edges[t]]['weight']
        
        # Find common elements and indices
        common_ff, common_sf = find_common_element_and_index(source_label, target_label)
        source_phasings = list(source_weights.keys())
        target_phasings = list(target_weights.keys())
        transitions_mtx = np.zeros((len(source_phasings), len(target_phasings)))
        
        # Recompute the transition matrix at time t
        for i, ffstr in enumerate(source_phasings):
            for j, sfstr in enumerate(target_phasings):
                matched_phasings = find_phasings_matches(
                    str_2_phas_1(ffstr, 3),
                    str_2_phas_1(sfstr, 3),
                    common_ff, common_sf
                )
                sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)
                
                matched_phasings_str = [phas_2_str(pm) for pm in sorted_phasings]
                this_weight = np.sum(
                    [e_weights[pm] for pm in matched_phasings_str if pm in e_weights.keys()]
                )
                transitions_mtx[i, j] = this_weight
        
        # Normalize the transition matrix
        transitions_mtx = transitions_mtx / transitions_mtx.sum(axis=1, keepdims=True)
        
        # Get the index of the next sampled state
        x_t1 = sampled_states[t + 1]
        x_t1_index = target_phasings.index(x_t1)
        
        # Compute the probabilities for the current state
        alpha_t = alpha[t].numpy()
        probs = alpha_t * transitions_mtx[:, x_t1_index]
        probs = probs / probs.sum()
        
        # Sample the current state
        x_t_index = np.random.choice(len(source_phasings), p=probs)
        x_t = source_phasings[x_t_index]
        sampled_states[t] = x_t
    return sampled_states


alpha = forward_sum(quotient_graph, spath_vertices, spath_edges)
seq_len = len(spath_vertices)

sampled_states = ffbs(alpha, seq_len, spath_vertices, spath_edges, quotient_graph)
revsered_sampled_states = sampled_states[::-1]





cliques = gt.max_cliques(quotient_graph)
lenss = []
for cli in cliques:
    lenss.append(len(cli))

def chordal_graph_tools(quotient_graph):

    new_graph = quotient_graph.copy()
    # self.v_label = self.graph.new_vertex_property("string")
    # self.v_weights = self.graph.new_vertex_property("object")
    # self.e_label = self.graph.new_edge_property("string")
    # self.e_weights = self.graph.new_edge_property("object")
    # self.v_label_reversed = {}
    # self.e_label_reversed = {}

    mst_graph, non_mst_graph, tree = get_minimum_spanning_tree(new_graph)
    
    cycle_basis_edges, index_to_cycles, edge_to_cycles, forbidden_edges = get_cycles_basis_info_graph_tool(non_mst_graph, mst_graph, new_graph)
    
    # e_weights = 
    # quotient_graph.edge_properties["e_weights"]
    
    for c_id, ccyc in index_to_cycles.items():
        for ccc in ccyc:
            print(new_graph.edge_properties["e_label"][ccc])

        if len(ccyc) > 3:
            sorted_edges = sorted(ccyc, key=lambda e: e_weights[e]['entropy'])
            contracting_edges = sorted_edges[:-3]
            for ce in contracting_edges:
                # print(e_weights[ce]['entropy'])
                source_vertex = ce.source()
                target_vertex = ce.target()

                # Get all edges incident to the source and target vertices
                incident_edges = set(source_vertex.all_edges()).union(set(target_vertex.all_edges()))

                # Optionally, remove the original edge from the list
                incident_edges.discard(ce)

                for eee in incident_edges:
                    print(new_graph.edge_properties["e_label"][eee])
                    print(new_graph.vertex_properties["v_label"][target_vertex])



                new_node_name = new_graph.edge_properties["e_label"][ce]

                # new_node_name = '-'.join([str(nnn) for nnn in sorted(list(set([int(nn) for nn in \
                # new_graph.vertex_properties["v_label"][source_vertex].split('-')] + \
                # [int(nn) for nn in new_graph.vertex_properties["v_label"][target_vertex].split('-')])))])

                ce_wei = new_graph.edge_properties["e_weights"][ce]
                
                v1 = new_graph.add_vertex()
                new_graph.vertex_properties["v_label"][v1] = new_node_name
                new_graph.vertex_properties["e_weights"][v1] = ce_wei
                
                for ie in list(incident_edges):
                    
                    neighbor = list(set([ie.source(), ie.target()]) - set([source_vertex, target_vertex]))[0]
                    neighbor_label = new_graph.vertex_properties["v_label"][neighbor]
                    poss = sorted(set([int(nn) for nn in new_node_name.split('-')] + [int(nn) for nn in neighbor_label.split('-')]))
                    

