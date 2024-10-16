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
from chordal_contraction import *
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


    qg = chordal_contraction_cycle_base(quotient_g, fragment_list, input_handler, config)

    # qg = chordal_contraction(quotient_g, fragment_list, input_handler, config)
    plot_graph(qg)    
    print('Chordal Graph constructed.')






    error_rate = 0.001
    # state_names, transition_matrix, emission_prob_matrix, emission_index_map = generate_hmm_with_weights_and_emissions2(qg, error_rate)


    # for nn in qg.nodes(data=True):
    #     print(nn)
    # Save the graph in graph-tool binary format (preserves all properties)
#     quotient_g.graph.save("/home/mok23003/BML/HaplOrbit/example/Contig1_k3/graphs/c2/ART_0_quotient.gt")

#     # Save the graph in GraphML format (to preserve properties for external tools)
#     quotient_g.graph.save("/home/mok23003/BML/HaplOrbit/example/Contig1_k3/graphs/c2/ART_0_quotient.graphml", fmt="xml")
    
#     v_label_loaded = g_loaded.vertex_properties["v_label"]  # Corrected name for vertex label
#     e_weights_loaded = g_loaded.edge_properties["e_weights"]  # Corrected name for edge weights
#     e_label_loaded = g_loaded.edge_properties["e_label"]  # Corrected name for edge label

#     # Print the vertex labels
#     print("Vertex Labels:")
#     for v in g_loaded.vertices():
#         print(f"Vertex {int(v)}: Label = {v_label_loaded[v]}")

#     # Print the edge labels and weights
#     print("\nEdge Labels and Weights:")
#     for e in g_loaded.edges():
#         print(f"Edge {int(e.source())} -> {int(e.target())}: Label = {e_label_loaded[e]}, Weight = {e_weights_loaded[e]}")


def is_chordless(g, cycle_vertices):
    n = len(cycle_vertices)
    
    # A cycle is chordless if there are no edges between non-adjacent vertices in the cycle
    for i in range(n):
        for j in range(i + 2, n):
            # print(i, j)
            if j != (i + 1) % n:  # Ensure we're not comparing consecutive vertices in the cycle
                # Check if there is an edge between non-adjacent vertices in the cycle
                if g.edge(cycle_vertices[i], cycle_vertices[j]) is not None:
                    return False  # Found a chord, so it's not chordless
    return True


main_path = '/home/mok23003/BML/HaplOrbit/old_simulated_data_graphs'

contig = 'Contig1_k3'
coverage = 'c2'
frag_file = 'ART_0.frag.txt'


quotient_graph, v_label_reversed, e_label_reversed = read_quotient_graph(main_path, contig, coverage, frag_file)

new_graph = quotient_graph.copy()
# simple_cycles = gt.all_circuits(new_graph)
e_weights = new_graph.edge_properties["e_weights"]
mst_graph, non_mst_graph, tree = get_minimum_spanning_tree(new_graph)
# cycle_basis_edges = []
for e in non_mst_graph.edges():

    ep = new_graph.new_edge_property("bool")
    ep.a = True
    ep[e] = False
    new_graph.set_edge_filter(ep)
    half_cyc = gt.shortest_path(new_graph, e.source(), e.target())
    cycle_e = half_cyc[1] + [e]
    cycle_v = half_cyc[0]



    cycle_e, cycle_v = find_cycle(mst_graph, new_graph, e)
    # cycle_basis_edges.append(cycle)
    print(len(cycle_e)) 

    if is_chordless(new_graph, cycle_v) or len(cycle_e) < 4:
        print('Chordless')




edge_to_exclude = g.edge(1, 3)

# Create an edge property map to filter edges
ep = g.new_edge_property("bool")  # A boolean property map for edges

# Set all edges to True (keep them)
ep.a = True

# Set the specific edge to False (exclude it)
ep[edge_to_exclude] = False

# Apply the edge filter to the graph
g.set_edge_filter(ep)

# Now the edge (1, 3) is excluded (but not deleted)

# Example: Check the edges that remain in the graph
g.set_edge_filter(None)












adjacency = gt.adjacency(new_graph)
dim = adjacency.shape[0]

def chordless_cycles(adjacency, dim):
    # Iterate over pairs of vertices (i, j)
    for i in range(dim - 2):
        for j in range(i + 1, dim - 1):
            # Check if there is an edge between i and j
            if adjacency[i][j] == 0:
                continue

            candidates = []  # Store potential cycles

            # Explore vertices k > j to form cycles with (i, j)
            for k in range(j + 1, dim):
                if adjacency[i][k] == 0:
                    continue

                if adjacency[j][k] == 1:
                    # Print the 3-cycle (i, j, k)
                    print(f"{i + 1} {j + 1} {k + 1}")
                    continue

                # Store candidates for further exploration
                candidates.append([j, i, k])

            # Explore candidates for larger chordless cycles
            while candidates:
                v = candidates.pop(0)  # Get the current candidate cycle
                k = v[-1]  # Last vertex in the candidate cycle

                # Explore new vertices m to extend the cycle
                for m in range(i + 1, dim):
                    if m in v:
                        continue  # Avoid duplicate vertices in the cycle
                    if adjacency[m][k] == 0:
                        continue  # Ensure m and k are connected

                    # Check for a chord (an edge between any non-adjacent vertices)
                    chord = False
                    for n in range(1, len(v) - 1):
                        if adjacency[m][v[n]] == 1:
                            chord = True
                            break

                    if chord:
                        continue  # Skip if a chord is found

                    if adjacency[m][j] == 1:
                        # Print the complete cycle
                        print(" ".join(str(x + 1) for x in v) + f" {m + 1}")
                        continue

                    # Extend the current candidate cycle with vertex m
                    new_cycle = v + [m]
                    candidates.append(new_cycle)






























for cyc in simple_cycles:
    print(len(cyc))    
    if is_chordless(new_graph, cyc):
        stop
        sorted_edges = sorted(cyc, key=lambda e: e_weights[e]['entropy'])
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
                




# quotient_graph, v_label_reversed, e_label_reversed = quotient_g, quotient_g.v_label_reversed, edges_map_quotient
# e_label = quotient_graph.edge_properties["e_label"]  # Edge property,
# e_weights = quotient_graph.edge_properties["e_weights"]  # Edge property, 

# e_weights[e_label_reversed['5-6-7']]


# quotient_graph.vertex_properties["v_label"]  # Vertex property,


# for v in quotient_graph.vertices():
#     print(quotient_graph.vertex_properties["v_label"][v])




# quotient_graph.edge_properties["e_label"][spath_edges[0]]
# quotient_graph.edge_properties["e_label"][spath_edges[1]]

# len_seq = len(spath_edges) + 1

# for v in spath_vertices:
#     print(quotient_graph.vertex_properties["v_label"][v])

# for ee in spath_edges:
#     print(quotient_graph.edge_properties["e_label"][ee])

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
                    
























e_label = quotient_graph.edge_properties["e_label"]
v_label = quotient_graph.vertex_properties["v_label"]


vertex_filter = quotient_graph.new_vertex_property("bool")
for v in quotient_graph.vertices():
    vertex_filter[v] = False

# Step 2: Selectively enable only certain vertices

# For example, enable only v1, v2, and v3
# vertex_filter[source_vertex] = True
# vertex_filter[target_vertex] = True
# vertex_filter[neighbor] = True
# for cc in ccyc:
#     vertex_filter[cc.source()] = True
#     vertex_filter[cc.target()] = True

for vv in mst_tree_vertices:
    vertex_filter[vv] = True


quotient_graph.set_vertex_filter(vertex_filter)

# Step 3: Draw the filtered graph
gt.graph_draw(quotient_graph, vertex_text=v_label, vertex_size=10, edge_font_size=10, edge_text=e_label, output_size=(500, 500))
quotient_graph.clear_filters()
