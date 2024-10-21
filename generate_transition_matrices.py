from utils.utils import *
import os
import numpy as np
import gzip
import _pickle as cPickle
from multiprocessing import Pool


def generate_transition_matrix_make_input():
    inputs = []
    contigs = ['Contig1_k3', 'Contig2_k3', 'Contig3_k3']
    # contigs = ['Contig2_k3']
    # simulated_data_path = '/home/mok23003/BML/HaplOrbit/simulated_data'
    main_path = '/home/mok23003/BML/HaplOrbit/simulated_data_graphs/'
    # save_path = '/home/mok23003/BML/HaplOrbit/simulated_data_graphs/transition_matrices'
    # '/home/mok23003/BML/HaplOrbit/simulated_data_graphs/quotient_graphs/Contig1_k3/c2'
    for cont in contigs:
        cont_path = os.path.join(main_path, 'quotient_graphs', cont)
        coverges = sorted([d for d in os.listdir(cont_path) if os.path.isdir(os.path.join(cont_path, d))])
        for coverage in coverges:
            
            this_cont_coverage_path = os.path.join(cont_path, coverage)
            this_save_path = os.path.join(main_path, 'transition_matrices', cont, coverage)

            if not os.path.exists(this_save_path):
                os.makedirs(this_save_path)
            existing_files = [ff for ff in os.listdir(this_save_path) if '.pkl.gz' in ff]

            graph_files = [f for f in os.listdir(this_cont_coverage_path) if '.gt.gz' in f]

            frag_files = [f.split('.gt.gz')[0] + '.frag.txt' for f in graph_files if not f.split('.gt.gz')[0] + '.pkl.gz' in existing_files]
            for frag_file in frag_files:
                this_input = [main_path, cont, coverage, frag_file, this_save_path]
                inputs.append(this_input)
    return inputs


def generate_transition_matrices(inp):
    main_path, contig, coverage, frag_file, this_save_path = inp
    transition_path = os.path.join(this_save_path, frag_file.split('.frag.txt')[0] + '.pkl.gz')
    print('Working on', transition_path)
    quotient_graph, v_label_reversed, e_label_reversed = read_quotient_graph(main_path, contig, coverage, frag_file)
    edge_weights_dict = {}
    edges = list(quotient_graph.edges())

    for edge in edges:
        source = edge.source()
        target = edge.target()
        source_weights = quotient_graph.vertex_properties["v_weights"][source]['weight']
        target_weights = quotient_graph.vertex_properties["v_weights"][target]['weight']
        source_label = quotient_graph.vertex_properties["v_label"][source]
        target_label = quotient_graph.vertex_properties["v_label"][target]
        # e_label = quotient_graph.edge_properties["e_label"][edge]
        e_label = '--'.join(sorted([source_label, target_label]))
        e_weights = quotient_graph.edge_properties["e_weights"][edge]['weight']
        common_ff, common_sf = find_common_element_and_index(source_label, target_label)
        source_phasings = list(source_weights.keys())
        target_phasings = list(target_weights.keys())
        # transitions_dict = {'source': source_phasings, 'target': target_phasings}
        transitions_mtx = np.zeros((len(source_phasings), len(target_phasings)))
        for i, ffstr in enumerate(source_phasings):
            for j, sfstr in enumerate(target_phasings):
                
                matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, 3), str_2_phas_1(sfstr, 3), common_ff, common_sf, source_label, target_label)
                sorted_phasings = []
                for mtx in matched_phasings:
                    sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                    sorted_phasings.append(sorted_matrix)
                matched_phasings_str = [phas_2_str(pm) for pm in sorted_phasings] 
                this_weight = np.sum([e_weights[pm] for pm in matched_phasings_str if pm in e_weights.keys()])
                transitions_mtx[i, j] = this_weight
        transitions_mtx = transitions_mtx / transitions_mtx.sum(axis=1, keepdims=True)
        edge_weights_dict[e_label] = {'transitions': transitions_mtx, 'source_label': source_label, 'target_label': target_label}

    # with open(transition_path, "wb") as f:
    #     pickle.dump(edge_weights_dict, f)
    with gzip.open(transition_path, 'wb') as f:
        cPickle.dump(edge_weights_dict, f)
    print('[Done]', transition_path)


def generate_transitions():
    inputs = generate_transition_matrix_make_input()
    print('Number of inputs:', len(inputs))
    # for inp in inputs:
    #     generate_quotient_graph(inp)
    pool = Pool(10)
    pool.map(generate_transition_matrices, inputs)



# path = '/home/mok23003/BML/HaplOrbit/simulated_data_graphs/transition_matrices/Contig1_k3/c2/ART_0.pkl.gz'
# with gzip.open(path, 'rb') as f:
#     transitions = cPickle.load(f)



if __name__ == '__main__':
    generate_transitions()
