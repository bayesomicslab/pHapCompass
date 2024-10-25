import numpy as np
import itertools
from utils.utils import *
from scipy.stats import entropy


class QuadNode:
    def __init__(self, label, depth=0):
        self.label = label
        self.depth = depth
        self.children = {}
    
    def add_child(self, label):
        self.children[label] = QuadNode(label, self.depth + 1)
        return self.children[label]


def create_pruned_quadtree_ploidy(max_depth, genotype, allel_set=[0, 1]):
    allels = [str(a) for a in allel_set]
    root = QuadNode("root")
    allels_product = [''.join(ll) for ll in list(itertools.product(allels, allels))]
    # label_to_index = {"00": 0, "01": 1, "10": 2, "11": 3}
    label_to_index = {k: allels_product.index(k) for k in allels_product}
    first_positions = [int(i[0]) for i in label_to_index]
    second_positions = [int(i[1]) for i in label_to_index]
    max_allel = int(max(allel_set))
    seen_sequences_at_depth = {}
    
    stack = [(root, [0] * len(label_to_index))]
    
    while stack:
        # print(' ------------ go to stack -------------- ')
        node, path_counts = stack.pop()
        # print('Path count', path_counts)
        remaining_levels = max_depth - node.depth
        # print('path_counts', path_counts, 'remaining_levels', remaining_levels)
        if node.depth not in seen_sequences_at_depth:
            seen_sequences_at_depth[node.depth] = set()
        
        for label, index in label_to_index.items():
            # print('   ', label)
            if node.depth >= max_depth:
                continue
            
            new_counts = path_counts.copy()
            new_counts[index] += 1
            sum_first_position = np.sum(np.array(new_counts) * np.array(first_positions))
            sum_second_position = np.sum(np.array(new_counts) * np.array(second_positions))
            
            # print(new_counts)
            # print('remaining_levels', remaining_levels, 'sum_first_position', sum_first_position, 'sum_second_position',
            #       sum_second_position)
            # if remaining_levels - 1 == 0 and ((sum_first_position!= int(genotype[0])) or (sum_second_position != int(genotype[1]))):
            #     print('inja 000')
            #     continue
            
            if (sum_first_position + max_allel * (remaining_levels - 1) < int(genotype[0])) or \
                    (sum_second_position + max_allel * (remaining_levels - 1) < int(genotype[1])):
                # print('never reaches')
                continue
            
            if sum_first_position > int(genotype[0]) or sum_second_position > int(genotype[1]):
                # print('larger')
                continue
            
            # print(new_counts)
            # print('  ', 'sum_first_position:', sum_first_position , ', max remaining', max_allel*(remaining_levels - 1))
            # print('  ', 'sum_second_position:', sum_second_position , ', max remaining', max_allel*(remaining_levels - 1))
            
            not_equality_condition = (sum_first_position != int(genotype[0]) or sum_second_position != int(genotype[1]))
            if (remaining_levels == 1) and not_equality_condition:
                # print('not equal at last level ==============================================================')
                continue
            
            if (tuple(new_counts) not in seen_sequences_at_depth[node.depth]):
                # print('oooooomad')
                seen_sequences_at_depth[node.depth].add(tuple(new_counts))
                child = node.add_child(label)
                stack.append((child, new_counts))
                # print('   ', 'stack len:', len(stack))
    
    return root


def get_max_length_paths_ploidy(root, max_depth, allel_set=[0, 1]):
    # label_to_index = {"00": 0, "01": 1, "10": 2, "11": 3}
    allels = [str(a) for a in allel_set]
    allels_product = [''.join(ll) for ll in list(itertools.product(allels, allels))]
    label_to_index = {k: allels_product.index(k) for k in allels_product}
    
    stack = [(root, [])]
    
    max_length_paths = []
    
    while stack:
        current, path = stack.pop()
        
        if not current.children and len(path) == max_depth:
            counts = [0] * len(label_to_index)
            for label in path:
                counts[label_to_index[label]] += 1
            max_length_paths.append(counts)
        
        for label, child in current.children.items():
            stack.append((child, path + [label]))
    
    return max_length_paths




def counts_to_phasing_ploidy(max_length_paths, allel_set=[0, 1]):
    allels = [str(a) for a in allel_set]
    allels_product = [''.join(ll) for ll in list(itertools.product(allels, allels))]
    label_to_index = {k: allels_product.index(k) for k in allels_product}
    
    # label_to_index = {"00": 0, "01": 1, "10": 2, "11": 3}
    index_to_label = {v: k for k, v in label_to_index.items()}
    aa = [[ss[i] * [index_to_label[i]] for i in range(len(ss))] for ss in max_length_paths]
    phasings = [list(itertools.chain.from_iterable(a)) for a in aa]
    phasing_np = [np.vstack([np.array([int(ph) for ph in phase]) for phase in phasing]) for phasing in phasings]
    return phasing_np


def generate_phasings_ploidy(ploidy, genotype, allel_set=[0, 1]):
    tree = create_pruned_quadtree_ploidy(ploidy, genotype, allel_set=allel_set)
    max_length_paths = get_max_length_paths_ploidy(tree, ploidy, allel_set=allel_set)
    phasing_np_list = counts_to_phasing_ploidy(max_length_paths, allel_set=allel_set)
    return phasing_np_list




def phas_2_str(phas):
    return ''.join([str(ph) for ph in list(np.ravel(phas))])


def find_matchings(nodes_part1, nodes_part2):
    # Sort both parts and remember the original indices.
    sorted_part1 = sorted(enumerate(nodes_part1), key=lambda x: x[1])
    sorted_part2 = sorted(enumerate(nodes_part2), key=lambda x: x[1])
    
    # Split nodes by type and collect their original indices.
    def split_by_type(sorted_nodes):
        grouped = {}
        for idx, t in sorted_nodes:
            if t not in grouped:
                grouped[t] = []
            grouped[t].append(idx)
        return grouped
    
    grouped_part1 = split_by_type(sorted_part1)
    grouped_part2 = split_by_type(sorted_part2)
    if grouped_part1.keys() != grouped_part2.keys():
        return []
    if any([len((grouped_part1[i])) != len((grouped_part2[i])) for i in grouped_part1.keys()]):
        return []
    # Start with a single empty matching.
    matchings = [[]]
    for node_type, indices1 in grouped_part1.items():
        indices2 = grouped_part2[node_type]
        
        # For each current matching, extend it with all possible permutations for the current type.
        new_matchings = []
        for perm in itertools.permutations(indices2, len(indices2)):
            for current_matching in matchings:
                # Add new matching to the results only if it doesn't conflict with the current matching.
                if all((i1, i2) not in current_matching for i1, i2 in zip(indices1, perm)):
                    new_matchings.append(current_matching + list(zip(indices1, perm)))
        matchings = new_matchings
    
    return matchings

# @profile
def generate_phasings_ploidy_long(ploidy, genotype, allel_set=[0, 1]):
    tree = create_pruned_quadtree_ploidy_long(ploidy, genotype, allel_set=allel_set)
    max_length_paths = get_max_length_paths_ploidy_long(tree, ploidy, genotype, allel_set=allel_set)
    phasing_np_list = counts_to_phasing_ploidy_long(max_length_paths, genotype, allel_set=allel_set)
    return phasing_np_list

# @profile
def create_pruned_quadtree_ploidy_long(max_depth, genotype, allel_set=[0, 1]):
    allels = [str(a) for a in allel_set]
    root = QuadNode("root")
    allels_product = [''.join(ll) for ll in list(itertools.product(allels, repeat=len(genotype)))]
    # label_to_index = {"00": 0, "01": 1, "10": 2, "11": 3}
    label_to_index = {k: allels_product.index(k) for k in allels_product}
    # first_positions = [int(i[0]) for i in label_to_index]
    # second_positions = [int(i[1]) for i in label_to_index]
    numbered_positions = {n: [int(i[n]) for i in label_to_index] for n in range(len(genotype))}
    
    max_allel = int(max(allel_set))
    seen_sequences_at_depth = {}
    
    stack = [(root, [0] * len(label_to_index))]
    
    while stack:
        # print(' ------------ go to stack -------------- ')
        node, path_counts = stack.pop()
        # print('Path count', path_counts)
        remaining_levels = max_depth - node.depth
        # print('path_counts', path_counts, 'remaining_levels', remaining_levels)
        if node.depth not in seen_sequences_at_depth:
            seen_sequences_at_depth[node.depth] = set()
        
        for label, index in label_to_index.items():
            # print('   ', label)
            if node.depth >= max_depth:
                continue
            
            new_counts = path_counts.copy()
            new_counts[index] += 1
            
            # sum_first_position = np.sum(np.array(new_counts) * np.array(first_positions))
            # sum_second_position = np.sum(np.array(new_counts) * np.array(second_positions))
            numbered_sum = {n: np.sum(np.array(new_counts) * np.array(numbered_positions[n])) for n in
                            range(len(genotype))}
            
            # print(new_counts)
            # print('remaining_levels', remaining_levels, 'sum_first_position', sum_first_position,
            # 'sum_second_position',
            #       sum_second_position)
            # if remaining_levels - 1 == 0 and ((sum_first_position!= int(genotype[0])) or
            # (sum_second_position != int(genotype[1]))):
            #     print('inja 000')
            #     continue
            
            never_reaches_condition = [numbered_sum[n] + max_allel * (remaining_levels - 1) < int(genotype[n])
                                       for n in range(len(genotype))]
            # if (sum_first_position + max_allel * (remaining_levels - 1) < int(genotype[0])) or \
            #         (sum_second_position + max_allel * (remaining_levels - 1) < int(genotype[1])):
            if any(never_reaches_condition):
                # print('never reaches')
                continue
            
            larger_condition = [numbered_sum[n] > int(genotype[n]) for n in range(len(genotype))]
            # if sum_first_position > int(genotype[0]) or sum_second_position > int(genotype[1]):
            if any(larger_condition):
                # print('larger')
                continue
            
            # print(new_counts)
            # print('  ', 'sum_first_position:', sum_first_position , ', max remaining', max_allel*(remaining_levels - 1))
            # print('  ', 'sum_second_position:', sum_second_position , ', max remaining', max_allel*(remaining_levels - 1))
            
            # not_equality_condition = (sum_first_position != int(genotype[0]) or sum_second_position != int(genotype[1]))
            not_equality_condition = any([numbered_sum[n] != int(genotype[n]) for n in range(len(genotype))])
            if (remaining_levels == 1) and not_equality_condition:
                # print('not equal at last level ==============================================================')
                continue
            
            if (tuple(new_counts) not in seen_sequences_at_depth[node.depth]):
                # print('oooooomad')
                seen_sequences_at_depth[node.depth].add(tuple(new_counts))
                child = node.add_child(label)
                stack.append((child, new_counts))
                # print('   ', 'stack len:', len(stack))
    
    return root

# @profile
def get_max_length_paths_ploidy_long(root, max_depth, genotype, allel_set=[0, 1]):
    # label_to_index = {"00": 0, "01": 1, "10": 2, "11": 3}
    allels = [str(a) for a in allel_set]
    allels_product = [''.join(ll) for ll in list(itertools.product(allels, repeat=len(genotype)))]
    label_to_index = {k: allels_product.index(k) for k in allels_product}
    
    stack = [(root, [])]
    
    max_length_paths = []
    
    while stack:
        current, path = stack.pop()
        
        if not current.children and len(path) == max_depth:
            counts = [0] * len(label_to_index)
            for label in path:
                counts[label_to_index[label]] += 1
            max_length_paths.append(counts)
        
        for label, child in current.children.items():
            stack.append((child, path + [label]))
    
    return max_length_paths

# @profile
def counts_to_phasing_ploidy_long(max_length_paths, genotype, allel_set=[0, 1]):
    allels = [str(a) for a in allel_set]
    allels_product = [''.join(ll) for ll in list(itertools.product(allels, repeat=len(genotype)))]
    label_to_index = {k: allels_product.index(k) for k in allels_product}
    
    # label_to_index = {"00": 0, "01": 1, "10": 2, "11": 3}
    index_to_label = {v: k for k, v in label_to_index.items()}
    aa = [[ss[i] * [index_to_label[i]] for i in range(len(ss))] for ss in max_length_paths]
    phasings = [list(itertools.chain.from_iterable(a)) for a in aa]
    phasing_np = [np.vstack([np.array([int(ph) for ph in phase]) for phase in phasing]) for phasing in phasings]
    return phasing_np


def compute_edge_weight(new_vertex_name, v_label, source_phasings, target_phasings, fragment_model, config):

    possitions = sorted(set([int(nn) for nn in new_vertex_name.split('-')] + [int(nn) for nn in v_label.split('-')]))
    common_ff, common_sf = find_common_element_and_index(new_vertex_name, v_label)
    all_phasings =[]
    for ffstr in source_phasings:
        for sfstr in target_phasings:
            
            matched_phasings = find_phasings_matches(str_2_phas_1(ffstr, config.ploidy), str_2_phas_1(sfstr, config.ploidy), 
                                                     common_ff, common_sf, new_vertex_name, v_label)
            
            sorted_phasings = []
            for mtx in matched_phasings:
                sorted_matrix = mtx[np.argsort([''.join(map(str, row)) for row in mtx])]
                sorted_phasings.append(sorted_matrix)
            matched_phasings_str = [phas_2_str(pm) for pm in sorted_phasings]
            all_phasings += matched_phasings_str
    matches = get_matching_reads_for_positions(possitions, fragment_model.fragment_list)
    weights = {phas_2_str(phas): 0 for phas in all_phasings}
    for phas in all_phasings:
        for indc, this_po, obs in matches:
            weights[phas_2_str(phas)] += compute_likelihood_generalized_plus(np.array(obs), phas, indc,
                                                                                list(range(len(indc))),
                                                                                config.error_rate)
    entr = entropy(list(weights.values()), base=10)
    final_weight = {"weight": weights, "entropy": entr}
    return final_weight