import numpy as np
import itertools
from algorithm.haplotype_assembly_helper import generate_phasings_ploidy_long, compute_likelihood_generalized_plus
from utils.utils import get_matching_reads_for_positions, phas_2_str, networkit_find_cliques, nx2nk
from scipy.stats import entropy
import random
import networkx as nx
import time


def chordal_contraction(quotient_g, fragment_list, inpt_handler, config):
    # plot_graph(quotient_g)
    qg = quotient_g.copy()
    start0 = time.time()
    while not nx.is_chordal(qg):
        s1 = time.time()
        # print('not chordal ---> time for checking:', s1 - start0)
        cliques_larger_than_2 = [cli for cli in nx.find_cliques(qg) if len(cli) > 2]
        s3 = time.time()
        # print('            ---> time for finding cliques:', s3 - s1)
        non_candicate_edges = []
        for cli in cliques_larger_than_2:
            non_candicate_edges += sorted(list(itertools.combinations(sorted(cli), 2)))
        non_candicate_edges = list(set(non_candicate_edges))
        # nx.is_chordal(quotient_g)
        # cycles = nx.chordless_cycles(quotient_g)
        candidate_edges = [edg for edg in list(qg.edges()) if edg not in non_candicate_edges]
        entropies = np.array([qg[edg[0]][edg[1]]['entropy'] for edg in candidate_edges])
        entropies[np.isnan(entropies)] = 0
        rev_entropies = [1 - e if 1 > e else 0 for e in entropies]
        if np.sum(rev_entropies) == 0:
            for ent_id in range(len(rev_entropies)):
                rev_entropies[ent_id] = 1 / len(rev_entropies)
        picked_edge = random.choices(candidate_edges, weights=rev_entropies, k=1)
        picked_edge = [picked_edge[0][0], picked_edge[0][1]]
        s4 = time.time()
        # print('            ---> time for pick a candidate edge:', s4 - s3)
        if picked_edge in [list(edg) for edg in list(qg.edges())]:
            qg = contract_one_edge(qg, picked_edge, inpt_handler, config, fragment_list)
        else:
            print(f"edge {picked_edge} not in graph")
        s5 = time.time()
        # print('            ---> time for edge contraction:', s5 - s4)
    # plot_graph(qg)
    
    return qg

# @profile
def chordal_contraction_networkit(quotient_g, fragment_list, inpt_handler, config):
    # quotient_g is a networkx graph
    qg_nx = quotient_g.copy()
    # start0 = time.time()
    while not nx.is_chordal(qg_nx):
        # s1 = time.time()
        # print('not chordal ---> time for checking:', s1 - start0)
        qg, reverse_map = nx2nk(qg_nx)
        # s2 = time.time()
        # print('            ---> time for converting nx2nk:', s2 - s1)
        # qg is neworkit
        cliques_larger_than_2 = [cli for cli in networkit_find_cliques(qg) if len(cli) > 2]
        # s3 = time.time()
        # print('            ---> time for finding cliques:', s3 - s2)
        non_candicate_edges = []
        for cli in cliques_larger_than_2:
            non_candicate_edges += sorted(list(itertools.combinations(sorted(cli), 2)))
        non_candicate_edges = list(set(non_candicate_edges))
        # nx.is_chordal(quotient_g)
        # cycles = nx.chordless_cycles(quotient_g)
        candidate_edges = [edg for edg in list(qg.iterEdges()) if edg not in non_candicate_edges]
        
        # entropies = np.array([qg[edg[0]][edg[1]]['entropy'] for edg in candidate_edges])
        entropies = np.array([qg_nx[reverse_map[edg[0]]][reverse_map[edg[1]]]['entropy'] for edg in candidate_edges])
        entropies[np.isnan(entropies)] = 0
        rev_entropies = [1 - e if 1 > e else 0 for e in entropies]
        if np.sum(rev_entropies) == 0:
            for ent_id in range(len(rev_entropies)):
                rev_entropies[ent_id] = 1 / len(rev_entropies)
        picked_edge = random.choices(candidate_edges, weights=rev_entropies, k=1)
        picked_edge = [reverse_map[picked_edge[0][0]], reverse_map[picked_edge[0][1]]]
        # s4 = time.time()
        # print('            ---> time for pick a candidate edge:', s4 - s3)
        qg_nx = contract_one_edge(qg_nx, picked_edge, inpt_handler, config, fragment_list)
        # s5 = time.time()
        # print('            ---> time for edge contraction:', s5 - s4)
    return qg_nx

# @profile
def contract_one_edge(quotient_g, picked_edge, inpt_handler, config, fragment_list):
    # s5 = time.time()
    new_graph = quotient_g.copy()
    node_neighbors = [nbr for nbr in list(quotient_g.adj[picked_edge[0]]) if nbr not in picked_edge] + \
                     [nbr for nbr in list(quotient_g.adj[picked_edge[1]]) if nbr not in picked_edge]
    
    removing_edge_attributes = quotient_g[picked_edge[0]][picked_edge[1]]
    # new_graph.remove_edge(picked_edge[0], picked_edge[1])
    new_graph.remove_node(picked_edge[0])
    new_graph.remove_node(picked_edge[1])

    new_node_name = '-'.join(str(n) for n in sorted(list(set([int(i) for i in picked_edge[0].split('-')] +
                                                             [int(i) for i in picked_edge[1].split('-')]))))
    
    new_graph.add_node(new_node_name, weight=removing_edge_attributes['weight'],
                       entropy=removing_edge_attributes['entropy'])
    adding_edges = [sorted([nod, new_node_name]) for nod in node_neighbors]
    # s6 = time.time()
    # print('                     ---> time removing 2 nodoes and adding 1 node:', s6 - s5)
    # print('                     ---> len of adding edges: ', len(adding_edges))
    # s7 = time.time()
    for ed in adding_edges:
        # print('                         ---> starting adding edge:', ed)
    
        poss = sorted(list(set(ed[0].split('-') + ed[1].split('-'))))
        # poss_genotype = ''.join(['1' for i in poss])
        poss_genotype = inpt_handler.get_genotype_positions([int(p) for p in poss])
        all_phasings = generate_phasings_ploidy_long(config.ploidy, poss_genotype, allel_set=config.alleles)
        
        # s8 = time.time()
        # print('                         ---> time for generate_phasings for edge:', ed, s8-s7)
        
        matches = get_matching_reads_for_positions([int(i) for i in poss], fragment_list)
        # s9 = time.time()
        # print('                         ---> time for get matching for edge:', ed, s9-s8)
        weights = {phas_2_str(phas): 0 for phas in all_phasings}
        for phas in all_phasings:
            for indc, this_po, obs in matches:
                # for obs in all_obs:
                #     obs_np = np.array([int(po) for po in obs])
                #     weights[phas_2_str(phas)] += compute_likelihood(obs_np, phas, error_rate)
                weights[phas_2_str(phas)] += compute_likelihood_generalized_plus(np.array(obs), phas, indc,
                                                                                 list(range(len(indc))),
                                                                                 config.error_rate)
                # s10 = time.time()
                # print('                             ---> time for compute likelihood:', ed, s10 - s9)
                # s9 = s10
                
        entr = entropy(list(weights.values()), base=10)
        new_graph.add_edge(ed[0], ed[1], weight=weights, entropy=entr)
        # s11 = time.time()
        # print('                         ---> finised adding edge:', ed, s11-s7)
        # s7 = s11
    return new_graph


def contract_one_edge2(quotient_g, picked_edge, inpt_handler, config, fragment_list):
    
    new_graph = quotient_g.copy()
    node_neighbors = [nbr for nbr in list(quotient_g.adj[picked_edge[0]]) if nbr not in picked_edge] + \
                     [nbr for nbr in list(quotient_g.adj[picked_edge[1]]) if nbr not in picked_edge]
    removing_edge_attributes = quotient_g[picked_edge[0]][picked_edge[1]]
    # new_graph.remove_edge(picked_edge[0], picked_edge[1])
    new_graph.remove_node(picked_edge[0])
    new_graph.remove_node(picked_edge[1])
    
    # plot_graph(new_graph)
    
    new_node_name = '-'.join(str(n) for n in sorted(list(set([int(i) for i in picked_edge[0].split('-')] +
                                                             [int(i) for i in picked_edge[1].split('-')]))))
    
    new_graph.add_node(new_node_name, weight=removing_edge_attributes['weight'],
                       entropy=removing_edge_attributes['entropy'])
    
    adding_edges = [sorted([nod, new_node_name]) for nod in node_neighbors]
    for ed in adding_edges:
        poss = sorted(list(set(ed[0].split('-') + ed[1].split('-'))))
        # poss_genotype = ''.join(['1' for i in poss])
        poss_genotype = inpt_handler.get_genotype_positions([int(p) for p in poss])
        all_phasings = generate_phasings_ploidy_long(config.ploidy, poss_genotype, allel_set=config.alleles)
        matches = get_matching_reads_for_positions([int(i) for i in poss], fragment_list)
        weights = {phas_2_str(phas): 0 for phas in all_phasings}
        for phas in all_phasings:
            for indc, this_po, obs in matches:
                # for obs in all_obs:
                #     obs_np = np.array([int(po) for po in obs])
                #     weights[phas_2_str(phas)] += compute_likelihood(obs_np, phas, error_rate)
                weights[phas_2_str(phas)] += compute_likelihood_generalized_plus(np.array(obs), phas, indc,
                                                                                 list(range(len(indc))),
                                                                                 config.error_rate)
        entr = entropy(list(weights.values()), base=10)
        new_graph.add_edge(ed[0], ed[1], weight=weights, entropy=entr)
    
    return new_graph

#
# def chordal_contraction(quotient_g, fragment_list, inpt_handler, config):
#     # plot_graph(quotient_g)
#     qg = quotient_g.copy()
#     while not nx.is_chordal(qg):
#         cliques_larger_than_2 = [cli for cli in nx.find_cliques(qg) if len(cli) > 2]
#
#         non_candicate_edges = []
#         for cli in cliques_larger_than_2:
#             non_candicate_edges += sorted(list(itertools.combinations(sorted(cli), 2)))
#         non_candicate_edges = list(set(non_candicate_edges))
#         # nx.is_chordal(quotient_g)
#         # cycles = nx.chordless_cycles(quotient_g)
#         candidate_edges = [edg for edg in list(qg.edges()) if edg not in non_candicate_edges]
#         entropies = np.array([qg[edg[0]][edg[1]]['entropy'] for edg in candidate_edges])
#         entropies[np.isnan(entropies)] = 0
#         rev_entropies = [1 - e if 1 > e else 0 for e in entropies]
#         if np.sum(rev_entropies) == 0:
#             for ent_id in range(len(rev_entropies)):
#                 rev_entropies[ent_id] = 1 / len(rev_entropies)
#         picked_edge = random.choices(candidate_edges, weights=rev_entropies, k=1)
#         picked_edge = [picked_edge[0][0], picked_edge[0][1]]
#         if picked_edge in [list(edg) for edg in list(qg.edges())]:
#             qg = contract_one_edge(qg, picked_edge, inpt_handler, config, fragment_list)
#         else:
#             print(f"edge {picked_edge} not in graph")
#     # plot_graph(qg)
#
#     return qg
#
