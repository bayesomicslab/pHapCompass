import numpy as np
import itertools
from haplotype_assembly_helper import generate_phasings_ploidy_long, compute_likelihood
from utils.utils import get_matching_reads_for_positions, phas_2_str
from scipy.stats import entropy
import random
import networkx as nx


def contract_one_edge(quotient_g, picked_edge, ploidy, fragment_list, error_rate):
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
        poss_genotype = ''.join(['1' for i in poss])
        all_phasings = generate_phasings_ploidy_long(ploidy, poss_genotype, allel_set=[0, 1])
        all_obs = get_matching_reads_for_positions([int(i) for i in poss], fragment_list)
        weights = {phas_2_str(phas): 0 for phas in all_phasings}
        for phas in all_phasings:
            for obs in all_obs:
                obs_np = np.array([int(po) for po in obs])
                weights[phas_2_str(phas)] += compute_likelihood(obs_np, phas, error_rate)
        entr = entropy(list(weights.values()), base=2)
        new_graph.add_edge(ed[0], ed[1], weight=weights, entropy=entr)
    return new_graph


def chordal_contraction(quotient_g, ploidy, fragment_list, error_rate):
    # plot_graph(quotient_g)
    qg = quotient_g.copy()
    while not nx.is_chordal(qg):
        cliques_larger_than_3 = [cli for cli in nx.find_cliques(qg) if len(cli) > 2]
        non_candicate_edges = []
        for cli in cliques_larger_than_3:
            non_candicate_edges += sorted(list(itertools.combinations(sorted(cli), 2)))
        non_candicate_edges = list(set(non_candicate_edges))
        # nx.is_chordal(quotient_g)
        # cycles = nx.chordless_cycles(quotient_g)
        candidate_edges = [edg for edg in list(qg.edges()) if edg not in non_candicate_edges]
        entropies = np.array([qg[edg[0]][edg[1]]['entropy'] for edg in candidate_edges])
        entropies[np.isnan(entropies)] = 0
        rev_entropies = [1 - e for e in entropies]
        picked_edge = random.choices(candidate_edges, weights=rev_entropies, k=1)
        picked_edge = [picked_edge[0][0], picked_edge[0][1]]
        if picked_edge in [list(edg) for edg in list(qg.edges())]:
            qg = contract_one_edge(qg, picked_edge, ploidy, fragment_list, error_rate)
        else:
            print(f"edge {picked_edge} not in graph")
    # plot_graph(qg)
    
    return qg
