import numpy as np
import itertools
from algorithm.haplotype_assembly_helper import generate_phasings_ploidy_long
from utils.utils import *
from scipy.stats import entropy
import random
import networkx as nx
import time
import graph_tool.all as gt


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


def contract_one_edge3(quotient_g, picked_edge, inpt_handler, config, fragment_list):
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


def get_cycles_basis_info(G):
    cycle_basis_orig = nx.minimum_cycle_basis(G)
    _ = [cyb.append(cyb[0]) for cyb in cycle_basis_orig]
    cycle_basic = [cyc for cyc in cycle_basis_orig]
    ready_cycles = [cyc for cyc in cycle_basis_orig if len(cyc) <= 4]
    ready_edges = [[G.edges[e]['original_order'] for e in zip(cyb, cyb[1:])] for cyb in ready_cycles]
    forbidden_edges = [x for xs in ready_edges for x in xs]
    # cycles_edges = [[get_edge_name(f, s) for f, s in zip(cyb, cyb[1:])] for cyb in cycle_basic]
    cycles_edges = [[G.edges[e]['original_order'] for e in zip(cyb, cyb[1:])] for cyb in cycle_basic]
    index_to_cycles = {i: cycle for i, cycle in enumerate(cycles_edges)}
    edge_to_cycles = {G.edges[e]['original_order']: [] for e in G.edges()}
    for cycle_index, cycle in index_to_cycles.items():
        [edge_to_cycles[ed].append(cycle_index) for ed in cycle]
    return cycle_basis_orig, index_to_cycles, edge_to_cycles, forbidden_edges


def get_cycles_basis_info_graph_tool(non_mst_graph, mst_graph, quotient_graph):
    cycle_basis_edges = []
    for e in non_mst_graph.edges():
        
        cycle = find_cycle(mst_graph, quotient_graph, e)
        cycle_basis_edges.append(cycle)
    
    ready_cycles = [cyc for cyc in cycle_basis_edges if len(cyc) <= 3]
    ready_edges = [[e for e in cyb] for cyb in ready_cycles]
    forbidden_edges = [x for xs in ready_edges for x in xs]

    cycles_edges = [[e for e in cyb] for cyb in cycle_basis_edges]
    index_to_cycles = {i: cycle for i, cycle in enumerate(cycles_edges)}
    edge_to_cycles = {e: [] for e in quotient_graph.edges()}
    for cycle_index, cycle in index_to_cycles.items():
        [edge_to_cycles[ed].append(cycle_index) for ed in cycle]
    return cycle_basis_edges, index_to_cycles, edge_to_cycles, forbidden_edges


def find_cycle(mst_graph, quotient_graph, e):
    """Finds a cycle in the graph after adding the given edge to the spanning tree"""
    # Add the edge to the MST
    # source, target = int(e.source()), int(e.target())
    source, target = e.source(),e.target()
    # mst_tree = gt.Graph(gt.GraphView(graph, efilt=mst_view))
    # Traverse the tree to find the path between source and target (forming a cycle)
    # mst_tree = gt.shortest_path(mst_view, graph.vertex(source), graph.vertex(target))[0]

    # tree_graph = gt.Graph(mst_graph)

    half_cyc = gt.shortest_path(mst_graph, source, target)
    
    # The cycle is the path from source to target in the MST + the edge being added
    cycle_e = half_cyc[1] + [e]
    return cycle_e, half_cyc[0]


def get_edge_name(a, b):
    mi = min(a, b)
    ma = max(a, b)
    return (mi, ma)


def pick_contracted_edge(cyc, new_graph):
    ents = [new_graph.edges[edg]['entropy'] for edg in cyc]
    sel = cyc[np.argmin(ents)]
    return sel


def chordal_contraction_cycle_base(quotient_g, fragment_list, inpt_handler, config):
    for edg in quotient_g.edges():
        quotient_g.edges[edg]['original_order'] = tuple(sorted(tuple(edg)))

    new_graph = quotient_g.copy()
    # plot_graph(new_graph)
    cycle_basic, index_to_cycles, edge_to_cycles, forbidden_edges = get_cycles_basis_info(new_graph)
    total_entropy = np.sum([new_graph.edges[edg]['entropy'] for edg in new_graph.edges()])
    print(total_entropy)
    for i in range(len(index_to_cycles.keys())):
        # print(i, len(index_to_cycles[i]))
        cyc = index_to_cycles[i]
        # print(i, len(cyc))
        # if len(cyc) > 3:
        #     print(i, len(cyc))
            
        while len(cyc) > 3:
            
            picked_edge = pick_contracted_edge(cyc, new_graph)
            print(picked_edge)
            removing_edge_attributes = new_graph[picked_edge[0]][picked_edge[1]]
            
            # 0: nodes neighbors
            ns = [nn for nn in new_graph.neighbors(picked_edge[0]) if nn != picked_edge[1]]
            ne = [nn for nn in new_graph.neighbors(picked_edge[1]) if nn != picked_edge[0]]
            
            # 1: figure out new node:
            new_node_name = '-'.join(str(n) for n in sorted(list(set([int(i) for i in picked_edge[0].split('-')] +
                                                                     [int(i) for i in picked_edge[1].split('-')]))))

            new_graph.add_node(new_node_name, weight=removing_edge_attributes['weight'],
                               entropy=removing_edge_attributes['entropy'])
            
            # 2: remove old edge and add new edge and update dictionaries
            incident_edges = list(set([(picked_edge[0], n) for n in ns] + [(picked_edge[1], n) for n in ne]))
            for ie in incident_edges:
                orig_order = new_graph.get_edge_data(*ie)['original_order']
                neighbor = list(set(ie) - set(picked_edge))[0]
                new_edge = tuple(sorted((neighbor, new_node_name)))

                poss = sorted(list(set(new_edge[0].split('-') + new_edge[1].split('-'))))
                # poss_genotype = ''.join(['1' for i in poss])
                poss_genotype = inpt_handler.get_genotype_positions([int(p) for p in poss])
                all_phasings = generate_phasings_ploidy_long(config.ploidy, poss_genotype, allel_set=config.alleles)

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
                entr = entropy(list(weights.values()), base=10)
                new_graph.add_edge(new_edge[0], new_edge[1], weight=weights, entropy=entr, original_order=new_edge)
                new_graph.remove_edges_from([orig_order])
                
                edge_to_cycles[new_edge] = []
                edge_cycles = edge_to_cycles[orig_order]
                for ecyc in edge_cycles:
                    index_to_cycles[ecyc].remove(orig_order)
                    index_to_cycles[ecyc].append(new_edge)
                    edge_to_cycles[new_edge].append(ecyc)
                
            # 3: update information about the contracted edge
            for cy in edge_to_cycles[picked_edge]:
                index_to_cycles[cy].remove(picked_edge)
            
            # 4: remove contracted edge
            new_graph.remove_edges_from([picked_edge])
            
            # 5: remove contracted nodes
            new_graph.remove_nodes_from([picked_edge[0], picked_edge[1]])
            
            # 6: update dictinary
            cyc = index_to_cycles[i]
            total_entropy = np.sum([new_graph.edges[edg]['entropy'] for edg in new_graph.edges()])
            print(total_entropy)
            # plot_graph(new_graph)
    return new_graph


def is_chordal(graph):
    """
    Check if the given graph is chordal using Lexicographical Breadth-First Search (Lex-BFS).

    Parameters:
    
    graph: A graph_tool Graph object.

    Returns:
    
    True if the graph is chordal, False otherwise."""
    labels = {}
    for v in graph.vertices():
        labels[v] = []

    ordering = []
    unnumbered = set(graph.vertices())

    n = graph.num_vertices()
    for i in range(n, 0, -1):
        # Select the unnumbered vertex with the largest label lexicographically
        max_label_vertex = None
        max_label = None
        for v in unnumbered:
            label = labels[v]
            if max_label is None or label > max_label:
                max_label = label
                max_label_vertex = v

        v = max_label_vertex
        ordering.append(v)
        unnumbered.remove(v)

        # Update labels of unnumbered neighbors
        for u in v.all_neighbors():
            if u in unnumbered:
                labels[u].append(i)

    # Reverse the ordering to get the perfect elimination ordering (PEO)
    ordering.reverse()
    position = {v: idx for idx, v in enumerate(ordering)}
    for idx, v in enumerate(ordering):
        nbrs = [u for u in v.all_neighbors() if position[u] > idx]
        if not nbrs:
            continue
        # Find the neighbor with the smallest position (earliest in PEO)
        w = min(nbrs, key=lambda u: position[u])
        for u in nbrs:
            if u != w and not graph.edge(u, w) and not graph.edge(w, u):
                # The neighbors do not form a clique
                return False
    return True

# @profile
def get_chordless_cycles(subgraph_copy):
    """
    Check if the given graph is chordal using Lexicographical Breadth-First Search (Lex-BFS).

    Parameters:
    
    graph: A graph_tool Graph object.

    Returns:
    
    True if the graph is chordal, False otherwise."""
    chordless_cycles = []
    labels = {}
    for v in subgraph_copy.vertices():
        labels[v] = []

    ordering = []
    unnumbered = set(subgraph_copy.vertices())

    n = subgraph_copy.num_vertices()
    for i in range(n, 0, -1):
        # print(i)
        # Select the unnumbered vertex with the largest label lexicographically
        max_label_vertex = None
        max_label = None
        for v in unnumbered:
            label = labels[v]
            if max_label is None or label > max_label:
                max_label = label
                max_label_vertex = v

        v = max_label_vertex
        ordering.append(v)
        unnumbered.remove(v)

        # Update labels of unnumbered neighbors
        for u in v.all_neighbors():
            if u in unnumbered:
                labels[u].append(i)

    # Reverse the ordering to get the perfect elimination ordering (PEO)
    ordering.reverse()
    position = {v: idx for idx, v in enumerate(ordering)}
    for idx, v in enumerate(ordering):
        # Check if the neighbors of v that appear later in the PEO form a clique
        
        nbrs = [u for u in v.all_neighbors() if position[u] > idx]
        if not nbrs:
            continue
        # Find the neighbor with the smallest position (earliest in PEO)
        w = min(nbrs, key=lambda u: position[u])
        # print('min neighbor:', w, 'neighbors:', nbrs)
        for u in nbrs:
            if u != w and not subgraph_copy.edge(u, w) and not subgraph_copy.edge(w, u):
                # The neighbors do not form a clique
                # print('ordering:', v, 'neighbor', u, 'min neighbor',w)
                # gt.shortest_path(subgraph_copy, u, w)
                v_filter = subgraph_copy.new_vertex_property("bool", val=True)
                v_filter[v] = False
                subgraph_copy.set_vertex_filter(v_filter)
                # Exclude the specified node by setting its filter to False
                # v_filter[subgraph_copy.vertex(v)] = False
                # subgraph_2 = gt.GraphView(subgraph_copy, vfilt=v_filter)
    
                # Find the shortest path in the filtered graph
                path_v = gt.shortest_path(subgraph_copy, source=w, target=u)[0]
                # path_edges = path_edges[::-1]
                path_v = path_v + [v]
                # print(path_v)
                chordless_cycles.append(path_v)
                # subgraph_copy.set_vertex_filter(None)

                # Reset the vertex filter after the computation
                subgraph_copy.set_vertex_filter(None)

    return chordless_cycles
   

def chordal_contraction_graph_tool(quotient_graph, config, fragment_model):
    # more efficient version of the chordal contraction
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
                sorted_labels = sort_nodes([new_vertex_name, v_label])
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

                sorted_labels = sort_nodes([new_vertex_name, v_label])
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
                sorted_labels = sort_nodes([new_vertex_name, v_label])
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
    
# @profile
def get_top_k_weights(mydict, k):
    # Extract the 'weight' dictionary
    weights = mydict['weight']
    
    # Sort the weights dictionary by value in descending order and take the top k items
    top_k_weights = dict(sorted(weights.items(), key=lambda item: item[1], reverse=True)[:k])
    
    # Return a new dictionary with the top k weights and the same entropy
    return {'weight': top_k_weights, 'entropy': mydict['entropy']}

# @profile
def chordal_contraction_graph_tool_approx(inp):
    save_path, subg_id, subg, config, fragment_model, k = inp
    this_path = os.path.join(save_path, 'chordal_sub_' + str(subg_id) + '.gt.gz')
    print('Working on', this_path)
    new_graph = subg.copy()
    e_weights = new_graph.edge_properties["e_weights"]
    # new_graph.clear_filters()
    e_entropy = new_graph.new_edge_property("double")

    # Loop over edges and assign entropy from the e_weights property
    for e in new_graph.edges():
        e_entropy[e] = e_weights[e]['entropy']

    new_graph.ep['e_entropy'] = e_entropy

    chordless_cycles = get_chordless_cycles(new_graph)

    to_be_removed_nodes = []

    for cyc_id, cyc in enumerate(chordless_cycles):
        
        # print(cyc_id)
        
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
            vertex_weights_appr = get_top_k_weights(vertex_weights, k)

            new_graph.vertex_properties["v_weights"][min_edge.source()] = vertex_weights_appr
            new_graph.vertex_properties["v_label"][min_edge.source()] = new_vertex_name

            source_nbrs = [n for n in min_edge.source().all_neighbors() if n != min_edge.target()]
            target_nbrs = [n for n in min_edge.target().all_neighbors() if n != min_edge.source()]
            common_nbrs = set(source_nbrs).intersection(set(target_nbrs))

            for n in common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]
                # e_poss = sorted(set([int(nn) for nn in v_label.split('-')] + poss))
                # print(len(e_poss))
                # new_edge_name = '-'.join([str(nnn) for nnn in e_poss])
                sorted_labels = sort_nodes([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)
                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)
                final_weights_appr = get_top_k_weights(final_weight, k)

                e1 = new_graph.edge(min_edge.source(), n)
                e2 = new_graph.edge(min_edge.target(), n)
                
                new_graph.edge_properties["e_weights"][e1] = final_weights_appr
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weights_appr['entropy']
                new_graph.remove_edge(e2)

            for n in set(source_nbrs)-common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]

                sorted_labels = sort_nodes([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)
                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)
                final_weights_appr = get_top_k_weights(final_weight, k)


                e1 = new_graph.edge(min_edge.source(), n)
                # e2 = new_graph.edge(min_edge.target(), n)
                new_graph.edge_properties["e_weights"][e1] = final_weights_appr
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weights_appr['entropy']
                # new_graph.edge_properties["e_weights"][e2]

            
            for n in set(target_nbrs)-common_nbrs:
                
                v_label = new_graph.vertex_properties["v_label"][n]
                sorted_labels = sort_nodes([new_vertex_name, v_label])
                new_edge_name = '--'.join(sorted_labels)

                (first_label, first_node), (second_label, second_node) = [(new_vertex_name, min_edge.source()),(v_label, n)]

                first_phasings = list(new_graph.vertex_properties["v_weights"][first_node]['weight'].keys())
                second_phasings = list(new_graph.vertex_properties["v_weights"][second_node]['weight'].keys())
                final_weight = compute_edge_weight(first_label, second_label, first_phasings, second_phasings, fragment_model, config)
                final_weights_appr = get_top_k_weights(final_weight, k)

                e2 = new_graph.edge(min_edge.target(), n)
                new_graph.remove_edge(e2)
                e1 = new_graph.add_edge(min_edge.source(), n)
                new_graph.edge_properties["e_weights"][e1] = final_weights_appr
                new_graph.edge_properties["e_label"][e1] = new_edge_name
                new_graph.edge_properties['e_entropy'][e1] = final_weights_appr['entropy']
            
            # to_be_removed_nodes += [min_edge.target()]
            to_be_removed_nodes.append(min_edge.target())
            new_graph.remove_edge(min_edge)
            edges.remove(min_edge)

    new_graph.remove_vertex(to_be_removed_nodes)
    new_graph.save(this_path)
    print('[Done]', this_path)

# @profile
def divide_graph_by_labels(graph, m):
    # Retrieve vertex labels
    vertex_labels = graph.vp["v_label"]

    # Create a sorted list of vertices based on their labels
    sorted_vertices = sorted(
        graph.vertices(), key=lambda v: tuple(map(int, vertex_labels[v].split('-')))
    )

    # Divide sorted vertices into chunks of approximately 'm' vertices
    subgraphs = []
    current_chunk = []

    for vertex in sorted_vertices:
        current_chunk.append(vertex)
        if len(current_chunk) >= m:
            # Create a subgraph for the current chunk
            subgraph = gt.GraphView(graph, vfilt=lambda v, chunk=current_chunk: v in chunk)
            subgraphs.append(subgraph)
            current_chunk = []  # Reset for the next chunk

    # Handle any remaining vertices in the last chunk
    if current_chunk:
        subgraph = gt.GraphView(graph, vfilt=lambda v, chunk=current_chunk: v in chunk)
        subgraphs.append(subgraph)

    return subgraphs