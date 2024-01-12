from pgmpy.inference import VariableElimination, BeliefPropagation
import random
from pgmpy.sampling import GibbsSampling
from pgmpy.models import MarkovNetwork
from pgmpy.factors.discrete import DiscreteFactor
from collections import Counter
from utils.utils import *


def factor_graph_inference(factor_graph):
    beliefs = BeliefPropagation(factor_graph)
    return beliefs


def path_2_markov_network(spath, beliefs):
    path_graph = MarkovNetwork()
    for i in range(len(spath) - 1):
        nod1, nod2 = spath[i:i + 2]
        path_graph.add_edge(nod1, nod2)
        
        joint = beliefs.query(variables=[nod1, nod2], show_progress=False).values
        
        edge_potential = DiscreteFactor(variables=[nod1, nod2], cardinality=[joint.shape[0], joint.shape[1]],
                                        values=joint)
        # Add edge potentials to graph
        path_graph.add_factors(edge_potential)
        
        # Add factor nodes and connect them to variable nodes
        # path_graph.add_edges_from([(nod1, edge_potential), (nod2, edge_potential)])
    
    return path_graph


def path_sample_gibbs(path_graph, n_samples, qg):
    gibbs_chain = GibbsSampling(path_graph)
    sample = gibbs_chain.sample(size=n_samples)
    
    sample_df_phases = pd.DataFrame(columns=sample.columns.values, index=range(len(sample)))
    for col in sample.columns.values:
        sample_df_phases[col] = sample[col].apply(lambda x: list(qg.nodes[col]['weight'].keys())[x])
    return sample_df_phases

def sample_gibbs(source, target, beliefs, qg, n_samples):
    paths = list(nx.all_shortest_paths(qg, source, target))
    if len(paths) > 0:
        spath = random.choices(paths, k=1)[0]
        
        path_graph = path_2_markov_network(spath, beliefs)
        gibbs_sample_df = path_sample_gibbs(path_graph, n_samples, qg)
    else:
        gibbs_sample_df = None
    return gibbs_sample_df


def query_paths_gibbs(fragment_list, qg, beliefs, n_samples):
    source, target, positions = get_first_last_positions(fragment_list, qg)
    samples = sample_gibbs(source, target, beliefs, qg, n_samples)
    
    sh_dict = make_shared_dict(samples)
    columns = list(samples.columns.values)
    
    samples['path_phasing'] = samples.apply(lambda row: merge_one_row(columns, row, sh_dict), axis=1)
    return samples['path_phasing']

def query_paths_gibbs_max(fragment_list, qg, beliefs, n_samples):
    source, target, _ = get_first_last_positions(fragment_list, qg)
    samples = sample_gibbs(source, target, beliefs, qg, n_samples)

    sh_dict = make_shared_dict(samples)
    columns = list(samples.columns.values)

    samples['path_phasing'] = samples.apply(lambda row: merge_one_row(columns, row, sh_dict), axis=1)
    
    counts = Counter(samples['path_phasing'])

    max_count = max(counts.values())

    max_phase = [num for num, count in counts.items() if count == max_count][0]
    
    obs_positions = [elem.split('-') for elem in list(samples.columns.values)]
    positions = [item for sublist in obs_positions for item in sublist]
    positions = sorted(set([int(pos) for pos in positions]))

    return max_phase, positions

def phas2hap(inp, ploidy):
    if isinstance(inp, str):
        phas = str_2_phas([inp], ploidy)[0]
    else:
        phas = inp
    haplotypes = [''.join([str(i) for i in (list(phas[r, :]))]) for r in range(ploidy)]
    return haplotypes


def get_first_last_positions(fragment_list, qg):
    obs_positions = fragment_list[::2]
    positions = [item for sublist in obs_positions for item in sublist]
    poss = sorted(set(positions))
    min_p = str(min(poss))
    max_p = str(max(poss))
    min_node = ''
    while min_node == '':
        for node in qg.nodes():
            # print(node, min_p in node.split('-'))
            if min_p in node.split('-'):
                min_node = node
                break
    
    max_node = ''
    while max_node == '':
        for node in list(qg.nodes())[::-1]:
            if max_p in node.split('-'):
                max_node = node
                break
    return min_node, max_node, poss