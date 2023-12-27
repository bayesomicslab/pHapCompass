from pgmpy.inference import BeliefPropagation
import networkx as nx
import random
import pandas as pd
from pgmpy.sampling import GibbsSampling
from pgmpy.models import MarkovNetwork
from pgmpy.factors.discrete import DiscreteFactor


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

def sample_gibbs(source, target, beliefs, qg):
    n_samples = 1000
    paths = list(nx.all_shortest_paths(qg, source, target))
    if len(paths) > 0:
        spath = random.choices(paths, k=1)[0]
        
        path_graph = path_2_markov_network(spath, beliefs)
        gibbs_sample_df = path_sample_gibbs(path_graph, n_samples, qg)
    else:
        gibbs_sample_df = None
    return gibbs_sample_df
