import numpy as np
from scipy.sparse import csr_matrix


def str_2_phas(phasings, ploidy):
    return [np.array([int(p) for p in [*phas]]).reshape(ploidy, -1) for phas in phasings]


def str_2_phas_1(phasing, ploidy):
    return np.array([int(p) for p in [*phasing]]).reshape(ploidy, -1)



def phas_2_str(phas):
    return ''.join([str(ph) for ph in list(np.ravel(phas))])


def load_test_data():
    # Load test data
    frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
    # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test_clique.frag'
    genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
    ploidy = 3
    class Args:
        def __init__(self):
            # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
            # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = frag_path
            # self.data_path = frag_path
            self.bam_path = 'example/example.bam'
            self.genotype_path = genotype_path
            self.ploidy = ploidy
            self.error_rate = 0.001
            self.epsilon = 0.0001
            self.output_path = 'output'
            self.root_dir = 'D:/UCONN/HaplOrbit'
            self.alleles = [0, 1]
    args = Args()
    return args


def load_test_data_long():
    # Load test data
    frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/sim_long_1based2.frag'
    genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/sim_long_haplotypes.csv'
    ploidy = 3
    class Args:
        def __init__(self):
            # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
            # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = frag_path
            # self.data_path = frag_path
            self.bam_path = 'example/example.bam'
            self.genotype_path = genotype_path
            self.ploidy = ploidy
            self.error_rate = 0.001
            self.epsilon = 0.0001
            self.output_path = 'output'
            self.root_dir = 'D:/UCONN/HaplOrbit'
            self.alleles = [0, 1]
    args = Args()
    return args


def load_sassafras():

    # Load test data
    frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/sassafras/CP142451.1_RagTag.frag'
    genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/sassafras/hapl.csv'
    ploidy = 4
    class Args:
        def __init__(self):
            # frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/test.frag'
            # genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/FFBS/haplotypes.csv'
            self.vcf_path = 'example/62_ID0.vcf'
            self.data_path = frag_path
            # self.data_path = frag_path
            self.bam_path = 'example/example.bam'
            self.genotype_path = genotype_path
            self.ploidy = ploidy
            self.error_rate = 0.001
            self.epsilon = 0.0001
            self.output_path = 'output'
            self.root_dir = 'D:/UCONN/HaplOrbit'
            self.alleles = [0, 1]
    args = Args()
    return args


def _check_edge_labels(edges):
    for u, v in edges:
        ui, uj = map(int, u.split("-"))
        vi, vj = map(int, v.split("-"))
        shared = {ui, uj} & {vi, vj}
        assert len(shared) == 1, f"edge does not share exactly one SNP: {u}--{v}"
        s = next(iter(shared))
        # ensure topological u={s,a}, v={s,b} with a<b
        a = (ui if ui != s else uj)
        b = (vi if vi != s else vj)
        assert a < b, f"edge not in topological order: {u}--{v}"


def sort_nodes(strings):
    # Sort the list of strings using the custom comparison logic
    return sorted(strings, key=lambda x: list(map(int, x.split('-'))))

def sort_edges(strings):
    # Sort the list of strings using the custom key
    return sorted(strings, key=lambda x: [
        list(map(int, part.split('-'))) for part in x.split('--')
    ])

def compute_likelihood(observed, phasing, error_rate):
    """This likelihood computation assumes the length of observation is the same as the length of phasing"""
    y = np.tile(observed, (phasing.shape[0], 1))
    diff = y - phasing
    diff[diff != 0] = 1
    comp_diff = 1 - diff
    term1 = diff * error_rate
    term2 = comp_diff * (1 - error_rate)
    terms = term1 + term2
    probs = np.prod(terms, axis=1)
    likelihood = np.sum(probs)
    return likelihood



def build_read_list_from_M(M: csr_matrix, one_based: bool = False) -> list:
    """
    Convert M (reads x SNPs; 0=no-call, 1=REF, 2=ALT) to read list format.
    
    Args:
        M: Sparse read matrix (n_reads x n_snps)
        one_based: If True, positions are 1-based; if False, 0-based
        
    Returns:
        read_list: [pos_array0, allele_array0, pos_array1, allele_array1, ...]
                   where alleles are 0=REF, 1=ALT
    """
    read_list = []
    n_reads = M.shape[0]
    
    for r in range(n_reads):
        row = M[r, :]
        cols = row.indices  # nonzero columns (0-based)
        vals = row.data     # {1, 2}
        
        if cols.size == 0:
            continue
        
        # Sort by column index (usually already sorted, but guarantee it)
        order = np.argsort(cols)
        cols_sorted = cols[order]
        vals_sorted = vals[order]
        
        # Convert to alleles: REF(1)→0, ALT(2)→1
        alleles = (vals_sorted == 2).astype(np.int32)
        
        # Convert to 1-based if needed
        if one_based:
            positions = cols_sorted + 1
        else:
            positions = cols_sorted
        
        read_list.append(positions.tolist())
        read_list.append(alleles.tolist())
    
    return read_list
