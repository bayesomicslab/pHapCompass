
from itertools import permutations
import numpy as np
import pandas as pd
import pickle
import os
import itertools
from itertools import combinations, permutations
from scipy.optimize import linear_sum_assignment
from collections import defaultdict
from utils.utils import sort_nodes, str_2_phas_1, compute_likelihood
import graph_tool.all as gt


def calculate_correct_phasing_rate(reconstructed_haplotypes, true_haplotypes):
    n = reconstructed_haplotypes.shape[1]  # Number of SNPs
    k = reconstructed_haplotypes.shape[0]  # Number of haplotypes

    similarity_score = 0
    for i in range(k):
        for j in range(n):
            if reconstructed_haplotypes[i, j] == true_haplotypes[i, j]:
                similarity_score += 1

    # Calculate the Correct Phasing Rate (Rc)
    Rc = similarity_score / (n * k)
    return Rc


def calculate_perfect_solution_rate(reconstructed_haplotypes, true_haplotypes):
    k = reconstructed_haplotypes.shape[0]  # Number of haplotypes
    perfect_count = 0

    # Check if each haplotype matches perfectly with the corresponding true haplotype
    for i in range(k):
        if np.array_equal(reconstructed_haplotypes[i], true_haplotypes[i]):
            perfect_count += 1

    # Calculate the Perfect Solution Rate (Rp)
    Rp = perfect_count / k
    return Rp


def find_matches(h_star_col, h_col):
  """
  input:
  h_col : list or array, one allele position across all haplotypes in true phase
  h_star_col : list or array, one allele position across all haplotypes in assembled phase
  returns:
  list of tuples, each tuple is a mapping of the indices 0, ..., k-1 in h_star_col to the indices 0, ..., k-1 in h_col
  """
  # find where all of the 1s and 0s are in H and H star
  one_indices_h = [i for i, val in enumerate(h_col) if val == 1]
  zero_indices_h = [i for i, val in enumerate(h_col) if val == 0]

  one_indices_h_star = [i for i, val in enumerate(h_star_col) if val == 1]
  zero_indices_h_star = [i for i, val in enumerate(h_star_col) if val == 0]

  # separately find all of the ways to permute the 1s and 0s in H star
  # this will give us all of the ways to pair elements in H to elements in H star
  all_one_permutations_h_star = list(permutations(one_indices_h_star))
  all_zero_permutations_h_star = list(permutations(zero_indices_h_star))

  # need to consider the cross product of
  # every way to match the 1s with every way to match the 0s
  # to find the full match which we store as a dictionary
  # dictionary maps index in H to index in H star
  matches = []
  for permutation_of_ones in all_one_permutations_h_star:
    for permutation_of_zeros in all_zero_permutations_h_star:

      temp_dict = {}

      for i, idx in enumerate(one_indices_h):
        temp_dict[idx] = permutation_of_ones[i]

      for i, idx in enumerate(zero_indices_h):
        temp_dict[idx] = permutation_of_zeros[i]

      matching = tuple(temp_dict[i] for i in range(len(temp_dict)))


      matches.append(matching)

  return matches


def find_partial_matches(h_star_col, h_col):
    """
    returns all of the ways to match the alleles
    in a partial column of H_star (length m) to a column of H (length k)
    """
    m = len(h_star_col)
    
    # Get indices for allele 1 and 0 in the reconstructed column (length m).
    one_indices_star = [i for i, val in enumerate(h_star_col) if val == 1]
    zero_indices_star = [i for i, val in enumerate(h_star_col) if val == 0]
    
    # Get indices for allele 1 and 0 in the ground truth column (length k).
    one_indices = [j for j, val in enumerate(h_col) if val == 1]
    zero_indices = [j for j, val in enumerate(h_col) if val == 0]
    
    matches = []
    # choose which of the 1s we will be using in the m-subset
    # then consider all orderings (matchings) of that chosen subset.
    for chosen_ones in combinations(one_indices, len(one_indices_star)):
        for perm_ones in permutations(chosen_ones):
            # Similarly for zeros.
            for chosen_zeros in combinations(zero_indices, len(zero_indices_star)):
                for perm_zeros in permutations(chosen_zeros):
                    # Build a mapping vector of length m.
                    mapping = [None] * m
                    for idx, gt_idx in zip(one_indices_star, perm_ones):
                        mapping[idx] = gt_idx
                    for idx, gt_idx in zip(zero_indices_star, perm_zeros):
                        mapping[idx] = gt_idx
                    matches.append(tuple(mapping))
    return matches


def compute_vector_error_rate(H, H_star):
  """
  input:
  H and H_star : numpy arrays of the true and assembled phasings
  retuns:
  the value of the vector error
  the backtracking steps - how we have to permute each column of H to make them match H_star
  """
  
  k, n = H.shape

  # i^th dictionary stores the ways to get through the i^th col
  # dictionary stores:
  # key = mapping from H to H star at that column
  # value = (vector error using that key, mapping form prev column that induces that vector error)

  dp_table = [{} for _ in range(n)]

  # vector error is 0 when we only consider the first column
  first_col_matchings = find_matches(H_star[:, 0], H[:, 0])
  first_col = H_star[:,0]
  m0 = np.count_nonzero(first_col == 0)
  m1 = np.count_nonzero(first_col == 1)
  for first_col_matching in first_col_matchings:
    dp_table[0][first_col_matching] = (int(m0 > 1)+int(m1 > 1), None)

  # Forward pass
  for col in range(1, n):
    all_matches = find_matches(H_star[:, col], H[:, col])
    # for every way to match that column of H to H star
    for matching in all_matches:
      best_ve = np.inf
      best_prev = None
      for prev_matching in dp_table[col-1].keys():
        prev_cost = dp_table[col-1][prev_matching][0]
        switches = np.count_nonzero(
          np.array(prev_matching) != np.array(matching)
        )
        # vector error = vector error up to the previous step + number of switches we take at that step
        temp_ve = prev_cost + switches
        if temp_ve < best_ve:
          best_ve = temp_ve
          best_prev = prev_matching
      if best_prev is not None:
          dp_table[col][matching] = (best_ve, best_prev)

  # Grab the minimal cost from the last column
  vector_error, _ = min(dp_table[n-1].values(), key=lambda x: x[0])

  vector_error_rate = vector_error/n
  # Identify the matching that gave that minimal cost
  last_matching = min(dp_table[n-1], key=lambda m: dp_table[n-1][m][0])

  # Backtracking
  backtracking_steps = []
  matching = last_matching

  # Loop through columns in reverse order
  for col in range(n-1, -1, -1):
      backtracking_steps.append(matching)
      cost, prev_match = dp_table[col][matching]
      if prev_match is None:
          break
      matching = prev_match


  backtracking_steps.reverse()

  return vector_error_rate, vector_error, backtracking_steps, dp_table


def compute_vector_error_rate_partial(H, H_star_full):
    '''
    for full ground truth matrix H and 
    partial reconstruction H_star_full (k by n but with only m rows phased)
      assuming that all rows are either fully phased or fully unphased 
    computes vector error rate for the phased rows of H_star 
    some of the rows in the ground truth matrix will not be matched with any reconstructed rows for a given SNP, 
      maybe it is matched for a part of the reconstructed row, 
      then the reconstructed row switches to another ground truth row and 
      incurs a vector error and that ground truth row is no longer used,
      but m haplotypes are used for any given SNP (not always the same m for different SNPs)
    '''
    k, n = H.shape
    H_star = H_star_full[~np.isnan(H_star_full).any(axis=1),:]
    m = H_star.shape[0]
    dp_table = [{} for _ in range(n)]
    
    # Use the new helper for the first column.
    first_col_matchings = find_partial_matches(H_star[:, 0], H[:, 0])
    first_col = H_star[:,0]
    m0 = np.count_nonzero(first_col == 0)
    m1 = np.count_nonzero(first_col == 1)
    for matching in first_col_matchings:
        dp_table[0][matching] = (int(m0 > 1)+int(m1 > 1), None)
    
    # Forward pass: iterate over columns.
    for col in range(1, n):
        all_matches = find_partial_matches(H_star[:, col], H[:, col])
        for matching in all_matches:
            best_ve = np.inf
            best_prev = None
            for prev_matching in dp_table[col-1]:
                prev_cost = dp_table[col-1][prev_matching][0]
                # Count the number of reconstructed rows that switch their ground truth match.
                switches = np.count_nonzero(np.array(prev_matching) != np.array(matching))
                temp_ve = prev_cost + switches
                if temp_ve < best_ve:
                    best_ve = temp_ve
                    best_prev = prev_matching
            if best_prev is not None:
                dp_table[col][matching] = (best_ve, best_prev)
    
    # Backtracking to obtain the matching sequence.
    vector_error, _ = min(dp_table[n-1].values(), key=lambda x: x[0])
    vector_error_rate = vector_error / n
    last_matching = min(dp_table[n-1], key=lambda m: dp_table[n-1][m][0])
    backtracking_steps = []
    matching = last_matching
    for col in range(n-1, -1, -1):
        backtracking_steps.append(matching)
        cost, prev_match = dp_table[col][matching]
        if prev_match is None:
            break
        matching = prev_match
    backtracking_steps.reverse()

    return vector_error_rate, vector_error, backtracking_steps, dp_table


def map_back_proportion(reconstructed_haplotypes, list_of_reads):
  """
  input:
  reconstucted haplotypes - numpy array assembled phasing
  list of reads - list of lists
    odd index lists are positions, even index reads are alleles
  returns:
  the proportion of reads that could be mapped to some row of H_star at their specified position
  """
  pos_allele_pairs = [(list_of_reads[i], list_of_reads[i + 1]) for i in range(0, len(list_of_reads), 2)]

  map_back = np.zeros(len(pos_allele_pairs))

  for idx, (position_seq, allele_seq) in enumerate(pos_allele_pairs):
      for reconstructed_hap in reconstructed_haplotypes:
          subarray = reconstructed_hap[position_seq[0] : position_seq[-1]+1]

          if np.array_equal(subarray, allele_seq):
              map_back[idx] = 1
              # once we find a match for this read, we can stop checking other rows
              break
  return np.sum(map_back)/len(pos_allele_pairs)

'''
map each read to the haplotype that it fits with best
returns the proportion of alleles that have to be coerced to make them match
** since we are taking the haplotype that matches best, it should never match with a nan row 
'''
def mec(reconstructed_haplotypes, list_of_reads):
  pos_allele_pairs = [(list_of_reads[i], list_of_reads[i + 1]) for i in range(0, len(list_of_reads), 2)]

  error_coercion_count = 0
  total_reads_size = 0
  
  for idx, (position_seq, allele_seq) in enumerate(pos_allele_pairs):
    subarray = reconstructed_haplotypes[position_seq[0] : position_seq[-1]+1]
    best_match_haplotype = np.argmax(np.sum(subarray==allele_seq, axis=1))
    error_coercion_count += np.sum(best_match_haplotype != allele_seq)
    total_reads_size += len(allele_seq)
  
  return error_coercion_count/total_reads_size

'''
helper for probabalistic MEC
likelihood that a read is truly from a given haplotype
works for unphased (nan)
potential edit - instead of having a penalty of 1 for unphased,
define penalty based on allele frequency
(penalty = prob that it is not the correct one)
'''
def likelihood_haplotype(read, haplotype, seq_error=0.001):
  difference = np.abs(read-haplotype)
  num_not_same = np.count_nonzero(difference)
  likelihood = (seq_error ** num_not_same) * ( (1 - seq_error) ** (len(read) - num_not_same) )
  return likelihood
  
'''
can be used when not all haplotypes are phased
'''
def mec_probabalistic(reconstructed_haplotypes, list_of_reads):
  pos_allele_pairs = [(list_of_reads[i], list_of_reads[i + 1]) for i in range(0, len(list_of_reads), 2)]

  expected_error = 0
  
  for idx, (position_seq, allele_seq) in enumerate(pos_allele_pairs):
    subarray = reconstructed_haplotypes[:, position_seq[0] : position_seq[-1]+1]
    for reconst_hap in subarray:
      likelihood = likelihood_haplotype(allele_seq, reconst_hap)
      error = np.sum(reconst_hap != allele_seq)
      expected_error += likelihood*error

  return expected_error

def mec_multiple_blocks(H_star, list_of_reads):
  '''
  H_star = reconstructed haplotypes, can be multiple blocks
  '''
  blocks = find_blocks(np.array(H_star, dtype=np.float64), np.array(H_star, dtype=np.float64))
  total_mec = 0
  for block in blocks:
    # print(block)
    block_H_star = H_star[:, block['start']:block['end']+1]
    if block['type'] == 1:
      total_mec += mec(block)
    elif block["type"] == 2 or block["type"] == 3:
      total_mec += mec_probabalistic(block, list_of_reads) 
    elif block["type"] ==4:
      continue
      # should i fill in the last row using the genotype?
  
def calculate_accuracy(reconstructed_haplotypes, true_haplptypes):
  """
  calculates how many SNPs line up perfectly (across all haplotypes)
  from the reconstructed haplotypes to the ground truth
  returns:
  the proportion of correct SNPs
  the permutation of the haplotypes that gives the best accuracy
  """

  n = reconstructed_haplotypes.shape[1]  # Number of SNPs
  k = reconstructed_haplotypes.shape[0]  # Number of haplotypes
  
  accuracy = 0
  
  for row_permutations in permutations(reconstructed_haplotypes):
    permuted_reconstructed_haplotypes = np.array(row_permutations)
    temp_accuracy = np.sum(np.all(true_haplptypes == permuted_reconstructed_haplotypes, axis=0))/n
    if temp_accuracy>accuracy:
      accuracy = temp_accuracy
      best_permutation = permuted_reconstructed_haplotypes
      
  return accuracy, best_permutation
  

def calculate_mismatch_error(reconstructed_haplotypes, true_haplotypes):
  """
  calculates the number of alleles that must be switched
  to make the reconstructed haplotypes line up perfectly
  with the ground truth
  returns:
  the numbers of alleles to be switched
  the permutation of haplotypes that gives the best mismatch error
  """
  k, n = reconstructed_haplotypes.shape
  mismatch_error = n*k
  
  mismatch_matrix = np.zeros((k,k))
  
  for reconstructed_row in range(k):
      for true_row in range(k):
          mismatch_matrix[reconstructed_row, true_row] = np.sum(
              reconstructed_haplotypes[reconstructed_row] != true_haplotypes[true_row]
          )
  # true indices give us the mapping from the rows 0,...k-1 in H to the rows sigma(0), ... in H*        
  recon_indices, true_indices = linear_sum_assignment(mismatch_matrix)
  mismatch_error = mismatch_matrix[recon_indices, true_indices].sum()
  return mismatch_error, true_indices


def calculate_fmpr(SNP_matrix, true_haplotypes):
  """
  fragment mapping phase relationship from hap compass
  takes all fragments and true haplotypes
  returns: counting all of the pairwise phase relationships defined 
  by the input set of fragments that do not exist in the solution
  """
  m, n = SNP_matrix.shape
  k = true_haplotypes.shape[0]
  
  fmpr_metric = 0
  
  for frag in SNP_matrix:
    # idx_j = snp location, f_ij = value read out of the fragment
    for idx_j, f_ij in enumerate(frag):
      for idx_k, f_ik in enumerate(frag):
        if idx_j != idx_k:
          current_min = np.inf
          for hap in true_haplotypes:
            temp_min = one_fmpr(f_ij, idx_j, f_ik, idx_k, hap)
            if temp_min<current_min:
              current_min = temp_min
          fmpr_metric += current_min
          
  return fmpr_metric


def one_fmpr(f_ij, j_idx, f_ik, k_idx, haplotype):
  """
  helper function for fragment mapping phase relationship
  np.nan might need to be replaced with
  however we represent no read for that location on that fragment
  """
  if ((f_ij != np.nan and f_ik != np.nan) and 
      (f_ij != haplotype[j_idx] or f_ik != haplotype[k_idx])):
    return 1
  return 0


def evaluate_ffbs_acc():
    results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    eval_df = pd.DataFrame(columns=['Contig', 'Ploidy', 'Coverage', 'Sample', 'FFBS_Acc'], index=range(2*3*5))
    counter = 0 
    for contig in [100]:
        for ploidy in [3, 4, 6, 8]:
            for cov in [10, 30, 50, 70, 100]:
                for sample in range(100):
                    sample_name = str(sample).zfill(2)
                    print('Contig:', contig, 'Ploidy:', ploidy, 'Coverage:', cov)
                    # ploidy = 4
                    sample_result = os.path.join(main_path, 'contig_{}/ploidy_{}/cov_{}/results_likelihood/FFBS_{}.pkl'.format(contig, ploidy, cov, sample_name)) 
                    with open(sample_result, 'rb') as f:
                        this_results = pickle.load(f)
                        f.close()
                    ffbs_acc = this_results['evaluation']['ffbs_acc']
                    eval_df.loc[counter, 'Contig'] = contig
                    eval_df.loc[counter, 'Ploidy'] = ploidy
                    eval_df.loc[counter, 'Coverage'] = cov
                    eval_df.loc[counter, 'Sample'] = sample_name
                    eval_df.loc[counter, 'FFBS_Acc'] = ffbs_acc
                    counter += 1
    # eval_groups = eval_df.groupby(['Contig', 'Ploidy', 'Coverage'])['FFBS_Acc'].mean().reset_index()
    eval_df.to_csv(os.path.join(results_path, 'ffbs_acc_NA12878.csv'), index=False)


def evaulate_ffbs_acc_sample(genotype_path, samples, ploidy):
    samples_brief = {}
    for t in samples.keys():
        for nn in samples[t].keys():
            if nn not in samples_brief.keys():
                samples_brief[nn] = samples[t][nn]
    # print(len(samples_brief))
    true_haplotypes = pd.read_csv(genotype_path).T
    sorted_nodes = sort_nodes(samples_brief.keys())
    eval_ffbs = {node: 0 for node in sorted_nodes}
    evals = []
    for node in sorted_nodes:
        positions = [int(i)-1 for i in node.split('-')]
        true_phasing = true_haplotypes.loc[:, positions].values
        true_phasing_permutations = np.array(list(itertools.permutations(true_phasing)))

        phasing = samples_brief[node]
        # print(phasing)
        phas_np = str_2_phas_1(phasing, ploidy)
        if not any(np.array_equal(phas_np, perm) for perm in true_phasing_permutations):
            evals.append(0)
            continue
        eval_ffbs[node] = 1
        evals.append(1)

    ffbs_acc = np.sum(evals)/len(evals)
    # print(ffbs_acc)
    return ffbs_acc


def calculate_pair_counts(fragment_list):
    """
    Parses the fragment list and creates a dictionary of dictionaries to count '00', '01', '10', and '11' 
    between every pair of positions across reads.
    
    Args:
        fragment_list: List of lists, where each pair of sublists represents positions and contents of a read.

    Returns:
        A dictionary of dictionaries with pair counts.
    """
    # pair_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    pair_counts = defaultdict(lambda: defaultdict(int))

    for i in range(0, len(fragment_list), 2):
        positions = fragment_list[i]
        contents = fragment_list[i + 1]

        for idx1 in range(len(positions)):
            for idx2 in range(idx1 + 1, len(positions)):
                pos1, pos2 = positions[idx1], positions[idx2]
                state1, state2 = contents[idx1], contents[idx2]

                # Sort positions to ensure consistency in key order
                # pos_pair = tuple(sorted((pos1, pos2)))
                pos_pair = '-'.join([str(s) for s in sorted((pos1, pos2))])
                state_pair = f"{state1}{state2}"

                pair_counts[pos_pair][state_pair] += 1

    return pair_counts

    
def get_block_info(quotient_g, predicted_haplotypes, true_haplotypes, fragment_model):
    component_labels, _ = gt.label_components(quotient_g.graph)
    components = {}
    for v in quotient_g.graph.vertices():
        comp_id = component_labels[v]  # Get component ID of the vertex
        if comp_id not in components:
            components[comp_id] = {'blocks': []}
        components[comp_id]['blocks'].append(quotient_g.graph.vertex_properties["v_label"][v])

    for key in components.keys():
        block = components[key]['blocks']
        positions = sorted(set(int(num) for r in block for num in r.split('-')))
        positions_ind = [p-1 for p in positions]
        block_pred_haplotype = predicted_haplotypes[positions_ind].to_numpy()
        block_true_haplotype = true_haplotypes[positions_ind].to_numpy()
        block_vector_error_rate, block_vector_error, _, _ = compute_vector_error_rate(block_pred_haplotype, block_true_haplotype)
        block_accuracy, _ = calculate_accuracy(block_pred_haplotype, block_true_haplotype)
        block_mismatch_error, _ = calculate_mismatch_error(block_pred_haplotype, block_true_haplotype)
        block_mec_ = mec(block_pred_haplotype, fragment_model.fragment_list)
        components[key]['evaluation'] = {'vector_error_rate': block_vector_error_rate, 'vector_error': block_vector_error, 
                                         'accuracy': block_accuracy, 'mismatch_error': block_mismatch_error, 'mec': block_mec_}
        components[key]['block_size'] = len(positions_ind)

    block_info = {'vector_error_rate': np.mean([components[key]['evaluation']['vector_error_rate'] for key in components.keys()]), 
                  'vector_error': np.mean([components[key]['evaluation']['vector_error'] for key in components.keys()]),
                  'accuracy': np.mean([components[key]['evaluation']['accuracy'] for key in components.keys()]),
                  'mismatch_error': np.mean([components[key]['evaluation']['mismatch_error'] for key in components.keys()]),
                  'mec': np.mean([components[key]['evaluation']['mec'] for key in components.keys()]), 
                  'average_block_size': np.mean([components[key]['block_size'] for key in components.keys()])}

    return block_info, components


def find_blocks(H, H_star):
    """
    Splits the columns of H (possibly containing NaNs) into consecutive blocks,
    based on how many rows (out of k) are non-NaN in each column.

    We define four block types:
      Type 1: completely filled (m == k)
      Type 2: completely empty (m == 0)
      Type 3: partially filled with m < k - 1
      Type 4: partially filled with m == k - 1

    Inputs:
      H      : (k x n) numpy array, can have NaNs
      H_star : (k x n) numpy array, assumed fully filled (though here
               we mainly just validate shape, not used in block splitting)

    Returns:
      blocks : a list of dictionaries, each having:
         {
           'start': start_col_index (0-based),
           'end'  : end_col_index (0-based),
           'type' : block_type (1, 2, 3, or 4),
           'm'    : number of non-NaN rows in these columns
         }
        The columns 'start' through 'end' (inclusive) all have the same m-value
        and thus fall into the same type of block.
    """

    # --- Validate shapes ---
    k, n = H.shape
    if H_star.shape != (k, n):
        raise ValueError("H and H_star must have the same shape.")

    # --- Helper to determine block type by 'm' ---
    def get_block_type(m, k):
        if m == 0:
            return 2  # completely empty
        elif m == k:
            return 1  # completely filled
        elif m == k - 1:
            return 4  # partially filled, m = k - 1
        else:
            return 3  # partially filled, m < k - 1

    # --- If there are no columns, return empty list ---
    if n == 0:
        return []

    # --- Identify blocks ---
    blocks = []

    # Count how many rows are non-NaN for the first column
    prev_m = np.sum(~np.isnan(H[:, 0]))
    start_col = 0

    # Iterate over columns 1..n-1
    for col in range(1, n):
        current_m = np.sum(~np.isnan(H[:, col]))
        # If current_m changes from prev_m, we finalize the previous block
        # not sure if this works, 
        # what if the unphased haplotype is a different one from prev to current?
        if current_m != prev_m:
            blocks.append({
                'start': start_col,
                'end'  : col - 1,
                'type' : get_block_type(prev_m, k),
                'm'    : prev_m
            })
            # Start a new block
            start_col = col
            prev_m = current_m

    # Close the last block (from start_col .. n-1)
    blocks.append({
        'start': start_col,
        'end'  : n - 1,
        'type' : get_block_type(prev_m, k),
        'm'    : prev_m
    })

    return blocks


def generate_example_H_Hstar():
    """
    Generates an example pair of H and H_star matrices for testing.

    H_star: Fully filled haplotype matrix (no NaNs).
    H     : Partially filled haplotype matrix with NaNs.
    
    Returns:
        H (numpy array): Partially filled haplotype matrix
        H_star (numpy array): Fully filled haplotype matrix (ground truth)
    """
    # Set random seed for reproducibility
    np.random.seed(42)

    # Define k (rows) and n (columns)
    k = 4  # Number of rows (haplotypes)
    n = 16  # Number of columns (positions)

    # Generate a fully filled matrix H_star with random values (0/1)
    H_star = np.random.randint(0, 2, size=(k, n))

    # Create a copy of H_star to modify for H (introducing NaNs)
    H = H_star.astype(float)  # Convert to float to allow NaNs

    # Introduce NaNs to match the given example blocks:
    # Block 1: columns 0-2 (fully filled, keep as is)
    # Block 2: columns 3-4 (completely empty)
    H[:, 3:5] = np.nan
    
    # Block 3: columns 5-6 (partially filled, 2 non-NaN per column)
    H[:, 5] = [np.nan, np.nan, 1, 0]  # Only last 2 rows filled
    H[:, 6] = [np.nan, np.nan, 0, 1]  # Only last 2 rows filled
    
    # Block 4: columns 7-8 (completely empty)
    H[:, 7:9] = np.nan
    
    # Block 5: columns 9-10 (fully filled, keep as is)
    
    # Block 6: columns 11-12 (partially filled, 2 non-NaN per column)
    H[:, 11] = [np.nan, np.nan, 1, 0]
    H[:, 12] = [np.nan, np.nan, 0, 1]
    
    # Block 7: column 13 (partially filled, only 1 row filled)
    H[:, 13] = [np.nan, np.nan, np.nan, 1]
    
    # Block 8: column 14 (partially filled, 3 out of 4 rows filled)
    H[:, 14] = [1, 0, 1, np.nan]
    H[:, 15] = [0, 0, 1, np.nan]
    
    return H, H_star


def compute_vector_error_rate_with_missing_positions(H_star, H):
    """
    Computes the vector error rate between two haplotype matrices, allowing for missing positions.
    H_star is assumed to be fully filled, while H may contain NaNs.
    H_star: True haplotype
    H: Assembled haplotype with NaNs for missing positions
    """
    # H, H_star = generate_example_H_Hstar()
    blocks = find_blocks(np.array(H, dtype=np.float64), np.array(H_star, dtype=np.float64))
    vector_error = 0
    for block in blocks:
        # print(block)
        block_H = H[:, block['start']:block['end']+1]
        block_H_star = H_star[:, block['start']:block['end']+1]
        if block['type'] == 1:
          this_vector_error = compute_vector_error_rate(block_H_star, block_H)[1]
          vector_error += this_vector_error
        elif block['type'] == 2:
          for col in range(block['start'], block['end']+1):
            H_star_col = H_star[:, col]
            m0 = np.count_nonzero(H_star_col == 0)
            m1 = np.count_nonzero(H_star_col == 1)
            s0 = int(m0 > 1)  # 1 if m0 > 1, else 0
            s1 = int(m1 > 1)  # 1 if m1 > 1, else 0
            vector_error += s0 + s1
        elif block['type'] == 3:
          # completed_rows = block['m']
          block_H = H[:, block['start']:block['end']+1]
          block_H_star = H_star[:, block['start']:block['end']+1]
          this_vector_error1 = compute_vector_error_rate_partial(block_H_star, block_H)[1]
          vector_error += this_vector_error1
          for col in range(block['start'], block['end']+1):
            H_star_col = H_star[:, col]
            H_col = H[:, col]
            H_star_0_count = np.count_nonzero(H_star_col == 0)
            H_star_1_count = np.count_nonzero(H_star_col == 1)
            H_0_count = np.count_nonzero(H_col == 0)
            H_1_count = np.count_nonzero(H_col == 1)
            m0 = H_star_0_count - H_0_count
            m1 = H_star_1_count - H_1_count
            s0 = int(m0 > 1)
            s1 = int(m1 > 1)
            vector_error += s0 + s1
        elif block['type'] == 4:
          block_H = H[:, block['start']:block['end']+1]
          block_H_star = H_star[:, block['start']:block['end']+1]
          for col_id, col in enumerate(range(block['start'], block['end']+1)):
            # print(col_id, col)
            H_star_col = H_star[:, col]
            H_col = H[:, col]
            unphased_val = np.sum(H_star_col) - np.nansum(H_col)
            if unphased_val > 1:
              unphased_val = 1
              print('Unphased value is greater than 1')
            nan_index = np.where(np.isnan(H_col))[0]
            block_H[nan_index, col_id] = unphased_val
          this_vector_error = compute_vector_error_rate(block_H_star, block_H)[1]
          vector_error += this_vector_error
        # print(block, vector_error)
    vector_error_rate = vector_error/H.shape[1]
    return vector_error_rate, vector_error, blocks

