
from itertools import permutations
import numpy as np

'''
input:
h_col : list or array, one allele position across all haplotypes in true phase
h_star_col : list or array, one allele position across all haplotypes in assembled phase
returns:
list of tuples, each tuple is a mapping of the indices 0, ..., k-1 in h_star_col to the indices 0, ..., k-1 in h_col
'''
def find_matches(h_star_col, h_col):

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


'''
input:
H and H_star : numpy arrays of the true and assembled phasings
retuns:
the value of the vector error
the backtracking steps - how we have to permute each column of H to make them match H_star
'''
def compute_vector_error_rate(H, H_star):
  k, n = H.shape

  # i^th dictionary stores the ways to get through the i^th col
  # dictionary stores:
  # key = mapping from H to H star at that column
  # value = (vector error using that key, mapping form prev column that induces that vector error)

  dp_table = [{} for _ in range(n)]

  # vector error is 0 when we only consider the first column
  first_col_matchings = find_matches(H_star[:, 0], H[:, 0])
  for first_col_matching in first_col_matchings:
    dp_table[0][first_col_matching] = (0, None)

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



'''
input:
H_star - numpy array assembled phasing
list of reads - list of tuples
  each tuple is of form (list representing the read, integer representing the starting index)
returns:
the proportion of reads that could be mapped to some row of H_star at their specified position
'''
def mec(H_star, list_of_reads):

    map_back = np.zeros(len(list_of_reads))

    for idx, (read_seq, start_pos) in enumerate(list_of_reads):
        for row in H_star:
            subarray = row[start_pos : start_pos + len(read_seq)]

            if np.array_equal(subarray, read_seq):
                map_back[idx] = 1
                # once we find a match for this read, we can stop checking other rows
                break
    return np.sum(map_back)/len(list_of_reads)
  
  
'''
calculates how many SNPs line up perfectly (across all haplotypes)
from the reconstructed haplotypes to the ground truth
returns:
the proportion of correct SNPs
the permutation of the haplotypes that gives the best accuracy
'''
def calculate_accuracy(reconstructed_haplotypes, true_haplptypes):
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
  
  
'''
calculates the number of alleles that must be switched
to make the reconstructed haplotypes line up perfectly
with the ground truth
returns:
the numbers of alleles to be switched
the permutation of haplotypes that gives the best mismatch error
'''  
def calculate_mismatch_error(reconstructed_haplotypes, true_haplotypes):
  n = reconstructed_haplotypes.shape[1]  # Number of SNPs
  k = reconstructed_haplotypes.shape[0]  # Number of haplotypes
    
  mismatch_error = n*k
  
  for row_permutations in permutations(reconstructed_haplotypes):
    permuted_reconstructed_haplotypes = np.array(row_permutations)
    temp_mismatch_error = np.sum(permuted_reconstructed_haplotypes != true_haplotypes)
    if temp_mismatch_error<mismatch_error:
      mismatch_error = temp_mismatch_error
      best_permutation = permuted_reconstructed_haplotypes
      
  return mismatch_error, best_permutation
  

'''
fragment mapping phase relationship from hap compass
takes all fragments and true haplotypes
returns: counting all of the pairwise phase relationships defined 
by the input set of fragments that do not exist in the solution
'''
def calculate_fmpr(SNP_matrix, true_haplotypes):
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


'''
helper function for fragment mapping phase relationship
np.nan might need to be replaced with
however we represent no read for that location on that fragment
'''
def one_fmpr(f_ij, j_idx, f_ik, k_idx, haplotype):
  
  if ((f_ij != np.nan and f_ik != np.nan) and 
      (f_ij != haplotype[j_idx] or f_ik != haplotype[k_idx])):
    return 1
  return 0

# correct phasing rate = 
# num alleles that are in the same place 
# in the true & reconstructed

# perfect solution rate = num haplotypes 
# w no switches (identical in true & reconst)

########
# examples
########

'''
h = [0, 0, 1]
h_star = [0, 1, 0]
find_matches(h, h_star)

H = np.array([
    [1,1,1,0,0,0,1],
    [1,0,1,0,0,1,1],
    [0,0,0,1,1,0,0]
])

H_star = np.array([
    [1,1,1,1,1,0,1],
    [1,0,1,0,0,0,1],
    [0,0,0,0,0,1,0]
])

compute_vector_error(H, H_star)

list_of_reads = [([1,1,1], 0), ([1,1], 2 ), ([0,0],5)]

mec(H_star, list_of_reads)
'''
