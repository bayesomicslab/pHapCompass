from itertools import permutations
import numpy as np
import pandas as pd
import pickle
import os
import itertools
from itertools import combinations, permutations
from scipy.optimize import linear_sum_assignment
from collections import defaultdict
from utils import *


########################################################################################################
# BLOCK PARSERS, ETC
########################################################################################################

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


def find_blocks_from_ids(H_star, block_ids):
    """
    Split columns of H_star (reconstructed phase) into blocks using the provided block_ids vector.

    Inputs:
      H_star    : (k x n) numpy array; may contain NaNs
      block_ids : (n,) array-like; entries equal within a block,
                  change at block boundaries; NaN runs are separate blocks.

    Returns:
      blocks : list of dicts with keys:
        {
          'start': start_col_index (0-based),
          'end'  : end_col_index (0-based),
          'type' : 1|2|3|4,
          'm'    : representative m for the block
        }
    """
    k, n = H_star.shape
    block_ids = np.asarray(block_ids)
    if block_ids.shape[0] != n:
        raise ValueError("block_ids must have length equal to number of columns in H.")
    if n == 0:
        return []

    def _isnan(x):
        try:
            return np.isnan(x)
        except TypeError:
            return False

    def _id_equal(a, b):
        # treat NaN == NaN as equal
        if _isnan(a) and _isnan(b):
            return True
        return a == b

    not_nan = ~np.isnan(H_star)
    blocks = []

    def _finalize_block(i0, i1):
        m_cols = np.sum(not_nan[:, i0:i1+1], axis=0)
        if np.all(m_cols == k):
            btype, m_val = 1, k
        elif np.all(m_cols == 0):
            btype, m_val = 2, 0
        elif np.all(m_cols == (k - 1)):
            btype, m_val = 4, k - 1
        else:
            vals, counts = np.unique(m_cols, return_counts=True)
            m_val = int(vals[np.argmax(counts)])
            btype = 3
        blocks.append({'start': i0, 'end': i1, 'type': btype, 'm': int(m_val)})

    start = 0
    current_id = block_ids[0]
    for col in range(1, n):
        if not _id_equal(block_ids[col], current_id):
            _finalize_block(start, col - 1)
            start = col
            current_id = block_ids[col]
    _finalize_block(start, n - 1)

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


########################################################################################################
# VECTOR ERROR 
########################################################################################################

# VECTOR ERROR HELPERS

def find_matches(h_col, h_star_col):
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

def find_partial_matches(h_col, h_star_col):
    """
    returns all of the ways to match the alleles
    in a partial column of H_star of length tilde{K}
    to a column of H of length K
    where 1 <= tilde{K} <= K-2
    """
    k_tilde = len(h_star_col)
    
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
                    mapping = [None] * k_tilde
                    for idx, gt_idx in zip(one_indices_star, perm_ones):
                        mapping[idx] = gt_idx
                    for idx, gt_idx in zip(zero_indices_star, perm_zeros):
                        mapping[idx] = gt_idx
                    matches.append(tuple(mapping))
    return matches

# VECTOR ERROR FOR EACH BLOCK TYPES

def vector_error_type_1(H_block, H_star_block):
    """
    H_block = true
    H_star_block = reconstructed
    H star is fully phased (K haplotypes)
    
    returns vector error, not vector error rate
    no first position penalty
    """
  
    K, L = H_block.shape
    
    # Create a copy once to avoid modifying the original
    H_star_block_copy = H_star_block.copy()
    
    # Counter for genotype disagreements
    genotype_disagreement_count = 0

    # i^th dictionary stores the ways to get through the i^th col
    # dictionary stores:
    # key = mapping from H to H star at that column
    # value = (vector error using that key, mapping form prev column that induces that vector error)

    dp_table = [{} for _ in range(L)]

    # Handle first column
    try:
        first_col_matchings = find_matches(H_block[:, 0], H_star_block_copy[:, 0])
    except IndexError:
        # If we're here, genotypes definitely don't match
        h_col = H_block[:, 0]
        h_star_col = H_star_block_copy[:, 0]
        
        genotype_h = np.sum(h_col)
        genotype_h_star = np.sum(h_star_col)
        
        # No need for if statement - we know they differ
        if genotype_h_star > genotype_h:
            num_flips = int(genotype_h_star - genotype_h)
            one_indices = [i for i, val in enumerate(h_star_col) if val == 1]
            flip_indices = np.random.choice(one_indices, size=num_flips, replace=False)
            H_star_block_copy[flip_indices, 0] = 0
        else:
            num_flips = int(genotype_h - genotype_h_star)
            zero_indices = [i for i, val in enumerate(h_star_col) if val == 0]
            flip_indices = np.random.choice(zero_indices, size=num_flips, replace=False)
            H_star_block_copy[flip_indices, 0] = 1
        
        genotype_disagreement_count += num_flips
        first_col_matchings = find_matches(H_block[:, 0], H_star_block_copy[:, 0])    
    for first_col_matching in first_col_matchings:
        dp_table[0][first_col_matching] = (0, None)

    # Forward pass
    for col in range(1, L):
        try:
            all_matches = find_matches(H_block[:, col], H_star_block_copy[:, col])
        except IndexError:
            # Handle genotype disagreement
            h_col = H_block[:, col]
            h_star_col = H_star_block_copy[:, col]
            
            genotype_h = np.sum(h_col)
            genotype_h_star = np.sum(h_star_col)
            
            if genotype_h_star > genotype_h:
                # Too many 1s in h_star_col, need to flip some 1s to 0s
                num_flips = int(genotype_h_star - genotype_h)
                one_indices = [i for i, val in enumerate(h_star_col) if val == 1]
                flip_indices = np.random.choice(one_indices, size=num_flips, replace=False)
                H_star_block_copy[flip_indices, col] = 0
            else:
                # Too many 0s in h_star_col, need to flip some 0s to 1s
                num_flips = int(genotype_h - genotype_h_star)
                zero_indices = [i for i, val in enumerate(h_star_col) if val == 0]
                flip_indices = np.random.choice(zero_indices, size=num_flips, replace=False)
                H_star_block_copy[flip_indices, col] = 1
            
            genotype_disagreement_count += num_flips
            all_matches = find_matches(H_block[:, col], H_star_block_copy[:, col])
        
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
    vector_error, _ = min(dp_table[L-1].values(), key=lambda x: x[0])

    return int(vector_error + genotype_disagreement_count)

def vector_error_type_2(H_block, H_star_block):
    '''
    H star fully empty, all nan
    applies expected number of vector errors to all cols EXCEPT FIRST
    vectorized over columns in the block
    '''
    
    K, L = H_block.shape
    
    if L <= 1:
        return 0
    
    # Get columns 1 onwards (second column onwards)
    H_subsequent_cols = H_block[:, 1:]  # Shape: (K, L-1)
    
    # Count 0s and 1s in each column
    m0 = np.sum(H_subsequent_cols == 0, axis=0)  # Shape: (L-1,)
    m1 = np.sum(H_subsequent_cols == 1, axis=0)  # Shape: (L-1,)
    
    # Compute K - I{m0>1} - I{m1>1} for each column
    vector_errors = K - (m0 > 1).astype(int) - (m1 > 1).astype(int)
    
    # Sum across all columns
    return int(np.sum(vector_errors))

def vector_error_type_3(H_block, H_star_block):
    '''
    H star partially phased, 
    tilde{K} haplotypes phased where 1 <= tilde{K} <= K-2
    
    returns sum of actual vector error from phased part 
    with expected vector error from unphased part
    (except first column)
    handles cases where H_star is not extendable to H
    '''
    
    K, L = H_block.shape
    H_star_block_phased = H_star_block[~np.isnan(H_star_block).any(axis=1),:]
    
    # Counter for genotype disagreements
    genotype_disagreement_count = 0
    
    # Create a copy to modify if needed
    H_star_block_phased_corrected = H_star_block_phased.copy()
    
    dp_table = [{} for _ in range(L)]
    
    # Handle first column
    first_col_matchings = find_partial_matches(H_block[:, 0], H_star_block_phased_corrected[:, 0])
    if len(first_col_matchings) == 0:
        # Not extendable, need to fix
        h_col = H_block[:, 0]
        h_star_col = H_star_block_phased_corrected[:, 0]
        
        # Count 1s and 0s
        ones_h = np.sum(h_col == 1)
        zeros_h = np.sum(h_col == 0)
        ones_h_star = np.sum(h_star_col == 1)
        zeros_h_star = np.sum(h_star_col == 0)
        
        # Determine what needs to be flipped
        if ones_h_star > ones_h:
            # Too many 1s in h_star, flip some 1s to 0s
            num_flips = ones_h_star - ones_h
            one_indices = np.where(h_star_col == 1)[0]
            flip_indices = np.random.choice(one_indices, size=num_flips, replace=False)
            H_star_block_phased_corrected[flip_indices, 0] = 0
            genotype_disagreement_count += num_flips
        
        if zeros_h_star > zeros_h:
            # Too many 0s in h_star, flip some 0s to 1s
            num_flips = zeros_h_star - zeros_h
            zero_indices = np.where(h_star_col == 0)[0]
            flip_indices = np.random.choice(zero_indices, size=num_flips, replace=False)
            H_star_block_phased_corrected[flip_indices, 0] = 1
            genotype_disagreement_count += num_flips

        print(H_star_block_phased_corrected)
        
        # Now try again
        first_col_matchings = find_partial_matches(H_block[:, 0], H_star_block_phased_corrected[:, 0])
    
    for matching in first_col_matchings:
        dp_table[0][matching] = (0, None)
        
    # Forward pass: iterate over columns.
    for col in range(1, L):
        all_matches = find_partial_matches(H_block[:, col], H_star_block_phased_corrected[:, col])
        
        if len(all_matches) == 0:
            # Not extendable, need to fix
            h_col = H_block[:, col]
            h_star_col = H_star_block_phased_corrected[:, col]
            
            # Count 1s and 0s
            ones_h = np.sum(h_col == 1)
            zeros_h = np.sum(h_col == 0)
            ones_h_star = np.sum(h_star_col == 1)
            zeros_h_star = np.sum(h_star_col == 0)
            
            # Determine what needs to be flipped
            if ones_h_star > ones_h:
                # Too many 1s in h_star, flip some 1s to 0s
                num_flips = ones_h_star - ones_h
                one_indices = np.where(h_star_col == 1)[0]
                flip_indices = np.random.choice(one_indices, size=num_flips, replace=False)
                H_star_block_phased_corrected[flip_indices, col] = 0
                genotype_disagreement_count += num_flips
            
            if zeros_h_star > zeros_h:
                # Too many 0s in h_star, flip some 0s to 1s
                num_flips = zeros_h_star - zeros_h
                zero_indices = np.where(h_star_col == 0)[0]
                flip_indices = np.random.choice(zero_indices, size=num_flips, replace=False)
                H_star_block_phased_corrected[flip_indices, col] = 1
                genotype_disagreement_count += num_flips
            
            # Now try again
            all_matches = find_partial_matches(H_block[:, col], H_star_block_phased_corrected[:, col])
        
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
    
    # best final alignment
    actual_vector_error, _ = min(dp_table[L-1].values(), key=lambda x: x[0])
    
    # to compute expected vector error from unphased part, 
    # consider "remaining" genotype implied by H and H star
    
    # Count 0s in each column of H_block
    H_zeros = np.sum(H_block == 0, axis=0)  # Shape: (L,)
    
    # Count 1s in each column of H_block
    H_ones = np.sum(H_block == 1, axis=0)   # Shape: (L,)
    
    # Use the corrected H_star_block_phased for counting
    H_star_phased_zeros = np.sum(H_star_block_phased_corrected == 0, axis=0)  # Shape: (L,)
    H_star_phased_ones = np.sum(H_star_block_phased_corrected == 1, axis=0)   # Shape: (L,)

    # Compute differences (exclude first column)
    m_0 = H_zeros[1:] - H_star_phased_zeros[1:]  # Shape: (L-1,)
    m_1 = H_ones[1:] - H_star_phased_ones[1:]    # Shape: (L-1,)
    # Compute K - I{m0>1} - I{m1>1} for each column
    expected_vector_error = np.sum(K - (m_0 > 1).astype(int) - (m_1 > 1).astype(int))
    
    return int(actual_vector_error + expected_vector_error)
  
  
def vector_error_type_4(H_block, H_star_block):
    """
    H star mostly phased
    just missing one row which can be added to meet genotype
    we assume that the genotype implied by each column of H block is
    either exactly the reported genotype (so we augment with a 0)
    or one less (so we augment with a 1)
    handles genotype disagreements by choosing closest value (0 or 1)
    and counting the disagreement
    """
    
    K, L = H_block.shape
    
    # Create a copy to avoid modifying the input
    H_star_block_augmented = H_star_block.copy()
    
    # Find the NaN row (should be the same for all columns in type 4, so take first)
    nan_mask = np.isnan(H_star_block)  # Shape: (K, L)
    nan_row_idx = np.where(np.any(nan_mask, axis=1))[0][0]
    
    # Compute column-wise sums
    H_star_sums = np.nansum(H_star_block, axis=0)  # Shape: (L,) - sum ignoring NaNs
    H_sums = np.sum(H_block, axis=0)  # Shape: (L,) 
    
    # Compute the difference for each column
    unphased_vals = H_sums - H_star_sums  # Shape: (L,)
    
    # instead of throwing an error if not extendable to genotype
    # just use the closest and count disagreement
    corrected_vals = np.where(unphased_vals >= 1, 1, 0)
    
    # Vectorized disagreement counting
    genotype_disagreement_count = (
        np.sum((unphased_vals - 1) * (unphased_vals >= 2)) +
        np.sum(np.abs(unphased_vals) * (unphased_vals <= -1))
    )
    
    # Fill the NaN row with the corrected values
    H_star_block_augmented[nan_row_idx, :] = corrected_vals
    
    return int(genotype_disagreement_count + vector_error_type_1(H_block, H_star_block_augmented))   

# VECTOR ERROR OVER ENTIRE PHASE
    
def vector_error_wrapper(H, H_star, block_ids):
    '''
    the vector error function we will actually call for evaluation
    vector error for each block together with block cutoff penalties
    '''
    K, L = H.shape
    blocks = find_blocks_from_ids(H_star, block_ids)
    total_vector_error = 0
    
    block_type_to_ve_function = {
        1: vector_error_type_1,
        2: vector_error_type_2,
        3: vector_error_type_3,
        4: vector_error_type_4
    }
    
    # First block (no boundary penalty)
    block_H = H[:, blocks[0]["start"]:blocks[0]["end"]+1]
    block_H_star = H_star[:, blocks[0]["start"]:blocks[0]["end"]+1]
    total_vector_error += block_type_to_ve_function[blocks[0]["type"]](block_H, block_H_star)
    
    # Subsequent blocks (add boundary penalty for first column of each block)
    for block in blocks[1:]:
        block_H = H[:, block["start"]:block["end"]+1]
        block_H_star = H_star[:, block["start"]:block["end"]+1]
        
        # Add vector error within the block
        total_vector_error += block_type_to_ve_function[block["type"]](block_H, block_H_star)
        
        # Add penalty for crossing block boundary 
        # based on genotype of first column of new block
        first_col_H = block_H[:, 0]
        m0 = np.count_nonzero(first_col_H == 0)
        m1 = np.count_nonzero(first_col_H == 1)
        # total_vector_error += K - (m0 > 1).astype(int) - (m1 > 1).astype(int)
        total_vector_error += K - int(m0 > 1) - int(m1 > 1)
    
    return total_vector_error / L

########################################################################################################
# MEC
########################################################################################################

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
helpers for mec, to be used
when the read and the reconstructed array fit together perfectly 
(same num of cols)
'''

def mec_helper(read, reconstruction):
  best_match_haplotype = reconstruction[np.argmax(np.sum(reconstruction==read, axis=1))]
  error_correction_count = np.sum(best_match_haplotype != read)
  return error_correction_count/len(read)
  
def mec_probabalistic_helper(read, reconstruction):
  expected_error = 0
  standardizing_constant = sum((likelihood_haplotype(read, reconst_hap)) for reconst_hap in reconstruction) 
  for reconst_hap in reconstruction:
    likelihood = likelihood_haplotype(read, reconst_hap) / standardizing_constant
    error = np.sum(reconst_hap != read)
    expected_error += likelihood*error
  return expected_error


def mec_full(H_star, block_ids, list_of_reads, probabalistic=False):
  '''
  given:
  full reconstruction (one big array, not split into blocks) ,
  list of reads (list of position/sequence pairs)
  whether or not you want to use the probabalistic version of error correction
  '''
  mec_running_total = 0
  k, n = H_star.shape
  pos_allele_pairs = [(list_of_reads[i], list_of_reads[i + 1]) for i in range(0, len(list_of_reads), 2)]
  
  # Get block definitions using new function
  block_defns = find_blocks_from_ids(H_star, block_ids)
  
  # Extract block slices from H_star
  blocks = [H_star[:, b['start']:b['end'] + 1] for b in block_defns]
  
  # Precompute block widths and start indices
  widths = [b['end'] - b['start'] + 1 for b in block_defns]
  starts = [b['start'] for b in block_defns]

  # make sure that sum(widths)==n
  # col2block maps cols to blocks
  # col2offset maps cols to index within the respective block
  col2block = np.zeros(n, dtype=int)
  col2offset = np.zeros(n, dtype=int)
  for i, (start, w) in enumerate(zip(starts, widths)):
      col2block[start:start + w] = i
      col2offset[start:start + w] = np.arange(w)
  # one read at a time, calculate mec
  for idx, (position_seq, allele_seq) in enumerate(pos_allele_pairs):
    blocks_covering_read = col2block[position_seq]
    indices_within_covering_blocks = col2offset[position_seq]
    uniq, first = np.unique(blocks_covering_read, return_index=True)
    ordered_blocks = uniq[np.argsort(first)]
    # directly gather exact columns within each block (no contiguity assumption)
    parts = [
        blocks[b][:, indices_within_covering_blocks[blocks_covering_read == b]]
        for b in ordered_blocks
    ]
    for (idx2,b) in enumerate(ordered_blocks):
      mask = blocks_covering_read == b
      idxs = np.where(mask)[0]
      read_slice = allele_seq[idxs[0] : idxs[-1] + 1]
      if probabalistic:
        mec_running_total += mec_probabalistic_helper(read_slice, parts[idx2])
      else: 
        mec_running_total += mec_helper(read_slice, parts[idx2])
    mec_running_total += (len(ordered_blocks)-1)*(k-1)/k
    return mec_running_total


def mec_full_geometric_penalty(H_star, block_ids, list_of_reads, probabalistic=False):
  '''
  given:
  full reconstruction (one big array, not split into blocks) ,
  list of reads (list of position/sequence pairs)
  whether or not you want to use the probabalistic version of error correction
  '''
  mec_running_total = 0
  k, n = H_star.shape
  pos_allele_pairs = [(list_of_reads[i], list_of_reads[i + 1]) for i in range(0, len(list_of_reads), 2)]
  
  # Get block definitions using new function
  block_defns = find_blocks_from_ids(H_star, block_ids)
  
  # Extract block slices from H_star
  blocks = [H_star[:, b['start']:b['end'] + 1] for b in block_defns]
  
  # Precompute block widths and start indices
  widths = [b['end'] - b['start'] + 1 for b in block_defns]
  starts = [b['start'] for b in block_defns]

  # make sure that sum(widths)==n
  # col2block maps cols to blocks
  # col2offset maps cols to index within the respective block
  col2block = np.zeros(n, dtype=int)
  col2offset = np.zeros(n, dtype=int)
  for i, (start, w) in enumerate(zip(starts, widths)):
      col2block[start:start + w] = i
      col2offset[start:start + w] = np.arange(w)
  # one read at a time, calculate mec
  for idx, (position_seq, allele_seq) in enumerate(pos_allele_pairs):
    blocks_covering_read = col2block[position_seq]
    indices_within_covering_blocks = col2offset[position_seq]
    uniq, first = np.unique(blocks_covering_read, return_index=True)
    ordered_blocks = uniq[np.argsort(first)]
    # directly gather exact columns within each block (no contiguity assumption)
    parts = [
        blocks[b][:, indices_within_covering_blocks[blocks_covering_read == b]]
        for b in ordered_blocks
    ]
    for (idx2,b) in enumerate(ordered_blocks):
      mask = blocks_covering_read == b
      idxs = np.where(mask)[0]
      read_slice = allele_seq[idxs[0] : idxs[-1] + 1]
      if probabalistic:
        mec_running_total += mec_probabalistic_helper(read_slice, parts[idx2])
      else: 
        mec_running_total += mec_helper(read_slice, parts[idx2])
    mec_running_total += 1 - (1/k)**(len(ordered_blocks)-1)
    return mec_running_total


########################################################################################################
# FFBS evaluations
########################################################################################################

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
  
  
 ########################################################################################################
 # UNCERTAINTY QUANTIFICATION
 ########################################################################################################
 
def sort_matrix_rows(matrix):
    """Sort matrix rows lexicographically to create canonical form."""
    sorted_indices = np.lexsort(matrix.T[::-1])
    return matrix[sorted_indices]


def optimal_hamming_distance_across_samples(u_ref, v_samples, cache):
    """
    Compute optimal mean Hamming distance between reference matrix and multiple sample matrices
    using a single optimal permutation that minimizes total distance across all samples.
    Uses cache for individual row pair distances.
    
    Args:
        u_ref: Reference binary matrix of shape (K, 2) - UNSORTED
        v_samples: List of m sample binary matrices, each of shape (K, 2) - UNSORTED
        cache: Dict mapping (sorted_row_u, sorted_row_v) -> distance
        
    Returns:
        total_distance: Total Hamming distance under optimal row matching
    """
    K = u_ref.shape[0]
    
    # Build cost matrix - sum of Hamming distances across all samples
    cost = np.zeros((K, K))
    
    for i in range(K):
        for j in range(K):
            # Sum Hamming distance for this (i,j) pair across all samples
            for v_sample in v_samples:
                # Create cache key from sorted rows
                row_u = u_ref[i]
                row_v = v_sample[j]
                
                # Sort the rows for cache key
                row_u_sorted = tuple(sorted(row_u))
                row_v_sorted = tuple(sorted(row_v))
                
                cache_key = (row_u_sorted, row_v_sorted)
                reverse_key = (row_v_sorted, row_u_sorted)
                
                # Check cache
                if cache_key in cache:
                    distance = cache[cache_key]
                elif reverse_key in cache:
                    distance = cache[reverse_key]
                else:
                    # Compute and cache
                    distance = np.sum(row_u != row_v)
                    cache[cache_key] = distance
                    cache[reverse_key] = distance
                
                cost[i, j] += distance
    
    # Solve assignment problem once for all samples
    row_ind, col_ind = linear_sum_assignment(cost)
    
    # Return total distance
    total_distance = cost[row_ind, col_ind].sum()
    
    return total_distance


def pair_uncertainty_quantification(H, samples):
    """
    Compute uncertainty quantification with caching of individual row pairs.
    Finds a single optimal permutation per column pair that minimizes
    total Hamming distance across all samples.
    
    Args:
        H: Reference matrix (K, L) - binary values
        samples: List of m sample matrices, each (K, L) - binary values
        
    Returns:
        results: Dict mapping (col1, col2) -> total distance across all samples
        cache: Dict mapping (sorted_row_u, sorted_row_v) -> distance
    """
    K, L = H.shape
    m = len(samples)
    
    # Validate inputs
    for i, sample in enumerate(samples):
        assert sample.shape == H.shape, f"Sample {i} shape mismatch"
    
    cache = {}
    results = {}
    
    # Process each column pair
    for h1, h2 in combinations(range(L), 2):
        if abs(h1 - h2) > 150:
          continue
        # Extract reference column pair (unsorted for computation)
        u_ref = H[:, [h1, h2]]
        
        # Extract all sample column pairs (unsorted for computation)
        v_samples = [sample[:, [h1, h2]] for sample in samples]
        
        # Compute optimal distance using unsorted matrices and cache
        total_distance = optimal_hamming_distance_across_samples(u_ref, v_samples, cache)
        
        results[(h1, h2)] = total_distance / m
    
    return results, cache
  

########################################################################################################
# miscellaneous old evaluations
########################################################################################################


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
################################


# old versions of MEC for reference
'''
def read_filtering(list_of_reads, block_start_index, block_end_index):
    pos_allele_pairs = [(list_of_reads[i], list_of_reads[i + 1]) for i in range(0, len(list_of_reads), 2)]
    # gets rid of the reads that end before the block starts or starts after the block ends
    pos_allele_pairs_filtered = [elem for elem in pos_allele_pairs if (elem[0][-1]>= block_start_index and elem[0][0]<=block_end_index)]
        
    return pos_allele_pairs_filtered
    
map each read to the haplotype that it fits with best
returns the proportion of alleles that have to be coerced to make them match
** since we are taking the haplotype that matches best, it should never match with a nan row 

def mec(reconstructed_haplotypes, pos_allele_pairs):
  # gets rid of the reads that end before the block starts or starts after the block ends
  error_coercion_count = 0
  total_reads_size = 0
  
  for idx, (position_seq, allele_seq) in enumerate(pos_allele_pairs):
    subarray = reconstructed_haplotypes[position_seq[0] : position_seq[-1]+1]
    best_match_haplotype = np.argmax(np.sum(subarray==allele_seq, axis=1))
    error_coercion_count += np.sum(best_match_haplotype != allele_seq)
    total_reads_size += len(allele_seq)
  
  return error_coercion_count/total_reads_size

MEC function for reads in sparse matrix
for now, assuming sparse matrix is numpy array with lots of nans

def mec_sparse(reconstructed_haplotypes, reads):
  error_coercion_count = 0
  total_reads_size = 0
  
  for read in reads:
    subarray = reconstructed_haplotypes[~np.isnan(read)]
    observed_read = read[~np.isnan(read)]
    best_match_haplotype = np.argmax(np.sum(subarray==observed_read, axis=1))
    error_coercion_count += np.sum(best_match_haplotype != observed_read)
    total_reads_size += len(observed_read)
    
    return error_coercion_count/total_reads_size

can be used when not all haplotypes are phased
original indices refers to where in the full reconstruction reconstructed_haplotypes sits

def mec_probabalistic(reconstructed_haplotypes, pos_allele_pairs):
  expected_error = 0
  
  for idx, (position_seq, allele_seq) in enumerate(pos_allele_pairs):
    subarray = reconstructed_haplotypes[:, position_seq[0] : position_seq[-1]+1]
    for reconst_hap in subarray:
      likelihood = likelihood_haplotype(allele_seq, reconst_hap)
      error = np.sum(reconst_hap != allele_seq)
      expected_error += likelihood*error

  return expected_error
  
# positions that the read covers in the full reconstruction
# the read itself
# the reconstructed matrix we map the read to
# the slices in the full reconstruction that give H_star
# pos seq = list of positions that the read covers
# allele seq = the read itself
# H_star = a piece of the full reconstruction
# because of how this function is called, H_star will fit within read
# H_star start index = position in reconstruction corresponding to first col of H_star
# H_star start index = position in reconstruction corresponding to last col of H_star 
# (so, if we want to slice out the same SNPs as H_star, we have to add 1 to end index)
def mec_singleton(pos_seq, allele_seq, H_star, H_star_start_index, H_star_end_index):
  if H_star_end_index - H_star_start_index +1 != H_star.shape[1]:
    raise Exception("dimensions dont match")
  
  start_pos = pos_seq[0]
  end_pos = pos_seq[-1]
  # because of preprocess slicing, haplotype array H_star_within_read should fit fully in read
  # so, only need to cut read to fit within haplotype array
  if end_pos > H_star_end_index:
    sub_read = allele_seq[max(H_star_start_index-start_pos, 0):H_star_end_index+1]
  else:
    sub_read = allele_seq[max(H_star_start_index-start_pos, 0):]
  best_match_haplotype = H_star[np.argmax(np.sum(H_star==sub_read, axis=1))]
  error_correction_count = np.sum(best_match_haplotype != sub_read)
  
  return error_correction_count/len(sub_read)



# computes the "expected" mec 
#by weighing error correction with likelihood that that read is from that haplotype
def mec_probabalistic_singleton(pos_seq, allele_seq, H_star, H_star_start_index, H_star_end_index):
  if H_star_end_index - H_star_start_index +1 != H_star.shape[1]:
    raise Exception("dimensions dont match")
  
  start_pos = pos_seq[0]
  end_pos = pos_seq[-1]
  
  expected_error = 0
  
  # because of preprocess slicing, haplotype array H_star_within_read should fit fully in read
  # so, only need to cut read to fit within haplotype array
  if end_pos > H_star_end_index:
    sub_read = allele_seq[max(H_star_start_index-start_pos, 0):H_star_end_index+1]
  else:
    sub_read = allele_seq[max(H_star_start_index-start_pos, 0):] 
  standardizing_constant = sum((likelihood_haplotype(sub_read, reconst_hap)) for reconst_hap in H_star) 
  for reconst_hap in H_star:
    likelihood = likelihood_haplotype(sub_read, reconst_hap) / standardizing_constant
    error = np.sum(reconst_hap != sub_read)
    expected_error += likelihood*error
    
  return expected_error

def mec_multi_block_over_reads(H_star, list_of_reads, genotype):
  
  k = H_star.shape[0]
  
  mec_running_sum = 0
  
  pos_allele_pairs = [(list_of_reads[i], list_of_reads[i + 1]) for i in range(0, len(list_of_reads), 2)]

  blocks = find_blocks(np.array(H_star, dtype=np.float64), np.array(H_star, dtype=np.float64))
  
  # go through every read
  for idx, (position_seq, allele_seq) in enumerate(pos_allele_pairs):
    
    read_start_pos = position_seq[0]
    read_end_pos = position_seq[-1]
    
    blocks_within_read = [block for block in blocks if (block["start"]<= read_end_pos and block["end"]>=read_start_pos)]
    
    for block in blocks_within_read:
      # penalty for error correction
      
      H_star_block = H_star[:, block["start"]:block["end"]+1]
      
      if block["type"]==1:
        mec_running_sum += mec_singleton(position_seq, allele_seq, H_star_block, block["start"], block["end"])
      elif block["type"]==2:
        mec_running_sum += min(len(allele_seq), block["end"]-block["start"]+1)
      elif block["type"]==3:
        mec_running_sum += mec_singleton(position_seq, allele_seq, H_star_block, block["start"], block["end"])
      elif block["type"]==4:
        genotype_subarray = genotype[block['start']:block['end']+1]
        remaining_genotype = np.sum(H_star_block, axis=1)-genotype_subarray
        nan_row_idx = np.where(np.isnan(H_star_block))[0][0]
        H_star_block[nan_row_idx] = remaining_genotype
        mec_running_sum += mec_singleton(position_seq, allele_seq, H_star_block, block["start"], block["end"])
      
    # penalty for requiring len(blocks) many blocks to cover one read
    # where (k-1)/k is the probability of a switch error if the blocks had been one
    # for discrete uniform transition probabilities
    mec_running_sum += (k-1)/k * (len(blocks_within_read)-1)
    
  return mec_running_sum
    
def mec_probabalistic_multi_block_over_reads(H_star, list_of_reads, genotype):
  
  k = H_star.shape[0]
  
  mec_running_sum = 0
  
  pos_allele_pairs = [(list_of_reads[i], list_of_reads[i + 1]) for i in range(0, len(list_of_reads), 2)]
  
  # go through every read
  for idx, (position_seq, allele_seq) in enumerate(pos_allele_pairs):
    
    start_pos = position_seq[0]
    end_pos = position_seq[-1]
    # look at part of reconstruction that covers that read
    H_star_within_read = H_star[:, start_pos:end_pos+1]
    # go through each block within that portion of the reconstruction
    blocks = find_blocks(np.array(H_star_within_read, dtype=np.float64), np.array(H_star_within_read, dtype=np.float64))
    
    for block in blocks:
      # penalty for error correction
      mec_running_sum += mec_probabalistic_singleton(position_seq, allele_seq, H_star_within_read, block["start"], block["end"])
    
    # penalty for requiring len(blocks) many blocks to cover one read
    # where (k-1)/k is the probability of a switch error if the blocks had been one
    # for discrete uniform transition probabilities
    mec_running_sum += (k-1)/k * (len(blocks)-1)
    
  return mec_running_sum
  
def mec_multiple_blocks(H_star, list_of_reads, genotype):
  #H_star = reconstructed haplotypes, can be multiple blocks
  #genotype = numpy array of length = num snps
  blocks = find_blocks(np.array(H_star, dtype=np.float64), np.array(H_star, dtype=np.float64))
  total_mec = 0
  for block in blocks:
    pos_allele_pairs_filtered = read_filtering(list_of_reads, block["start"], block["end"])
    # print(block)
    block_H_star = H_star[:, block['start']:block['end']+1]
    if block['type'] == 1:
      total_mec += mec(block_H_star, pos_allele_pairs_filtered)
    elif block["type"] == 2 or block["type"] == 3:
      total_mec += mec_probabalistic(block_H_star, pos_allele_pairs_filtered) 
    elif block["type"] ==4:
      # fill in last row using the genotype
      genotype_subarray = genotype[block['start']:block['end']+1]
      remaining_genotype = np.sum(block_H_star, axis=1)-genotype_subarray
      nan_row_idx = np.where(np.isnan(block_H_star))[0][0]
      block_H_star[nan_row_idx] = remaining_genotype
      total_mec += mec(block_H_star, pos_allele_pairs_filtered)
      
  return total_mec
'''

    
