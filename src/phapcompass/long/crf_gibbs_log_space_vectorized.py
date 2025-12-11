import numpy as np
from typing import List, Dict, Tuple
from scipy.special import logsumexp

class HaplotypeGibbsSampler:
    """
    Vectorized Gibbs sampler for haplotype phasing using mixture of CRFs
    All computations done in log space to prevent underflow
    
    New algorithm:
    1. Updates read assignments using LOO vs standard potentials
    2. Constructs joint CRF over phases with hard genotype constraints
    3. Uses Viterbi on joint CRF to decode haplotypes
    4. Updates individual CRF potentials via gradient descent towards Viterbi solution
    """
    
    def __init__(self, reads, genotype, K, epsilon=0.00001, delta=5, learning_rate=0.01):
        """
        Initialize the Gibbs sampler
        
        Args:
            reads: numpy array of shape (n_reads, L) with values in {0, 1, np.nan}
                   np.nan indicates missing data at that position
            genotype: numpy array of shape (L,) with observed genotype values
            K: Ploidy (number of haplotypes)
            epsilon: Sequencing error rate
            delta: Smoothing parameter for transitions
            learning_rate: How much to gradient descent towards Viterbi solution (0 to 1)
                          0 = no update, 1 = full replacement with Viterbi
        """
        # Store input data
        self.reads = reads
        self.genotype = genotype
        self.K = K
        self.n_reads, self.L = reads.shape
        self.epsilon = epsilon
        self.delta = delta
        self.learning_rate = learning_rate
        
        # Log emission probabilities for each state to emit a match/mismatch
        self.log_emit_match = np.log(1 - epsilon) if epsilon < 1 else -np.inf
        self.log_emit_mismatch = np.log(epsilon) if epsilon > 0 else -np.inf
        
        # Precompute coverage masks (vectorized)
        # reads covering a certain position
        self.coverage_mask = ~np.isnan(reads)  # Shape: (n_reads, L)
        # reads covering a certain position with allele 0 or 1
        self.allele_0_mask = (reads == 0)  # Shape: (n_reads, L)
        self.allele_1_mask = (reads == 1)  # Shape: (n_reads, L)
        
        # Precompute config <-> index mappings for joint CRF
        # Number of joint states: 2^K configurations
        self.n_joint_states = 2 ** self.K

        # idx_to_config: array of shape (2^K, K) where row i is the binary config for state i
        self.idx_to_config = np.zeros((self.n_joint_states, self.K), dtype=int)
        for idx in range(self.n_joint_states):
            temp_idx = idx
            for k in range(self.K):
                self.idx_to_config[idx, k] = temp_idx % 2
                temp_idx //= 2

        # config_to_idx: use array indexing trick
        # A K-bit binary number maps directly to an integer
        # We can compute idx = sum(config[k] * 2^k for k in range(K))
        # But for fast lookup, we'll just use the array computation inline when needed
        
        # Storage for tracking quantities over iterations
        self.history = {
            'log_likelihood': [],
            'cluster_sizes': [],
            'n_switches': [],
            "marginals": [],
            "transitions": [],
            "phase": []
        }
    
    
    def fit(self, n_iterations=100, burn_in=50, verbose=False):
        """
        Run Gibbs sampling
        
        Args:
            n_iterations: Total number of iterations
            burn_in: Number of burn-in iterations to discard
            verbose: Whether to print progress
        
        Returns:
            self (for method chaining)
        """
        
        if verbose:
            print("Starting Gibbs sampling...")
            print(f"  Reads: {self.n_reads}, SNPs: {self.L}, Haplotypes: {self.K}")
            print(f"  Iterations: {n_iterations}, Burn-in: {burn_in}")
            print(f"  Learning rate: {self.learning_rate}")
        
        # ===== INITIALIZATION =====
        
        # Initialize cluster assignments randomly
        # according to π ~ Cat(1/K, ..., 1/K)
        self.z = np.random.randint(0, self.K, size=self.n_reads)
        
        # Build cluster membership indicator matrix (|R| x K)
        # true/false whether read j is in cluster k
        self.cluster_indicators = np.zeros((self.n_reads, self.K), dtype=bool)
        self.cluster_indicators[np.arange(self.n_reads), self.z] = True
        
        # Initialize log potential arrays for individual CRFs
        # psi_marginal[k, ℓ, a]
        # marginal potential for haplotype k at position ℓ with allele a
        self.log_psi_marginal = np.zeros((self.K, self.L, 2))
        self.log_psi_transition = np.zeros((self.K, self.L - 1, 2, 2))
        
        # Initialize log message arrays
        # fwd_msg[k, ℓ, a] and bwd_msg[k, ℓ, a]
        self.log_fwd_msg = np.zeros((self.K, self.L, 2))
        self.log_bwd_msg = np.zeros((self.K, self.L, 2))
        
        # Storage for samples after burn-in
        self.samples = []
        
        # Initial potential computation
        self._update_all_transition_potentials_vectorized()
        self._update_all_marginal_potentials_vectorized()
        self._compute_all_forward_messages_vectorized()
        self._compute_all_backward_messages_vectorized()
        
        
        # ===== MAIN GIBBS LOOP =====
        
        for iteration in range(n_iterations):
            
            if verbose and (iteration + 1) % 100 == 0:
                print(f"Iteration {iteration + 1}/{n_iterations}")
            
            # ----- STEP 1: Update read assignments (sequential) -----
            
            read_order = np.random.permutation(self.n_reads)
            n_switches = 0
            
            for j in read_order:
                switched = self._update_single_read(j)
                if switched:
                    n_switches += 1
            
            # ----- STEP 2 & 3: Construct joint CRF and run Viterbi -----
            # ----- STEP 4: Update potentials via gradient descent -----
            if iteration < 20:
                # Get MAP estimate of haplotypes via joint Viterbi
                ffbs_haplotypes = self._ffbs_decode_joint_phase()
                
                self._gradient_descent_update_potentials(ffbs_haplotypes)
            
                # Recompute messages with updated potentials
                self._compute_all_forward_messages_vectorized()
                self._compute_all_backward_messages_vectorized()
                
            viterbi_haplotypes = self._viterbi_decode_joint_phase()
            
            # ----- Track quantities -----
            
            log_likelihood = self._compute_log_likelihood()
            cluster_sizes = np.sum(self.cluster_indicators, axis=0).tolist()
            
            self.history['log_likelihood'].append(log_likelihood)
            self.history['cluster_sizes'].append(cluster_sizes)
            self.history['n_switches'].append(n_switches)
            self.history["marginals"].append(np.exp(self.log_psi_marginal.copy()))
            self.history["transitions"].append(np.exp(self.log_psi_transition.copy()))
            self.history["phase"].append(viterbi_haplotypes.copy())
            
            # ----- Store sample -----
            
            if iteration >= burn_in:
                self.samples.append(self.z.copy())
                
            if  iteration > 20:
                moving_avg = np.mean(self.history["log_likelihood"][-6:-1])
                if np.abs(moving_avg - log_likelihood)/np.abs(moving_avg) < 0.001:
                    print(f"terminated on iteration {iteration}")
                    return self
            
            # ----- Print diagnostics -----
            
            if verbose and (iteration + 1) % 10 == 0:
                print(f"  Cluster sizes: {cluster_sizes}")
                print(f"  Log-likelihood: {log_likelihood:.2f}")
                print(f"  Switches: {n_switches}/{self.n_reads}")            
        
        if verbose:
            print("Gibbs sampling complete!")
        
        return self
    
    
    def _update_all_transition_potentials_vectorized(self):
        """
        Vectorized update of log ψ^k_{ℓ,ℓ+1}(a₁, a₂) for all k, ℓ, a₁, a₂
        based on vote counts from reads in each cluster
        vectorized over positions
        not vectorized over clusters or allele pairs
        """
        # Precompute masks for consecutive positions
        # Shape: (n_reads, L-1)
        # covers_both[j, ell] = read j covers both position ell and ell+1
        covers_both = self.coverage_mask[:, :-1] & self.coverage_mask[:, 1:]
        
        # For each cluster k
        for k in range(self.K):
            cluster_mask = self.cluster_indicators[:, k]  # Shape: (n_reads,)
            
            # Vectorized computation over all ℓ
            # Shape: (n_reads, L-1)
            in_cluster_covers_both = cluster_mask[:, None] & covers_both
            
            # Count denominators: reads in cluster k covering both positions
            # Shape: (L-1,)
            denominators = np.sum(in_cluster_covers_both, axis=0)
            
            # For each allele pair (a1, a2)
            for a1 in [0, 1]:
                for a2 in [0, 1]:
                    # Get masks for alleles at consecutive positions
                    if a1 == 0:
                        mask_ell = self.allele_0_mask[:, :-1]
                    else:
                        mask_ell = self.allele_1_mask[:, :-1]
                    
                    if a2 == 0:
                        mask_ell_next = self.allele_0_mask[:, 1:]
                    else:
                        mask_ell_next = self.allele_1_mask[:, 1:]
                    
                    # Combined mask: in cluster, covers both, has correct alleles
                    # ie, reads in cluster k that exactly match (a1, a2)
                    # Shape: (n_reads, L-1)
                    full_mask = (cluster_mask[:, None] & 
                                mask_ell & mask_ell_next)
                    
                    # Count numerators
                    # going from yes/no tallies of matching reads to counts
                    # Shape: (L-1,)
                    numerators = np.sum(full_mask, axis=0)
                    
                    # Compute log potentials with smoothing
                    # Avoid division by zero
                    valid = denominators > 0
                    
                    # set to just delta in case not valid
                    log_psi = np.full(self.L - 1, np.log(self.delta))
                    # When we have data: log(num/denom + delta)
                    log_psi[valid] = np.log(numerators[valid] + self.delta * denominators[valid]) - np.log(denominators[valid])
                    
                    self.log_psi_transition[k, :, a1, a2] = log_psi
            
            # Normalize: convert to log probabilities
            # Shape: (L-1, 2, 2) -> normalize over last dimension for each (ℓ, a1)
            for a1 in [0, 1]:
                log_vals = self.log_psi_transition[k, :, a1, :]  # Shape: (L-1, 2)
                log_sums = logsumexp(log_vals, axis=1)  # Shape: (L-1,)
                
                # Handle non-finite sums
                finite_mask = np.isfinite(log_sums)
                self.log_psi_transition[k, finite_mask, a1, :] -= log_sums[finite_mask, None]
                self.log_psi_transition[k, ~finite_mask, a1, :] = np.log(0.5)
    
    
    def _update_all_marginal_potentials_vectorized(self):
        """
        Vectorized update of log ψ^k_ℓ(a) for all k, ℓ, a
        based on emission probabilities from reads in each cluster
        vectorized over positions
        
        NOTE: No longer uses genotype term here - that's enforced in joint CRF
        """
        for k in range(self.K):
            cluster_mask = self.cluster_indicators[:, k]  # Shape: (n_reads,)
            
            for a in [0, 1]:
                # Get reads in cluster k with allele a at each position
                if a == 0:
                    allele_mask = self.allele_0_mask  # Shape: (n_reads, L)
                else:
                    allele_mask = self.allele_1_mask
                
                # Reads in cluster k covering each position with allele a
                # Shape: (n_reads, L)
                # covers with the correct allele
                match_mask = cluster_mask[:, None] & allele_mask
                # covers but with the wrong allele
                mismatch_mask = cluster_mask[:, None] & self.coverage_mask & ~allele_mask
                
                # Count matches and mismatches per position
                # Shape: (L,)
                n_matches = np.sum(match_mask, axis=0)
                n_mismatches = np.sum(mismatch_mask, axis=0)
                
                # Emission product over all reads in cluster k covering position ℓ
                # how well do the reads in cluster k support allele a at position ℓ
                # Compute emission sum in log space
                # Handle epsilon=0 case where log_emit_mismatch=-inf
                log_emission_sum = n_matches * self.log_emit_match
                if self.epsilon > 0:
                    log_emission_sum += n_mismatches * self.log_emit_mismatch
                else:
                    # When epsilon=0, any mismatch makes potential -inf
                    log_emission_sum = np.where(n_mismatches > 0, -np.inf, log_emission_sum)
                
                self.log_psi_marginal[k, :, a] = log_emission_sum
            
            # Normalize to log probabilities (vectorized over ℓ)
            log_sums = logsumexp(self.log_psi_marginal[k, :, :], axis=1)  # Shape: (L,)
            finite_mask = np.isfinite(log_sums)
            # if finite, normalize as usual
            self.log_psi_marginal[k, finite_mask, :] -= log_sums[finite_mask, None]
            # If both are zero (-inf), use uniform
            self.log_psi_marginal[k, ~finite_mask, :] = np.log(0.5)
    
    
    def _compute_all_forward_messages_vectorized(self):
        """
        Compute forward messages in log space (vectorized over positions)
        vectorized over alleles
        """
        for k in range(self.K):
            # Base case: ℓ = 0
            self.log_fwd_msg[k, 0, :] = self.log_psi_marginal[k, 0, :]
            
            # Vectorized recurrence over positions
            for ell in range(1, self.L):
                # Compute transition contributions
                # For each a1 (next state), sum over a2 (prev state)
                # log_fwd_msg[k, ell-1, :] has shape (2,)
                # log_psi_transition[k, ell-1, :, :] has shape (2, 2)
                
                # Broadcast: (2,) + (2, 2) -> (2, 2) via broadcasting
                # log_terms[a2, a1] = log_fwd_msg[a2] + log_psi_transition[a2, a1]
                log_terms = (self.log_fwd_msg[k, ell - 1, :, None] + 
                            self.log_psi_transition[k, ell - 1, :, :])
                
                # Sum over previous states (axis 0)
                log_transition_contrib = logsumexp(log_terms, axis=0)  # Shape: (2,)
                
                # Add marginal potential
                self.log_fwd_msg[k, ell, :] = (
                    self.log_psi_marginal[k, ell, :] + log_transition_contrib
                )
    
    
    def _compute_all_backward_messages_vectorized(self):
        """
        Compute backward messages in log space (vectorized over positions)
        """
        for k in range(self.K):
            # Base case: ℓ = L-1
            self.log_bwd_msg[k, self.L - 1, :] = self.log_psi_marginal[k, self.L - 1, :]
            
            # Vectorized recurrence over positions (backward)
            for ell in range(self.L - 2, -1, -1):
                # For each a1 (current state), sum over a2 (next state)
                # log_bwd_msg[k, ell+1, :] has shape (2,)
                # log_psi_transition[k, ell, :, :] has shape (2, 2)
                
                # Broadcast: (2,) + (2, 2) -> (2, 2)
                # log_terms[a1, a2] = log_bwd_msg[a2] + log_psi_transition[a1, a2]
                log_terms = (self.log_bwd_msg[k, ell + 1, None, :] + 
                            self.log_psi_transition[k, ell, :, :])
                
                # Sum over next states (axis 1)
                log_transition_contrib = logsumexp(log_terms, axis=1)  # Shape: (2,)
                
                # Add marginal potential
                self.log_bwd_msg[k, ell, :] = (
                    self.log_psi_marginal[k, ell, :] + log_transition_contrib
                )
    
    
    def _update_single_read(self, j):
        """
        Update assignment z[j] for read j (must remain sequential)
        
        Uses LOO potentials when computing likelihood for current cluster,
        standard potentials for other clusters.
        
        Args:
            j: Read index
        
        Returns:
            bool: Whether the read switched clusters
        """
        k_old = self.z[j]
        
        positions = self._get_read_positions(j)
        if len(positions) == 0:
            return False
        
        ell_start = positions[0]
        ell_end = positions[-1]
        
        # Compute log likelihoods for all clusters
        log_likelihoods = np.zeros(self.K)
        for k in range(self.K):
            if k == k_old:
                # Use leave-one-out potentials for current cluster
                log_likelihoods[k] = self._compute_read_log_likelihood_LOO(j, k)
            else:
                # Use standard potentials for other clusters
                log_likelihoods[k] = self._compute_read_log_likelihood_standard(j, k)
        
        # Handle case where all log likelihoods are -inf or contain NaN
        if np.all(np.isinf(log_likelihoods)) or np.any(np.isnan(log_likelihoods)):
            probabilities = np.ones(self.K) / self.K
        else:
            log_probs = log_likelihoods - logsumexp(log_likelihoods)
            probabilities = np.exp(log_probs)
            
            if np.any(np.isnan(probabilities)):
                probabilities = np.ones(self.K) / self.K
        
        # Sample new assignment
        k_new = np.random.choice(self.K, p=probabilities)
        
        # Update if changed
        if k_new != k_old:
            self.cluster_indicators[j, k_old] = False
            self.z[j] = k_new
            self.cluster_indicators[j, k_new] = True
            
            # Update potentials for both clusters in affected region
            self._update_potentials_for_cluster(k_old, ell_start, ell_end)
            self._update_potentials_for_cluster(k_new, ell_start, ell_end)
            
            # Recompute forward and backward messages for affected clusters
            self._recompute_messages_for_cluster(k_old)
            self._recompute_messages_for_cluster(k_new)
            
            return True
        else:
            return False


    def _recompute_messages_for_cluster(self, k):
        """
        Recompute forward and backward messages for a single cluster k
        """
        # Forward messages
        self.log_fwd_msg[k, 0, :] = self.log_psi_marginal[k, 0, :]
        
        for ell in range(1, self.L):
            log_terms = (self.log_fwd_msg[k, ell - 1, :, None] + 
                        self.log_psi_transition[k, ell - 1, :, :])
            log_transition_contrib = logsumexp(log_terms, axis=0)
            self.log_fwd_msg[k, ell, :] = (
                self.log_psi_marginal[k, ell, :] + log_transition_contrib
            )
        
        # Backward messages
        self.log_bwd_msg[k, self.L - 1, :] = self.log_psi_marginal[k, self.L - 1, :]
        
        for ell in range(self.L - 2, -1, -1):
            log_terms = (self.log_bwd_msg[k, ell + 1, None, :] + 
                        self.log_psi_transition[k, ell, :, :])
            log_transition_contrib = logsumexp(log_terms, axis=1)
            self.log_bwd_msg[k, ell, :] = (
                self.log_psi_marginal[k, ell, :] + log_transition_contrib
            )
        
        
    def _update_potentials_for_cluster(self, k, ell_start=None, ell_end=None):
        """
        Update both transition and marginal potentials for a single cluster k
        in a specific region (used after read reassignment)
        vectorized over positions
        """
        if ell_start is None:
            ell_start = 0
        if ell_end is None:
            ell_end = self.L - 1
        
        cluster_mask = self.cluster_indicators[:, k]
        
        # Update transition potentials
        # Need to update transitions from (ell_start - 1) → ell_start up to ell_end → (ell_end + 1)
        trans_start = max(0, ell_start - 1)
        # L-2 because transitions are indexed up to L-2
        # the last transition is from L-2 to L-1, transitions indexed by source position
        trans_end = min(self.L - 2, ell_end)
        
        if trans_end >= trans_start:
            # Extract relevant slice
            ell_range = np.arange(trans_start, trans_end + 1)
            
            covers_both = (self.coverage_mask[:, trans_start:trans_end+1] & 
                          self.coverage_mask[:, trans_start+1:trans_end+2])
            
            in_cluster_covers_both = cluster_mask[:, None] & covers_both
            denominators = np.sum(in_cluster_covers_both, axis=0)
            
            for a1 in [0, 1]:
                for a2 in [0, 1]:
                    if a1 == 0:
                        mask_ell = self.allele_0_mask[:, trans_start:trans_end+1]
                    else:
                        mask_ell = self.allele_1_mask[:, trans_start:trans_end+1]
                    
                    if a2 == 0:
                        mask_ell_next = self.allele_0_mask[:, trans_start+1:trans_end+2]
                    else:
                        mask_ell_next = self.allele_1_mask[:, trans_start+1:trans_end+2]
                    
                    full_mask = cluster_mask[:, None] & mask_ell & mask_ell_next
                    numerators = np.sum(full_mask, axis=0)
                    
                    valid = denominators > 0
                    log_psi = np.full(len(ell_range), np.log(self.delta))
                    log_psi[valid] = np.log(numerators[valid] + self.delta * denominators[valid]) - np.log(denominators[valid])
                    
                    self.log_psi_transition[k, ell_range, a1, a2] = log_psi
            
            # Normalize
            for a1 in [0, 1]:
                log_vals = self.log_psi_transition[k, ell_range, a1, :]
                log_sums = logsumexp(log_vals, axis=1)
                finite_mask = np.isfinite(log_sums)
                
                self.log_psi_transition[k, ell_range[finite_mask], a1, :] -= log_sums[finite_mask, None]
                self.log_psi_transition[k, ell_range[~finite_mask], a1, :] = np.log(0.5)
        
        # Update marginal potentials
        # Update from ell_start to ell_end (inclusive)
        ell_range = np.arange(ell_start, ell_end + 1)
        
        for a in [0, 1]:
            if a == 0:
                allele_mask = self.allele_0_mask[:, ell_range]
            else:
                allele_mask = self.allele_1_mask[:, ell_range]
            
            match_mask = cluster_mask[:, None] & allele_mask
            mismatch_mask = cluster_mask[:, None] & self.coverage_mask[:, ell_range] & ~allele_mask
            
            n_matches = np.sum(match_mask, axis=0)
            n_mismatches = np.sum(mismatch_mask, axis=0)
            
            # Compute emission sum in log space
            # Handle epsilon=0 case where log_emit_mismatch=-inf
            log_emission_sum = n_matches * self.log_emit_match
            if self.epsilon > 0:
                log_emission_sum += n_mismatches * self.log_emit_mismatch
            else:
                # When epsilon=0, any mismatch makes potential -inf
                log_emission_sum = np.where(n_mismatches > 0, -np.inf, log_emission_sum)
            
            self.log_psi_marginal[k, ell_range, a] = log_emission_sum
        
        # Normalize
        log_sums = logsumexp(self.log_psi_marginal[k, ell_range, :], axis=1)
        finite_mask = np.isfinite(log_sums)
        
        self.log_psi_marginal[k, ell_range[finite_mask], :] -= log_sums[finite_mask, None]
        self.log_psi_marginal[k, ell_range[~finite_mask], :] = np.log(0.5)
    
    
    def _get_read_positions(self, j):
        """
        Get the positions (SNP indices) covered by read j
        """
        return np.where(self.coverage_mask[j, :])[0]
    
    
    def _compute_read_log_likelihood_LOO(self, j, k):
        """
        Compute log Pr(r_j ~ k) using leave-one-out potentials
        
        Divides out read j's contribution from marginals and transitions
        """
        positions = self._get_read_positions(j)
        if len(positions) == 0:
            return 0.0
        
        ell_start = positions[0]
        ell_end = positions[-1]
        
        # Get LOO marginal potentials (divide out read j's emission probabilities)
        log_psi_marginal_loo = self.log_psi_marginal[k, :, :].copy()
        
        # Convert positions to array for vectorization
        pos_array = np.array(positions)
        read_allele = self.reads[j, pos_array]
        
        # Divide out read j's contribution to marginals
        for i, ell in enumerate(pos_array):
            for a in [0, 1]:
                if read_allele[i] == a:
                    log_psi_marginal_loo[ell, a] -= self.log_emit_match
                else:
                    if self.epsilon > 0:
                        log_psi_marginal_loo[ell, a] -= self.log_emit_mismatch
            # Re-normalize
            log_sum = logsumexp(log_psi_marginal_loo[ell, :])
            log_psi_marginal_loo[ell, :] -= log_sum
                
        # Get LOO transition potentials
        log_psi_transition_loo = self._get_LOO_transition_potentials(j, k)
        
        # Get forward message from before read
        # the forward message is supposed to summarize 
        # everything before and up to the read
        # if the read starts at position 0, there is nothing before it, 
        # so nothing to summarize => uniform
        if ell_start == 0:
            log_fwd_msg_in = np.zeros(2)
        else:
            log_fwd_msg_in = self.log_fwd_msg[k, ell_start, :]
        # Get backward message from after read
        if ell_end == self.L - 1:
            log_bwd_msg_out = np.zeros(2)
        else:
            log_bwd_msg_out = self.log_bwd_msg[k, ell_end, :]
        
        # Compute log likelihood
        log_likelihood = self._compute_log_likelihood_matrix_product(
            j, k, positions, log_fwd_msg_in, log_bwd_msg_out,
            log_psi_marginal_loo, log_psi_transition_loo
        )
        
        return log_likelihood
    
    
    def _compute_read_log_likelihood_standard(self, j, k):
        """
        Compute log Pr(r_j ~ k) using standard potentials (not LOO)
        """
        positions = self._get_read_positions(j)
        if len(positions) == 0:
            return 0.0
        
        ell_start = positions[0]
        ell_end = positions[-1]
        # Use standard potentials
        log_psi_marginal_std = self.log_psi_marginal[k, :, :]
        log_psi_transition_std = self.log_psi_transition[k, :, :, :]
        # Get forward/backward messages
        if ell_start == 0:
            log_fwd_msg_in = np.zeros(2)
        else:
            log_fwd_msg_in = self.log_fwd_msg[k, ell_start, :]
        
        if ell_end == self.L - 1:
            log_bwd_msg_out = np.zeros(2)
        else:
            log_bwd_msg_out = self.log_bwd_msg[k, ell_end, :]
        
        log_likelihood = self._compute_log_likelihood_matrix_product(
            j, k, positions, log_fwd_msg_in, log_bwd_msg_out,
            log_psi_marginal_std, log_psi_transition_std
        )
        
        return log_likelihood
    
    
    def _get_LOO_transition_potentials(self, j, k):
        """
        Compute LOO transition potentials in log space for read j in cluster k
        
        Excludes read j from vote counts
        """
        log_psi_transition_loo = self.log_psi_transition[k, :, :, :].copy()
        
        positions = self._get_read_positions(j)
        # Only update transitions where read j covers both positions
        for i in range(len(positions) - 1):
            ell = positions[i]
            ell_next = positions[i + 1]
            # Only handle adjacent positions
            if ell_next == ell + 1:
                # Exclude read j from counts
                # cluster mask defined locally to safely exclude read j
                # so don't count j towards any of the 4 transitions
                cluster_mask = self.cluster_indicators[:, k].copy()
                cluster_mask[j] = False
                
                for a1 in [0, 1]:
                    for a2 in [0, 1]:
                        if a1 == 0:
                            mask_ell = self.allele_0_mask[:, ell]
                        else:
                            mask_ell = self.allele_1_mask[:, ell]
                        
                        if a2 == 0:
                            mask_ell_next = self.allele_0_mask[:, ell + 1]
                        else:
                            mask_ell_next = self.allele_1_mask[:, ell + 1]
                        
                        numerator = np.sum(cluster_mask & mask_ell & mask_ell_next)
                        denominator = np.sum(cluster_mask & 
                                            self.coverage_mask[:, ell] & 
                                            self.coverage_mask[:, ell + 1])
                        
                        if denominator > 0:
                            log_psi_transition_loo[ell, a1, a2] = (
                                np.log(numerator + self.delta * denominator) - 
                                np.log(denominator)
                            )
                        else:
                            log_psi_transition_loo[ell, a1, a2] = np.log(self.delta)
        
        return log_psi_transition_loo
    
    
    def _compute_log_likelihood_matrix_product(self, j, k, positions, 
                                               log_fwd_msg_in, log_bwd_msg_out,
                                               log_psi_marginal, log_psi_transition):
        """
        Compute read log likelihood using log-space matrix multiplication
        
        Chains together marginal×emission terms and transitions across the read
        note - not vectorized over positions
        """
        # Start with forward message (2D vector)
        log_result = log_fwd_msg_in.copy()
        
        # Process each position in the read one by one
        for i, ell in enumerate(positions):
            # Create log diagonal matrix: marginal × emission
            log_emission_diag = np.full((2, 2), -np.inf)
            read_allele = self.reads[j, ell]
            
            for a in [0, 1]:
                if read_allele == a:
                    log_emission = self.log_emit_match
                else:
                    log_emission = self.log_emit_mismatch
                
                log_emission_diag[a, a] = log_psi_marginal[ell, a] + log_emission
            
            # Log-space matrix-vector multiplication of likelihood we've
            # accumulated so far by the marginal × emission diagonal matrix
            log_result = self._log_matmul_vec(log_emission_diag, log_result)
            
            # If not the last position, handle transition to next position
            if i < len(positions) - 1:
                ell_next = positions[i + 1]
                
                if ell_next == ell + 1:
                    # Adjacent positions: use single transition matrix between them
                    log_trans_matrix = log_psi_transition[ell, :, :]
                    log_result = self._log_matmul_vec(log_trans_matrix, log_result)
                else:
                    # Non-adjacent: chain all transition matrices over gap
                    log_result = self._bridge_gap_with_transitions(
                        k, ell, ell_next, log_result, log_psi_transition
                    )
        
        # Finally multiply by backward message
        log_likelihood = logsumexp(log_bwd_msg_out + log_result)
        
        return log_likelihood


    def _bridge_gap_with_transitions(self, k, ell_start, ell_end, 
                                        log_result, log_psi_transition):
        """
        Bridge a gap in read coverage using forward message products
        
        Chains transition matrices across positions not covered by the read
        
        Args:
            k: Cluster index
            ell_start: Last covered position before gap
            ell_end: First covered position after gap
            log_result: Current log probability vector (shape: (2,))
            log_psi_transition: Transition potentials to use
        
        Returns:
            Updated log probability vector after bridging the gap
        """
        # Start from position after ell_start and propagate forward to ell_end
        current_log_vec = log_result.copy()
        
        # Multiply by transitions through the gap
        for ell in range(ell_start, ell_end):
            log_trans_matrix = log_psi_transition[ell, :, :]
            current_log_vec = self._log_matmul_vec(log_trans_matrix, current_log_vec)
        
        return current_log_vec    
    
    def _log_matmul_vec(self, log_matrix, log_vec):
        """
        Compute matrix-vector product in log space (vectorized)
        
        result[i] = logsumexp_j(log_matrix[i, j] + log_vec[j])
        """
        # Broadcast: (n, m) + (m,) -> (n, m)
        log_terms = log_matrix + log_vec[None, :]
        # Sum over columns (axis 1)
        return logsumexp(log_terms, axis=1)
    
    
    def _ffbs_decode_joint_phase(self):
        """
        Sample a joint phase path using forward-backward algorithm
        
        1. Forward pass: compute forward messages
        2. Backward sampling: sample states sequentially from end to start
        
        Returns:
            haplotypes: Array of shape (K, L) with sampled haplotypes
        """
        n_states = self.n_joint_states
        
        # ===== FORWARD PASS =====
        # Initialize log forward messages
        log_alpha = np.full((self.L, n_states), -np.inf)
        
        # Base case: ℓ = 0
        log_alpha[0, :] = self._compute_joint_marginal_potentials(0)
        
        # Forward recurrence
        for ell in range(1, self.L):
            log_marginals_ell = self._compute_joint_marginal_potentials(ell)
            
            for next_state_idx in range(n_states):
                log_terms = []
                for prev_state_idx in range(n_states):
                    log_transitions = self._compute_joint_transition_potentials(ell - 1, prev_state_idx)
                    log_term = log_alpha[ell - 1, prev_state_idx] + log_transitions[next_state_idx]
                    log_terms.append(log_term)
                
                log_alpha[ell, next_state_idx] = log_marginals_ell[next_state_idx] + logsumexp(log_terms)
        
        # ===== BACKWARD SAMPLING =====
        path = np.zeros(self.L, dtype=int)
        
        # Sample final state from forward messages
        log_probs_final = log_alpha[-1, :]
        log_probs_final -= logsumexp(log_probs_final)  # Normalize
        probs_final = np.exp(log_probs_final)
        
        # Handle numerical issues
        if np.any(np.isnan(probs_final)) or np.sum(probs_final) == 0:
            probs_final = np.ones(n_states) / n_states
        else:
            probs_final = probs_final / np.sum(probs_final)
        
        path[-1] = np.random.choice(n_states, p=probs_final)
        
        # Sample backwards through the sequence
        for ell in range(self.L - 2, -1, -1):
            next_state_idx = path[ell + 1]
            
            # Compute conditional probabilities p(state_ell | state_ell+1, observations)
            log_probs = np.full(n_states, -np.inf)
            
            for prev_state_idx in range(n_states):
                log_transitions = self._compute_joint_transition_potentials(ell, prev_state_idx)
                log_probs[prev_state_idx] = (
                    log_alpha[ell, prev_state_idx] + 
                    log_transitions[next_state_idx]
                )
            
            # Normalize
            log_probs -= logsumexp(log_probs)
            probs = np.exp(log_probs)
            
            # Handle numerical issues
            if np.any(np.isnan(probs)) or np.sum(probs) == 0:
                probs = np.ones(n_states) / n_states
            else:
                probs = probs / np.sum(probs)
            
            path[ell] = np.random.choice(n_states, p=probs)
        
        # Convert path to haplotype matrix
        haplotypes = np.zeros((self.K, self.L), dtype=int)
        for ell in range(self.L):
            config = self._get_config(path[ell])
            haplotypes[:, ell] = config
        
        return haplotypes
    
    def _ffbs_decode_joint_phase_multiple(self, m):
        """
        Sample m independent joint phase paths using forward-backward algorithm
        
        1. Forward pass: compute forward messages (done once)
        2. Backward sampling: sample m independent state sequences from end to start
        3. For each path, compute log likelihood as sum of marginal and transition potentials
        4. NORMALIZE log likelihoods 
            first, normalized in log space to be a valid log density
            then, rescaled so that the log densities sum to 1 for interpretability
        
        Args:
            m: Number of independent samples to generate
            
        Returns:
            haplotypes_samples: List of m tuples (haplotypes, normalized_log_likelihood) where:
                - haplotypes: Array of shape (K, L) with sampled haplotypes
                - normalized_log_likelihood: Normalized log probability such that 
                sum over m samples of exp(normalized_log_likelihood) = 1
        """
        n_states = self.n_joint_states
        
        # ===== FORWARD PASS (computed once) =====
        # Initialize log forward messages
        log_alpha = np.full((self.L, n_states), -np.inf)
        
        # Base case: ℓ = 0
        log_alpha[0, :] = self._compute_joint_marginal_potentials(0)
        
        # Forward recurrence
        for ell in range(1, self.L):
            log_marginals_ell = self._compute_joint_marginal_potentials(ell)
            
            for next_state_idx in range(n_states):
                log_terms = []
                for prev_state_idx in range(n_states):
                    log_transitions = self._compute_joint_transition_potentials(ell - 1, prev_state_idx)
                    log_term = log_alpha[ell - 1, prev_state_idx] + log_transitions[next_state_idx]
                    log_terms.append(log_term)
                
                log_alpha[ell, next_state_idx] = log_marginals_ell[next_state_idx] + logsumexp(log_terms)
        
        # ===== BACKWARD SAMPLING (repeated m times) =====
        haplotypes_list = []
        path_log_likelihoods = []  # Collect all path scores for normalization
        
        for sample_idx in range(m):
            path = np.zeros(self.L, dtype=int)
            path_log_likelihood = 0.0  # Initialize path score accumulator
            
            # Sample final state from forward messages
            log_probs_final = log_alpha[-1, :]
            log_probs_final -= logsumexp(log_probs_final)  # Normalize
            probs_final = np.exp(log_probs_final)
            
            # Handle numerical issues
            if np.any(np.isnan(probs_final)) or np.sum(probs_final) == 0:
                probs_final = np.ones(n_states) / n_states
            else:
                probs_final = probs_final / np.sum(probs_final)
            
            path[-1] = np.random.choice(n_states, p=probs_final)
            
            # Add marginal potential for final position to path score
            log_marginals_final = self._compute_joint_marginal_potentials(self.L - 1)
            path_log_likelihood += log_marginals_final[path[-1]]
            
            # Sample backwards through the sequence
            for ell in range(self.L - 2, -1, -1):
                next_state_idx = path[ell + 1]
                
                # Compute conditional probabilities p(state_ell | state_ell+1, observations)
                log_probs = np.full(n_states, -np.inf)
                
                for prev_state_idx in range(n_states):
                    log_transitions = self._compute_joint_transition_potentials(ell, prev_state_idx)
                    log_probs[prev_state_idx] = (
                        log_alpha[ell, prev_state_idx] + 
                        log_transitions[next_state_idx]
                    )
                
                # Normalize
                log_probs -= logsumexp(log_probs)
                probs = np.exp(log_probs)
                
                # Handle numerical issues
                if np.any(np.isnan(probs)) or np.sum(probs) == 0:
                    probs = np.ones(n_states) / n_states
                else:
                    probs = probs / np.sum(probs)
                
                path[ell] = np.random.choice(n_states, p=probs)
                
                # Add marginal potential for current position to path score
                log_marginals_ell = self._compute_joint_marginal_potentials(ell)
                path_log_likelihood += log_marginals_ell[path[ell]]
                
                # Add transition potential from ell to ell+1 to path score
                log_transitions_ell = self._compute_joint_transition_potentials(ell, path[ell])
                path_log_likelihood += log_transitions_ell[next_state_idx]
            
            # Convert path to haplotype matrix
            haplotypes = np.zeros((self.K, self.L), dtype=int)
            for ell in range(self.L):
                config = self._get_config(path[ell])
                haplotypes[:, ell] = config
            
            # Store haplotypes and unnormalized log likelihood
            haplotypes_list.append(haplotypes)
            path_log_likelihoods.append(path_log_likelihood)
            
        # ===== NORMALIZE PATH LOG LIKELIHOODS =====
        # Convert to numpy array for vectorized operations
        path_log_likelihoods = np.array(path_log_likelihoods)
        
        # Step 1: Normalize into valid log density (probabilities sum to 1 in linear space)
        log_Z = logsumexp(path_log_likelihoods)
        log_density = path_log_likelihoods - log_Z
        
        # Step 2: Make log densities sum to 1 (rescale in log space)
        log_density_sum = np.sum(log_density)
        normalized_log_likelihoods = log_density / log_density_sum
            
        # Create final list of tuples with normalized scores
        haplotypes_samples = [
            (haplotypes_list[i], normalized_log_likelihoods[i]) 
            for i in range(m)
        ]

        return haplotypes_samples

    def _viterbi_decode_joint_phase(self):
        """
        Run Viterbi algorithm over joint phase CRF
        
        Returns:
            haplotypes: Array of shape (K, L) with inferred haplotypes
        """
        n_states = self.n_joint_states
        
        # Initialize log Viterbi arrays
        log_V = _np.full((self.L, n_states), -np.inf)
        backpointer = np.zeros((self.L, n_states), dtype=int)
        
        # Base case: ℓ = 0
        log_V[0, :] = self._compute_joint_marginal_potentials(0)
        
        # Recurrence
        for ell in range(1, self.L):
            # Get marginals for current position
            log_marginals_ell = self._compute_joint_marginal_potentials(ell)
            
            for next_state_idx in range(n_states):
                best_log_val = -np.inf
                best_prev_idx = 0
                
                for prev_state_idx in range(n_states):
                    # Get normalized transitions from prev to next
                    log_transitions = self._compute_joint_transition_potentials(ell - 1, prev_state_idx)
                    
                    log_val = log_V[ell - 1, prev_state_idx] + log_transitions[next_state_idx]
                    
                    if log_val > best_log_val:
                        best_log_val = log_val
                        best_prev_idx = prev_state_idx
                
                backpointer[ell, next_state_idx] = best_prev_idx
                log_V[ell, next_state_idx] = log_marginals_ell[next_state_idx] + best_log_val
        
        # Backtracking
        path = np.zeros(self.L, dtype=int)
        path[-1] = np.argmax(log_V[-1, :])
        
        for ell in range(self.L - 2, -1, -1):
            path[ell] = backpointer[ell + 1, path[ell + 1]]
        
        # Convert path to haplotype matrix
        haplotypes = np.zeros((self.K, self.L), dtype=int)
        for ell in range(self.L):
            config = self._get_config(path[ell])
            haplotypes[:, ell] = config
        
        return haplotypes
    
    def _get_config(self, idx):
        """
        Get configuration for a given state index (fast array lookup)
        
        Args:
            idx: State index (0 to 2^K - 1)
        
        Returns:
            config: Binary vector of length K
        """
        return self.idx_to_config[idx]


    def _get_idx(self, config):
        """
        Get state index for a given configuration (fast computation)
        
        Args:
            config: Binary vector of length K
        
        Returns:
            idx: State index (0 to 2^K - 1)
        """
        idx = 0
        for k in range(self.K):
            idx += config[k] * (2 ** k)
        return idx


    def _compute_joint_marginal_potentials(self, ell):
        """
        Compute NORMALIZED marginal potentials for all configs at position ℓ
        
        For each config: potential = 𝟙[sum(config) == genotype] · ∏_k ψ^k_ℓ(config[k])
        Then normalize over all 2^K configurations
        
        Args:
            ell: Position index
        
        Returns:
            log_potentials: array of shape (2^K,) with normalized log probabilities
        """
        n_states = self.n_joint_states
        log_potentials = np.full(n_states, -np.inf)
        
        # Compute unnormalized potentials for all configs
        for state_idx in range(n_states):
            config = self._get_config(state_idx)
            config_sum = np.sum(config)
            
            if config_sum == self.genotype[ell]:
                # Product of marginal potentials
                log_pot = 0.0
                for k in range(self.K):
                    a = config[k]
                    log_pot += self.log_psi_marginal[k, ell, a]
                log_potentials[state_idx] = log_pot
            # else: stays -np.inf (genotype constraint violated)
        
        # Normalize
        log_Z = logsumexp(log_potentials)
        if np.isfinite(log_Z):
            log_potentials -= log_Z
        
        return log_potentials


    def _compute_joint_transition_potentials(self, ell, prev_state_idx):
        """
        Compute NORMALIZED transition potentials from prev_state_idx at ℓ to all next configs at ℓ+1
        
        For each next config: potential = 𝟙[sum(next_config) == genotype[ℓ+1]] · ∏_k ψ^k_{ℓ,ℓ+1}(prev_config[k], next_config[k])
        Then normalize over all 2^K next configurations
        
        Args:
            ell: Source position index
            prev_state_idx: Previous state index
        
        Returns:
            log_potentials: array of shape (2^K,) with normalized log probabilities
        """
        n_states = self.n_joint_states
        log_potentials = np.full(n_states, -np.inf)
        
        prev_config = self._get_config(prev_state_idx)
        
        # Compute unnormalized potentials for all next configs
        for next_state_idx in range(n_states):
            next_config = self._get_config(next_state_idx)
            next_sum = np.sum(next_config)
            
            if next_sum == self.genotype[ell + 1]:
                # Product of individual transition potentials
                log_pot = 0.0
                for k in range(self.K):
                    a_prev = prev_config[k]
                    a_next = next_config[k]
                    log_pot += self.log_psi_transition[k, ell, a_prev, a_next]
                log_potentials[next_state_idx] = log_pot
            # else: stays -np.inf (genotype constraint violated)
        
        # Normalize
        log_Z = logsumexp(log_potentials)
        if np.isfinite(log_Z):
            log_potentials -= log_Z
        
        return log_potentials


    def _gradient_descent_update_potentials(self, haplotypes):
        """
        Update individual CRF potentials via gradient descent towards Viterbi solution
        
        Treats Viterbi solution as a Dirac distribution (point mass) and moves
        current potentials towards it using learning rate.
        
        For marginals: ψ^k_ℓ(a) ← (1 - α) · ψ^k_ℓ(a) + α · 𝟙[h^k_ℓ == a]
        For transitions: ψ^k_{ℓ,ℓ+1}(a₁,a₂) ← (1 - α) · ψ^k_{ℓ,ℓ+1}(a₁,a₂) + α · 𝟙[h^k_ℓ == a₁ and h^k_{ℓ+1} == a₂]
        
        Args:
            haplotypes: Array of shape (K, L) from Viterbi decoding
        """
        alpha = self.learning_rate
        
        # Update marginal potentials
        for k in range(self.K):
            for ell in range(self.L):
                # Get Viterbi allele at this position for this haplotype
                viterbi_allele = haplotypes[k, ell]
                
                # Create target Dirac distribution in log space
                log_target = np.full(2, -np.inf)
                log_target[viterbi_allele] = 0.0  # log(1) = 0
                
                # Gradient descent update in log space
                # Convert to probability space, interpolate, convert back
                current_probs = np.exp(self.log_psi_marginal[k, ell, :])
                target_probs = np.exp(log_target)
                
                # Weighted combination
                new_probs = (1 - alpha) * current_probs + alpha * target_probs
                
                # Ensure non-zero probabilities for numerical stability
                new_probs = np.maximum(new_probs, 1e-10)
                new_probs = new_probs / np.sum(new_probs)
                
                # Convert back to log space
                self.log_psi_marginal[k, ell, :] = np.log(new_probs)
        
        # Update transition potentials
        for k in range(self.K):
            for ell in range(self.L - 1):
                # Get Viterbi alleles at consecutive positions
                viterbi_a1 = haplotypes[k, ell]
                viterbi_a2 = haplotypes[k, ell + 1]
                
                # For each previous state, create target distribution
                for a1 in [0, 1]:
                    # Create target Dirac distribution
                    log_target = np.full(2, -np.inf)
                    if a1 == viterbi_a1:
                        log_target[viterbi_a2] = 0.0
                    # If a1 != viterbi_a1, target is uniform -inf (impossible)
                    
                    # Gradient descent update
                    current_probs = np.exp(self.log_psi_transition[k, ell, a1, :])
                    target_probs = np.exp(log_target)
                    
                    # Weighted combination
                    new_probs = (1 - alpha) * current_probs + alpha * target_probs
                    
                    # Ensure non-zero probabilities and normalize
                    new_probs = np.maximum(new_probs, 1e-10)
                    new_probs = new_probs / np.sum(new_probs)
                    
                    # Convert back to log space
                    self.log_psi_transition[k, ell, a1, :] = np.log(new_probs)
    
    
    def _compute_log_likelihood(self):
        """
        Compute total log-likelihood of current configuration
        
        Includes:
        1. Read likelihoods (how well reads match their assigned clusters)
        2. Genotype compatibility (how well inferred haplotypes sum to genotype)
        
        vectorized over positions
        not vectorized over reads
        """
        total_log_likelihood = 0.0
        #reads_likelihood = 0.0
        #genotype_likelihood = 0.0
        
        # Read likelihoods
        for j in range(self.n_reads):
            k = self.z[j]
            log_likelihood = self._compute_read_log_likelihood_standard(j, k)
            total_log_likelihood += log_likelihood
            #reads_likelihood += log_likelihood
        
        # Genotype compatibility term
        # Get inferred haplotypes via Viterbi decoding
        # Shape: (K, L)
        # haplotypes = self.get_haplotypes()
        
        # Sum of all haplotype alleles at each position
        # inferred_sums = np.sum(haplotypes, axis=0)  # Shape: (L,)
        
        # For each position, check compatibility with observed genotype
        # Hard constraint: log(1) if match, log(0)=-inf if mismatch
        # log_genotype_compat = np.where(
        #     inferred_sums == self.genotype,
        #     0.0,  # log(1) = 0
        #     -np.inf  # log(0) = -inf
        # )
        # total_log_likelihood += np.sum(log_genotype_compat)
        #genotype_likelihood += np.sum(log_genotype_compat)
        
        return total_log_likelihood#, reads_likelihood, genotype_likelihood
    
    
    def get_haplotypes(self):
        """
        Extract consensus haplotypes using joint phase Viterbi decoding
        
        Returns:
            haplotypes: Array of shape (K, L) with inferred haplotypes
        """
        return self._viterbi_decode_joint_phase()
    
    
    def get_posterior_summary(self):
        """
        Summarize posterior samples
        
        Returns:
            Dictionary with posterior statistics
        """
        if len(self.samples) == 0:
            return {}
        
        samples_array = np.array(self.samples)
        
        # Vectorized computation of assignment probabilities
        assignment_probs = np.zeros((self.n_reads, self.K))
        for k in range(self.K):
            assignment_probs[:, k] = np.mean(samples_array == k, axis=0)
        
        return {
            'assignment_probabilities': assignment_probs,
            'mean_cluster_sizes': np.mean(self.history['cluster_sizes'], axis=0),
            'final_cluster_sizes': self.history['cluster_sizes'][-1],
            'log_likelihood_trace': self.history['log_likelihood'],
            'n_switches_trace': self.history['n_switches']
        }



