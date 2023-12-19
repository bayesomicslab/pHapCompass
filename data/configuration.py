from utils.utils import phas_2_str
from algorithm.haplotype_assembly_helper import *


class Configuration:
    def __init__(self, ploidy, error_rate, epsilon, alleles):
        self.ploidy = ploidy
        self.error_rate = error_rate
        self.alleles = alleles
        self.epsilon = epsilon
        self.dynamic_config = {}
        self.global_phasings, self.global_likelihoods = self.initialization()

    def set_algorithm_type(self, algorithm_type):
        self.algorithm_type = algorithm_type

    def get_algorithm_type(self):
        return self.algorithm_type

    def set_setting(self, key, value):
        self.settings[key] = value

    def get_setting(self, key):
        return self.settings.get(key)

    def initialization2(self):
        # ploidy = 3
        # species_alleles = {0, 1}
        max_genotype = self.ploidy * max(list(self.alleles))
        all_genotypes = [str(gg) for gg in list(range(max_genotype + 1))]
        all_alleles = [str(gg) for gg in list(self.alleles)]
        global_phasings = {}
        # global_likelihoods = {}
        global_likelihoods2 = {}
        max_snp = 2
        # error_rate = 0.01
    
        for n_positions in range(2, max_snp + 1):
            global_phasings[n_positions] = {}
            # global_likelihoods[n_positions] = {}
            global_likelihoods2[n_positions] = {}
        
            this_genotypes = [''.join(list(com)) for com in list(itertools.product(all_genotypes, repeat=n_positions))]
            this_observations = [''.join(list(com)) for com in list(itertools.product(all_alleles, repeat=n_positions))]
            for gen in this_genotypes:
                global_phasings[n_positions][gen] = [phas_2_str(phas) for phas in
                                                     generate_phasings_ploidy(self.ploidy, gen)]
                pos_phasings = generate_phasings_ploidy(self.ploidy, gen)
                # global_likelihoods[n_positions][gen] = {}
                global_likelihoods2[n_positions][gen] = {}
            
                for phas in pos_phasings:
                    global_likelihoods2[n_positions][gen][phas_2_str(phas)] = {}
                    for obs in this_observations:
                        obs_np = np.array([int(po) for po in obs])
                        global_likelihoods2[n_positions][gen][phas_2_str(phas)][obs] = compute_likelihood(obs_np, phas,
                                                                                                          self.error_rate)
    
        return global_phasings, global_likelihoods2

    def initialization(self):
        max_genotype = self.ploidy * max(list(self.alleles))
        all_genotypes = [str(gg) for gg in list(range(max_genotype + 1))]
        all_alleles = [str(gg) for gg in list(self.alleles)]
        global_phasings = {}
        global_likelihoods = {}
        max_snp = 3
        # error_rate = 0.01
    
        for n_positions in range(2, max_snp + 1):
            global_phasings[n_positions] = {}
            # global_likelihoods[n_positions] = {}
            global_likelihoods[n_positions] = {}
        
            this_genotypes = [''.join(list(com)) for com in list(itertools.product(all_genotypes, repeat=n_positions))]
            this_observations = [''.join(list(com)) for com in list(itertools.product(all_alleles, repeat=n_positions))]
            for gen in this_genotypes:
                global_phasings[n_positions][gen] = [phas_2_str(phas) for phas in
                                                     generate_phasings_ploidy_long(self.ploidy, gen)]
                pos_phasings = generate_phasings_ploidy_long(self.ploidy, gen)
                # global_likelihoods[n_positions][gen] = {}
                global_likelihoods[n_positions][gen] = {}
            
                for phas in pos_phasings:
                    global_likelihoods[n_positions][gen][phas_2_str(phas)] = {}
                    for obs in this_observations:
                        obs_np = np.array([int(po) for po in obs])
                        global_likelihoods[n_positions][gen][phas_2_str(phas)][obs] = compute_likelihood(obs_np, phas,
                                                                                                          self.error_rate)
    
        return global_phasings, global_likelihoods
