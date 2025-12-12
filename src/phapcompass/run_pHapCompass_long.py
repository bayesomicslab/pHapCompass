from operator import index
from .evaluations import *
from .utils import *
from .read_input import *
import pandas as pd
import numpy as np
from .long.crf_gibbs_log_space_vectorized import *



def run_pHapCompass_long(args):

    frag = args.frag_path
    vcf  = args.vcf_path
    ploidy   = args.ploidy
    delta, epsilon, learning_rate = args.delta, args.epsilon, args.learning_rate
    gen_df = args.genotype

    g = gen_df.sum(axis=1).to_numpy(dtype=np.int16)   # length = #SNP positions

    n_snps = len(gen_df)

    # ------------------- 2) Read fragments (.frag) into sparse CSR ------------------------
    cfg  = InputConfigSparse(data_path=frag, genotype_path=gen_df, ploidy=ploidy)
    frag = SparseFragment(cfg, positions_from_genotype=list(range(n_snps)))
    data_matrix = frag.csr  
    reads = data_matrix.toarray().astype(float)
    reads[reads==0] = np.nan
    reads[reads==1] = 0
    reads[reads==2] = 1

    sampler = HaplotypeGibbsSampler(reads=reads, genotype=g, K=ploidy, epsilon=epsilon, delta=delta, learning_rate=learning_rate)
    if args.uncertainty is not None: 
        block_ids = [np.zeros(n_snps) for _ in range(args.uncertainty)]
        probabilities = [] # fix it
        sampler.fit(n_iterations=500, burn_in=0, verbose=False)
        solutions = sampler._ffbs_decode_joint_phase_multiple(args.uncertainty)
        predicted_haplotypes, probabilities = map(list, zip(*solutions))

        predicted_haplotypes = [pd.DataFrame(data=sol, index=['haplotype_' + str(i + 1) for i in range(ploidy)]) for sol in predicted_haplotypes]


    else:

        sampler.fit(n_iterations=500, burn_in=0, verbose=False)
        predicted_haplotype = sampler.history["phase"][-1]
        block_id = np.zeros(n_snps)
        probability = 1

        predicted_haplotypes = [pd.DataFrame(data=predicted_haplotype, index=['haplotype_' + str(i + 1) for i in range(ploidy)])]
        # print(predicted_haplotypes)
        block_ids = [block_id]
        probabilities = [probability]  # or [likelihood] if you later want LK even without uncertainty

    norm_scores = list(probabilities/np.sum(probabilities))

    write_phased_vcf(vcf, args.result_path, predicted_haplotypes, block_ids, norm_scores, ploidy)

