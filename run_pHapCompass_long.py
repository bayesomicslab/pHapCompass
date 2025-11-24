import os, sys, time
from evaluations import *
from utils import *
from read_input import *
import pandas as pd
import numpy as np
import pickle
from pHapCompass_long.crf_gibbs_log_space_vectorized import *


def collect_mix_model_grid():
    columns = ['ploidy', 'Mutation Rate', 'coverage', 'sample', 'delta', 'epsilon', 'learning_rate', 'VER', 'MEC', 'Log-Likelihood', 'time']
    df = pd.DataFrame(columns=columns, index=range(700))
    # results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results_long_auto'
    results_path = '/archive/labs/aguiar/pHapCompass/results/results_long_auto'
    coverages = [3, 5, 10, 20, 40] # 14400
    ploidies = [2, 3, 4, 6] # 2, 4, 6, 8
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(20)
    deltas = [1, 2, 5]
    epsilons = [0.00001, 0.0001]
    learning_rates = [0.01]
    counter = 0
    for ploidy in ploidies:
        for mr in mut_rates:
            for coverage in coverages:
                for sample in samples:
                    for delta in deltas:
                        for epsilon in epsilons:
                            for lr in learning_rates:
                                this_results = os.path.join(results_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(2))
                                result_file = f'delta{delta}_ep{epsilon}_lr{lr}.pkl.gz'
                                res = os.path.join(this_results, result_file)
                                
                                if os.path.exists(res):
                                    # all_results.append(os.path.join(this_results, result_file))
                                    with gzip.open(res, "rb") as f:
                                        result_pkl = pickle.load(f)
                                    this_delta = result_pkl['delta']
                                    this_epsilon = result_pkl['epsilon']
                                    this_lr = result_pkl['learning_rate']
                                    time = result_pkl['time']
                                    sampler = result_pkl['sampler']
                                    vec = result_pkl['vector_error_rate']
                                    mec = result_pkl['mec']
                                    # predicted_haplotypes = result_pkl['predicted_haplotypes']
                                    # true_haplotypes = result_pkl['true_haplotypes']
                                    # ll_mean = np.mean([l[0] for l in sampler.history['log_likelihood'][-10:]])
                                    ll_mean = np.mean([l for l in sampler.history['log_likelihood'][-10:]])
                                    # stop
                                    df.loc[counter, :] = ploidy, mr, coverage, str(sample).zfill(2), this_delta, this_epsilon, this_lr, vec, mec, ll_mean, time
                                    counter += 1
    df = df.dropna().reset_index(drop=True)
    # df.to_csv()
    df['Method'] = 'pHapCompass-long'
    df.to_csv('/labs/Aguiar/pHapCompass/results/pHapCompass_long_auto.csv')


def collect_mix_model_auto():
    columns = ['Method', 'Data', 'Ploidy', 'Subgenome', 'Mutation Rate', 'Coverage', 'Sample', 'Metric', 'Value']
    df = pd.DataFrame(columns=columns, index=range(700))
    # results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results_long_auto'
    results_path = '/archive/labs/aguiar/pHapCompass/results/results_long_auto'
    coverages = [3, 5, 10, 20, 40] # 14400
    ploidies = [2, 3, 4, 6] # 2, 4, 6, 8
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(20)
    deltas = [1, 2, 5]
    epsilons = [0.00001, 0.0001]
    learning_rates = [0.01]
    counter = 0
    for ploidy in ploidies:
        for mr in mut_rates:
            for coverage in coverages:
                for sample in samples:
                    for delta in deltas:
                        for epsilon in epsilons:
                            for lr in learning_rates:
                                this_results = os.path.join(results_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(2))
                                result_file = f'delta{delta}_ep{epsilon}_lr{lr}.pkl.gz'
                                res = os.path.join(this_results, result_file)
                                
                                if os.path.exists(res):
                                    # all_results.append(os.path.join(this_results, result_file))
                                    with gzip.open(res, "rb") as f:
                                        result_pkl = pickle.load(f)
                                    this_delta = result_pkl['delta']
                                    this_epsilon = result_pkl['epsilon']
                                    this_lr = result_pkl['learning_rate']
                                    time = result_pkl['time']
                                    sampler = result_pkl['sampler']
                                    vec = result_pkl['vector_error_rate']
                                    mec = result_pkl['mec']
                                    # predicted_haplotypes = result_pkl['predicted_haplotypes']
                                    # true_haplotypes = result_pkl['true_haplotypes']
                                    # ll_mean = np.mean([l[0] for l in sampler.history['log_likelihood'][-10:]])
                                    ll_mean = np.mean([l for l in sampler.history['log_likelihood'][-10:]])
                                    # stop
                                    df.loc[counter, :] = ploidy, mr, coverage, str(sample).zfill(2), this_delta, this_epsilon, this_lr, vec, mec, ll_mean, time
                                    counter += 1
    df = df.dropna().reset_index(drop=True)
    # df.to_csv()
    df['Method'] = 'pHapCompass-long'
    df.to_csv('/labs/Aguiar/pHapCompass/results/pHapCompass_long_auto.csv')


def make_input_for_mix_model_auto():
    # data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/dataset/simulated_data_auto_long'
    # output_dir = '/mnt/research/aguiarlab/proj/HaplOrbit/inputs_paper_long'
    # results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results_long_auto'
    # reference_path = '/mnt/research/aguiarlab/proj/HaplOrbit/reference/simulated_haplotypes/auto'

    data_path = '/labs/Aguiar/pHapCompass/dataset/simulated_data_auto_long'
    output_dir = '/labs/Aguiar/pHapCompass/scripts/long_auto/input'
    results_path = '/archive/labs/aguiar/pHapCompass/results/results_long_auto'
    reference_path = '/labs/Aguiar/pHapCompass/reference/simulated_haplotypes/auto'


    coverages = [3, 5, 10, 20, 40] # 1200
    ploidies = [6] #[2, 3, 4, 6] #, 6] # 2, 4, 6, 8
    mut_rates = [0.001, 0.005, 0.01]
    samples = range(6)
    inputs = []
    # n_runs = 0
    deltas = [5]
    epsilons = [0.00001]
    learning_rates = [0.02]
    # expected = [f'delta{delta}_ep{epsilon}_gen_coef{gen_coef}.pkl.gz' for delta in deltas for epsilon in epsilons for gen_coef in gen_coefficients]
    for ploidy in ploidies:
        for mr in mut_rates:
            for coverage in coverages:
                for sample in samples:
                    for delta in deltas:
                        for epsilon in epsilons:
                            for lr in learning_rates:
                                this_results = os.path.join(results_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(2))
                                if not os.path.exists(this_results):
                                    os.makedirs(this_results, exist_ok=True)
                                frag_path = os.path.join(data_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(2) + '_mbq4.frag')
                                bam_path = os.path.join(data_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(2) + '.bam')
                                vcf_path = os.path.join(reference_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2), 'Chr1.vcf')
                                genotype_path = os.path.join(reference_path, 'ploidy_' + str(ploidy), 'mut_' + str(mr), str(sample).zfill(2), 'genotype.csv')
                                # ress = [f for f in os.listdir(this_results) if 'pkl.gz' in f]
                                # not_done = [e for e in expected if e not in ress]
                                # n_runs += len(not_done)
                                this_inp = [frag_path, bam_path, vcf_path, genotype_path, ploidy, this_results, delta, epsilon, lr]
                                # if len(ress) < 16:
                                #    inputs.append(this_inp)
                                result_file = f'delta{delta}_ep{epsilon}_lr{lr}.pkl.gz'
                                if not os.path.exists(os.path.join(this_results, result_file)):
                                    inputs.append(this_inp)
                                    print(os.path.join(this_results, result_file))
                                # else:
                                #     print(os.path.join(this_results, result_file))
    print(len(inputs))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, inp in enumerate(inputs):
        input_file = os.path.join(output_dir, f"input_{i}.pkl")
        with open(input_file, "wb") as f:
            pickle.dump(inp, f)
    print(f"Saved {len(inputs)} inputs to {output_dir}")


def make_input_for_mix_model_allo():
    data_path = '/labs/Aguiar/pHapCompass/dataset/simulated_data_allo_long_new'
    output_dir = '/labs/Aguiar/pHapCompass/scripts/long_allo/input'
    results_path = '/labs/Aguiar/pHapCompass/results/long_allo_new/'
    reference_path = '/labs/Aguiar/pHapCompass/reference/simulated_haplotypes/allo'
    inputs = []
    # n_runs = 0
    deltas = [5]
    epsilons = [0.00001]
    learning_rates = [0.02]
    # subgenome_configs = [0.01, 0.05]
    # within_mutation_rates = [0.001, 0.005, 0.0075, 0.01, 0.05]
    subgenome_configs = ['0.0005', '0.0001']
    within_mutation_rates = ['0.00005', '0.0001']
    ploidies = [3, 4, 6]
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C', 8:'8ploid_4A+4B'}
    coverages = [5, 10, 20, 40]
    samples = range(10)
    for sgc in subgenome_configs:
        for mr in within_mutation_rates:
            for ploidy in ploidies:
                for coverage in coverages:
                    for sample in samples:
                        for delta in deltas:
                            for epsilon in epsilons:
                                for lr in learning_rates:
                                    this_results = os.path.join(results_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(1))
                                    # stop
                                    if not os.path.exists(this_results):
                                        os.makedirs(this_results, exist_ok=True)
                                    
                                    frag_path = os.path.join(data_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(1) + '_mbq4.frag')
                                    bam_path = os.path.join(data_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(1) + '.bam')                        
                                    vcf_path = os.path.join(reference_path, 'subgenome_config_mut' + str(sgc), ploidy_folder_name[ploidy], 'within_mut_' + str(mr), str(sample).zfill(1), 'Chr1.vcf')
                                    genotype_path = os.path.join(reference_path, 'subgenome_config_mut' + str(sgc), ploidy_folder_name[ploidy], 'within_mut_' + str(mr), str(sample).zfill(1), 'genotype.csv')
                                    result_file = os.path.join(this_results, f'delta{delta}_ep{epsilon}_lr{lr}.pkl.gz')
                                    this_inp = [frag_path, bam_path, vcf_path, genotype_path, ploidy, this_results, delta, epsilon, lr]
                                    if not os.path.exists(result_file):
                                        inputs.append(this_inp)
    print(len(inputs))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, inp in enumerate(inputs):
        input_file = os.path.join(output_dir, f"input_{i}.pkl")
        with open(input_file, "wb") as f:
            pickle.dump(inp, f)
    print(f"Saved {len(inputs)} inputs to {output_dir}")


def make_input_for_mix_model_allo_short_reads():
    data_path = '/labs/Aguiar/pHapCompass/dataset/simulated_data_allo_short_new'
    output_dir = '/labs/Aguiar/pHapCompass/scripts/long_model_short_allo/input'
    results_path = '/labs/Aguiar/pHapCompass/results/long_model_short_allo/'
    reference_path = '/labs/Aguiar/pHapCompass/reference/simulated_haplotypes/allo'
    inputs = []
    # n_runs = 0
    deltas = [5]
    epsilons = [0.00001]
    learning_rates = [0.02]
    # subgenome_configs = [0.01, 0.05]
    # within_mutation_rates = [0.001, 0.005, 0.0075, 0.01, 0.05]
    subgenome_configs = ['0.0005', '0.0001']
    within_mutation_rates = ['0.00005', '0.0001']
    ploidies = [3, 4, 6]
    ploidy_folder_name = {3:'3ploid_2A+1B', 4: '4ploid_2A+2B', 6: '6ploid_2A+2B+2C', 8:'8ploid_4A+4B'}
    coverages = [5, 10, 20, 40]
    samples = range(10)
    for sgc in subgenome_configs:
        for mr in within_mutation_rates:
            for ploidy in ploidies:
                for coverage in coverages:
                    for sample in samples:
                        for delta in deltas:
                            for epsilon in epsilons:
                                for lr in learning_rates:
                                    this_results = os.path.join(results_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(1))
                                    # stop
                                    if not os.path.exists(this_results):
                                        os.makedirs(this_results, exist_ok=True)
                                    
                                    frag_path = os.path.join(data_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(1) + '.frag')
                                    bam_path = os.path.join(data_path, 'ploidy_' + str(ploidy), 'subgenome_' + str(sgc), 'mut_' + str(mr), 'cov_' + str(coverage), str(sample).zfill(1) + '.bam')                        
                                    vcf_path = os.path.join(reference_path, 'subgenome_config_mut' + str(sgc), ploidy_folder_name[ploidy], 'within_mut_' + str(mr), str(sample).zfill(1), 'Chr1.vcf')
                                    genotype_path = os.path.join(reference_path, 'subgenome_config_mut' + str(sgc), ploidy_folder_name[ploidy], 'within_mut_' + str(mr), str(sample).zfill(1), 'genotype.csv')
                                    result_file = os.path.join(this_results, f'delta{delta}_ep{epsilon}_lr{lr}.pkl.gz')
                                    this_inp = [frag_path, bam_path, vcf_path, genotype_path, ploidy, this_results, delta, epsilon, lr]
                                    if not os.path.exists(result_file):
                                        inputs.append(this_inp)
    print(len(inputs))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, inp in enumerate(inputs):
        input_file = os.path.join(output_dir, f"input_{i}.pkl")
        with open(input_file, "wb") as f:
            pickle.dump(inp, f)
    print(f"Saved {len(inputs)} inputs to {output_dir}")


def run_pHapCompass_long(inp):

    this_frag_path, bam_file_path, vcf_file_path, genotype_path, ploidy, results_path, delta, epsilon, learning_rate = inp

    if not os.path.exists(results_path):
        os.makedirs(results_path, exist_ok=True)

    result_name = f'delta{delta}_ep{epsilon}_lr{learning_rate}.pkl.gz'
    print(f'Working on {os.path.join(results_path, result_name)} ...')


    if genotype_path is None:
        os.makedirs(results_path, exist_ok=True)
        genotype_path = os.path.join(results_path, "genotype.csv")
        # writes per-variant GT columns haplotype_0..haplotype_{K-1}
        vcf_gt_to_csv(vcf_file_path, genotype_path, ploidy=ploidy)

    gen_df = pd.read_csv(genotype_path)
    g = gen_df.sum(axis=1).to_numpy(dtype=np.int16)   # length = #SNP positions

    N = len(gen_df)

    start_time = time.time()

    # ------------------- 2) Read fragments (.frag) into sparse CSR ------------------------
    cfg  = InputConfigSparse(data_path=this_frag_path, genotype_path=genotype_path, ploidy=ploidy)
    frag = SparseFragment(cfg, positions_from_genotype=list(range(N)))
    M = frag.csr  
    reads = M.toarray().astype(float)
    reads[reads==0] = np.nan
    reads[reads==1] = 0
    reads[reads==2] = 1

    sampler = HaplotypeGibbsSampler(reads=reads, genotype=g, K=ploidy, epsilon=epsilon, delta=delta, learning_rate=learning_rate)

    sampler.fit(n_iterations=500, burn_in=0, verbose=False)

    end_time = time.time()

    elapsed_time = round(end_time - start_time, 2)

    predicted_haplotypes = sampler.history["phase"][-1]
    true_haplotypes = pd.read_csv(genotype_path).T.to_numpy()
    block_ids = np.zeros(N)
    ver = vector_error_wrapper(true_haplotypes, predicted_haplotypes, block_ids)
    list_of_reads_M = build_read_list_from_M(M)
    mec = mec_full(predicted_haplotypes, block_ids, list_of_reads_M, probabalistic=False)
    # mecpt = mec_full(predicted_haplotypes, block_ids, list_of_reads_M, probabalistic=True)
    geo_mec = mec_full_geometric_penalty(predicted_haplotypes, block_ids, list_of_reads_M)
    results = {"sampler": sampler, 'time': elapsed_time, 'vector_error_rate': ver, 'mec': mec, 'geo_mec':geo_mec, 
    'predicted_haplotypes': predicted_haplotypes, 'true_haplotypes': true_haplotypes, 'delta': delta, 
    'epsilon': epsilon, 'learning_rate': learning_rate}
    with gzip.open(os.path.join(results_path, result_name), "wb") as f:
        pickle.dump(results, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f'[Done] {os.path.join(results_path, result_name)}.')

    # # later (must have the class definition importable in the same module path)
    # with gzip.open(os.path.join(results_path, result_name), "rb") as f:
    #     res = pickle.load(f)


def inspect_mix_model_results():
    results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results_long/ploidy_3/mut_0.001/cov_10'
    this_frag_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_auto_long/ploidy_3/mut_0.001/cov_10/00.frag'

    ploidy = 3

    results = sorted([os.path.join(results_path, f) for f in os.listdir(results_path) if '.pkl' in f])
    for res in results:
        with gzip.open(res, "rb") as f:
            result_pkl = pickle.load(f)
              


        delta = res.split('/')[-1].split('_')[0].split('delta')[1]
        epsilon = res.split('/')[-1].split('_')[1].split('ep')[1]
        gen_coef = res.split('/')[-1].split('_')[3].split('.pkl')[0].split('coef')[1]
        print('-----------------------------------------------------------------')
        print('delta:', delta, ', epsilon:', epsilon, ', gen coef:', gen_coef)
        sampler = result_pkl["sampler"]
        time = result_pkl["time"]

        likelihood_trace = sampler.history["log_likelihood"]

        iterations = range(len(likelihood_trace))

        predicted_haplotypes = sampler.get_haplotypes()
        genotype_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results_long/ploidy_3/mut_0.001/cov_10/genotype.csv'
        true_haplotypes = pd.read_csv(genotype_path).to_numpy().T
        acc = 0

        num_rows = predicted_haplotypes.shape[0]
        perms = itertools.permutations(range(num_rows))
        for hap1 in  [true_haplotypes[list(p), :] for p in perms]:
            temp_acc = np.sum(hap1 == predicted_haplotypes)/(hap1.shape[1] * hap1.shape[0])
            if temp_acc > acc:
                acc = temp_acc

        N = true_haplotypes.shape[1]
        cfg  = InputConfigSparse(data_path=this_frag_path, genotype_path=genotype_path, ploidy=ploidy)
        frag = SparseFragment(cfg, positions_from_genotype=list(range(N)))
        M = frag.csr  
        # reads = M.toarray().astype(float)
        
        list_of_reads_M = build_read_list_from_M(M)

        mec = mec_full(predicted_haplotypes, true_haplotypes, list_of_reads_M)


        print('accuracy:', acc)
        print('log likelihood:', sampler.history['log_likelihood'][-1])
        print('MEC:', mec)

        plt.figure(figsize=(10, 6))
        plt.plot(iterations, likelihood_trace)
        plt.xlabel("Iteration")
        plt.ylabel("log_likelihood")
        plt.title("Log Likelihood Trace over Gibbs Sampling Iterations")
        plt.grid(True)
        plt.show()


def unpack_results_mix_model(res):
    res = '/mnt/research/aguiarlab/proj/HaplOrbit/results_long_auto/ploidy_2/mut_0.001/cov_3/00/delta10_ep0.0001_gen_coef2.pkl.gz'
    with gzip.open(res, "rb") as f:
        result_pkl = pickle.load(f)
        
    delta = result_pkl['delta']
    epsilon = result_pkl['epsilon']
    gen_coef = result_pkl['gen_coef']
    time = result_pkl['time']
    sampler = result_pkl['sampler']
    vector_error_rate = result_pkl['vector_error_rate']
    mec = result_pkl['mec']
    predicted_haplotypes = result_pkl['predicted_haplotypes']
    true_haplotypes = result_pkl['true_haplotypes']

        
def run_pHapcompass_from_input(input_file):
    with open(input_file, "rb") as f:
        inp = pickle.load(f)

    mix_model(inp)



if __name__ == '__main__':

    # mix_model()

    if len(sys.argv) != 2:
        print("Usage: python3 simulator_paper.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    run_pHapcompass_from_input(input_file)