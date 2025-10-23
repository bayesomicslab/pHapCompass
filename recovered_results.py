import os
import sys
import pickle
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from utils.utils import sort_nodes, str_2_phas_1
import numpy as np
import itertools
from collections import defaultdict
import pysam
from evaluation.evaluation import *


def prepare_results_ismb():
    agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    sim_data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test'
    metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    contigs = ['100'] 
    ploidies = ['6']
    coverages = ['10', '50', '100']
    results_dfs = []
    for contig in contigs:
        for ploidy in ploidies:
            for coverage in coverages:
                # results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results')
                results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results_likelihood')
                if os.path.exists(results_path):
                    samples = [f for f in os.listdir(results_path) if 'FFBS' in f]
                    this_result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value', 'length_phased'], index=range(len(samples)*len(metrics)))
                    this_result_df['Contig'] = contig
                    this_result_df['Ploidy'] = ploidy
                    this_result_df['Coverage'] = coverage
                    this_result_df['Method'] = 'pHapCompass'                
                    counter = 0
                    for sample in samples:
                        sample_result = os.path.join(results_path, sample)
                        with open(sample_result, 'rb') as f:
                            this_results = pickle.load(f)
                        sample_name = sample.split('.pkl')[0].split('_')[-1]
                        evals = this_results['evaluation']
                        phased_snp = this_results['predicted_haplotypes'].shape[1]
                        for metric in metrics:
                            this_result_df.loc[counter, 'Sample'] = sample_name
                            this_result_df.loc[counter, 'Metric'] = metric
                            this_result_df.loc[counter, 'Value'] = evals[metric]
                            this_result_df.loc[counter, 'length_phased'] = phased_snp
                            counter += 1
                    results_dfs.append(this_result_df)

    results_df = pd.concat(results_dfs, ignore_index=True)
    results_df.to_csv(os.path.join(agg_results_path, 'pHapCompass_results_simulated_data_test.csv'), index=False)


    # server3_results = pd.read_csv(os.path.join(agg_results_path, 'sim_awri_results_likelihood10_100_server3.csv'))
    beagle_results = pd.read_csv(os.path.join(agg_results_path, 'pHapCompass_results_simulated_data_test.csv'))
    # ser3_beagle = pd.concat([server3_results, beagle_results], ignore_index=True)
    # ser3_beagle_non_redundant = ser3_beagle.groupby(["Contig", "Ploidy", "Coverage"], group_keys=False).apply(lambda group: group.drop_duplicates(subset="Sample")).reset_index(drop=True)
    # ser3_beagle_non_redundant = ser3_beagle_non_redundant[ser3_beagle_non_redundant['Ploidy'] != 8].reset_index(drop=True)
    # non_redundant_sorted_df = (
    #     ser3_beagle_non_redundant.groupby(["Contig", "Ploidy", "Coverage"], group_keys=False)
    #     .apply(lambda group: group.drop_duplicates(subset="Sample").sort_values(by="Sample"))
    #     .reset_index(drop=True)
    # )
    
    group_sizes = beagle_results.groupby(["Contig", "Ploidy", "Coverage"]).size()
    small_groups = group_sizes[group_sizes < 100]


    capital = {'method': 'Method', 'contig': 'Contig', 'ploidy': 'Ploidy', 'coverage': 'Coverage', 'sample': 'Sample', 'metric': 'Metric', 'value': 'Value'}
    hpopg_sim10 = pd.read_csv(os.path.join(agg_results_path, 'hpopg_sim_awri_results10_346.csv'))
    hpopg_sim100 = pd.read_csv(os.path.join(agg_results_path, 'hpopg_sim_awri_results100_3_4_6.csv'))
    hpopg_sim10 = hpopg_sim10.rename(columns=capital)
    hpopg_sim100 = hpopg_sim100.rename(columns=capital)
    whatshapp_10 = pd.read_csv(os.path.join(agg_results_path, 'whatshap_sim_awri_results_10_3_4_6.csv'))
    whatshapp_100 = pd.read_csv(os.path.join(agg_results_path, 'whatshap_sim_awri_results_100_3_4_6.csv'))
    whatshapp_10 = whatshapp_10.rename(columns=capital)
    whatshapp_100 = whatshapp_100.rename(columns=capital)

    all_results = pd.concat([beagle_results, hpopg_sim10, hpopg_sim100, whatshapp_10, whatshapp_100], ignore_index=True)
    all_results['Contig'] = all_results['Contig'].astype(int)
    all_results['Ploidy'] = all_results['Ploidy'].astype(int)
    all_results['Coverage'] = all_results['Coverage'].astype(int)
    all_results['length_phased'] = all_results['length_phased'].astype(int)

    all_results = all_results[all_results['Ploidy'] != 8].reset_index(drop=True)
    all_results = all_results[all_results['Metric'].isin(['vector_error_rate', 'mec', 'mismatch_error'])].reset_index(drop=True)
    all_results = all_results.sort_values(by=['Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric']).reset_index(drop=True)
    methods_order = ['pHapCompass', 'HPoP-G', 'WhatsHap']
    all_results = all_results.sort_values(by=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric'],
        key=lambda col: col.map({method: i for i, method in enumerate(methods_order)})).reset_index(drop=True)

    all_results.to_csv(os.path.join(agg_results_path, 'all_methods_results.csv'), index=False)

    # results_df_100 = results_df[results_df['contig'] == '100'].reset_index(drop=True)
    results_df_10 = all_results[all_results['Contig'] == 10].reset_index(drop=True)
    results_df_10.to_csv(os.path.join(agg_results_path, 'all_methods_results_10.csv'), index=False)

    results_df_100 = all_results[all_results['Contig'] == 100].reset_index(drop=True)
    results_df_100.to_csv(os.path.join(agg_results_path, 'all_methods_results_100.csv'), index=False)


def prepare_results():
    agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    sim_data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    contigs = ['100'] 
    ploidies = ['3', '4', '6', '8']
    coverages = ['10', '30', '50', '70', '100']
    results_dfs = []
    for contig in contigs:
        for ploidy in ploidies:
            for coverage in coverages:
                # results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results')
                # results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results_likelihood')
                results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results_FFBS_Multiple')
                if os.path.exists(results_path):
                    samples = [f for f in os.listdir(results_path) if 'FFBS' in f]
                    this_result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value', 'length_phased'], index=range(len(samples)*len(metrics)))
                    this_result_df['Contig'] = contig
                    this_result_df['Ploidy'] = ploidy
                    this_result_df['Coverage'] = coverage
                    this_result_df['Method'] = 'pHapCompass'                
                    counter = 0
                    for sample in samples:
                        sample_result = os.path.join(results_path, sample)
                        with open(sample_result, 'rb') as f:
                            this_results = pickle.load(f)
                        sample_name = sample.split('.pkl')[0].split('_')[-1]
                        evals = this_results['evaluation']
                        phased_snp = this_results['predicted_haplotypes'].shape[1]
                        for metric in metrics:
                            this_result_df.loc[counter, 'Sample'] = sample_name
                            this_result_df.loc[counter, 'Metric'] = metric
                            this_result_df.loc[counter, 'Value'] = evals[metric]
                            this_result_df.loc[counter, 'length_phased'] = phased_snp
                            counter += 1
                    results_dfs.append(this_result_df)

    results_df = pd.concat(results_dfs, ignore_index=True)
    results_df.to_csv(os.path.join(agg_results_path, 'pHapCompass_results_simulated_data_NA12878.csv'), index=False)


    # server3_results = pd.read_csv(os.path.join(agg_results_path, 'sim_awri_results_likelihood10_100_server3.csv'))
    phapcompass_results = pd.read_csv(os.path.join(agg_results_path, 'pHapCompass_results_simulated_data_NA12878.csv'))
    # whatshapp_results = pd.read_csv(os.path.join(agg_results_path, 'whatshap_results_simulated_data_test.csv'))
    # hpopg_results = pd.read_csv(os.path.join(agg_results_path, 'hpop_results_simulated_data_test.csv'))
    whatshapp_results346 = pd.read_csv(os.path.join(agg_results_path, 'whatshap_results_simulated_data_NA12878_346.csv'))
    hpopg_results346 = pd.read_csv(os.path.join(agg_results_path, 'hpop_results_simulated_data_NA12878_346.csv'))
    whatshapp_results8 = pd.read_csv(os.path.join(agg_results_path, 'whatshap_results_simulated_data_NA12878_ploidy8.csv'))
    hpopg_results8 = pd.read_csv(os.path.join(agg_results_path, 'hpop_results_simulated_data_NA12878_8.csv'))


    # all_results = pd.concat([phapcompass_results, whatshapp_results, hpopg_results], ignore_index=True)
    all_results = pd.concat([phapcompass_results, whatshapp_results346, whatshapp_results8, hpopg_results346, hpopg_results8], ignore_index=True)
    all_results['Contig'] = all_results['Contig'].astype(int)
    all_results['Ploidy'] = all_results['Ploidy'].astype(int)
    all_results['Coverage'] = all_results['Coverage'].astype(int)
    all_results['length_phased'] = all_results['length_phased'].astype(int)

    all_results = all_results.sort_values(by=['Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric']).reset_index(drop=True)
    methods_order = ['pHapCompass', 'WhatsHap' , 'HPoP-G']
    all_results = all_results.sort_values(by=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric'],
        key=lambda col: col.map({method: i for i, method in enumerate(methods_order)})).reset_index(drop=True)

    all_results.to_csv(os.path.join(agg_results_path, 'all_methods_results_simulated_data_NA12878.csv'), index=False)


def prepare_results_partial_haplotypes_pHapcompass346():
    agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    sim_data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    contigs = ['100'] 
    ploidies = ['3', '4', '6', '8']
    coverages = ['10', '30', '50', '70', '100']
    results_dfs = []
    for contig in contigs:
        for ploidy in ploidies:
            for coverage in coverages:
                # results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results')
                # results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results_LBP_FFBS_Single')
                results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results_FFBS_Multiple')
                if os.path.exists(results_path):
                    samples = [f for f in os.listdir(results_path) if 'FFBS' in f]
                    this_result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(samples)*2))
                    this_result_df['Contig'] = contig
                    this_result_df['Ploidy'] = ploidy
                    this_result_df['Coverage'] = coverage
                    this_result_df['Method'] = 'pHapCompass + M.V.'                
                    counter = 0
                    for sample in samples:
                        # print(contig, ploidy, coverage, sample) 
                        sample_result = os.path.join(results_path, sample)
                        with open(sample_result, 'rb') as f:
                            this_results = pickle.load(f)
                        
                        sample_name = sample.split('.pkl')[0].split('_')[-1]
                        
                        # H = this_results['predicted_haplotypes'].to_numpy()
                        # H_star = this_results['true_haplotypes'].to_numpy()
                        # vector_error_rate, _, _ = compute_vector_error_rate_with_missing_positions(H_star, H)
                        # length_phased = np.sum(~np.isnan(np.array(H, dtype=np.float64)).any(axis=0))
                        vector_error_rate = this_results['evaluation']['vector_error_rate']
                        length_phased = this_results['length_phased']
                        this_result_df.loc[counter, 'Sample'] = sample_name
                        this_result_df.loc[counter, 'Metric'] = 'Vector Error Rate'
                        this_result_df.loc[counter, 'Value'] = vector_error_rate
                        counter += 1
                        this_result_df.loc[counter, 'Sample'] = sample_name
                        this_result_df.loc[counter, 'Metric'] = '# Phased Variants'
                        this_result_df.loc[counter, 'Value'] = length_phased
                        counter += 1
                    results_dfs.append(this_result_df)

    results_df = pd.concat(results_dfs, ignore_index=True)
    results_df.to_csv(os.path.join(agg_results_path, 'pHapCompass_LBP_results_simulated_data_NA12878.csv'), index=False)


def prepare_results_partial_haplotypes_pHapcompass8():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapCompass_results_simulated_data_NA12878_8.csv'
    contig_lens = [100]
    ploidies = [8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    # metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*2))
    counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                results_path = os.path.join(this_cov_path, 'results_likelihood')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, 'dict_' + str(rd).zfill(2) + '.pkl')
                    with open(result_file, "rb") as f:
                        evals = pickle.load(f)
                    result_df.loc[counter, 'Sample'] = evals['Sample']
                    result_df.loc[counter, 'Metric'] = 'Vector Error Rate'
                    result_df.loc[counter, 'Value'] = evals['vector_error_rate']
                    result_df.loc[counter, 'Contig'] = contig_len
                    result_df.loc[counter, 'Ploidy'] = ploidy
                    result_df.loc[counter, 'Coverage'] = coverage
                    counter += 1
                    result_df.loc[counter, 'Contig'] = contig_len
                    result_df.loc[counter, 'Ploidy'] = ploidy
                    result_df.loc[counter, 'Coverage'] = coverage
                    result_df.loc[counter, 'Sample'] = evals['Sample']
                    result_df.loc[counter, 'Metric'] = '# Phased Variants'
                    result_df.loc[counter, 'Value'] = evals['length_phased']
                    counter += 1

    result_df['Method'] = 'pHapCompass'
    result_df.to_csv(output_path, index=False)
    print("pHapcompass results collected")


def collect_ground_truth_results():
    agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    sim_data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    # metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    contigs = ['100'] 
    ploidies = ['3', '4', '6', '8']
    coverages = ['10', '30', '50', '70', '100']
    results_dfs = []
    for contig in contigs:
        for ploidy in ploidies:
            genotype_path = os.path.join(sim_data_path, 'contig_{}/ploidy_{}/haplotypes.csv'.format(contig, ploidy))
            for coverage in coverages:
                # results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results')
                results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results_likelihood')
                if os.path.exists(results_path):
                    samples = [f for f in os.listdir(results_path) if 'FFBS' in f]
                    this_result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Vector Error Rate', 'FFBS Accuracy'], index=range(len(samples)))
                    this_result_df['Contig'] = contig
                    this_result_df['Ploidy'] = ploidy
                    this_result_df['Coverage'] = coverage
                    this_result_df['Method'] = 'pHapCompass'                
                    counter = 0
                    for sample in samples:
                        sample_result = os.path.join(results_path, sample)
                        with open(sample_result, 'rb') as f:
                            this_results = pickle.load(f)
                        sample_name = sample.split('.pkl')[0].split('_')[-1]
                        # graph_samples = this_results['samples']
                        evals = this_results['evaluation']
                        # ffbs_acc = evaulate_ffbs_acc_sample(genotype_path, graph_samples, int(ploidy))
                        this_result_df.loc[counter, 'Vector Error Rate'] = evals['vector_error_rate']
                        this_result_df.loc[counter, 'FFBS Accuracy'] = evals['ffbs_acc']
                        this_result_df.loc[counter, 'Sample'] = sample_name
                        counter += 1

                    results_dfs.append(this_result_df)

    results_df = pd.concat(results_dfs, ignore_index=True)
    results_df.to_csv(os.path.join(agg_results_path, 'simulated_data_NA12878_ground_truth_ffbs_vector_error_rate.csv'), index=False)


def collect_pairwise_samples_accuracy():
    save_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pairwise_samples_accuracy.csv'
    sim_data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    # metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    contigs = ['100'] 
    ploidies = ['3', '4', '6', '8']
    coverages = ['10', '30', '50', '70', '100']
    results_df = pd.DataFrame(columns=['Contig', 'Ploidy', 'Coverage', 'Sample', 'Sampling', 'Value'], index=range(100*len(contigs)*len(ploidies)*len(coverages)*5))
    counter = 0
    for contig in contigs:
        for ploidy in ploidies:
            for coverage in coverages:
                results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'pairwise_results')
                if os.path.exists(results_path):
                    samples = [f for f in os.listdir(results_path) if 'samples_accuracy' in f]
                    for sample in samples:
                        sample_result = os.path.join(results_path, sample)
                        with open(sample_result, 'rb') as f:
                            this_results = pickle.load(f)
                        sample_name = sample.split('.pkl')[0].split('_')[-1]
                        for key in this_results.keys():
                            results_df.loc[counter, 'Contig'] = contig
                            results_df.loc[counter, 'Ploidy'] = ploidy
                            results_df.loc[counter, 'Coverage'] = coverage
                            results_df.loc[counter, 'Sample'] = sample_name
                            results_df.loc[counter, 'Sampling'] = key
                            results_df.loc[counter, 'Value'] = this_results[key]
                            counter += 1
    results_df.to_csv(save_path, index=False)


def collect_results_blocks_pHapCompass():
    agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    sim_data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec', 'length_phased', 'n_blocks', 'average_block_size', 'ffbs_acc']
    contigs = ['100'] 
    ploidies = ['3', '4', '6', '8']
    coverages = ['10', '30', '50', '70', '100']
    results_dfs = []
    for contig in contigs:
        for ploidy in ploidies:
            for coverage in coverages:
                # results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results')
                results_path = os.path.join(sim_data_path, 'contig_' + contig, 'ploidy_' + ploidy, 'cov_' + coverage, 'results_likelihood')
                if os.path.exists(results_path):
                    samples = [f for f in os.listdir(results_path) if 'FFBS' in f]
                    this_result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(samples)*len(metrics)))
                    this_result_df['Contig'] = contig
                    this_result_df['Ploidy'] = ploidy
                    this_result_df['Coverage'] = coverage
                    this_result_df['Method'] = 'pHapCompass'                
                    counter = 0
                    for sample in samples:
                        sample_result = os.path.join(results_path, sample)
                        with open(sample_result, 'rb') as f:
                            this_results = pickle.load(f)

                        sample_name = sample.split('.pkl')[0].split('_')[-1]
                        vector_error_rate = this_results['block_evaluation']['vector_error_rate']
                        vector_error = this_results['block_evaluation']['vector_error']
                        accuracy = this_results['block_evaluation']['accuracy']
                        mismatch_error = this_results['block_evaluation']['mismatch_error']
                        mec_ = this_results['block_evaluation']['mec']
                        ffbs_acc = this_results['evaluation']['ffbs_acc']
                        average_block_size = this_results['average_block_size']
                        # components = this_results['components']
                        # average_block_size = np.mean([components[key]['block_size'] for key in components.keys()])
                        n_blocks = this_results['n_blocks']
                        length_phased = this_results['length_phased']
                        evals = {'vector_error_rate': vector_error_rate, 'vector_error': vector_error, 'accuracy': accuracy, 
                                 'mismatch_error': mismatch_error, 'mec': mec_, 'length_phased': length_phased, 
                                 'n_blocks': n_blocks, 'average_block_size': average_block_size, 'ffbs_acc': ffbs_acc}

                        for metric in metrics:
                            this_result_df.loc[counter, 'Sample'] = sample_name
                            this_result_df.loc[counter, 'Metric'] = metric
                            this_result_df.loc[counter, 'Value'] = evals[metric]
                            counter += 1
                    results_dfs.append(this_result_df)

    results_df = pd.concat(results_dfs, ignore_index=True)
    results_df.to_csv(os.path.join(agg_results_path, 'pHapCompass_results_simulated_data_NA12878_block.csv'), index=False)


def collect_results_blocks_all():
    agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    sim_data_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
 

    # server3_results = pd.read_csv(os.path.join(agg_results_path, 'sim_awri_results_likelihood10_100_server3.csv'))
    phapcompass_results = pd.read_csv(os.path.join(agg_results_path, 'pHapCompass_results_simulated_data_NA12878_block.csv'))
    # whatshapp_results = pd.read_csv(os.path.join(agg_results_path, 'whatshap_results_simulated_data_test.csv'))
    # hpopg_results = pd.read_csv(os.path.join(agg_results_path, 'hpop_results_simulated_data_test.csv'))
    whatshapp_results346 = pd.read_csv(os.path.join(agg_results_path, 'whatshap_results_simulated_data_NA12878_346_block.csv'))
    hpopg_results346 = pd.read_csv(os.path.join(agg_results_path, 'hpop_results_simulated_data_NA12878_346_block.csv'))
    whatshapp_results8 = pd.read_csv(os.path.join(agg_results_path, 'whatshap_results_simulated_data_NA12878_ploidy8_block.csv'))
    hpopg_results8 = pd.read_csv(os.path.join(agg_results_path, 'hpop_results_simulated_data_NA12878_8_block.csv'))


    # all_results = pd.concat([phapcompass_results, whatshapp_results, hpopg_results], ignore_index=True)
    all_results = pd.concat([phapcompass_results, whatshapp_results346, whatshapp_results8, hpopg_results346, hpopg_results8], ignore_index=True)
    all_results['Contig'] = all_results['Contig'].astype(int)
    all_results['Ploidy'] = all_results['Ploidy'].astype(int)
    all_results['Coverage'] = all_results['Coverage'].astype(int)
    # all_results['length_phased'] = all_results['length_phased'].astype(int)

    all_results = all_results.sort_values(by=['Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric']).reset_index(drop=True)
    methods_order = ['pHapCompass', 'WhatsHap' , 'HPoP-G']
    all_results = all_results.sort_values(by=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric'],
        key=lambda col: col.map({method: i for i, method in enumerate(methods_order)})).reset_index(drop=True)

    all_results.to_csv(os.path.join(agg_results_path, 'all_methods_results_simulated_data_NA12878_block.csv'), index=False)


def prepare_results_partial_haplotypes_pHapcompass_MV():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_controlled'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapCompass_simulated_controlled.csv'
    contig_lens = [100]
    ploidies = [2, 3, 4, 6] #, 8]
    coverages = [5, 10, 30, 50, 70, 100]
    n_samples = 100
    # metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*2))
    counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                # results_path = os.path.join(this_cov_path, 'results_FFBS_Multiple')
                results_path = os.path.join(this_cov_path, 'results_FFBS_single')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, 'FFBS_' + str(rd).zfill(2) + '.pkl')
                    with open(result_file, "rb") as f:
                        evals = pickle.load(f)
                
                    result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                    result_df.loc[counter, 'Metric'] = 'Vector Error Rate'
                    result_df.loc[counter, 'Value'] = evals['evaluation']['vector_error_rate']
                    result_df.loc[counter, 'Contig'] = contig_len
                    result_df.loc[counter, 'Ploidy'] = ploidy
                    result_df.loc[counter, 'Coverage'] = coverage
                    counter += 1
                    result_df.loc[counter, 'Contig'] = contig_len
                    result_df.loc[counter, 'Ploidy'] = ploidy
                    result_df.loc[counter, 'Coverage'] = coverage
                    result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                    result_df.loc[counter, 'Metric'] = '# Phased Variants'
                    result_df.loc[counter, 'Value'] = evals['length_phased']
                    counter += 1

    result_df['Method'] = 'pHapCompass'
    result_df.to_csv(output_path, index=False)
    print("pHapcompass results collected")


def FFBS_entropy_confidence():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/ffbs_simulated_NA12878.csv'
    contig_lens = [100]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    metrics = ['accuracy', 'entropy_correct_mean', 'entropy_wrong_mean', 'confidence_correct_mean', 'confidence_wrong_mean']
    result_df = pd.DataFrame(columns=['Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*len(metrics)))
    counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                # results_path = os.path.join(this_cov_path, 'results_FFBS_Multiple')
                results_path = os.path.join(this_cov_path, 'results_FFBS_single_evidence')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, 'FFBS_' + str(rd).zfill(2) + '.pkl')
                    with open(result_file, "rb") as f:
                        evals = pickle.load(f)
                    for metric in metrics:
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = metric
                        # result_df.loc[counter, 'Value'] = evals['evaluation']['ffbs_stats_dict'][metric]
                        result_df.loc[counter, 'Value'] = evals['ffbs_stats_dict'][metric]
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        counter += 1
    result_df.to_csv(output_path, index=False)
    print("pHapcompass results collected")

def collect_single_mv_results():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapcompass_simulated_NA12878_count_evidence_mv.csv'
    contig_lens = [100]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    methods = ['Edge count','Evidence count', 'M.V.']
    result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*2*3))
    counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                # results_path = os.path.join(this_cov_path, 'results_FFBS_Multiple')
                results_path = os.path.join(this_cov_path, 'results_FFBS_single_evidence')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, 'FFBS_' + str(rd).zfill(2) + '.pkl')
                    with open(result_file, "rb") as f:
                        evals = pickle.load(f)
                    for method in methods:
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = 'Vector Error Rate'
                        result_df.loc[counter, 'Value'] = evals[method]['vector_error_rate']
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Method'] = method
                        counter += 1
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = '# Phased Variants'
                        result_df.loc[counter, 'Value'] = evals[method]['length_phased']
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Method'] = method

                        counter += 1
    result_df.to_csv(output_path, index=False)
    print("pHapcompass results collected")


def collect_single_mv_sample_results():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapcompass_simulated_NA12878_solutions_sampled.csv'
    contig_lens = [100]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    methods = ['Single Variant','Multi Variant']
    # metrics = ['vector_error_rate', 'length_phased']
    result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*4*2))
    counter = 0

    for contig_len in contig_lens:
        for ploidy in ploidies:
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                # results_path = os.path.join(this_cov_path, 'results_FFBS_Multiple')
                results_path = os.path.join(this_cov_path, 'results_FFBS_single_mv_samples')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, 'FFBS_' + str(rd).zfill(2) + '.pkl')
                    if result_file == '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878/contig_100/ploidy_8/cov_70/results_FFBS_single_mv_samples/FFBS_01.pkl':
                        continue
                    with open(result_file, "rb") as f:
                        evals = pickle.load(f)
                    for method in methods:
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = 'Vector Error Rate'
                        result_df.loc[counter, 'Value'] = evals[method]['vector_error_rate']
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Method'] = method
                        counter += 1
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = '# Phased Variants'
                        result_df.loc[counter, 'Value'] = evals[method]['length_phased']
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Method'] = method
                        counter += 1
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = 'Confidence'
                        result_df.loc[counter, 'Value'] = evals[method]['solution_confidence']
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Method'] = method
                        counter += 1
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = 'Entropy'
                        result_df.loc[counter, 'Value'] = evals[method]['solution_entropy']
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Method'] = method
                        counter += 1
    result_df.to_csv(output_path, index=False)
    print("pHapcompass results collected")


def collect_single_viterbi_results():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapcompass_simulated_NA12878_single_viterbi_true22.csv'
    contig_lens = [100]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    methods = ['Single Variant', 'Multi Variant']
    # metrics = ['vector_error_rate', 'length_phased']
    result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*2*2))
    counter = 0

    for contig_len in contig_lens:
        for ploidy in ploidies:
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                # results_path = os.path.join(this_cov_path, 'results_FFBS_Multiple')
                results_path = os.path.join(this_cov_path, 'results_viterbi')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, 'FFBS_' + str(rd).zfill(2) + '.pkl')
                    
                    with open(result_file, "rb") as f:
                        evals = pickle.load(f)
                    for method in methods:
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = 'Vector Error Rate'
                        result_df.loc[counter, 'Value'] = evals[method]['vector_error_rate']
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Method'] = 'Viterbi + ' + method.split(' ')[0]
                        counter += 1
                        result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                        result_df.loc[counter, 'Metric'] = '# Phased Variants'
                        result_df.loc[counter, 'Value'] = evals[method]['length_phased']
                        result_df.loc[counter, 'Contig'] = contig_len
                        result_df.loc[counter, 'Ploidy'] = ploidy
                        result_df.loc[counter, 'Coverage'] = coverage
                        result_df.loc[counter, 'Method'] = 'Viterbi + ' + method.split(' ')[0]
                        counter += 1

    result_df.to_csv(output_path, index=False)
    print("pHapcompass results collected")


def collect_results_all():
    agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    phapcompass_results346 = pd.read_csv(os.path.join(agg_results_path, 'pHapCompass_simulated_controlled.csv'))
    # phapcompass_results8 = pd.read_csv(os.path.join(agg_results_path, 'pHapCompass_results_simulated_data_NA12878_8.csv'))
    whatshapp_results346 = pd.read_csv(os.path.join(agg_results_path, 'whatshap_simulated_controlled.csv'))
    hpopg_results346 = pd.read_csv(os.path.join(agg_results_path, 'hpop_simulated_controlled.csv'))
    # whatshapp_results8 = pd.read_csv(os.path.join(agg_results_path, 'whatshap_results_simulated_data_NA12878_ploidy8.csv'))
    # hpopg_results8 = pd.read_csv(os.path.join(agg_results_path, 'hpop_results_simulated_data_NA12878_8.csv'))
    # all_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/all_methods_results_simulated_data_NA12878_LBP_included.csv'
    # all_results = pd.read_csv(all_results_path)
    # mv_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapCompass_MV_results_simulated_data_NA12878.csv'
    # mv_results = pd.read_csv(mv_results_path)


    # all_results = pd.concat([all_results, mv_results], ignore_index=True)
    all_results = pd.concat([phapcompass_results346, whatshapp_results346, hpopg_results346], ignore_index=True)

    # all_results = pd.concat([phapcompass_results346, phapcompass_results8, whatshapp_results346, whatshapp_results8, hpopg_results346, hpopg_results8], ignore_index=True)
    all_results['Contig'] = all_results['Contig'].astype(int)
    all_results['Ploidy'] = all_results['Ploidy'].astype(int)
    all_results['Coverage'] = all_results['Coverage'].astype(int)
    # all_results['length_phased'] = all_results['length_phased'].astype(int)

    all_results = all_results.sort_values(by=['Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric']).reset_index(drop=True)
    
    # methods_order = ['pHapCompass', 'pHapCompass + LBP' , 'pHapCompass + M.V.', 'H-PoPG', 'WhatsHap']
    methods_order = ['pHapCompass', 'H-PoPG', 'WhatsHap']
    all_results = all_results.sort_values(by=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric'],
        key=lambda col: col.map({method: i for i, method in enumerate(methods_order)})).reset_index(drop=True)

    all_results.to_csv(os.path.join(agg_results_path, 'all_methods_controlled.csv'), index=False)


def compare_metric_methods(all_results, agg_results_path):
    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/plots/'

    all_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/all_methods_results_simulated_data_NA12878_missing.csv'
    method_palette = {
    "pHapCompass": "tab:blue",  # Blue
    "WhatsHap": "tab:orange",  # Orange
    "H-PoPG": "tab:green",  # Green  
}
    
    # compare metrics for different methods:
    for contig in [100]:
        # for metric in ['vector_error_rate', 'mismatch_error', 'mec', 'accuracy', 'vector_error' ,'length_phased', 'n_blocks', 'average_block_size',]:
        for metric in ['# Phased Variants', 'Vector Error Rate']:
            for ploidy in [3, 4, 6, 8]:
                print('Contig:', contig, 'Metric:', metric, 'Ploidy:', ploidy)
                metric_df = all_results[(all_results['Metric'] == metric) & (all_results['Ploidy'] == ploidy)].reset_index(drop=True)
                
                g = sns.catplot(x="Coverage", y="Value", hue="Method", data=metric_df, kind="box", height=6, aspect=1.5, palette=method_palette)

                # Add the title and labels
                g.fig.suptitle(f"Contig: {str(contig)}, Metric: {metric.capitalize()}, Ploidy: {str(ploidy).capitalize()}", y=1.05)
                g.set_axis_labels("Coverage", "Value")

                # Move the legend to the top-right, inside the plot area
                g._legend.set_bbox_to_anchor((0.95, 0.9))  # Adjust the position
                g._legend.set_frame_on(True)  # Optional: Add a frame around the legend
                g._legend.set_title("Method")  # Optional: Customize legend title

                # Adjust the layout to ensure everything fits within the figure
                g.fig.subplots_adjust(top=0.85, right=0.9)

                # # Save the figure
                # plt.show()
                plt.savefig(os.path.join(plot_path, f"compare_{contig}_{metric}_{ploidy}.png"), bbox_inches="tight", dpi=300)
                plt.close()


def compare_metric_methods_line(all_results, agg_results_path):
    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/plots/'

    all_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/all_methods_results_simulated_data_NA12878_missing.csv'
    all_results = pd.read_csv(all_results_path)
    method_palette = {
    "pHapCompass": "tab:blue",  # Blue
    "WhatsHap": "tab:orange",  # Orange
    "H-PoPG": "tab:green"}
    
    # compare metrics for different methods:
    for contig in [100]:
        for metric in ['# Phased Variants', 'Vector Error Rate']:
            for ploidy in [3, 4, 6, 8]:
                print('Contig:', contig, 'Metric:', metric, 'Ploidy:', ploidy)
                metric_df = all_results[(all_results['Metric'] == metric) & (all_results['Ploidy'] == ploidy)].reset_index(drop=True)
                
                # Create figure and axis
                plt.figure(figsize=(6, 6))
                ax = sns.lineplot(x="Coverage", y="Value", hue="Method", data=metric_df, marker="o", linewidth=4, palette=method_palette)

                # Add the title and labels
                # plt.title(f"Contig: {str(contig)}, Metric: {metric.capitalize()}, Ploidy: {str(ploidy).capitalize()}", y=1.05, fontdict={"size": 20})
                plt.xlabel("Coverage", fontdict={"size": 20})
                # plt.ylabel("Value", fontdict={"size": 20})
                plt.ylabel(metric, fontdict={"size": 20})

                # Ensure x-ticks match only the unique Coverage values
                unique_coverage = sorted(metric_df["Coverage"].unique())
                plt.xticks(unique_coverage)
                plt.xticks(fontsize=16)
                plt.yticks(fontsize=16)
                # Move the legend to the top-right, inside the plot area
                plt.legend(title="Method", loc="upper right", bbox_to_anchor=(0.95, 0.9), frameon=True, fontsize=16, title_fontsize=20)

                # Adjust layout
                plt.tight_layout()

                # Save the figure
                plt.savefig(os.path.join(plot_path, f"line_compare_{contig}_{metric}_{ploidy}.png"), bbox_inches="tight", dpi=300)
                plt.close()


def plot_phased_length_box(all_results, agg_results_path):
    for contig in [100]:
        for metric in ['vector_error_rate']:
            for ploidy in [3, 4, 6, 8]:
                # stop
                print('Contig:', contig, 'Ploidy:', ploidy)
                metric_df = all_results[(all_results['Metric'] == metric) & (all_results['Ploidy'] == ploidy) & (all_results['Contig'] == contig)].reset_index(drop=True)

                metric_df['length_phased'] = metric_df['length_phased'].astype(int)
                metric_df['length_phased'] = metric_df['length_phased']/contig * 100

                g = sns.catplot(
                    x="Coverage", y="length_phased", hue="Method",
                    data=metric_df, kind="box", height=6, aspect=1.5
                )

                # Add the title and labels
                g.fig.suptitle("Percentage of phased variants (%)," + f"Contig: {str(contig)}" + f"Ploidy: {str(ploidy)}", y=1.05)
                g.set_axis_labels("Coverage", "phased variants (%)")

                # Move the legend to the top-right, inside the plot area
                g._legend.set_bbox_to_anchor((0.95, 0.9))  # Adjust the position
                g._legend.set_frame_on(True)  # Optional: Add a frame around the legend
                g._legend.set_title("Method")  # Optional: Customize legend title

                # Adjust the layout to ensure everything fits within the figure
                # g.fig.subplots_adjust(top=0.85, right=0.9)
                g.fig.subplots_adjust(top=1, right=0.9)

                # Save the figure
                # plt.show()
                plt.savefig(os.path.join(agg_results_path, f"phased_length_{ploidy}_{contig}.png"), bbox_inches="tight", dpi=300)
                plt.close()


def plot_phased_length_line(all_results, agg_results_path):
    for contig in [100]:
        for metric in ['vector_error_rate']:
            for ploidy in [3, 4, 6, 8]:
                # stop
                print('Contig:', contig, 'Ploidy:', ploidy)
                metric_df = all_results[(all_results['Metric'] == metric) & (all_results['Ploidy'] == ploidy) & (all_results['Contig'] == contig)].reset_index(drop=True)

                metric_df['length_phased'] = metric_df['length_phased'].astype(int)
                metric_df['length_phased'] = metric_df['length_phased']/contig * 100

                sns.lineplot(
                    x="Coverage", y="length_phased", hue="Method",
                    data=metric_df, marker="o", linewidth=2
                )

                # Title and labels
                plt.title(f"Percentage of Phased Variants (%), Contig: {contig}, Ploidy: {ploidy}", fontsize=14)
                plt.xlabel("Coverage", fontsize=12)
                plt.ylabel("Phased Variants (%)", fontsize=12)

                # Move the legend to the top-right inside the plot area
                plt.legend(title="Method", loc="upper right", bbox_to_anchor=(0.95, 0.8), frameon=True)
                # plt.show()
                # Save and close the figure
                plt.savefig(os.path.join(agg_results_path, f"phased_length_{ploidy}_{contig}_line.png"), bbox_inches="tight", dpi=300)
                plt.close()


def plot_ffbs_accuracy():
    agg_results_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    ffbs_df_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/ffbs_acc_NA12878.csv'

    ffbs_df = pd.read_csv(ffbs_df_path)
    ffbs_df['Ploidy'] = ffbs_df['Ploidy'].astype(int)
    ffbs_df['Coverage'] = ffbs_df['Coverage'].astype(int)
    g = sns.catplot(x="Ploidy", y="FFBS_Acc", hue="Coverage",  data=ffbs_df, kind="box", height=6, aspect=1.5)

    # Add the title and labels
    g.fig.suptitle("FFBS accuracy" , y=1.05)
    g.set_axis_labels("Ploidy", "")

    # Move the legend to the top-right, inside the plot area
    g._legend.set_bbox_to_anchor((0.95, 0.25))  # Adjust the position
    g._legend.set_frame_on(True)  # Optional: Add a frame around the legend
    g._legend.set_title("Coverage")  # Optional: Customize legend title

    # Adjust the layout to ensure everything fits within the figure
    # g.fig.subplots_adjust(top=0.85, right=0.9)
    g.fig.subplots_adjust(top=1, right=1)

    # Save the figure
    plt.show()
    plt.savefig(os.path.join(agg_results_path, f"FFBS_acc_NA12878.png"), bbox_inches="tight", dpi=300)
    plt.close()

    # linegraph pallet:
    ploidy_palette = {
    3: "tab:blue",  # Blue
    4: "tab:orange",  # Orange
    6: "tab:green", 
    8: "tab:red"}

    plt.figure(figsize=(6, 6))
    ax = sns.lineplot(x="Coverage", y="FFBS_Acc", hue="Ploidy", data=ffbs_df, marker="o", linewidth=4, palette=ploidy_palette)
    
    plt.xlabel("Coverage", fontdict={"size": 20})
    # plt.ylabel("Value", fontdict={"size": 20})
    plt.ylabel('FFBS Acuuracy', fontdict={"size": 20})
    unique_coverage = sorted(ffbs_df["Coverage"].unique())
    plt.xticks(unique_coverage)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    # plt.legend(title="Method", loc="upper right", bbox_to_anchor=(0.95, 0.9), frameon=True, fontsize=16, title_fontsize=20)
    plt.legend(title="Ploidy", loc="lower right", bbox_to_anchor=(0.9, 0.1), frameon=True, fontsize=16, title_fontsize=20)

    plt.tight_layout()
    plt.show()
    # Adjust layout


def generate_fragmentfile_test():
    vcf_file_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/AWRI_ploidy6_contig10.vcf'
    modified_vcf_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/modified_AWRI_ploidy6_contig10.vcf'
    bam_file_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/cov_50/bam/00.bam'
    frag_file_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/cov_50/frag/00.frag'

    vcf_file = pysam.VariantFile(vcf_file_path)
    original_snps = [record.pos for record in vcf_file.fetch()] # 1-based positions
    orig_base = [record.ref for record in vcf_file.fetch()]
    alt_base = [record.alts[0] for record in vcf_file.fetch()]

    # Open the BAM file
    read_cnt = 0
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for read in bam_file.fetch():
            # Get the read name
            read_name = read.query_name

            # Get the aligned positions (list of genomic coordinates)
            aligned_positions = read.get_reference_positions(full_length=True)

            # Extract the base at each position of interest
            for pos in original_snps:
                # Check if the position is covered by the read
                if pos in aligned_positions:
                    base_index = aligned_positions.index(pos)
                    base = read.query_sequence[base_index]
                    # print(f"Read: {read_name}, Position: {pos}, Base: {base}")
                    if base == orig_base[original_snps.index(pos)]:
                        # print(0)
                        print(f"Read: {read_name}, Position: {pos}, Base: {base}: 0")
                    elif base == alt_base[original_snps.index(pos)]:
                        # print(1)
                        print(f"Read: {read_name}, Position: {pos}, Base: {base}: 1")
                    else:
                        print(f"Read: {read_name}, Position: {pos}, Base: {base}: -")
            print('------------------------------------------------')
            read_cnt += 1    # else:
                    # print(f"Read: {read_name}, Position: {pos} not covered.")

    fasta_file_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/contig_10_ploidy_6.fa'
    with open(fasta_file_path, 'r') as f:
        haplotypes = f.readlines()
        f.close()

    haplotypes = [h.strip() for h in haplotypes if h.startswith('>') == False]
    for ii, haplotype in enumerate(haplotypes):
        print([haplotype[pos-1] for pos in original_snps])


    readlist = []
    fragment_list = []

    for fragment in open(frag_file_path, 'r'):
        # print(fragment)
        # read index when there are pairs, readname, start pos (in variants), allele sequence, alleles
        parts = fragment.split()
        readlist.append(parts[1])
        positions = []
        alleles = []
        for iii in range(int(parts[0])):
            # process i+2,i+3.... i+4,i+5...
            start_idx_of_read = iii * 2 + 3
            seq_len = len(parts[start_idx_of_read])
            positions.extend(list(range(int(parts[start_idx_of_read - 1]), int(parts[start_idx_of_read - 1]) + seq_len)))
            [alleles.append(int(a)) for a in parts[start_idx_of_read]]
            fragment_list.append(positions)
            fragment_list.append(alleles)


    # List of missing chromosomes to add
    missing_chromosomes = ["haplotype_1", "haplotype_2", "haplotype_3", "haplotype_4", "haplotype_5", "haplotype_6"]

    # Read the original VCF file
    with open(vcf_file_path, "r") as infile:
        lines = infile.readlines()

    # Open the output VCF file for writing
    with open(modified_vcf_path, "w") as outfile:
        header_written = False

        # Track if contigs are already in the VCF
        existing_contigs = {line.split("=")[1].split(",")[0] for line in lines if line.startswith("##contig")}

        for line in lines:
            # Write the original lines
            if not header_written and not line.startswith("#CHROM"):
                outfile.write(line)
                # If the line is a contig line, mark it
                if line.startswith("##contig"):
                    last_contig_line = line
            elif not header_written and line.startswith("#CHROM"):
                # Add missing chromosomes just before the data starts
                for chrom in missing_chromosomes:
                    if chrom not in existing_contigs:
                        outfile.write(f"##contig=<ID={chrom}>\n")
                header_written = True
                outfile.write(line)  # Write the column headers
            else:
                # Write the rest of the data lines
                outfile.write(line)

    print(f"Updated VCF file saved to: {modified_vcf_path}")




    input_vcf = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/modified_AWRI_ploidy6_contig10.vcf"
    output_vcf = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/1_modified_AWRI_ploidy6_contig10.vcf"

    # List of haplotypes to add
    haplotypes = ["haplotype_1", "haplotype_2", "haplotype_3", "haplotype_4", "haplotype_5", "haplotype_6"]

    # Open the input and output files
    with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            # Write header lines as-is
            if line.startswith("#"):
                outfile.write(line)
            else:
                # Split the original line into columns
                columns = line.strip().split("\t")
                original_chrom, pos = columns[0], columns[1]
                # Write the original line for AHIQ01000001.1
                outfile.write(line)
                # Add corresponding lines for haplotypes
                for hap in haplotypes:
                    columns[0] = hap  # Replace chromosome name
                    outfile.write("\t".join(columns) + "\n")

    print(f"Expanded VCF file saved to: {output_vcf}")





    # File paths
    bam_file_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/test_out.bam"
    vcf_file_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/1_modified_AWRI_ploidy6_contig10.vcf"

    # Extract contig names from BAM
    bam_header = pysam.AlignmentFile(bam_file_path, "rb").header
    bam_contigs = [ref['SN'] for ref in bam_header['SQ']]

    # Extract contig names from VCF
    vcf_contigs = []
    with open(vcf_file_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith("##contig"):
                contig_name = line.split("ID=")[1].split(",")[0].strip(">")
                vcf_contigs.append(contig_name)
            elif not line.startswith("#"):
                break  # Stop after header

    # Check for contig mismatches
    mismatched_contigs = set(bam_contigs) - set(vcf_contigs)

    print("BAM Contigs:", bam_contigs)
    print("VCF Contigs:", vcf_contigs)
    print("Mismatched Contigs:", mismatched_contigs)

    # Updated BAM file path
    bam_file_path = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/test_out.bam"

    # Gather summary statistics

    # Open the BAM file
    bam_file = pysam.AlignmentFile(bam_file_path, "rb")

    # Gather summary statistics
    bam_stats = {
        "total_reads": bam_file.mapped + bam_file.unmapped,
        "mapped_reads": bam_file.mapped,
        "unmapped_reads": bam_file.unmapped,
        "references": bam_file.references,  # Contig names
    }

    print("BAM File Statistics:", bam_stats)


    input_vcf = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/1_modified_AWRI_ploidy6_contig10.vcf"
    output_vcf = "/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/2_modified_AWRI_ploidy6_contig10.vcf" # Output VCF file

    with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                # Write header lines as is
                outfile.write(line)
            else:
                # Split the line into columns
                columns = line.strip().split("\t")
                format_field = columns[8]  # FORMAT column
                sample_data = columns[9]  # Sample data
                # stop
                # Modify the GT field if it exists in the FORMAT column
                if "GT" in format_field:
                    gt_index = format_field.split(":").index("GT")
                    sample_data_fields = sample_data.split(":")
                    sample_data_fields[gt_index] = sample_data_fields[gt_index].replace("|", "/")
                    columns[9] = ":".join(sample_data_fields)

                # Write the modified line to the output
                outfile.write("\t".join(columns) + "\n")

    print(f"Fixed VCF file saved as {output_vcf}")


def extract_fragments_0(vcf_file, bam_file, output_file):
    # Parse VCF file and collect variant information
    vcf_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/1_modified_AWRI_ploidy6_contig10.vcf'
    bam_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/test_out.bam'
    output_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/test_out7.frag'
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        bam_contigs = set(bam.references)

    # Parse VCF file and collect variant information
    variants = []
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip header lines
            columns = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, _, format_field, sample_data = columns[:10]
            
            # Skip variants on the excluded chromosome
            if chrom == "AHIQ01000001.1":
                continue

            # Skip variants on contigs not in the BAM file
            if chrom not in bam_contigs:
                continue

            gt_index = format_field.split(":").index("GT")
            genotype = sample_data.split(":")[gt_index]
            if "|" in genotype or "/" in genotype:  # Keep heterozygous positions
                variants.append({
                    "chrom": chrom,
                    "pos": int(pos),  # Keep 1-based position as is
                    "ref": ref,
                    "alt": alt.split(",")
                })

    # Open BAM file and process reads for each variant
    fragments = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for variant in variants:
            chrom, pos, ref, alt_alleles = variant["chrom"], variant["pos"], variant["ref"], variant["alt"]
            
            # Fetch reads overlapping the variant position
            for read in bam.fetch(chrom, pos - 1, pos):  # BAM fetch uses 0-based start, 1-based end
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue  # Skip unmapped or non-primary alignments

                aligned_positions = read.get_reference_positions(full_length=True)
                if pos in aligned_positions:  # Directly use 1-based position for checking
                    # Get the base at the variant position
                    base_index = aligned_positions.index(pos)
                    base = read.query_sequence[base_index]

                    # Determine haplotype informative base (0 for REF, 1 for ALT)
                    if base == ref:
                        haplotype_base = "0"
                    elif base in alt_alleles:
                        haplotype_base = "1"
                    else:
                        haplotype_base = "-"  # Non-informative or mismatched base

                    # Create fragment entry
                    fragment_entry = f"{read.query_name}\t{chrom}-{pos}\t{haplotype_base}\t{base}\n"
                    fragments.append(fragment_entry)

    # Write fragments to output file
    with open(output_file, "w") as output:
        output.writelines(fragments)

    print(f"Fragments saved to {output_file}")


def extract_fragments_old(vcf_file, bam_file, output_file):
    # Parse VCF file and collect variant information into a dictionary
    variants = {}
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip header lines
            columns = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, _, _, _ = columns[:10]
            variants.setdefault(chrom, {})[int(pos)] = {
                "ref": ref,
                "alt": alt.split(","),
            }
    variants = {chrom: variants[chrom] for chrom in variants if chrom != "AHIQ01000001.1"}
    # Open BAM file and process reads
    fragments = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue  # Skip unmapped or non-primary alignments

            read_name = read.query_name
            chrom = bam.get_reference_name(read.reference_id)  # Chromosome for this read
            if chrom not in variants:
                continue  # Skip reads from chromosomes without variants

            aligned_positions = read.get_reference_positions(full_length=True)
            read_sequence = read.query_sequence

            # Process only positions in the current chromosome's variants
            fragment_positions = []
            fragment_bases = []
            for pos in variants[chrom]:
                if pos - 1 in aligned_positions:  # BAM positions are 0-based
                    base_index = aligned_positions.index(pos - 1)
                    base = read_sequence[base_index]
                    if base == variants[chrom][pos]["ref"]:
                        fragment_bases.append("0")
                    elif base in variants[chrom][pos]["alt"]:
                        fragment_bases.append("1")
                    else:
                        fragment_bases.append("-")
                    fragment_positions.append(pos)  # Keep as 1-based for output

            # Write a fragment if the read spans multiple variants
            if len(fragment_positions) > 1:
                fragments.append(
                    f"{read_name}\t{chrom}\t{len(fragment_positions)}\t{''.join(fragment_bases)}\n"
                )


    # Write fragments to output file
    with open(output_file, "w") as output:
        output.writelines(fragments)

    print(f"Fragments saved to {output_file}")



    import pysam


def extract_custom_fragments(vcf_file, bam_file, output_file, distance_threshold=10):
    # Parse VCF file and collect variant information into a dictionary
    # Parse VCF file and collect variant information into a dictionary
    variants = {}
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip header lines
            columns = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, _, _, _ = columns[:10]
            variants.setdefault(chrom, {})[int(pos)] = {
                "ref": ref,
                "alt": alt.split(","),
            }

    # Open BAM file and process reads
    fragments = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue  # Skip unmapped or non-primary alignments

            read_name = read.query_name
            chrom = bam.get_reference_name(read.reference_id)  # Chromosome for this read
            if chrom not in variants:
                continue  # Skip reads from chromosomes without variants

            aligned_positions = read.get_reference_positions(full_length=True)
            read_sequence = read.query_sequence

            # Process only positions in the current chromosome's variants
            variant_positions = []
            haplotype_bases = []
            for pos in variants[chrom]:
                if pos - 1 in aligned_positions:  # BAM positions are 0-based
                    base_index = aligned_positions.index(pos - 1)
                    base = read_sequence[base_index]
                    if base == variants[chrom][pos]["ref"]:
                        haplotype_bases.append("0")
                    elif base in variants[chrom][pos]["alt"]:
                        haplotype_bases.append("1")
                    else:
                        haplotype_bases.append("-")
                    variant_positions.append(pos)  # Keep as 1-based for output

            # Split into parts based on gap_threshold = 1
            parts = []
            current_part = []
            current_bases = []
            for i in range(len(variant_positions)):
                if current_part and variant_positions[i] != current_part[-1] + 1:
                    parts.append((current_part, current_bases))
                    current_part = []
                    current_bases = []
                current_part.append(variant_positions[i])
                current_bases.append(haplotype_bases[i])
            if current_part:
                parts.append((current_part, current_bases))

            # Write the fragment parts
            if parts:
                parts_count = len(parts)
                fragment_line = f"{parts_count}\t{chrom}-{variant_positions[0]}_MP\t"
                for part_positions, part_bases in parts:
                    fragment_line += f"{part_positions[0]}\t{len(part_bases)}\t{''.join(part_bases)}\t"
                fragment_line += "GGCGGGG\n"  # Arbitrary sequence
                fragments.append(fragment_line)

    # Write fragments to output file
    with open(output_file, "w") as output:
        output.writelines(fragments)

    print(f"Fragments saved to {output_file}")


def change_haplotype_orientation():
    hap_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_test/contig_100/ploidy_6/haplotypes.csv'
    hap_df = pd.read_csv(hap_path).T
    columns = ['haplotype_' + str(i + 1) for i in range(len(hap_df.columns))]
    hap_df.columns = columns
    hap_df.to_csv(hap_path, index=False)


def extract_fragments(vcf_file, bam_file, output_file):

    # Parse VCF file and collect variant information into a dictionary
    vcf_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/AWRI_ploidy6_contig10.vcf'
    bam_file = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_awri/contig_10/ploidy_6/cov_50/bam/00.bam'
    variants = {}
    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue  # Skip header lines
            columns = line.strip().split("\t")
            chrom, pos, _, ref, alt, _, _, _, _, _ = columns[:10]
            variants.setdefault(chrom, {})[int(pos)] = {
                "ref": ref,
                "alt": alt.split(","),
            }

    # Open BAM file and process reads
    fragments = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue  # Skip unmapped or non-primary alignments
            
            read_name = read.query_name
            chrom = bam.get_reference_name(read.reference_id)  # Chromosome for this read
            
            if chrom not in variants:
                continue  # Skip reads from chromosomes without variants

            aligned_positions = read.get_reference_positions(full_length=True)
            read_sequence = read.query_sequence

            # Process only positions in the current chromosome's variants
            fragment_positions = []
            fragment_bases = []
            for pos in variants[chrom]:
                if pos - 1 in aligned_positions:  # BAM positions are 0-based
                    base_index = aligned_positions.index(pos - 1)
                    base = read_sequence[base_index]
                    if base == variants[chrom][pos]["ref"]:
                        fragment_bases.append("0")
                    elif base in variants[chrom][pos]["alt"]:
                        fragment_bases.append("1")
                    else:
                        fragment_bases.append("-")
                    fragment_positions.append(pos)  # Keep as 1-based for output
            
            # Write a fragment if the read spans multiple variants
            if len(fragment_positions) > 1:
                # stop
                chrom_variants = sorted(variants[chrom].keys())
                fragment_indices = [chrom_variants.index(pos) + 1 for pos in fragment_positions]

                # Segment into consecutive parts
                parts = []
                current_part_indices = []
                current_part_bases = []

                for idx, base in zip(fragment_indices, fragment_bases):
                    if current_part_indices and idx != current_part_indices[-1] + 1:
                        # Append the current part when a gap is found
                        parts.append((current_part_indices, current_part_bases))
                        current_part_indices = []
                        current_part_bases = []
                    # Add the current index and base to the current part
                    current_part_indices.append(idx)
                    current_part_bases.append(base)
                            # Append the last part
                if current_part_indices:
                    parts.append((current_part_indices, current_part_bases))

            # Construct the fragment line
                fragment_line = f"{len(parts)}\t{chrom}-{read_name}\t"
                for part_indices, part_bases in parts:
                    print(part_indices, part_bases)
                    fragment_line += f"{part_indices[0]}\t{''.join(part_bases)}\t"
                fragment_line += "GG\n"


                # fragments.append(
                #     f"{read_name}\t{chrom}-{'-'.join(map(str, fragment_positions))}\t{len(fragment_positions)}\t{''.join(fragment_bases)}\n"
                # )
                fragments.append(fragment_line)

    # Write fragments to output file
    with open(output_file, "w") as output:
        output.writelines(fragments)

    print(f"Fragments saved to {output_file}")


def make_ffbs_vs_vector_error_rate_plot():
    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    ffbs_df_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/ffbs_acc_NA12878.csv'
    results_df_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapCompass_results_simulated_data_NA12878.csv'

    # Load data
    ffbs_df = pd.read_csv(ffbs_df_path)
    results_df = pd.read_csv(results_df_path)

    # Convert relevant columns to integer type
    ffbs_df['Ploidy'] = ffbs_df['Ploidy'].astype(int)
    ffbs_df['Coverage'] = ffbs_df['Coverage'].astype(int)
    results_df['Ploidy'] = results_df['Ploidy'].astype(int)
    results_df['Coverage'] = results_df['Coverage'].astype(int)

    # Filter results_df to keep only 'vector_error_rate' rows
    results_df = results_df[results_df['Metric'] == 'vector_error_rate'].reset_index(drop=True)

    # Merge dataframes
    merged_df = results_df.merge(ffbs_df, on=['Contig', 'Ploidy', 'Coverage', 'Sample'], how='inner')
    merged_df = merged_df.drop(columns=['Method', 'Metric', 'length_phased'])
    merged_df = merged_df.rename(columns={'Value': 'Vector Error Rate', 'FFBS_Acc': 'FFBS Accuracy'})

    # Define color palette
    palette = sns.color_palette("tab10")

    # Loop through specific contigs and ploidies
    for contig in [100]:
        for ploidy in [3, 4, 6, 8]:
            print('Contig:', contig, 'Ploidy:', ploidy)
            
            # Subset data for this ploidy
            ploidy_df = merged_df[(merged_df['Ploidy'] == ploidy) & (merged_df['Contig'] == contig)].reset_index(drop=True)

            # Create scatter plot with regression lines
            plt.figure(figsize=(8, 6))

            # Scatter plot (swapped axes)
            scatter = sns.scatterplot(data=ploidy_df, x="Vector Error Rate", y="FFBS Accuracy",
                                    hue="Coverage", palette=palette, s=10, edgecolor="black", alpha=0.5)

            # Regression lines (swapped axes) - NO LEGEND
            coverages = sorted(ploidy_df['Coverage'].unique())  # Ensure same order for colors
            for i, coverage in enumerate(coverages):
                subset = ploidy_df[ploidy_df["Coverage"] == coverage]
                sns.regplot(data=subset, x="Vector Error Rate", y="FFBS Accuracy",
                            scatter=False, color=palette[i], ci=95)  # Removed label to avoid extra legend

            # Title and labels (swapped)
            plt.title(f"Ploidy {ploidy}, Contig {contig}: Vector Error Rate vs. FFBS Accuracy")
            plt.xlabel("Vector Error Rate")  
            plt.ylabel("FFBS Accuracy")  

            # Keep only scatter plot legend
            handles, labels = scatter.get_legend_handles_labels()
            plt.legend(handles, labels, title="Coverage")

            # Save plot (optional)
            plt.savefig(os.path.join(plot_path, f"ffbs_vs_vector_error_rate_ploidy_{ploidy}_contig_{contig}.png"), dpi=300, bbox_inches="tight")

            # Show plot
            plt.show()


def make_vector_error_rate_vs_FFBS_plot():
    plot_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results'
    ffbs_df_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/ffbs_acc_NA12878.csv'
    results_df_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/pHapCompass_results_simulated_data_NA12878.csv'

    # Load data
    ffbs_df = pd.read_csv(ffbs_df_path)
    results_df = pd.read_csv(results_df_path)

    # Convert relevant columns to integer type
    ffbs_df['Ploidy'] = ffbs_df['Ploidy'].astype(int)
    ffbs_df['Coverage'] = ffbs_df['Coverage'].astype(int)
    results_df['Ploidy'] = results_df['Ploidy'].astype(int)
    results_df['Coverage'] = results_df['Coverage'].astype(int)

    # Filter results_df to keep only 'vector_error_rate' rows
    results_df = results_df[results_df['Metric'] == 'vector_error_rate'].reset_index(drop=True)

    # Merge dataframes
    merged_df = results_df.merge(ffbs_df, on=['Contig', 'Ploidy', 'Coverage', 'Sample'], how='inner')
    merged_df = merged_df.drop(columns=['Method', 'Metric', 'length_phased'])
    merged_df = merged_df.rename(columns={'Value': 'Vector Error Rate', 'FFBS_Acc': 'FFBS Accuracy'})

    # Define color palette
    palette = sns.color_palette("tab10")

    # Loop through specific contigs and ploidies
    for contig in [100]:
        for ploidy in [3, 4, 6, 8]:
            print('Contig:', contig, 'Ploidy:', ploidy)
            
            # Subset data for this ploidy
            ploidy_df = merged_df[(merged_df['Ploidy'] == ploidy) & (merged_df['Contig'] == contig)].reset_index(drop=True)

            # Create scatter plot with regression lines
            plt.figure(figsize=(8, 6))

            # Scatter plot (swapped axes)
            scatter = sns.scatterplot(data=ploidy_df, x="FFBS Accuracy", y="Vector Error Rate",
                                    hue="Coverage", palette=palette, s=10, edgecolor="black", alpha=0.5)

            # Regression lines (swapped axes) - NO LEGEND
            coverages = sorted(ploidy_df['Coverage'].unique())  # Ensure same order for colors
            for i, coverage in enumerate(coverages):
                subset = ploidy_df[ploidy_df["Coverage"] == coverage]
                sns.regplot(data=subset, x="FFBS Accuracy", y="Vector Error Rate",
                            scatter=False, color=palette[i], ci=95)  # Removed label to avoid extra legend

            # Title and labels (swapped)
            plt.title(f"Ploidy {ploidy}, Contig {contig}: FFBS Accuracy vs. Vector Error Rate")
            plt.xlabel("FFBS Accuracy")  # Swapped
            plt.ylabel("Vector Error Rate")  # Swapped

            # Keep only scatter plot legend
            handles, labels = scatter.get_legend_handles_labels()
            plt.legend(handles, labels, title="Coverage")

            # Save plot (optional)
            plt.savefig(os.path.join(plot_path, f"vector_error_rate_vs_ffbs_ploidy_{ploidy}_contig_{contig}.png"), dpi=300, bbox_inches="tight")

            # Show plot
            plt.show()


def make_inputs_for_pHapcompass_eval():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/simulated_data_NA12878'
    output_dir = '/mnt/research/aguiarlab/proj/HaplOrbit/pHapcompass_evals'
    # output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/hpop_results_simulated_data_NA12878.csv'
    contig_lens = [100]
    ploidies = [8]
    coverages = [10, 30, 50, 70, 100]
    n_samples = 100
    inputs = []
    # metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    # result_df = pd.DataFrame(columns=['Method', 'Contig', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value', 'length_phased'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*len(metrics)))
    # counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            true_haplotypes_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                results_path = os.path.join(this_cov_path, 'results_likelihood')
                frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, 'FFBS_' + str(rd).zfill(2) + '.pkl')
                    inputs.append([result_file, true_haplotypes_path, frag_path, contig_len, ploidy, coverage, rd])

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for i, inp in enumerate(inputs):
        input_file = os.path.join(output_dir, f"input_{i}.pkl")
        with open(input_file, "wb") as f:
            pickle.dump(inp, f)
    print(f"Saved {len(inputs)} inputs to {output_dir}")


def eval_one_input_pHapcompass_missing(inp):
    result_file, true_haplotypes_path, frag_path, contig_len, ploidy, coverage, rd = inp
    # true_haplotypes = pd.read_csv(true_haplotypes_path).T.to_numpy()

    with open(result_file, 'rb') as f:
        this_results = pickle.load(f)
    
    # sample_name = result_file.split('.pkl')[0].split('_')[-1]
    H = this_results['predicted_haplotypes'].to_numpy()
    H_star = this_results['true_haplotypes'].to_numpy()
    vector_error_rate, _, _ = compute_vector_error_rate_with_missing_positions(H_star, H)
    length_phased = np.sum(~np.isnan(np.array(H, dtype=np.float64)).any(axis=0))
    evals = {'Contig': contig_len, 'Ploidy': ploidy, 'Coverage': coverage , 'Sample': str(rd).zfill(2), 
             'vector_error_rate': vector_error_rate, 'length_phased': length_phased}
    pkl_name = result_file.replace('FFBS', 'dict')

    with open(pkl_name, "wb") as f:
        pickle.dump(evals, f)
    print(f"Done with {pkl_name}")


def prepare_results_recomb_auto_short():
    main_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results_short'
    output_path = '/mnt/research/aguiarlab/proj/HaplOrbit/results/recomb/pHapCompass_auto_short_2346.csv'
    ploidies = [2, 3, 4, 6] #, 8]
    mut_rates = [0.001, 0.005, 0.01]
    coverages = [3, 5, 10, 20, 40, 70]
    n_samples = 20
    # metrics = ['vector_error_rate', 'vector_error', 'accuracy', 'mismatch_error', 'mec']
    result_df = pd.DataFrame(columns=['Method', 'Data', 'Ploidy', 'Coverage', 'Sample', 'Metric', 'Value'], index=range(len(contig_lens)*len(ploidies)*len(coverages)*n_samples*2))
    counter = 0
    for contig_len in contig_lens:
        for ploidy in ploidies:
            # true_haplotypes = pd.read_csv(os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'haplotypes.csv')).T.to_numpy()
            for coverage in coverages:
                this_cov_path = os.path.join(main_path, 'contig_{}'.format(contig_len), 'ploidy_{}'.format(ploidy), 'cov_{}'.format(coverage))
                # results_path = os.path.join(this_cov_path, 'results_FFBS_Multiple')
                results_path = os.path.join(this_cov_path, 'results_FFBS_single')
                # frag_path = os.path.join(this_cov_path, 'frag')
                for rd in range(n_samples):
                    print(f"Collecting results for contig {contig_len}, ploidy {ploidy}, coverage {coverage}, sample {rd}")
                    result_file = os.path.join(results_path, 'FFBS_' + str(rd).zfill(2) + '.pkl')
                    with open(result_file, "rb") as f:
                        evals = pickle.load(f)
                
                    result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                    result_df.loc[counter, 'Metric'] = 'Vector Error Rate'
                    result_df.loc[counter, 'Value'] = evals['evaluation']['vector_error_rate']
                    result_df.loc[counter, 'Contig'] = contig_len
                    result_df.loc[counter, 'Ploidy'] = ploidy
                    result_df.loc[counter, 'Coverage'] = coverage
                    counter += 1
                    result_df.loc[counter, 'Contig'] = contig_len
                    result_df.loc[counter, 'Ploidy'] = ploidy
                    result_df.loc[counter, 'Coverage'] = coverage
                    result_df.loc[counter, 'Sample'] = str(rd).zfill(2)
                    result_df.loc[counter, 'Metric'] = '# Phased Variants'
                    result_df.loc[counter, 'Value'] = evals['length_phased']
                    counter += 1

    result_df['Method'] = 'pHapCompass'
    result_df.to_csv(output_path, index=False)
    print("pHapcompass results collected")


    with open(result_file, 'rb') as f:
        this_results = pickle.load(f)
    
    # sample_name = result_file.split('.pkl')[0].split('_')[-1]
    H = this_results['predicted_haplotypes'].to_numpy()
    H_star = this_results['true_haplotypes'].to_numpy()
    vector_error_rate, _, _ = compute_vector_error_rate_with_missing_positions(H_star, H)
    length_phased = np.sum(~np.isnan(np.array(H, dtype=np.float64)).any(axis=0))








def eval_one_input_from_input_pHapcompass(input_file):
    with open(input_file, "rb") as f:
        inp = pickle.load(f)

    # eval_one_input_hpop_block(inp)
    eval_one_input_pHapcompass_missing(inp)


if __name__ == '__main__':
    # # results from 3 4 6 ploidies
    # prepare_results_partial_haplotypes_pHapcompass346()
    # make_inputs_for_pHapcompass_eval()
    # prepare_results_partial_haplotypes_pHapcompass8()

    if len(sys.argv) != 2:
        print("Usage: python3 results.py <input_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    eval_one_input_from_input_pHapcompass(input_file)