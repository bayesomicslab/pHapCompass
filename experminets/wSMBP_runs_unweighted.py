import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from algorithm.wSMBP import pipeline
import pickle
import pandas as pd
from multiprocessing import Pool

def run_pipeline(args):
    contig, ploidy, coverage, sample, data, divide = args
    try:
        pipeline(ploidy=ploidy, coverage=coverage, contig=contig, sample=sample, data=data, is_weighted=False, divide=divide)
    except Exception as e:
        print(f"Error running pipeline for contig {contig}, ploidy {ploidy}, coverage {coverage}, sample {sample}, divide {divide}")

def run_pipelines(data, contigs, ploidies, coverages, samples, divides, num_processes=15):
    args_list = [(contig, ploidy, coverage, sample, data, divide) for contig in contigs for ploidy in ploidies for coverage in coverages for sample in samples for divide in divides]
    processes_pool = Pool(num_processes)
    processes_pool.map(run_pipeline, args_list)
    processes_pool.close()
    processes_pool.join()
    print("All pipelines completed.")

if __name__ == '__main__':
    data = 'NA12878'
    contigs = [100, 1000]
    ploidies = [3, 4, 6, 8]
    coverages = [10, 30, 50, 70, 100]
    divides = ['none']
    samples = [str(i).zfill(2) for i in range(100)]

    run_pipelines(data, contigs, ploidies, coverages, samples, divides, num_processes=15)