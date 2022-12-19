import argparse
import os
import pandas as pd
import numpy as np
from multiprocessing import Process

import time


def args_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="input file in genepop (xlsx) format with raw sequencing")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Name of output directory. Program will generate a new directory with that name,"
                             " and all output files will be there")
    parser.add_argument("-w", "--weighted", dest="weighted_metric", default=False, action="store_true",
                        help="If used, compute weighted metric, from Greenbaum et al., 2016. If not use, compute "
                             "unweighted metric from  Li and Horvitz, 1953")
    parser.add_argument("--max_memo", dest="max_mb", default=2, type=float,
                        help="Max number of cells (individuals lultiple by sites) to use in a single matrix "
                             "(in millions). if doesn't know, don't touch. If there are memory failures, reduce it.")
    parser.add_argument("--max_threads", dest="max_threads", default=8, type=int,
                        help="Maximum number of threads to compute small matrices simultaneously")
    parser.add_argument("--max_sites", dest="max_sites", default=0, type=int,
                        help="If assigned, compute similarity only based on the first n sites")
    parser.add_argument("--args", dest="args", help="Any additional args")
    options = parser.parse_args()
    if options.output[-1] != '/':
        options.output += '/'
    os.makedirs(options.output, exist_ok=True)
    return options


def wait_for_jobs_to_be_done(jobs, max_number_of_threads):
    finished_jobs = []
    while len(jobs) - len(finished_jobs) >= max_number_of_threads:
        time.sleep(1)
        for job in jobs:
            job.join(timeout=0)
            if not job.is_alive():
                finished_jobs.append(job)
    running_jobs = [j for j in jobs if j not in finished_jobs]
    return running_jobs


def dfs2_3d_numpy(list_of_dfs):
    num_of_matrices = len(list_of_dfs)
    res = np.empty(shape=(num_of_matrices, list_of_dfs[0].shape[0], list_of_dfs[0].shape[1]))
    for idx, matrix in enumerate(list_of_dfs):
        res[idx] = matrix
    return res


def parse_input_file_format(input):
    if input.endswith(".xlsx") or input.endswith(".csv"):
        return "DF"
    if input.endswith("vcf.gz"):
        return "VCF-GZ"
    if input.endswith(".vcf") or input.endswith(".txt"):
        return "VCF"
    raise IOError("File format is not supported")


def read_df_file(f_path):
    if f_path.endswith(".xlsx"):
        df = pd.read_excel(f_path)
    elif f_path.endswith(".csv"):
        df = pd.read_csv(f_path)
    elif f_path.endswith('.txt'):
        with open(f_path, "r") as f:
            data = f.readlines()
            if data[-1] == '\n':
                data = data[:-1]
            data_split = [e[:-1].split('\t') for e in data]
            line_indices_to_remove = []
            for idx, line in enumerate(data_split):
                non_fail = [e for e in line if e != "FAIL"]
                if len(non_fail) == 0:
                    line_indices_to_remove.append(idx)
                alleles = non_fail[0].split('/')
                if all([non_fail[0] == e for e in non_fail]) and alleles[0] == alleles[1]:
                    line_indices_to_remove.append(idx)
            data_split = [j for i, j in enumerate(data_split) if i not in line_indices_to_remove]
            df = pd.DataFrame(data_split).T
    else:
        assert False, "ERROR Parsing GENEpop format file"
    return df


