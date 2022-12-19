import argparse
import os

import pandas as pd
import numpy as np


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
    parser.add_argument("--max_memo", dest="max_mb", default=20, type=int,
                        help="Max number of cells (individuals multiple by sites) to use in a single matrix "
                             "(in millions). if doesn't know, don't touch. If there are memory failures, reduce it.")
    parser.add_argument("--args", dest="args", help="Any additional args")
    options = parser.parse_args()
    if options.output[-1] != '/':
        options.output += '/'
    os.makedirs(options.output, exist_ok=True)
    return options


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
    if input.endswith(".vcf"):
        return "VCF"
    raise IOError("File format is not supported")


def read_df_file(f_path):
    if f_path.endswith(".xlsx"):
        data = pd.read_excel(f_path)
    elif f_path.endswith(".csv"):
        data = pd.read_csv(f_path)
    else:
        assert False, "ERROR Parsing GENEpop format file"
    return data

