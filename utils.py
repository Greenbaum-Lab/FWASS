import argparse
import pandas as pd
import numpy as np


def args_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", required=True,
                        help="input file in genepop (xlsx) format with raw sequencing")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Name of output directory. Program will generate a new directory with that name,"
                             " and all output files will be there")

    parser.add_argument("--args", dest="args", help="Any additional args")
    options = parser.parse_args()
    if options.output[-1] != '/':
        options.output += '/'
    return options


def dfs2_3d_numpy(list_of_dfs):
    num_of_matrices = len(list_of_dfs)
    res = np.empty(shape=(num_of_matrices, list_of_dfs[0].shape[0], list_of_dfs[0].shape[1]))
    for idx, matrix in enumerate(list_of_dfs):
        res[idx] = matrix
    return res


def read_df_file(f_path):
    if f_path.endswith(".xlsx"):
        data = pd.read_excel(f_path)
    elif f_path.endswith(".csv"):
        data = pd.read_csv(f_path)
    else:
        assert False, "ERROR Parsing GENEpop format file"
    return data

