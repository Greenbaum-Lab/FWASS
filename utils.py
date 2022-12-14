import argparse
import os
import pandas as pd
import numpy as n
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
    parser.add_argument("--max_memo", dest="max_mb", default=10, type=float,
                        help="Max number of cells (individuals multiply by sites) to use in a single matrix "
                             "(in millions). If you don't know, don't touch. If there are memory failures, reduce it."
                             " Default is 10")
    parser.add_argument("--max_threads", dest="max_threads", default=8, type=int,
                        help="Maximum number of threads to compute small matrices simultaneously")
    parser.add_argument("--max_sites", dest="max_sites", default=0, type=int,
                        help="If assigned, compute similarity only based on the first n sites")
    options = parser.parse_args()
    if options.output[-1] != '/':
        options.output += '/'
    os.makedirs(options.output, exist_ok=False)
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


def parse_input_file_format(input):
    if input.endswith(".xlsx") or input.endswith(".csv"):
        return "DF"
    if input.endswith("vcf.gz"):
        return "VCF-GZ"
    if input.endswith(".vcf"):
        return "VCF"
    raise IOError("File format is not supported")


def filter_out_bad_sites(df):
    loci_to_remove = []
    for locus_name, col in df.items():
        if locus_name == 'ID':
            continue
        non_fail = col[col != 'FAIL']
        if non_fail.size == 0:
            loci_to_remove.append(locus_name)
            continue
        first_individual = [e.strip() for e in non_fail.iloc[0].split('/')]
        if first_individual[0] == first_individual[1]:
            same_as_first = non_fail[non_fail == non_fail.iloc[0]]
            if same_as_first.size == non_fail.size:
                loci_to_remove.append(locus_name)
    if loci_to_remove:
        df = df.drop(loci_to_remove, axis=1)
        num_of_sites = len(df.columns) - int('ID' in df.columns)
        print(f"Removed {len(loci_to_remove)} invalid loci. Compute similarity based on {num_of_sites} sites.")
    return df


def read_df_file(f_path):
    if f_path.endswith(".xlsx"):
        df = pd.read_excel(f_path)
        df = filter_out_bad_sites(df)
    elif f_path.endswith(".csv"):
        df = pd.read_csv(f_path)
        df = filter_out_bad_sites(df)
    elif f_path.endswith('.tmp'):
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


