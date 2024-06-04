import argparse
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
import time

METHODS = ['asd', 'similarity']


def args_parser(parser=None):
    if parser is None:
        parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input",
                        help="input file in genepop (xlsx) format with raw sequencing")
    parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Name of output directory. Program will generate a new directory with that name,"
                             " and all output files will be there")
    parser = add_method_args(parser)
    return handle_args(parser)


def is_valid_method_arg(options):
    if options.method == 'asd' and not isinstance(options.weight, bool):
        raise Exception("ASD metric doesn't support a weighted method. Did you mean to use --method similarity?")
    if options.method not in METHODS:
        raise Exception(f"Metric {options.method} is unknown. Available metrics to compute are {METHODS}")


def handle_args(parser, check_method=True):
    options = parser.parse_args()
    if isinstance(options.weight, str):
        options.weight = float(options.weight)
    if check_method:
        is_valid_method_arg(options)
    os.makedirs(options.output, exist_ok=True)
    return options


def add_method_args(parser):
    parser.add_argument("-w", "--weighted", nargs='?', default=False, const=1, dest="weight",
                        help="If used, compute weighted metric, from Greenbaum et al., 2016. If not use, compute "
                             "unweighted metric from  Li and Horvitz, 1953. If used with a number, it will be used as "
                             "the power of the frequency vector. The default (if used) is linear (same as setting 1) "
                             "and the default if not used is same as setting 0. This flag is available with similarity"
                             " method only. With other method the program will fail.")
    parser.add_argument("--method",
                        help="Method to use. Use 'asd' for allele sharing distance (ASD)"
                             " as in Xiaoyi Gaoa and Eden R. Martin, 2009."
                             " The other option is 'similarity'. If used similarity see the 'weighted' flag for info "
                             "between the 3 options to use.")
    parser.add_argument("--max_memo", dest="max_mb", default=1, type=float,
                        help="Max number of cells (individuals multiply by sites) to use in a single matrix "
                             "(in millions). If you don't know, don't touch. If there are memory failures, reduce it,"
                             "if you have problem of writing too many files increase it."
                             "The default is 1")
    parser.add_argument("--max_threads", dest="max_threads", default=8, type=int,
                        help="Maximum number of threads to compute small matrices simultaneously")
    parser.add_argument("--max_sites", dest="max_sites", default=0, type=int,
                        help="If assigned, compute similarity only based on the first 'max_site' sites")
    parser.add_argument("--format", dest="output_format", default='csv',
                        help='output metric matrix in that format. Default is csv, other option is `ns` or `netstruct`'
                             'for text file in the NetStruct format.')
    return parser


def wait_for_jobs_to_be_done(jobs, max_number_of_threads):
    finished_jobs = []
    while len(jobs) - len(finished_jobs) >= max_number_of_threads:
        time.sleep(.1)
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
    elif f_path.endswith(".csv"):
        df = pd.read_csv(f_path)
    else:
        assert False, "ERROR Parsing GENEpop format file"
    df = filter_out_bad_sites(df)
    return df


def read_vcf_tmp_file(f_path):
    with open(f_path, "r") as f:
        data = f.readlines()
    if data[-1] == '\n':
        data = data[:-1]
    data_split = [e[:-1].split('\t') for e in data]
    line_indices_to_remove = []
    num_of_not_diploid_sites = 0
    for idx, line in enumerate(data_split):
        non_fail = [e for e in line if e != "FAIL"]
        if len(non_fail) == 0:
            line_indices_to_remove.append(idx)
            continue
        alleles = non_fail[0].split('/')
        if len(alleles) != 2:
            num_of_not_diploid_sites += 1
            line_indices_to_remove.append(idx)
            continue
        if all([non_fail[0] == e for e in non_fail]) and alleles[0] == alleles[1]:
            line_indices_to_remove.append(idx)
    data_split = [j for i, j in enumerate(data_split) if i not in line_indices_to_remove]
    df = pd.DataFrame(data_split).T
    return df, num_of_not_diploid_sites


def df2ns_format(similarity_data_frame):
    new_data = ""
    number_of_cells_to_skip = 1
    for col_name, line in similarity_data_frame.iterrows():
        for num in line[number_of_cells_to_skip:]:
            new_data += f'{num} '
        number_of_cells_to_skip += 1
        new_data = new_data[:-1] + '\n'
    return new_data[:-1]

def write_random_vcf(num_indv, num_snps, output_path, num_alleles=2):
    if os.path.exists(output_path):
        os.remove(output_path)
    txt = f"""##fileformat=VCFv4.2
##source=tskit 0.4.1
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1,length={num_snps}>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""
    indv_names = '\t'.join([f'indv_{i}' for i in range(num_indv)])
    txt += f"#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  {indv_names}\n"
    for snp in range(num_snps):
        rnd_samples = np.random.randint(num_alleles, size=num_indv * 2).reshape(-1, 2)
        line = f"1\t{snp}\t.\t0\t1\t.\tPASS\t.\tGT\t"
        lst_str_samples = [f'{e[0]}|{e[1]}' for e in rnd_samples]
        num_fails = round(np.random.exponential(0.05) * num_indv)
        for _ in range(num_fails):
            lst_str_samples[np.random.randint(num_indv)] = '.'
        line += '\t'.join(lst_str_samples)
        txt += line + '\n'
        if snp % 1000 == 0:
            with open(output_path, "a+") as f:
                f.write(txt)
                txt = ""
    with open(output_path, "a+") as f:
        f.write(txt)
