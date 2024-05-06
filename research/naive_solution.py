import os
import gzip
from tqdm import tqdm
import numpy as np
import time
import pandas as pd

from utils import args_parser, parse_input_file_format, write_random_vcf
from numba import jit
import warnings
warnings.filterwarnings("ignore")

# Slower if compiled! Don't use numba
# @jit(nopython=True)
def single_site_asd_compute(indv_gt, asd_mat, count_mat, weight=None):
    len2_sort = lambda x: [x[1], x[0]] if x[1] < x[0] else x
    allele_split = [e.split('/') for e in indv_gt]
    for idx1, gt1 in enumerate(allele_split):
        if gt1 == ["FAIL"]:
            continue
        count_mat[idx1, idx1] += 0.5   # based on the fact its a diploid
        a = gt1[0]
        b = gt1[1]
        idx2 = idx1
        for gt2 in allele_split[idx1:]:
            if gt2 == ["FAIL"]:
                idx2 += 1
                continue
            count_mat[idx1, idx2] += 1
            if len2_sort([a, b]) == len2_sort(gt2):
                pass  # sim_mat[idx1, idx2] += 0
            elif a not in gt2 and b not in gt2:
                asd_mat[idx1, idx2] += 2
            else:
                asd_mat[idx1, idx2] += 1
            idx2 += 1
    return


# @jit(nopython=True)
def single_site_sim_weighted_compute(indv_gt, sim_mat, count_mat, freq_weight=1):
    allele_split = [e.split('/') for e in indv_gt]
    non_fails = [e for e in allele_split if e != ["FAIL"]]
    if len(non_fails) == 0:
        return
    all_alleles = [j for sub in non_fails for j in sub]
    alleles = set(all_alleles)
    allele_counts = {a: all_alleles.count(a) for a in alleles}
    freqs = {k: v /sum(allele_counts.values()) for k,v in allele_counts.items()}
    nice_freq = {k: (1 - v) ** freq_weight for k, v in freqs.items()}

    for idx1, gt1 in enumerate(allele_split):
        if gt1 == ["FAIL"]:
            continue
        a = gt1[0]
        b = gt1[1]
        idx2 = idx1
        for gt2 in allele_split[idx1:]:
            if gt2 == ["FAIL"]:
                idx2 += 1
                continue
            sim_mat[idx1, idx2] += gt2.count(a) * nice_freq[a] + gt2.count(b) * nice_freq[b]
            count_mat[idx1, idx2] += 1
            idx2 += 1
    return

@jit(nopython=True)
def single_site_sim_compute(indv_gt, sim_mat, count_mat, weight=None):
    allele_split = [e.split('/') for e in indv_gt]
    for idx1, gt1 in enumerate(allele_split):
        if gt1 == ["FAIL"]:
            continue
        a = gt1[0]
        b = gt1[1]
        idx2 = idx1
        for gt2 in allele_split[idx1:]:
            if gt2 == ["FAIL"]:
                idx2 += 1
                continue
            sim_mat[idx1, idx2] += int(a == gt2[0]) + int(a == gt2[1]) + int(b == gt2[0]) + int(b == gt2[1])
            count_mat[idx1, idx2] += 1
            idx2 += 1
    return


class Naive:
    method_name2func = {'similarity': single_site_sim_compute,
                        'weighted': single_site_sim_weighted_compute,
                        'asd': single_site_asd_compute}

    method_name2constant_multiplication = {'similarity': .5,
                                           'weighted': .25,
                                           'asd': 1}

    def __init__(self, arguments, output):
        self.arguments = arguments
        os.makedirs(output, exist_ok=True)
        self.output_dir = output
        method_name = 'weighted' if isinstance(arguments.weight, float) else arguments.method
        self.single_site_method = self.method_name2func[method_name]
        self.single_site_method(["FAIL"], np.zeros((2, 2)), np.zeros((2, 2)))  # Compile out of time measure
        self.multiplication_constant = self.method_name2constant_multiplication[method_name]
        self.times = []

    def vcf_to_metric(self, input_format, file_name, iter_num):
        """
        Parse the VCF file, and save small matrices to compute on different threads.
        :param input_format: "VCF-GZ" if it compressed as gzip, "VCF" if it's a text file in a VCF format.
        :param file_name: vcf file name
        :param mid_outputs_path: path of directory to throw mid-results for computations.
        :return: None, The function parse the file and save the small matrices in a mid_res directory.
        """
        start_time = time.time()
        read_file_func = gzip.open if "GZ" in input_format else open
        with read_file_func(file_name, "rb") as f:
            last_line = f.readline().decode()
            while last_line.startswith("##"):
                last_line = f.readline().decode()
            line = last_line.split()
            individuals = line[9:]
            format_dict = {x: line.index(x) for x in ["ID", "FORMAT", "#CHROM", "POS"]}
            num_of_indv = len(individuals)
            metric_mat = np.zeros(shape=(num_of_indv, num_of_indv))
            count_mat = np.zeros(shape=(num_of_indv, num_of_indv))
            sites_counter = 1
            last_line = f.readline().decode()
            while last_line:
                line = last_line.split()
                assert(len(line[9:]) == len(individuals) and line[format_dict["FORMAT"]].startswith("GT"))
                indv_gt = ['FAIL' if '.' in e[:3] else e[:3].replace('|', '/') for e in line[9:]]
                if not all([e == indv_gt[0] for e in indv_gt]):
                    self.single_site_method(indv_gt, metric_mat, count_mat, self.arguments.weight)
                sites_counter += 1
                last_line = f.readline().decode()
        metric_mat += metric_mat.T
        count_mat += count_mat.T
        self.metric = np.true_divide(metric_mat, count_mat)
        self.metric *= self.multiplication_constant
        df = pd.DataFrame(self.metric, index=individuals, columns=individuals)
        np.fill_diagonal(count_mat, count_mat.diagonal() / 2)
        count_df = pd.DataFrame(count_mat, index=individuals, columns=individuals)

        file_name = 'similarity_weighted' if isinstance(self.arguments.weight, float) else self.arguments.method
        df.to_csv(os.path.join(self.output_dir, f'{file_name}.csv'))
        count_df.to_csv(os.path.join(self.output_dir, 'count.csv'))
        self.times.append(time.time() - start_time)
