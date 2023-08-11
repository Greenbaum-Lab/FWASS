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
@jit(nopython=True)
def single_site_asd_compute(indv_gt, sim_mat, count_mat):
    len2_sort = lambda x: [x[1], x[0]] if x[1] < x[0] else x
    allele_split = [e.split('/') for e in indv_gt]
    for idx1, gt1 in enumerate(allele_split):
        if gt1 == ["FAIL"]:
            continue
        a = gt1[0]
        b = gt1[1]
        idx2 = idx1 + 1
        for gt2 in allele_split[idx1 + 1:]:
            if gt2 == ["FAIL"]:
                idx2 += 1
                continue
            count_mat[idx1, idx2] += 1
            if len2_sort([a, b]) == len2_sort(gt2):
                pass  # sim_mat[idx1, idx2] += 0
            elif a not in gt2 and b not in gt2:
                sim_mat[idx1, idx2] += 2
            else:
                sim_mat[idx1, idx2] += 1
            idx2 += 1
    return


@jit(nopython=True)
def single_site_sim_weighted_compute(indv_gt, freq_weight, sim_mat, count_mat):
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
        idx2 = idx1 + 1
        for gt2 in allele_split[idx1 + 1:]:
            if gt2 == ["FAIL"]:
                idx2 += 1
                continue
            sim_mat[idx1, idx2] += gt2.count(a) * nice_freq[a] + gt2.count(b) * nice_freq[b]
            count_mat[idx1, idx2] += 1
            idx2 += 1
    return

@jit(nopython=True)
def single_site_sim_compute(indv_gt, sim_mat, count_mat):
    allele_split = [e.split('/') for e in indv_gt]
    for idx1, gt1 in enumerate(allele_split):
        if gt1 == ["FAIL"]:
            continue
        a = gt1[0]
        b = gt1[1]
        idx2 = idx1 + 1
        for gt2 in allele_split[idx1 + 1:]:
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

    def __init__(self, arguments, input, output):
        self.arguments = arguments
        self.input_data_dir = input
        os.makedirs(output, exist_ok=True)
        self.output_dir = output
        self.single_site_method = self.method_name2func[arguments.method]
        self.single_site_method(["FAIL"], np.zeros((2, 2)), np.zeros((2, 2)))  # Compile out of time measure
        self.multiplication_constant = self.method_name2constant_multiplication[arguments.method]
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
        pbar = tqdm(desc="Run over sites")
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
            pbar.update(1)
            last_line = f.readline().decode()
            while last_line:
                line = last_line.split()
                assert(len(line[9:]) == len(individuals) and line[format_dict["FORMAT"]].startswith("GT"))
                indv_gt = ['FAIL' if '.' in e[:3] else e[:3].replace('|', '/') for e in line[9:]]
                if not all([e == indv_gt[0] for e in indv_gt]):
                    self.single_site_method(indv_gt, metric_mat, count_mat)
                sites_counter += 1
                last_line = f.readline().decode()
                pbar.update(1)
        metric_mat += metric_mat.T
        count_mat += count_mat.T
        np.fill_diagonal(count_mat, val=1)
        metric = np.true_divide(metric_mat, count_mat)
        metric *= self.multiplication_constant
        np.fill_diagonal(metric, val=0)
        df = pd.DataFrame(metric, index=individuals, columns=individuals)
        self.save_output(file_name, df, iter_num)
        self.times.append(time.time() - start_time)

    def save_output(self, file_name, df, iter_num):
        base_name = os.path.join(self.output_dir, os.path.basename(file_name)[:-4])
        output_path = base_name + f'_{iter_num}.csv'
        df.to_csv(output_path)


def main():
    arguments = args_parser()
    naive = Naive(arguments, input=arguments.input, output=arguments.output)
    file_path = "C:\\Users\\shaharma\\LAB\\test.vcf"
    if not os.path.exists(file_path):
        write_random_vcf(100, 10000, file_path, num_alleles=2)
    naive.vcf_to_metric('VCF', file_path, 0)
    print("Done computation")


if __name__ == '__main__':
    main()
