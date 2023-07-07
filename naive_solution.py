import os
import gzip
from tqdm import tqdm
import numpy as np
import time

from utils import args_parser, parse_input_file_format


def vcf_to_metric(input_format, options):
    """
    Parse the VCF file, and save small matrices to compute on different threads.
    :param input_format: "VCF-GZ" if it compressed as gzip, "VCF" if it's a text file in a VCF format.
    :param options: running arguments
    :param mid_outputs_path: path of directory to throw mid-results for computations.
    :return: None, The function parse the file and save the small matrices in a mid_res directory.
    """
    read_file_func = gzip.open if "GZ" in input_format else open
    print(f"Analyzing VCF file {options.input}")
    pbar = tqdm(desc="Run over sites")
    with read_file_func(options.input, "rb") as f:
        last_line = f.readline().decode()
        while last_line.startswith("##"):
            last_line = f.readline().decode()
        line = last_line.split()
        individuals = line[9:]
        format_dict = {x: line.index(x) for x in ["ID", "FORMAT", "#CHROM", "POS"]}
        num_of_indv = len(individuals)
        similarity_mat = np.zeros(shape=(num_of_indv, num_of_indv))
        count_mat = np.zeros(shape=(num_of_indv, num_of_indv))
        sites_counter = 1
        pbar.update(1)
        last_line = f.readline().decode()
        while last_line:
            line = last_line.split()
            assert(len(line[9:]) == len(individuals) and line[format_dict["FORMAT"]].startswith("GT"))
            indv_gt = ['FAIL' if '.' in e[:3] else e[:3].replace('|', '/') for e in line[9:]]
            single_site_compute(indv_gt, similarity_mat, count_mat)
            sites_counter += 1
            last_line = f.readline().decode()
            pbar.update(1)
    similarity_mat += similarity_mat.T
    count_mat += count_mat.T
    np.fill_diagonal(count_mat, val=1)
    return np.true_divide(similarity_mat, count_mat * 4)


def single_site_freq_compute(indv_gt, freq_weight, sim_mat, count_mat):
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
        for idx2, gt2 in enumerate(allele_split[idx1:]):
            if gt2 == ["FAIL"]:
                continue
            sim_mat[idx1, idx2] += gt2.count(a) * nice_freq[a] + gt2.count(b) * nice_freq[b]
            count_mat[idx1, idx2] += 1

    return


def single_site_compute(indv_gt, sim_mat, count_mat):
    allele_split = [e.split('/') for e in indv_gt]
    for idx1, gt1 in enumerate(allele_split):
        if gt1 == ["FAIL"]:
            continue
        a = gt1[0]
        b = gt1[1]
        for idx2, gt2 in enumerate(allele_split[idx1:]):
            if gt2 == ["FAIL"]:
                continue
            sim_mat[idx1, idx2] += int(a == gt2[0]) + int(a == gt2[1]) + int(b == gt2[0]) + int(b == gt2[1])
            count_mat[idx1, idx2] += 1
    return

def main():
    s_time = time.time()
    arguments = args_parser()
    input_format = parse_input_file_format(arguments.input)
    if 'VCF' in input_format:
        vcf_to_metric(input_format, arguments)
    print("Done computation")
    print(f"It all took {(time.time() - s_time) / 60} minutes!")


if __name__ == '__main__':
    main()
