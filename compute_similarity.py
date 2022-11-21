from utils import args_parser, read_df_file, dfs2_3d_numpy
from time import time
import numpy as np
import matplotlib.pyplot as plt


def assign_allele_numbers(data_df):
    snps_names = list(data_df.columns)
    snps_names.remove("ID")
    max_num_of_alleles = 1
    name_to_num = {s: {} for s in snps_names}
    num_to_name = {s: {} for s in snps_names}
    for snp_name in snps_names:
        snp_data = data_df[snp_name]
        for single_snp_data in snp_data.items():
            if single_snp_data[1] == 'FAIL':
                continue

            alleles = single_snp_data[1].split('/')
            for allele in [e.strip() for e in alleles]:
                if allele in name_to_num[snp_name]:
                    continue
                else:
                    if len(name_to_num[snp_name]) == 0:
                        name_to_num[snp_name][allele] = 1
                        num_to_name[snp_name][1] = allele
                    else:
                        next_idx = max(name_to_num[snp_name].values()) + 1
                        if next_idx > max_num_of_alleles:
                            max_num_of_alleles = next_idx
                        name_to_num[snp_name][allele] = next_idx
                        num_to_name[snp_name][next_idx] = allele

    return num_to_name, name_to_num, max_num_of_alleles
#
# def rows_assign_allele_numbers(data_df):
#     snps_names = list(data_df.columns)
#     snps_names.remove("ID")
#     name_to_num = {s: {} for s in snps_names}
#     num_to_name = {s: {} for s in snps_names}
#     for idx, sample in data_df.iterrows():
#         for snp_name in snps_names:
#             alleles = sample[snp_name]
#             if alleles == 'FAIL':
#                 continue
#             alleles = [e.strip() for e in alleles.split('/')]
#             for allele in alleles:
#                 if allele in name_to_num[snp_name]:
#                     continue
#                 else:
#                     if len(name_to_num[snp_name]) == 0:
#                         name_to_num[snp_name][allele] = 1
#                         num_to_name[snp_name][1] = allele
#                     else:
#                         next_idx = max(name_to_num[snp_name].values()) + 1
#                         name_to_num[snp_name][allele] = next_idx
#                         num_to_name[snp_name][next_idx] = allele
#
#     return num_to_name, name_to_num
#
#


def genepop_to_012matrix(df, num_to_name, max_num_of_alleles):
    dfs = [df.copy() for _ in range(max_num_of_alleles)]
    snps_names = list(df.columns)
    snps_names.remove("ID")
    for allele_num, allele_matrix in enumerate(dfs, start=1):
        del allele_matrix['ID']
        for snp_name in snps_names:
            if allele_num in num_to_name[snp_name]:
                allele_name = num_to_name[snp_name][allele_num]
                allele_matrix[snp_name] = allele_matrix[snp_name].apply(lambda x: x.count(allele_name) if x != "FAIL" else 0)
            else:
                allele_matrix[snp_name] = 0
    return dfs


def matrices012_to_similarity_matrix(input_matrices, weighted=False):
    """
    Get numpy array of input (count of ones of the alleles in 012 format). row per individual, column per site.
    Return similarity matrix (row for individual, column for individual)
    """
    is_valid_matrix = np.sum(input_matrices, axis=0) / 2
    pairwise_count_valid_sites = is_valid_matrix @ is_valid_matrix.T
    similarity = np.zeros_like(pairwise_count_valid_sites)
    for matrix in input_matrices:
        if weighted:
            num_valid_genotypes = np.sum(is_valid_matrix, axis=0)
            allele_count = np.sum(matrix, axis=0)
            freq = allele_count / (2 * num_valid_genotypes)
            similarity += ((1 - freq) * matrix) @ matrix.T
        else:
            similarity += matrix @ matrix.T * 2

    similarity = 1/4 * similarity / pairwise_count_valid_sites
    np.fill_diagonal(similarity, -1)
    return similarity


if __name__ == '__main__':
    arguments = args_parser()
    df = read_df_file(arguments.input)
    num_to_name, name_to_num, max_num_of_alleles = assign_allele_numbers(df)
    dfs_list = genepop_to_012matrix(df, num_to_name, max_num_of_alleles)
    similarity = matrices012_to_similarity_matrix(dfs2_3d_numpy(dfs_list))
    print()