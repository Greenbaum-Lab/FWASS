import gzip
import time
from multiprocessing import Process
import pandas as pd
from tqdm import tqdm
from utils import args_parser, read_df_file, dfs2_3d_numpy, parse_input_file_format
import numpy as np
import os


def assign_allele_numbers(data_df, snps_names):
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


def txtdf_to_012matrix(df, num_to_name, max_num_of_alleles):
    dfs = [df.copy() for _ in range(max_num_of_alleles)]
    snps_names = list(df.columns)
    for allele_num, allele_matrix in enumerate(dfs, start=1):
        for snp_name in snps_names:
            if allele_num in num_to_name[snp_name]:
                allele_name = num_to_name[snp_name][allele_num]
                allele_matrix[snp_name] = allele_matrix[snp_name].apply(lambda x: x.count(allele_name) if x != "FAIL" else 0)
            else:
                allele_matrix[snp_name] = 0
    return dfs


def vcf_to_small_matrices(input_format, options, mid_outputs_path):
    read_file_func = gzip.open if "GZ" in input_format else open
    max_num_of_cells = options.max_mb * 10**6
    pbar = tqdm(desc="Run over sites: ")
    with read_file_func(options.input, "rb") as f:
        last_line = f.readline().decode()
        while last_line.startswith("##"):
            last_line = f.readline().decode()
        line = last_line.split()
        individuals = line[9:]
        format_dict = {x: line.index(x) for x in ["ID", "FORMAT", "#CHROM", "POS"]}
        num_of_indv = len(individuals)
        num_sites_to_read = max_num_of_cells // num_of_indv
        matrices_counter = 0
        sites_names = []
        current_matrix = ""
        sites_counter = 1
        pbar.update(1)
        last_line = f.readline().decode()
        while last_line:
            if sites_counter % num_sites_to_read == 0:
                # df = pd.DataFrame(list(zip(*current_matrix)), index=individuals, columns=sites_names)
                # df.index.name = "ID"
                # df.to_csv(os.path.join(mid_outputs_path, f'mat_{matrices_counter}.csv'))
                with open(os.path.join(mid_outputs_path, f'mat_{matrices_counter}.txt'), "w") as g:
                    g.write(current_matrix)
                sites_names = []
                current_matrix = ""
                matrices_counter += 1
            line = last_line.split()
            assert(len(line[9:]) == len(individuals) and line[format_dict["FORMAT"]].startswith("GT"))
            indv_gt = ['FAIL' if '.' in e[:3] else e[:3].replace('|', '/') for e in line[9:]]
            current_matrix += "\t".join(indv_gt) + '\n'
            site_uniq_name = f"{line[format_dict['#CHROM']]}_{line[format_dict['POS']]}_{line[format_dict['ID']]}"
            sites_names.append(site_uniq_name)
            last_line = f.readline().decode()
            sites_counter += 1
            pbar.update(1)

        if sites_names:
            # df = pd.DataFrame(list(zip(*current_matrix)), index=individuals, columns=sites_names)
            # df.index.name = "ID"
            # df.to_csv(os.path.join(mid_outputs_path, f'mat_{matrices_counter}.csv'))
            with open(os.path.join(mid_outputs_path, f'mat_{matrices_counter}.txt'), "w") as g:
                g.write(current_matrix)

        with open(os.path.join(mid_outputs_path, "done.txt"), "w") as f:
            f.write("\t".join(individuals))
        print(f"DONE!\n{sites_counter} sites in total, divided over {matrices_counter} matrices."
              f" {num_sites_to_read} in each matrix.")


def matrices012_to_similarity_matrix(input_matrices, weighted):
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

    return 1/4 * similarity, pairwise_count_valid_sites


def save_numpy_outputs(options, sim_mat, count_mat, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    raw_sim_out_path = os.path.join(output_dir, "raw_similarity" + '_weighted' * options.weighted_metric + '.npy')
    count_out_path = os.path.join(output_dir, "count.npy")
    np.save(raw_sim_out_path, sim_mat)
    np.save(count_out_path, count_mat)


def save_df_outputs(options, sim_mat, count_mat, individual_names, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    raw_sim_out_path = os.path.join(output_dir, "raw_similarity" + '_weighted' * options.weighted_metric + '.csv')
    sim_out_path = os.path.join(output_dir, "similarity" + '_weighted' * options.weighted_metric + '.csv')
    count_out_path = os.path.join(output_dir, "count.csv")
    count_df = pd.DataFrame(count_mat, columns=individual_names, index=individual_names)
    raw_similarity_df = pd.DataFrame(sim_mat, columns=individual_names, index=individual_names)
    similarity_df = pd.DataFrame(sim_mat / count_mat, columns=individual_names, index=individual_names)
    similarity_df.to_csv(sim_out_path)
    raw_similarity_df.to_csv(raw_sim_out_path)
    count_df.to_csv(count_out_path)


def submit_all_matrices_ready(options, last_computed_matrix, mid_outputs_path, procs):
    files = [e for e in os.listdir(mid_outputs_path) if os.path.isfile(os.path.join(mid_outputs_path, e)) and e != "done.txt"]
    files_numbers = [int(e.split('.')[0][4:]) for e in files]
    files_to_process = sorted([i for i in files_numbers if i > last_computed_matrix])
    time.sleep(2)
    if files_to_process:
        last_computed_matrix = max(files_to_process) - 1
    for idx in files_to_process:
        if idx == max(files_to_process):
            continue
        input_name = os.path.join(mid_outputs_path, f"mat_{idx}.txt")
        output_name = os.path.join(mid_outputs_path, f"similarity_{idx}")
        run_similarity_job = Process(target=vcf_single_job, args=(options, input_name, output_name))
        procs.append(run_similarity_job)
        run_similarity_job.start()
    return last_computed_matrix, procs


def get_sim_and_count_from_directory(options, directory):
    raw_sim_path = os.path.join(directory, "raw_similarity" + '_weighted' * options.weighted_metric + '.npy')
    count_path = os.path.join(directory, "count.npy")
    sim_np = np.load(raw_sim_path)
    count_np = np.load(count_path)
    return sim_np, count_np


def merge_matrices(options, individual_names):
    mid_outputs_path = os.path.join(options.output, "mid_res")
    similarity_dirs = [os.path.join(mid_outputs_path, e) for e in os.listdir(mid_outputs_path) if 'similarity' in e]
    sim_np, count_np = get_sim_and_count_from_directory(options, similarity_dirs[0])
    for directory in similarity_dirs[1:]:
        new_sim_np, new_count_np = get_sim_and_count_from_directory(options, directory)
        sim_np += new_sim_np
        count_np += new_count_np
    similarity = sim_np / count_np
    similarity_df = pd.DataFrame(similarity, columns=individual_names, index=individual_names)
    similarity_df.to_csv(os.path.join(options.output, "similarity" + '_weighted' * options.weighted_metric + '.csv'))


def analyze_by_vcf(input_format, options):
    mid_outputs_path = os.path.join(options.output, "mid_res")
    os.makedirs(mid_outputs_path, exist_ok=True)
    reader_job = Process(target=vcf_to_small_matrices, args=(input_format, options, mid_outputs_path))
    procs = [reader_job]
    reader_job.start()
    last_computed_matrix = -1
    done_reading_file = os.path.join(mid_outputs_path, "done.txt")
    while not os.path.exists(done_reading_file):
        last_computed_matrix, jobs = submit_all_matrices_ready(options, last_computed_matrix, mid_outputs_path, procs)
    last_computed_matrix, jobs = submit_all_matrices_ready(options, last_computed_matrix, mid_outputs_path, procs)
    for proc in procs:
        proc.join()
    idx = last_computed_matrix + 1
    input_name = os.path.join(mid_outputs_path, f"mat_{idx}.txt")
    output_name = os.path.join(mid_outputs_path, f"similarity_{idx}")
    vcf_single_job(options, input_name, output_name)
    with open(done_reading_file, "r") as f:
        individual_names = f.read().strip().split('\t')
    merge_matrices(options, individual_names)


def vcf_single_job(options, input_name, output_name):
    s_time = time.time()
    df = read_df_file(input_name)
    num_to_name, name_to_num, max_num_of_alleles = assign_allele_numbers(df, df.columns)
    dfs_list = txtdf_to_012matrix(df, num_to_name, max_num_of_alleles)
    similarity, counts = matrices012_to_similarity_matrix(dfs2_3d_numpy(dfs_list), options.weighted_metric)
    save_numpy_outputs(options, similarity, counts, output_name)
    with open(os.path.join(output_name, "time.txt"), "w") as f:
        f.write(str(time.time() - s_time))


def analyze_by_df(options, input_file, output_dir):
    s_time = time.time()
    df = read_df_file(input_file)
    snps_names = list(df.columns); snps_names.remove("ID")
    num_to_name, name_to_num, max_num_of_alleles = assign_allele_numbers(df, snps_names)
    dfs_list = genepop_to_012matrix(df, num_to_name, max_num_of_alleles)
    similarity, counts = matrices012_to_similarity_matrix(dfs2_3d_numpy(dfs_list), options.weighted_metric)
    individual_names = list(df['ID'])
    save_df_outputs(options, similarity, counts, individual_names, output_dir)
    with open(os.path.join(output_dir, "time.txt"), "w") as f:
        f.write(str(time.time() - s_time))


if __name__ == '__main__':
    s_time = time.time()
    arguments = args_parser()
    input_format = parse_input_file_format(arguments.input)
    if 'VCF' in input_format:
        analyze_by_vcf(input_format, arguments)
    elif 'DF' in input_format:
        analyze_by_df(arguments, arguments.input, arguments.output)
    print(f"It all took {(time.time() - s_time) / 60} minutes!")
