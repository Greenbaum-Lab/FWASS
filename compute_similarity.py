import gzip
import time
from multiprocessing import Process
import pandas as pd
from tqdm import tqdm
from utils import args_parser, read_df_file, parse_input_file_format, wait_for_jobs_to_be_done, df2ns_format, \
    read_vcf_tmp_file
import numpy as np
import shutil
import os
from numba import jit


def assign_allele_numbers(data_df, loci_names):
    """
    Give a natural number (including 0) to each allele, without loss of generality.
    The function meanwhile build matrices as number of the maximum number of alleles per locus.
    It also builds the 3D matrix on the way.
    :param data_df: Raw data frame. Every diploid sample is mark as "x/y" where x and y are some alleles,
     not necessarily a single letter, marked as any way. Every invalid sample is filled with "FAIL"
     The data_df param is in shape(N*M), every row is individual, and every column is a locus.
    :param loci_names: list of loci names
    :return:  a 3D matrix of shape (T*N*M) where T is maximum number of alleles per locus,
     N is number of individuals, and M is number of loci.
    """
    new_arrs = [np.zeros(shape=(data_df.shape[0], len(loci_names)))]
    max_num_of_alleles = 0
    name_to_num = {s: {} for s in loci_names}
    for locus_idx, locus_name in enumerate(loci_names):
        locus_data = data_df[locus_name]
        for single_locus_data in locus_data.items():
            if single_locus_data[1] == 'FAIL':
                continue

            alleles = single_locus_data[1].split('/')
            for allele in [e.strip() for e in alleles]:
                if allele in name_to_num[locus_name]:
                    new_arrs[name_to_num[locus_name][allele]][single_locus_data[0]][locus_idx] += 1
                    continue
                else:
                    if len(name_to_num[locus_name]) == 0:
                        name_to_num[locus_name][allele] = 0
                        new_arrs[0][single_locus_data[0]][locus_idx] += 1
                    else:
                        next_idx = max(name_to_num[locus_name].values()) + 1
                        if next_idx > max_num_of_alleles:
                            max_num_of_alleles = next_idx
                            new_arrs.append(np.zeros(shape=(data_df.shape[0], len(loci_names))))
                        name_to_num[locus_name][allele] = next_idx
                        new_arrs[next_idx][single_locus_data[0]][locus_idx] += 1
    np_arr = np.array(new_arrs)

    return np_arr, max_num_of_alleles

def vcf_to_small_matrices(input_format, options, mid_outputs_path, input_file_path, verbose):
    """
    Parse the VCF file, and save small matrices to compute on different threads.
    :param input_format: "VCF-GZ" if it compressed as gzip, "VCF" if it's a text file in a VCF format.
    :param options: running arguments
    :param mid_outputs_path: path of directory to throw mid-results for computations.
    :return: None, The function parse the file and save the small matrices in a mid_res directory.
    """
    read_file_func = gzip.open if "GZ" in input_format else open
    max_num_of_cells = options.max_mb * 10**6
    # pbar = tqdm(desc="Run over sites")
    with read_file_func(input_file_path, "rb") as f:
        last_line = f.readline().decode()
        while last_line.startswith("##"):
            last_line = f.readline().decode()
        line = last_line.split()
        individuals = line[9:]
        format_dict = {x: line.index(x) for x in ["ID", "FORMAT", "#CHROM", "POS"]}
        num_of_indv = len(individuals)
        num_sites_to_read = int(max_num_of_cells // num_of_indv)
        matrices_counter = 0
        is_matrix_empty = True
        current_matrix = ""
        sites_counter = 1
        # pbar.update(1)
        last_line = f.readline().decode()
        while last_line:
            if sites_counter % num_sites_to_read == 0:
                with open(os.path.join(mid_outputs_path, f'mat_{matrices_counter}.tmp'), "w") as g:
                    g.write(current_matrix)
                current_matrix = ""
                matrices_counter += 1
            line = last_line.split()
            assert(len(line[9:]) == len(individuals) and line[format_dict["FORMAT"]].startswith("GT"))
            indv_gt = ['FAIL' if '.' in e.split(':')[0] else e.split(':')[0].replace('|', '/') for e in line[9:]]
            current_matrix += "\t".join(indv_gt) + '\n'
            is_matrix_empty = False
            sites_counter += 1
            last_line = f.readline().decode()
            # pbar.update(1)
            if options.max_sites and sites_counter > options.max_sites:
                break

        if not is_matrix_empty:
            with open(os.path.join(mid_outputs_path, f'mat_{matrices_counter}.tmp'), "w") as g:
                g.write(current_matrix)

        with open(os.path.join(mid_outputs_path, "done.txt"), "w") as f:
            f.write("\t".join(individuals))
        if verbose:
            print(f"Done reading VCF file!\n{sites_counter -1} sites in total, divided over {matrices_counter + 1} matrices."
              f" Max of {num_sites_to_read} in each matrix.\n Currently computing the last few similarity matrices.")


def matrices012_to_similarity_matrix(input_matrices, weighted, method):
    """
    Compute similarity (see options in README)
    :param method: Name of metric to use. Supports the list in utils.METHODS
    :param input_matrices: 3D of 012 matrix, first axis is per allele number,
     second axis is per individual, third axis is per locus. In cell[x,y,z] there is a count
     of how may times do individual y has the allele x in locus z. It can be up to 2.
    :param weighted: If False (default) use relatedness measure from (Li and Horvitz, 1953).
     If True, use frequency weighted metric from (Greenbaum et al., 2016).
    :return: Similarity matrix (rows are individuals, columns are individuals), and count matrix which
    is also individuals * individuals matrix, which tells us on how many sites each pair was computed.
    This verifies between pairs thanks to invalid data.
    """
    is_valid_matrix = np.sum(input_matrices, axis=0) / 2
    pairwise_count_valid_sites = is_valid_matrix @ is_valid_matrix.T
    similarity = np.zeros_like(pairwise_count_valid_sites)
    for matrix in input_matrices:
        if method == 'asd':
            bin_matrix = (matrix > 0).astype(int)
            similarity += bin_matrix @ bin_matrix.T
            similarity += np.clip(matrix - 1, 0, 1) @ np.clip(matrix - 1, 0, 1).T
        elif method == 'similarity':
            if weighted:
                num_valid_genotypes = np.sum(is_valid_matrix, axis=0)
                allele_count = np.sum(matrix, axis=0)
                freq = allele_count / (2 * num_valid_genotypes)
                similarity += (((1 - np.nan_to_num(freq)) ** weighted) * matrix) @ matrix.T
            else:
                similarity += matrix @ matrix.T * 2

    if method == 'asd':
        return 2 * pairwise_count_valid_sites - similarity, pairwise_count_valid_sites
    elif method == 'similarity':
        return 1/4 * similarity, pairwise_count_valid_sites


def save_numpy_outputs(options, sim_mat, count_mat, output_dir):
    """
    Save numpy arrays as npy files for later computations
    :param options: running arguments
    :param sim_mat: similarity matrix
    :param count_mat: counts matrix
    :param output_dir: directory to write outputs
    :return: None
    """
    os.makedirs(output_dir, exist_ok=True)
    raw_sim_out_path = os.path.join(output_dir, "raw_similarity" + '_weighted' * bool(options.weight) + '.npy')
    count_out_path = os.path.join(output_dir, "count.npy")
    np.save(raw_sim_out_path, sim_mat)
    np.save(count_out_path, count_mat)


def save_ns_format(sim_df, counts_df, sim_pref, count_pref):
    similarity = df2ns_format(sim_df)
    counts = df2ns_format(counts_df)
    with open(sim_pref + '.txt', "w") as f:
        f.write(similarity)
    with open(count_pref + '.txt', 'w') as f:
        f.write(counts)

def save_df_outputs(options, sim_mat, count_mat, individual_names, output_dir):
    """
    Save the numpy array as csv files with names of individuals.
    :param options: running arguments
    :param sim_mat: similarity matrix
    :param count_mat: counts matrix
    :param individual_names: list of individual names
    :param output_dir: directory to write outputs
    :return: None
    """
    os.makedirs(output_dir, exist_ok=True)
    sim_out_pref = os.path.join(output_dir, "similarity" + '_weighted' * bool(options.weight))
    count_out_pref = os.path.join(output_dir, "count")
    counts_df = pd.DataFrame(count_mat, columns=individual_names, index=individual_names)
    similarity_df = pd.DataFrame(sim_mat / count_mat, columns=individual_names, index=individual_names)
    if options.ns_format:
        save_ns_format(similarity_df, counts_df, sim_out_pref, count_out_pref)
    else:
        similarity_df.to_csv(sim_out_pref + '.csv')
        counts_df.to_csv(count_out_pref + '.csv')


def submit_all_matrices_ready(options, last_computed_matrix, mid_outputs_path, procs):
    """
    Check if the parser job prepare new small matrices outputs, if so,
     submit jobs to compute their similarities.
    :param options: running arguments
    :param last_computed_matrix: last index of computed matrix
    :param mid_outputs_path: mid-results output directory
    :param procs: list of running processes
    :return: new last computed index matrix and new list of processes
    """
    files = [e for e in os.listdir(mid_outputs_path) if os.path.isfile(os.path.join(mid_outputs_path, e)) and e != "done.txt"]
    files_numbers = [int(e.split('.')[0][4:]) for e in files]
    files_to_process = sorted([i for i in files_numbers if i > last_computed_matrix])
    if files_to_process:
        last_computed_matrix = max(files_to_process) - 1
    for idx in files_to_process:
        if idx == max(files_to_process):
            continue
        procs = wait_for_jobs_to_be_done(procs, options.max_threads)
        input_name = os.path.join(mid_outputs_path, f"mat_{idx}.tmp")
        output_name = os.path.join(mid_outputs_path, f"similarity_{idx}")
        run_similarity_job = Process(target=vcf_single_job, args=(options, input_name, output_name))
        procs.append(run_similarity_job)
        run_similarity_job.start()
    return last_computed_matrix, procs


def get_sim_and_count_from_directory(options, directory):
    """
    Read raw similarity (not divided by counts yet) and counts array from a certain directory.
    Used from mid-computations.
    :param options: running arguments
    :param directory: path to directory to read from.
    :return: numpy array of raw similarity matrix, and numpy array of count matrix
    """
    raw_sim_path = os.path.join(directory, "raw_similarity" + '_weighted' * bool(options.weight) + '.npy')
    count_path = os.path.join(directory, "count.npy")
    sim_np = np.load(raw_sim_path)
    count_np = np.load(count_path)
    return sim_np, count_np


def merge_matrices(options, individual_names):
    """
    Sum all raw similarity matrices (not divided by counts),
    sum all counts matrices, and then divided the sums.
    Resulting a true average similarity computations of all sites.
    :param options: running arguments
    :param individual_names: list of individual names
    :return: None, it saves the similarity and counts matrices in the output directory as csv files.
    """
    sum_of_non_diploid = 0
    mid_outputs_path = os.path.join(options.output, "mid_res")
    similarity_dirs = [os.path.join(mid_outputs_path, e) for e in os.listdir(mid_outputs_path) if 'similarity' in e]
    sim_np, count_np = get_sim_and_count_from_directory(options, similarity_dirs[0])
    for directory in similarity_dirs[1:]:
        with open(os.path.join(directory, 'monoploid.txt'), "r") as f:
            sum_of_non_diploid += int(f.read())
        new_sim_np, new_count_np = get_sim_and_count_from_directory(options, directory)
        sim_np += new_sim_np
        count_np += new_count_np
    similarity = sim_np / count_np
    np.fill_diagonal(similarity, val=0)
    similarity_df = pd.DataFrame(similarity, columns=individual_names, index=individual_names)
    counts_df = pd.DataFrame(count_np, columns=individual_names, index=individual_names)
    simi_path_prefix = os.path.join(options.output, "similarity" + '_weighted' * bool(options.weight))
    counts_path_pref = os.path.join(options.output, 'counts')
    if options.save_outputs:
        if options.ns_format:
            save_ns_format(similarity_df, counts_df, simi_path_prefix, counts_path_pref)
        else:
            similarity_df.to_csv(simi_path_prefix + '.csv')
            counts_df.to_csv(counts_path_pref + '.csv')
    with open(os.path.join(similarity_dirs[0], 'monoploid.txt'), "r") as f:
        sum_of_non_diploid += int(f.read())
    if sum_of_non_diploid > 0:
        print(f"WARNING! : There are {sum_of_non_diploid} sites with number of alleles per individual different than 2!"
              f"These sites were ignored")
    shutil.rmtree(mid_outputs_path)
    return similarity_df


def analyze_by_vcf(input_format, options, input_file_name=None, verbose=True):
    """
    compute similarity from a VCF file
    :param input_format: VCF-GZ for a VCF file compressed with gzip, VCF for uncompressed text file in VCF format.
    :param options: running arguments
    :return: None, it saves the similarity and counts matrices in the output directory.
    """
    if input_file_name is None:
        input_file_name = options.input
    mid_outputs_path = os.path.join(options.output, "mid_res")
    os.makedirs(mid_outputs_path, exist_ok=True)
    reader_job = Process(target=vcf_to_small_matrices, args=(input_format, options, mid_outputs_path, input_file_name,
                                                             verbose))
    procs = []
    reader_job.start()
    last_computed_matrix = -1
    done_reading_file = os.path.join(mid_outputs_path, "done.txt")
    while not os.path.exists(done_reading_file):
        time.sleep(.1)
        last_computed_matrix, procs = submit_all_matrices_ready(options, last_computed_matrix, mid_outputs_path, procs)
    reader_job.join()
    last_computed_matrix, procs = submit_all_matrices_ready(options, last_computed_matrix, mid_outputs_path, procs)
    idx = last_computed_matrix + 1
    input_name = os.path.join(mid_outputs_path, f"mat_{idx}.tmp")
    output_name = os.path.join(mid_outputs_path, f"similarity_{idx}")
    vcf_single_job(options, input_name, output_name)
    for proc in procs:
        proc.join()
    with open(done_reading_file, "r") as f:
        individual_names = f.read().strip().split('\t')
    return merge_matrices(options, individual_names)


def vcf_single_job(options, input_name, output_dir):
    """
    Compute similarity on a mid-results small matrix
    :param options: running arguments
    :param input_name: path of data file. The data is a matrix in a text format made in the VCF parser job.
    :param output_dir: directory to save raw similarity matrix and counts matrix
    :return: None, it saves the outputs in output_dir.
    """
    s_time = time.time()
    df, num_of_not_diploid_sites = read_vcf_tmp_file(input_name)
    np_3d_arr, max_num_of_alleles = assign_allele_numbers(df, df.columns)
    similarity, counts = matrices012_to_similarity_matrix(np_3d_arr, options.weight, options.method)
    save_numpy_outputs(options, similarity, counts, output_dir)
    with open(os.path.join(output_dir, "time.txt"), "w") as f:
        f.write(str(time.time() - s_time))
    with open(os.path.join(output_dir, "monoploid.txt"), "w") as f:
        f.write(str(num_of_not_diploid_sites))

def analyze_by_df(options):
    """
    compute similarity from a genepop file.
    First column should be named "ID" with individual names under it.
    Rest columns names are the sites names.
    Every cell should be "x/y" where x and y are alleles marked any desired way.
    Missing data should be mark as "FAIL"
    :param options: running arguments
    :return: None, it saves the similarity and counts matrices in the output directory.
    """
    df = read_df_file(options.input)
    loci_names = list(df.columns)
    loci_names.remove("ID")
    np_3d_arr, max_num_of_alleles = assign_allele_numbers(df, loci_names)
    similarity, counts = matrices012_to_similarity_matrix(np_3d_arr, options.weight, options.method)
    individual_names = list(df['ID'])
    save_df_outputs(options, similarity, counts, individual_names, options.output)


def main():
    arguments = args_parser()
    input_format = parse_input_file_format(arguments.input)
    if 'VCF' in input_format:
        analyze_by_vcf(input_format, arguments)
    elif 'DF' in input_format:
        analyze_by_df(arguments)


if __name__ == '__main__':
    main()