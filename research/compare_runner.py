import json
import os
import time

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from research.fast_solution import Fast
from research.naive_solution import Naive
from utils import args_parser, write_random_vcf

MOCK = True


def get_file_name(num_indv, num_snps):
    if num_indv > 1e3:
        num_indv_name = f'{num_indv // 1e3}K'
    else:
        num_indv_name = str(num_indv)
    if num_snps > 1e6:
        num_snps_name = f'{int(num_snps / 1e6)}M'
    elif num_snps > 1e3:
        num_snps_name = f'{int(num_snps // 1e3)}K'
    else:
        num_snps_name = str(num_snps)
    return f'indv{num_indv_name}_snps{num_snps_name}.vcf'


def file_name_to_indv(file_name):
    num_snps = file_name[file_name.find('snps') + 4:]
    if num_snps.endswith('K'):
        return int(num_snps[:-1]) * 1e3
    elif num_snps.endswith('M'):
        return int(num_snps[:-1]) * 1e6
    else:
        return int(num_snps)

def compare_different_num_of_snps(arguments):
    max_alleles_num = [4, 10, 50]
    avg_times = {'naive': [], 'fast': []}; avg_times.update({f'fast {x}': [] for x in max_alleles_num})
    std_times = {'naive': [], 'fast': []}; std_times.update({f'fast {x}': [] for x in max_alleles_num})
    num_indv = 100
    num_snps_lst = [1_000, 10_000] if MOCK else [500, 1000, 2000, 5000, 10_000, 20_000, 50_000, 75_000,
                                                          100_000, 200_000, 500_000, 1_000_000]
    repetitions = 2 if MOCK else 10
    computers = {}
    colors = {}
    for max_alleles in max_alleles_num:
        computers[f'fast {max_alleles}'] = Fast(arguments, output=os.path.join(arguments.output, f'fast_{max_alleles}'))

    computers['fast'] = Fast(arguments, output=os.path.join(arguments.output, 'fast'))
    computers['naive'] = Naive(arguments, output=os.path.join(arguments.output, 'naive'))
    for idx, name in enumerate(computers.keys()):
        colors[name] = list(mcolors.TABLEAU_COLORS.keys())[idx]
    if not os.path.exists(os.path.join(arguments.output, "mean_times.json")):
        for num_snps in num_snps_lst:   # tqdm(num_snps_lst, desc='Number of SNPs'):
            file_path = os.path.join(arguments.output, get_file_name(num_indv, num_snps))
            for rep in range(repetitions):   # tqdm(range(repetitions), desc='repetition', leave=False):
                write_random_vcf(num_indv, num_snps, file_path, num_alleles=2)
                computers['naive'].vcf_to_metric('VCF', file_path, rep)
                computers['fast'].vcf_to_metric('VCF', file_path, rep)
                # Todo - compare output files
                os.remove(file_path)
                for num_alleles in max_alleles_num:
                    write_random_vcf(num_indv, num_snps, file_path, num_alleles=num_alleles)
                    computers[f'fast {num_alleles}'].vcf_to_metric('VCF', file_path, rep)
                    os.remove(file_path)
            for name, alg in computers.items():
                avg_times[name].append(np.mean(alg.times))
                std_times[name].append(np.std(alg.times))
                alg.times = []
        with open(os.path.join(arguments.output, "mean_times.json"), "w") as f:
            json.dump(avg_times, f)
        with open(os.path.join(arguments.output, "std_times.json"), "w") as f:
            json.dump(std_times, f)
    else:
        with open(os.path.join(arguments.output, "mean_times.json"), "r") as f:
            avg_times = json.load(f)
        with open(os.path.join(arguments.output, "std_times.json"), "r") as f:
            std_times = json.load(f)
    for name in computers.keys():
        avg_times[name] = np.array(avg_times[name])
        std_times[name] = np.array(std_times[name])
        plt.plot(num_snps_lst, avg_times[name], label=name, color=colors[name])
        plt.fill_between(num_snps_lst, avg_times[name] - std_times[name], avg_times[name] + std_times[name], alpha=0.3,
                         color=colors[name])
    plt.xlabel("Number of sites")
    plt.ylabel("Time (seconds)")
    plt.title(f"Similarity computation time with {num_indv} individuals")
    plt.legend()
    plt.savefig(os.path.join(arguments.output, "naive_running_time_comparison_linear.svg"))

def compare_different_num_of_individuals(arguments):
    avg_times = {'naive': [], 'fast': []}
    std_times = {'naive': [], 'fast': []}
    num_indv_lst = [10, 50, 150] if MOCK else [10, 20, 50, 100, 200, 500, 1000]
    num_snps = 2_000 if MOCK else 50_000
    repetitions = 5
    arguments.save_outputs = False
    computers = {}
    colors = {}

    computers['fast'] = Fast(arguments, output=os.path.join(arguments.output, 'fast'))
    computers['naive'] = Naive(arguments, output=os.path.join(arguments.output, 'naive'))
    for idx, name in enumerate(computers.keys()):
        colors[name] = list(mcolors.TABLEAU_COLORS.keys())[idx]
    if not os.path.exists(os.path.join(arguments.output, "mean_times.json")):
        for num_indv in tqdm(num_indv_lst, desc='Num individuals: '):
            file_path = os.path.join(arguments.output, get_file_name(num_indv, num_snps))
            for rep in range(repetitions):
                write_random_vcf(num_indv, num_snps, file_path, num_alleles=2)
                # computers['naive'].vcf_to_metric('VCF', file_path, rep)
                computers['fast'].vcf_to_metric('VCF', file_path, rep)
                # Todo - compare output files
                os.remove(file_path)
            for name, alg in computers.items():
                avg_times[name].append(np.mean(alg.times))
                std_times[name].append(np.std(alg.times))
                alg.times = []
        with open(os.path.join(arguments.output, "mean_times.json"), "w") as f:
            json.dump(avg_times, f)
        with open(os.path.join(arguments.output, "std_times.json"), "w") as f:
            json.dump(std_times, f)
    else:
        with open(os.path.join(arguments.output, "mean_times.json"), "r") as f:
            avg_times = json.load(f)
        with open(os.path.join(arguments.output, "std_times.json"), "r") as f:
            std_times = json.load(f)
    for name in computers.keys():
        avg_times[name] = np.array(avg_times[name])
        std_times[name] = np.array(std_times[name])
        plt.plot(num_indv_lst, avg_times[name], label=name, color=colors[name])
        plt.fill_between(num_indv_lst, avg_times[name] - std_times[name], avg_times[name] + std_times[name], alpha=0.3,
                         color=colors[name])

    plt.xlabel("Number of individuals")
    plt.ylabel("Time (seconds)")
    plt.title(f"Similarity computation time with {num_snps} SNPs")
    plt.legend()
    plt.savefig(os.path.join(arguments.output, "naive_running_time_comparison_linear.svg"))

if __name__ == '__main__':
    arguments = args_parser()
    os.makedirs(arguments.output, exist_ok=True)
    arguments.save_outputs = False
    if False:
        compare_different_num_of_snps(arguments)
    else:
        compare_different_num_of_individuals(arguments)
