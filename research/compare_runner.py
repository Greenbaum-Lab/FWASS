import argparse
import json
import os
import time
from copy import copy
import pandas as pd
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from research.fast_solution import Fast
from research.naive_solution import Naive
from utils import args_parser, write_random_vcf, add_method_args, handle_args


class Comparator:
    def __init__(self, arguments):
        self.options = arguments
        self.mock = arguments.mock

    def load_times(self):
        with open(os.path.join(self.options.output, "mean_times.json"), "r") as f:
            avg_times = json.load(f)
        with open(os.path.join(self.options.output, "std_times.json"), "r") as f:
            std_times = json.load(f)
        return avg_times, std_times

    @staticmethod
    def verify_same_output(dir1, dir2):
        dir1_files = sorted(os.listdir(dir1))
        dir2_files = sorted(os.listdir(dir2))
        assert dir1_files == dir2_files, "Different files in naive directory and fast directory"
        for file in dir1_files:
            d1_df = pd.read_csv(os.path.join(dir1, file))
            d2_df = pd.read_csv(os.path.join(dir2, file))
            assert all(d1_df.columns == d2_df.columns)
            np1 = d1_df.to_numpy()
            np2 = d2_df.to_numpy()
            assert np.allclose(np1[:, 1:].astype(float), np2[:, 1:].astype(float)), f"There are differences in file {os.path.join(dir1, file)}"
            os.remove(os.path.join(dir1, file))
            os.remove(os.path.join(dir2, file))

    def run_comparison(self):
        raise NotImplementedError


class ComparatorDiffNumIndividuals(Comparator):
    def run_comparison(self):
        arguments = self.options
        avg_times = {'naive': [], 'fast': []}
        std_times = {'naive': [], 'fast': []}
        num_indv_lst = [10, 50, 100] if self.mock else [10, 20, 50, 100, 200, 500, 1000]
        num_snps = 2_000 if self.mock else 50_000
        repetitions = 5
        computers = {}
        colors = {}

        computers['fast'] = Fast(arguments, output=os.path.join(arguments.output, 'fast'))
        computers['naive'] = Naive(arguments, output=os.path.join(arguments.output, 'naive'))
        for idx, name in enumerate(computers.keys()):
            colors[name] = list(mcolors.TABLEAU_COLORS.keys())[idx]
        if not os.path.exists(os.path.join(arguments.output, "mean_times.json")):
            for num_indv in tqdm(num_indv_lst, desc='Num individuals: '):
                file_path = os.path.join(arguments.output, get_file_name(num_indv, num_snps))
                arguments.input = file_path
                for rep in range(repetitions):
                    write_random_vcf(num_indv, num_snps, file_path, num_alleles=2)
                    computers['naive'].vcf_to_metric('VCF', file_path, rep)
                    computers['fast'].vcf_to_metric('VCF', file_path, rep)
                    self.verify_same_output(computers['fast'].output_dir, computers['naive'].output_dir)
                    os.remove(file_path)
                for name, alg in computers.items():
                    avg_times[name].append(np.mean(alg.times))
                    std_times[name].append(np.std(alg.times))
                    alg.times = []
            with open(os.path.join(arguments.output, "mean_times.json"), "w") as f:
                json.dump(avg_times, f)
            with open(os.path.join(arguments.output, "std_times.json"), "w") as f:
                json.dump(std_times, f)

        avg_times, std_times = self.load_times()
        for name in computers.keys():
            avg_times[name] = np.array(avg_times[name])
            std_times[name] = np.array(std_times[name])
            plt.plot(num_indv_lst, avg_times[name], label=name, color=colors[name])
            plt.fill_between(num_indv_lst, avg_times[name] - std_times[name], avg_times[name] + std_times[name],
                             alpha=0.3,
                             color=colors[name])

        plt.xlabel("Number of individuals")
        plt.ylabel("Time (seconds)")
        plt.title(f"Similarity computation time with {num_snps} SNPs")
        plt.legend()
        plt.savefig(os.path.join(arguments.output, "running_time_comparison_diff_individuals.pdf"))


class ComparatorDiffNumSNPs(Comparator):

    def run_comparison(self):
        arguments = self.options
        max_alleles_num = [4, 10, 50]
        avg_times = {'naive': [], 'fast': []}
        avg_times.update({f'fast {x}': [] for x in max_alleles_num})
        std_times = {'naive': [], 'fast': []}
        std_times.update({f'fast {x}': [] for x in max_alleles_num})
        num_indv = 100
        num_snps_lst = [1_000, 3_000, 10_000] if self.mock else [1000, 5000, 10_000, 50_000, 100_000, 200_000, 500_000,
                                                                 1_000_000]
        repetitions = 2 if self.mock else 10
        computers = {}
        colors = {}
        for max_alleles in max_alleles_num:
            computers[f'fast {max_alleles}'] = Fast(arguments,
                                                    output=os.path.join(arguments.output, f'fast_{max_alleles}'))

        computers['fast'] = Fast(arguments, output=os.path.join(arguments.output, 'fast'))
        computers['naive'] = Naive(arguments, output=os.path.join(arguments.output, 'naive'))
        for idx, name in enumerate(computers.keys()):
            colors[name] = list(mcolors.TABLEAU_COLORS.keys())[idx]
        if not os.path.exists(os.path.join(arguments.output, "mean_times.json")):
            for num_snps in num_snps_lst:  # tqdm(num_snps_lst, desc='Number of SNPs'):
                file_path = os.path.join(arguments.output, get_file_name(num_indv, num_snps))
                arguments.input = file_path
                for rep in range(repetitions):  # tqdm(range(repetitions), desc='repetition', leave=False):
                    write_random_vcf(num_indv, num_snps, file_path, num_alleles=2)
                    computers['naive'].vcf_to_metric('VCF', file_path, rep)
                    computers['fast'].vcf_to_metric('VCF', file_path, rep)
                    self.verify_same_output(computers['fast'].output_dir, computers['naive'].output_dir)
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
        avg_times, std_times = self.load_times()
        for name in computers.keys():
            avg_times[name] = np.array(avg_times[name])
            std_times[name] = np.array(std_times[name])
            plt.plot(num_snps_lst, avg_times[name], label=name, color=colors[name])
            plt.fill_between(num_snps_lst, avg_times[name] - std_times[name], avg_times[name] + std_times[name],
                             alpha=0.3,
                             color=colors[name])
        plt.xlabel("Number of sites")
        plt.ylabel("Time (seconds)")
        plt.title(f"Similarity computation time with {num_indv} individuals")
        plt.legend()
        plt.savefig(os.path.join(arguments.output, "running_time_comparison_diff_snps.pdf"))


class CompareDiffMemorySize(Comparator):
    def run_comparison(self):
        arguments = self.options
        max_mems = [1, 4, 16] if self.mock else [0.1, 0.5, 1, 2, 5, 10]
        num_indv_list = [1000, 2000, 4000] if self.mock else [10, 30, 50, 100, 250, 1_000, 2_500]
        avg_times = {f'indv_{j}': {f'max_{i}': [] for i in max_mems} for j in num_indv_list}
        std_times = {f'indv_{j}': {f'max_{i}': [] for i in max_mems} for j in num_indv_list}
        file_path = os.path.join(arguments.output, get_file_name(num_indv_list[0], num_indv_list[0]))
        arguments.input = file_path
        num_snps = 50_000 if self.mock else 100_000
        repetitions = 2 if self.mock else 5
        computers = {}
        for max_mem in max_mems:
            temp_args = copy(arguments)
            temp_args.max_mb = max_mem
            computers[max_mem] = Fast(temp_args, output=os.path.join(arguments.output, f'max_{max_mem}'))
        for num_indv in num_indv_list:
            for rep in range(repetitions):
                write_random_vcf(num_indv, num_snps, file_path, num_alleles=2)
                for comp in computers.values():
                    comp.vcf_to_metric('VCF', file_path, rep)
            for name, alg in computers.items():
                avg_times[f'indv_{num_indv}'][f'max_{name}'] = np.mean(alg.times)
                std_times[f'indv_{num_indv}'][f'max_{name}'] = np.std(alg.times)
        print(f"avg_times: {avg_times}")
        print(f"std_times: {std_times}")







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


def experiment_argument_parser(check_method=True):
        parser = argparse.ArgumentParser()
        parser.add_argument("--comparison_name", help="Name of the required comparison to run.")
        parser.add_argument("-o", "--output", required=True, help="Name of output directory")
        parser.add_argument("--mock", default=False, action='store_true')
        parser = add_method_args(parser)
        return handle_args(parser, check_method=check_method)


if __name__ == '__main__':
    print('hiiiiii')
    np.random.seed(0)
    arguments = experiment_argument_parser()
    if arguments.comparison_name == 'num_of_individuals':
        print("IN ComparatorDiffNumIndividuals")
        ComparatorDiffNumIndividuals(arguments).run_comparison()
    elif arguments.comparison_name == 'num_of_snps':
        print("IN ComparatorDiffNumSNPs")
        ComparatorDiffNumSNPs(arguments).run_comparison()
    elif arguments.comparison_name == 'memory_size':
        print("IN CompareDiffMemorySize")
        CompareDiffMemorySize(arguments).run_comparison()
    else:
        print("OOOOOPS")
        raise NameError(f"Comparison name {arguments.comparison_name} is not valid")
