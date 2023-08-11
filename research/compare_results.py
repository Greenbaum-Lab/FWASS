import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from research.compare_runner import file_name_to_indv
from utils import args_parser


class Compare:
    def __init__(self, data_dir, ref_res_dir, other_res_dir=None):
        self.data_dir = data_dir
        self.ref_res_dir = ref_res_dir
        self.other_res_dir = other_res_dir

    def compare_repetitions_equal(self, dfs_dir):
        for file in os.listdir(self.data_dir):
            base_name = file[:-4]
            df_files_to_compare = [f for f in os.listdir(self.ref_res_dir) if f.startswith(f'{base_name}_')
                                   and f.endswith('csv')]
            df = pd.read_csv(os.path.join(dfs_dir, df_files_to_compare[0]))
            for other_df in df_files_to_compare[1:]:
                df2compare = pd.read_csv(os.path.join(dfs_dir, other_df))
                assert df.equals(df2compare), f"Different results for the same data with same algorithm! on {base_name}"

    def compare_between_algorithms(self):
        for file in os.listdir(self.data_dir):
            base_name = file[:-4]
            ref_file = os.path.join(self.ref_res_dir, base_name) + '_0.csv'
            other_file = os.path.join(self.other_res_dir, base_name) + '_0.csv'
            ref_df = pd.read_csv(ref_file)
            other_df = pd.read_csv(other_file)
            assert ref_df.equals(other_df), f"Different results between different algorithms! on {base_name}"

    def plot_times(self):
        all_times = []
        for file in os.listdir(self.data_dir):
            base_name = file[:-4]
            all_times.append([file_name_to_indv(base_name)])
            with open(os.path.join(self.ref_res_dir, f'{base_name}_time.txt'), "r") as f:
                times = f.readlines()
            times = np.array([float(e[:-1]) for e in times])
            all_times[-1].append(np.mean(times))
            all_times[-1].append(np.std(times))
            with open(os.path.join(self.other_res_dir, f'{base_name}_time.txt'), "r") as f:
                times = f.readlines()
            times = np.array([float(e[:-1]) for e in times])
            all_times[-1].append(np.mean(times))
            all_times[-1].append(np.std(times))
        all_times = np.array(sorted(all_times))
        num_snps = all_times[:, 0]
        ref_mean = all_times[:, 1]
        ref_std = all_times[:, 2]
        dev_mean = all_times[:, 3]
        dev_std = all_times[:, 4]
        plt.plot(num_snps, ref_mean, color='tab:blue', label='naive')
        plt.fill_between(num_snps, ref_mean - ref_std, ref_mean + ref_std, alpha=0.3, color='tab:blue')
        plt.plot(num_snps, dev_mean, color='tab:red', label='matrix form')
        plt.fill_between(num_snps, dev_mean - dev_std, dev_mean + dev_std, alpha=0.3, color='tab:red')
        plt.xlabel("Number of sites")
        plt.ylabel("Time (seconds)")
        plt.title("Similarity computation time with 100 individuals")
        plt.legend()
        plt.xscale('log')
        plt.xscale('log')
        plt.show()



if __name__ == '__main__':
    args = args_parser()
    compare = Compare(args.input, args.output, other_res_dir='C:\\Users\\shaharma\\LAB\\fast_sim_mock/')
    compare.compare_repetitions_equal(compare.ref_res_dir)
    compare.compare_repetitions_equal(compare.other_res_dir)
    compare.compare_between_algorithms()
    compare.plot_times()
    print("Done")