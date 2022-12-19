from utils import args_parser
import os
import pandas as pd
import numpy as np


def arguments_check(options):
    if not os.path.isfile(options.input_path):
        assert False, f"Input file {options.input_path} doesnt exists!"
    if os.path.exists(options.output_path):
        assert False, f"output path: {options.output_path} exists! Please use a new directory name. We don't wish to" \
                      f" run over your existing files"
    if options.data_type not in ['012', 'raw_data']:
        assert False, f"{options.data_type} is not a familiar data type. Please run --help for more info."


def read_df_file(arguments):
    i_path = arguments.input_path
    if i_path.endswith(".xlsx"):
        return pd.read_excel(i_path)
    if i_path.endswith(".csv"):
        return pd.read_csv(i_path)
    input_file_size = os.path.getsize(i_path)
    if input_file_size <= arguments.max_single_array * 8:
        return [i_path]
    os.makedirs(f"{arguments.output}/split_files")
    files_list = []
    with open(arguments.input_path, "r") as f:
        line = f.readline()
        num_of_individuals = len(line.split('\t'))
        max_num_of_sites = arguments.max_single_array // num_of_individuals
        file_num = -1
        while line:
            file_num += 1
            site_num = 1
            file = [line]
            while site_num < max_num_of_sites and line:
                line = f.readline()
                site_num += 1
                file.append(line)
            file_output_path = f"{arguments.output}/split_files/f_{file_num}.npy"
            np.save(file_output_path, np.array('\n'.join(file)))
            files_list.append(file_output_path)



if __name__ == '__main__':
    arguments = args_parser()
    arguments_check(arguments)
    os.makedirs(arguments.output_path)
    if arguments.data_type == '012':
        allele_numbers_matrix = read_df_file(arguments.input_path)