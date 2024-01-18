import os
import gzip
from tqdm import tqdm
import numpy as np
import time
import pandas as pd

from compute_similarity import analyze_by_vcf
from research.naive_solution import Naive
from utils import args_parser, parse_input_file_format
from numba import jit

class Fast(Naive):
    def vcf_to_metric(self, input_format, file_name, iter_num):
        start_time = time.time()
        df = analyze_by_vcf(input_format, self.arguments, file_name, verbose=False)
        self.save_output(file_name, df, iter_num)
        self.times.append(time.time() - start_time)

# def main():
#     arguments = args_parser()
#     arguments.save_outputs = False
#     fast = Fast(arguments, input=arguments.input, output=arguments.output, repetitions=3)
#     for file in os.listdir(arguments.input):
#         file_path = os.path.join(arguments.input, file)
#         input_format = parse_input_file_format(file_path)
#         fast.compute_with_repetitions(input_format, file_path)
#     print("Done computation")


# if __name__ == '__main__':
#     main()
