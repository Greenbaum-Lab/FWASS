import argparse


def args_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input_path", help="path to input file", required=True)
    parser.add_argument("-o", "--output", dest="output_path", help="path to a new output directory", required=True)
    parser.add_argument("-t", "--type", dest="data_type", help="The format presentation of the input file", required=True)
    parser.add_argument("--max_alleles", dest="max_alleles", help="The maximum number of alleles for a single site")
    parser.add_argument("--max_single_array", dest="max_single_array", default=10 ** 8,
                        help="The maximum number of cells in a single martix")
    return parser.parse_args()