import doctest
import gzip
import os
import pandas as pd

from tqdm import tqdm


class InputReader:
    def __init__(self, file_path):
        self.file_path = file_path

    def get_object(self):
        if self.file_path.endswith(".csv"):
            return CSVInputReader(self.file_path)
        if self.file_path.endswith(".xlsx"):
            return XLSXInputReader(self.file_path)
        if self.file_path.endswith("vcf.gz"):
            return VcfGZInputReader(self.file_path)
        if self.file_path.endswith(".vcf"):
            return VcfInputReader(self.file_path)
        raise IOError("File format is not supported")

    def file2small_matrices(self, options, mid_outputs_path, verbose):
        raise NotImplementedError


class VcfInputReader(InputReader):
    def __init__(self, file_path):
        self.file_path = file_path
        self.reading_method = open

    def file2small_matrices(self, options, mid_outputs_path, verbose):
        """
        Parse the VCF file, and save small matrices to compute on different threads.
        :param input_format: "VCF-GZ" if it compressed as gzip, "VCF" if it's a text file in a VCF format.
        :param options: running arguments
        :param mid_outputs_path: path of directory to throw mid-results for computations.
        :return: None, The function parse the file and save the small matrices in a mid_res directory.
        """
        max_num_of_cells = options.max_mb * 10 ** 6
        if verbose:
            pbar = tqdm(desc="Run over sites")
        with self.reading_method(self.file_path, "rb") as f:
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
            if verbose:
                pbar.update(1)
            last_line = f.readline().decode()
            while last_line:
                if sites_counter % num_sites_to_read == 0:
                    with open(os.path.join(mid_outputs_path, f'mat_{matrices_counter}.tmp'), "w") as g:
                        g.write(current_matrix)
                    current_matrix = ""
                    matrices_counter += 1
                line = last_line.split()
                assert (len(line[9:]) == len(individuals) and line[format_dict["FORMAT"]].startswith("GT"))
                indv_gt = ['FAIL' if '.' in e.split(':')[0] else e.split(':')[0].replace('|', '/') for e in line[9:]]
                current_matrix += "\t".join(indv_gt) + '\n'
                is_matrix_empty = False
                sites_counter += 1
                last_line = f.readline().decode()
                if verbose:
                    pbar.update(1)
                if options.max_sites and sites_counter > options.max_sites:
                    break

        if not is_matrix_empty:
            with open(os.path.join(mid_outputs_path, f'mat_{matrices_counter}.tmp'), "w") as g:
                g.write(current_matrix)

        with open(os.path.join(mid_outputs_path, "done.txt"), "w") as f:
            f.write("\t".join(individuals))
        if verbose:
            if matrices_counter == 0:
                print(f"Done reading VCF file!\n{sites_counter - 1} sites in total, divided over a single matrix.")
            else:
                print(
                    f"Done reading VCF file!\n{sites_counter - 1} sites in total, divided over {matrices_counter + 1}"
                    f" matrices. Max of {num_sites_to_read} in each matrix.\n Currently computing the last few"
                    f" similarity matrices.")

class VcfGZInputReader(VcfInputReader):
    def __init__(self, file_path):
        self.file_path = file_path
        self.reading_method = gzip.open


class CSVInputReader(InputReader):
    def __init__(self, file_path):
        self.file_type = 'CSV'
        self.file_path = file_path
        self.reading_method = pd.read_csv

    def file2small_matrices(self, options, mid_outputs_path, verbose):
        """
        Parse the csv file, and save a single matrix.
        :param options: running arguments
        :param mid_outputs_path: path of directory to throw mid-results for computations.
        :return: None, The function parse the file and save the small matrices in a mid_res directory.
        """
        df = self.reading_method(self.file_path)
        df = self.filter_out_bad_sites(df, verbose)
        trans_df = df.T
        individuals = list(trans_df.iloc(0)[0])
        individuals = [str(e) for e in individuals]
        df = trans_df.iloc[1:]
        with open(os.path.join(mid_outputs_path, "done.txt"), "w") as f:
            f.write("\t".join(individuals))
        mini_mat_path = os.path.join(mid_outputs_path, f'mat_{0}.tmp')
        df.to_csv(mini_mat_path, sep='\t', index=False,  header=False)
        if verbose:
            print(f"Done reading {self.file_type} file!\n{df.shape[0]} sites in total, divided over a single matrix.")

    def filter_out_bad_sites(self, df,  verbose=True):
        loci_to_remove = []
        for locus_name, col in df.items():
            if locus_name == 'ID':
                continue
            non_fail = col[col != 'FAIL']
            if non_fail.size == 0:
                loci_to_remove.append(locus_name)
                continue
            first_individual = [e.strip() for e in non_fail.iloc[0].split('/')]
            if first_individual[0] == first_individual[1]:
                same_as_first = non_fail[non_fail == non_fail.iloc[0]]
                if same_as_first.size == non_fail.size:
                    loci_to_remove.append(locus_name)
        if loci_to_remove:
            df = df.drop(loci_to_remove, axis=1)
            num_of_sites = len(df.columns) - int('ID' in df.columns)
            if verbose:
                print(f"Removed {len(loci_to_remove)} invalid loci. Compute metric based on {num_of_sites} sites.")
        return df


class XLSXInputReader(CSVInputReader):

    def __init__(self, file_path):
        self.file_type = 'XLSX'
        self.file_path = file_path
        self.reading_method = pd.read_excel