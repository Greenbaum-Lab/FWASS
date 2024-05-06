import os
import pandas as pd


class OutputWriter:
    def __init__(self, options):
        self.options = options
        self.output_dir = options.output

    def save_outputs(self, metric_mat, count_mat, individual_names=None):
        raise NotImplementedError

    def get_object(self):
        format = self.options.output_format
        if format == 'csv':
            return CSVOutputWriter(self.options)
        elif format == 'ns' or format == 'netstruct':
            return NetStructOutputWriter(self.options)
        else:
            raise NameError(f"format {format} is not valid.")

class CSVOutputWriter(OutputWriter):

    def save_outputs(self, metric_mat, count_mat, individual_names=None):
        """
        Save the numpy array as csv files with names of individuals.
        :param metric_mat: metric matrix
        :param count_mat: counts matrix
        :param individual_names: list of individual names
        :return: None
        """
        os.makedirs(self.output_dir, exist_ok=True)
        metric_out_path = os.path.join(self.output_dir, f"{self.options.method}" +
                                       '_weighted' * bool(self.options.weight) + '.csv')
        count_out_path = os.path.join(self.output_dir, "count.csv")
        counts_df = pd.DataFrame(count_mat, columns=individual_names, index=individual_names)
        metric_df = pd.DataFrame(metric_mat, columns=individual_names, index=individual_names)
        metric_df.to_csv(metric_out_path)
        counts_df.to_csv(count_out_path)


class NetStructOutputWriter(OutputWriter):

    def save_outputs(self, metric_mat, count_mat, individual_names=None):
        os.makedirs(self.output_dir, exist_ok=True)
        metric_out_path = os.path.join(self.output_dir, f"{self.options.method}" +
                                       '_weighted' * bool(self.options.weight) + '.txt')
        count_out_path = os.path.join(self.output_dir, "count.txt")
        with open(metric_out_path, "w") as f:
            f.write(self.numpy2net_struct_text(metric_mat))
        with open(count_out_path, 'w') as f:
            f.write(self.numpy2net_struct_text(count_mat))

    @staticmethod
    def numpy2net_struct_text(matrix):
        new_data = ""
        number_of_cells_to_skip = 1
        for line in matrix:
            for num in line[number_of_cells_to_skip:]:
                new_data += f'{num} '
            number_of_cells_to_skip += 1
            new_data = new_data[:-1] + '\n'
        return new_data[:-1]
