import os
import time
from copy import copy

from metric_compute.runner import AbstractMetricCompute
from research.naive_solution import Naive


class Fast(Naive):
    def vcf_to_metric(self, input_format, file_name, iter_num):
        temp_args = copy(self.arguments)
        temp_args.output = self.output_dir
        start_time = time.time()
        single_runner = AbstractMetricCompute(temp_args).get_object()
        single_runner.analyze_file(verbose=False)
        self.times.append(time.time() - start_time)

