import os
import subprocess

from research.compare_runner import experiment_argument_parser

if __name__ == '__main__':
    arguments = experiment_argument_parser(check_method=False)
    submit_helper_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'slurm.sh')
    compare_runner_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'compare_runner.py')
    os.makedirs(os.path.join(arguments.output, 'logs'), exist_ok=True)
    for comp in ['num_of_individuals', 'num_of_snps']:
        os.makedirs(os.path.join(arguments.output, comp), exist_ok=True)
        for job_name, method in {'asd': 'asd', 'similarity': 'similarity', "weighted_similarity": "similarity -w 1"}.items():
            job_err = os.path.join(arguments.output, 'logs', f'{comp}_{job_name}_err')
            job_out = os.path.join(arguments.output, 'logs', f'{comp}_{job_name}_out')
            sbatch_settings = f'sbatch --time=72:00:00 --mem=8G --error="{job_err}" --output="{job_out}" --job-name="{job_name}"'
            command = f"python {compare_runner_path} --comparison_name {comp} -o " \
                      f"{os.path.join(arguments.output, comp, job_name)} --method {method}"
            if arguments.mock:
                command += " --mock"
            cmd2run = f"{sbatch_settings} {submit_helper_path} {command}"
            subprocess.run([submit_helper_path, cmd2run])
