#!/bin/bash
#

#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process.
module load python
echo "Start running jobs $comp to $outputs"
mkdir $output
mkdir $output/$comp/
srun python research/compare_runner.py -o $output/$comp/asd --method asd --mock --comparison_name $comp
srun python research/compare_runner.py -o $output/$comp/similarity --method similarity --mock --comparison_name $comp
srun python research/compare_runner.py -o $output/$comp/weighted_similarity --method similarity -w 1 --mock --comparison_name $comp
