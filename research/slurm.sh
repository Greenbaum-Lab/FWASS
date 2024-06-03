#!/bin/bash
#

#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH -n 3 # Request that n cpus
module load python
echo "Start running jobs $comp to $output"
mkdir $output
mkdir $output/$comp/
srun python research/compare_runner.py -o $output/$comp/asd --method asd --mock --comparison_name $comp &
srun python research/compare_runner.py -o $output/$comp/similarity --method similarity --mock --comparison_name $comp &
srun python research/compare_runner.py -o $output/$comp/weighted_similarity --method similarity -w 1 --mock --comparison_name $comp &
wait
