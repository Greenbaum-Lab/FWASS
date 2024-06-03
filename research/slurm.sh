#!/bin/bash
#

#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH --ntasks 3  # Request that n cpus
#SBATCH --cpus-per-task 1  # Request that n cpus
module load python
echo "Start running jobs $comp to $output"
mkdir -p $output/$comp/
srun python research/compare_runner.py --mock --comparison_name $comp -o $output/$comp/asd --method asd &
srun python research/compare_runner.py --mock --comparison_name $comp -o $output/$comp/similarity --method similarity &
srun python research/compare_runner.py --mock --comparison_name $comp -o $output/$comp/weighted_similarity --method similarity -w 1 &
wait
