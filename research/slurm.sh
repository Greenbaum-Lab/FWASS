#!/bin/bash
#

#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process.

module load python
source /sci/labs/gilig/shahar.mazie/icore-data/snpnmi_venv/bin/activate.csh

echo "Start running jobs $0"
mkdir /sci/labs/gilig/shahar.mazie/icore-data/FWASS_research/$0/
srun python /sci/labs/gilig/shahar.mazie/icore-data/code/FWASS/research/compare_runner.py -o /sci/labs/gilig/shahar.mazie/icore-data/FWASS_research/$0/asd --method asd --mock --comparison_name $0
srun python /sci/labs/gilig/shahar.mazie/icore-data/code/FWASS/research/compare_runner.py -o /sci/labs/gilig/shahar.mazie/icore-data/FWASS_research/$0/similarity --method similarity --mock --comparison_name $0
srun python /sci/labs/gilig/shahar.mazie/icore-data/code/FWASS/research/compare_runner.py -o /sci/labs/gilig/shahar.mazie/icore-data/FWASS_research/$0/weighted_similarity --method similarity -w 1 --mock --comparison_name $0
