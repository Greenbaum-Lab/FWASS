#!/bin/bash
#

#SBATCH -t 3-00:00 # Runtime in D-HH:MM
#SBATCH --cpus-per-task=1 # Request that ncpus be allocated per process.

module load python
source /sci/labs/gilig/shahar.mazie/icore-data/snpnmi_venv/bin/activate.csh

echo "Start running jobs $comp"
mkdir /sci/labs/gilig/shahar.mazie/icore-data/FWASS_research/$comp/
srun python /sci/labs/gilig/shahar.mazie/icore-data/code/FWASS/research/compare_runner.py -o /sci/labs/gilig/shahar.mazie/icore-data/FWASS_research/$comp/asd --method asd --mock --comparison_name $comp
srun python /sci/labs/gilig/shahar.mazie/icore-data/code/FWASS/research/compare_runner.py -o /sci/labs/gilig/shahar.mazie/icore-data/FWASS_research/$comp/similarity --method similarity --mock --comparison_name $comp
srun python /sci/labs/gilig/shahar.mazie/icore-data/code/FWASS/research/compare_runner.py -o /sci/labs/gilig/shahar.mazie/icore-data/FWASS_research/$comp/weighted_similarity --method similarity -w 1 --mock --comparison_name $comp
