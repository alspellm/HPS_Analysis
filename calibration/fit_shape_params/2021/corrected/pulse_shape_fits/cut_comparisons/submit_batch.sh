#!/usr/bin/scl enable devtoolset-8 -- /bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=25000M
#SBATCH --array=1-10
#SBATCH --partition=shared
#SBATCH --job-name=pulses

export FIRST_ID=0
export JOB_ID=$(($SLURM_ARRAY_TASK_ID+FIRST_ID-1))
echo "JOB ID IS"
echo $JOB_ID
source /sdf/home/a/alspellm/.bashrc

python3 /sdf/group/hps/users/alspellm/projects/fit_shape_params/2021/corrected/pulse_shape_fits/cut_comparisons/compareReco.py -n ${JOB_ID} -o hps_14191_data_events_calgroup7_cut_comparison_${JOB_ID}.root
