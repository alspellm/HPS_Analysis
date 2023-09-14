#!/usr/bin/scl enable devtoolset-8 -- /bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=15000M
#SBATCH --array=1-10
#SBATCH --partition=shared
#SBATCH --job-name=pulses

export FIRST_ID=0
export JOB_ID=$(($SLURM_ARRAY_TASK_ID+FIRST_ID-1))
echo "JOB ID IS"
echo $JOB_ID
source /sdf/home/a/alspellm/.bashrc

python3 /sdf/group/hps/users/alspellm/projects/fit_shape_params/2019/corrected/pulse_shape_fits/comparisons/compare_to_25ns_delay/compareReco.py -n ${JOB_ID} -o hps_10030_data_events_cut_comparison_${JOB_ID}.root
