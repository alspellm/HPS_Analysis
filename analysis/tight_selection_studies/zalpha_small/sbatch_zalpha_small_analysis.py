#!/bin/bash
#SBATCH --job-name=wide
#SBATCH --output=array_job.%A.%a.out
#SBATCH --error=array_job.%A.%a.err
#SBATCH --array=1-23
#SBATCH --partition=hps  # Specify the partition (queue)
#SBATCH --time=8:00:00            # Maximum runtime (1 hour)

# Check if the number of command-line arguments is correct
if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch script.sh <path_to_python_script> <run_directory>"
    exit 1
fi

# Assign the command-line arguments to variables
python_script="$1"
python_script=$(readlink -f "$python_script")

run_directory="$2"
run_directory=$(readlink -f "$run_directory")

source /sdf/home/a/alspellm/.bashrc
source /sdf/group/hps/users/alspellm/src/root_src/root_install/bin/thisroot.sh

# Set the working directory
cd "$run_directory"
echo "Run directory: $run_directory"
echo "Running script $python_script"

# Run your Python script
mass=$((35 + $SLURM_ARRAY_TASK_ID * 5))
python3 -u ${python_script} --mass $mass --logeps2 -6.0 > run_log_mass_${mass}.txt

