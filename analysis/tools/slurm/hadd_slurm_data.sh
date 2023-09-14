#!/bin/bash
#SBATCH --job-name=hadd_job_array
#SBATCH --output=hadd_job_array_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000M
#SBATCH --partition=hps
#SBATCH --array=0-19

echo "Job started at $(date)"

# Load any necessary modules or set up environment variables here
source /sdf/home/a/alspellm/.bashrc

# Define the input file containing the list of files
input_file=$(readlink -f blpass4_ana_files.txt)

# Set the working directory
working_dir=$(readlink -f /sdf/group/hps/users/alspellm/run/hadd_jobs/working_dir)
mkdir -p "$working_dir"
cd $working_dir

# Read the list of files into an array
mapfile -t input_files < "$input_file"

# Define the batch size
batch_size=100

# Calculate the number of batches
num_files=${#input_files[@]}
num_jobs=$(( (num_files + batch_size - 1) / batch_size ))

# Calculate the start and end indices for this job
start=$((SLURM_ARRAY_TASK_ID * batch_size))
end=$((start + batch_size - 1))

# Ensure that end does not exceed the number of files
if ((end >= num_files)); then
    end=$((num_files - 1))
fi

# Extract the files for this batch
batch_files=("${input_files[@]:start:end-start+1}")

# Run the hadd command
hadd "hadd_blpass4_$SLURM_ARRAY_TASK_ID.root" "${batch_files[@]}"
