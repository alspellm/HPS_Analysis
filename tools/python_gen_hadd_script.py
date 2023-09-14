#!/usr/bin/python3
import argparse
import math

# Define the SLURM template
SLURM_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=hadd_job_array
#SBATCH --output=hadd_job_array_%A_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000M
#SBATCH --partition=hps
#SBATCH --array=0-{num_jobs}

echo "Job started at $(date)"

# Load any necessary modules or set up environment variables here
source /sdf/home/a/alspellm/.bashrc

# Define the input file containing the list of files
input_file=$(readlink -f {input_textfile})

# Set the working directory
working_dir=$(readlink -f {working_dir})
mkdir -p "$working_dir"
cd $working_dir

# Read the list of files into an array
mapfile -t input_files < "$input_file"

# Define the batch size
batch_size={batch_size}

# Calculate the number of batches
num_files=${{#input_files[@]}}
num_jobs=$(( (num_files + batch_size - 1) / batch_size ))

# Calculate the start and end indices for this job
start=$((SLURM_ARRAY_TASK_ID * batch_size))
end=$((start + batch_size - 1))

# Ensure that end does not exceed the number of files
if ((end >= num_files)); then
    end=$((num_files - 1))
fi

# Extract the files for this batch
batch_files=("${{input_files[@]:start:end-start+1}}")

# Run the hadd command
hadd "{output_name}_$SLURM_ARRAY_TASK_ID.root" "${{batch_files[@]}}"
"""

def create_slurm_script(input_textfile, output_name, working_dir, output_script, batch_size):
    num_jobs = math.ceil(len(open(input_textfile).readlines()) / batch_size)
    slurm_script = SLURM_TEMPLATE.format(num_jobs=num_jobs, batch_size=batch_size, input_textfile=input_textfile,
        output_name=output_name, working_dir=working_dir)
    with open("%s.sh"%(output_script), 'w') as f:
        f.write(slurm_script)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate SLURM submission script for hadd')
    parser.add_argument('-i', '--input-textfile', required=True, type=str, help='list of input root files to hadd')
    parser.add_argument('-o', '--output-name', required=True, type=str, help='hadd output prefix', default='hadd_sample')
    parser.add_argument('-d', '--working-dir', required=True, type=str, help='working directory',
            default='/scratch/alspellm/hadd')
    parser.add_argument('-s', '--output-script', required=True, help='Output SLURM script filename', 
            default='hadd_script.sh')
    parser.add_argument('-n', '--batch-size', required=True, type=int, help='N files to hadd together', default='10')
    
    args = parser.parse_args()
    
    create_slurm_script(args.input_textfile, args.output_name, args.working_dir, args.output_script,
            args.batch_size)

