#!/bin/bash

# Load any necessary modules or set up environment variables here
bashrc=/sdf/home/a/alspellm/.bashrc

#config slurm args
user=alspellm
max_jobs_allowed=20
partition=hps
time_hrs=1
mem=4000

# Set the working directory
working_dir=$(readlink -f /sdf/group/hps/users/alspellm/run/hadd_jobs/working_dir)
logs_dir=${working_dir}/logs
output_dir=$(reradlink -f /sdf/group/hps/users/alspellm/run/hadd_jobs/working_dir/output)
mkdir -p "${logs_dir}"
mkdir -p "${output_dir}"

# Define the input file containing the list of files
list_of_job_files=$(readlink -f blpass4_ana_files.txt)

mkdir -p "$working_dir"
cd $working_dir

# Read the list of files into an array
mapfile -t jobs < "$list_of_job_files"

# Number of jobs to be run
num_jobs=${#jobs[@]}

jobn=1
njobs_running=0
for job in "${jobs[@]}"; do
    job_name=$(basename "$job")
    job_name_woext="${jobname%.*}"

    #Count number of your jobs running on partition
    njobs_running=$(squeue | grep "$partition" | grep "$user" | wc -l)
    echo " '$count' jobs running "
    room=$((max_jobs_allowed-count))

    if [ $njobs_running -lt $room ] ; then

        #Slurm job submission command
        sbatch <<EOT
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=${time_hrs}:00:00  # Adjust the time limit as needed
#SBATCH --mem=${mem}M
#SBATCH --job-name=${job_name_woext}
#SBATCH --output=${logs_dir}/${job_name}.out
#SBATCH --error=${logs_dir}/${job_name}.err
#SBATCH --partition=${partition}

echo "Loading environment (.bashrc)"
source ${bashrc}
echo "Running job ${job_name}"

sbatch ${job}
EOT

    echo "submitted job number ${jobn}"
    jobn=$((jobn + 1))
    else
        echo "No room for new jobs"
        while [ $room -lt 1 ]; do
            sleep 10
            count=$(squeue | grep "$partition" | grep "$uuser" | wc -l)
            room=$((max_jobs_allowed-count))
        done
    fi

done

