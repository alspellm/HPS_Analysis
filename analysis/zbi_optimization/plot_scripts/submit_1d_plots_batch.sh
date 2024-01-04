#!/bin/bash
# Initialize variables
OPTIND=1
input_directory=""
rundir=""
while getopts :r:i: flag
do
    case "${flag}" in
        r) rundir=${OPTARG};;
        i) input_directory=${OPTARG};;
    esac
done

# Check if the rundir and input_directory are provided
if [ -z "$rundir" ] || [ -z "$input_directory" ]; then
    echo "Both -r (rundir) and -i (input_directory) flags are required."
fi

# Check if the rundir and input_directory exist
if [ ! -d "$rundir" ]; then
    echo "rundir does not exist: $rundir"
fi

if [ ! -d "$input_directory" ]; then
    echo "input_directory does not exist: $input_directory"
fi

# Path to your Python script
python_script="/sdf/group/hps/users/alspellm/projects/THESIS/analysis/zbi_optimization/plot_scripts/make_1D_sig_bkg_plots.py"

#Make logs directory in run
rundir=$(readlink -f "$rundir")
input_directory=$(readlink -f "$input_directory")
mkdir -p ${rundir}/logs
mkdir -p ${rundir}/plots_1d

# Use find to generate a list of input files in the input_directory
input_files=()
while IFS= read -r -d '' file; do
    input_files+=("$file")
done < <(find "$input_directory" -name '*.root' -type f -print0)

# Check if any input files were found
if [ ${#input_files[@]} -eq 0 ]; then
    echo "No ROOT files found in $input_directory."
fi

# Print the list of input files
for file in "${input_files[@]}"; do
    echo "Input file: $file"
done

# Loop through input files and submit jobs
for input_file in "${input_files[@]}"; do
    # Extract filename without path and extension
    filename=$(basename "$input_file")
    filename_without_extension="${filename%.*}"
    job_name=${filename_without_extension}
    echo "JOB NAME is $job_name"

    # SLURM job submission command
    sbatch <<EOT
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=00:30:00  # Adjust the time limit as needed
#SBATCH --mem=10000M
#SBATCH --job-name=${rundir}/${job_name}
#SBATCH --output=${rundir}/logs/${job_name}.out
#SBATCH --error=${rundir}/logs/${job_name}.err
#SBATCH --partition=shared

echo "Running on input file ${input_file}"
source /sdf/home/a/alspellm/.bashrc

python3 ${python_script} --infileName ${input_file} --plotsDir ${rundir}/plots_1d/${filename_without_extension}

EOT
done
