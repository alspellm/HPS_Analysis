#!/bin/bash

# Reset OPTIND to its initial value (1)
OPTIND=1
while getopts r: flag
do
    case "${flag}" in
        r) rundir=${OPTARG};;
    esac
done

echo "LOOK $rundir"
export Rundir=$(readlink -f "$rundir")
echo "$Rundir"

# Path to your Python script
python_script="/sdf/group/hps/users/alspellm/projects/THESIS/analysis/zbi_optimization/new_html_plotting/make_2D_sig_bkg_plots.py"

# List of input files
input_files=(
     "/sdf/group/hps/users/alspellm/projects/THESIS/analysis/tight_selection_studies/v0_projection/zbi_opt/output/output_exp_tail/test_deltaZ_and_projsig_stepsize_2pct.root"
    "/sdf/group/hps/users/alspellm/projects/THESIS/analysis/tight_selection_studies/v0_projection/zbi_opt/output/output_exp_tail/test_deltaZ_stepsize_2pct.root"
    "/sdf/group/hps/users/alspellm/projects/THESIS/analysis/tight_selection_studies/v0_projection/zbi_opt/output/output_exp_tail/test_projsig_stepsize_2pct.root"
)

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
#SBATCH --job-name=${job_name}
#SBATCH --output=${Rundir}/${job_name}.out
#SBATCH --error=${Rundir}/${job_name}.err
#SBATCH --partition=hps

echo "Running on input file ${input_file}"
source /sdf/home/a/alspellm/.bashrc

python3 ${python_script} --infileName ${input_file} --plotsDir ${Rundir}/plots_2d/${filename_without_extension}

EOT
done
