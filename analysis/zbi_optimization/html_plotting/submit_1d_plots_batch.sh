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
python_script="/sdf/group/hps/users/alspellm/projects/THESIS/ana/zbi_html/make_1D_sig_bkg_plots.py"

# List of input files
input_files=(
    "/sdf/group/hps/users/alspellm/projects/THESIS/ana/isolation_cut_dev/20230724_100MeV/cut_tests/zalpha_gt_lt/signal_55_bkg20_sig5_zalpha_slope_0.05_gt_lt.root"
    "/sdf/group/hps/users/alspellm/projects/THESIS/ana/isolation_cut_dev/20230724_100MeV/cut_tests/zalpha_gt_lt/signal_55_bkg20_sig5_zalpha_slope_0.07_gt_lt.root"
    "/sdf/group/hps/users/alspellm/projects/THESIS/ana/isolation_cut_dev/20230724_100MeV/cut_tests/zalpha_gt_lt/signal_55_bkg20_sig5_zalpha_slope_0.09_gt_lt.root"
    "/sdf/group/hps/users/alspellm/projects/THESIS/ana/isolation_cut_dev/20230724_100MeV/cut_tests/zalpha_gt_lt/signal_55_bkg20_sig5_zalpha_slope_0.11_gt_lt.root"
    "/sdf/group/hps/users/alspellm/projects/THESIS/ana/isolation_cut_dev/20230724_100MeV/cut_tests/zalpha_gt_lt/signal_55_bkg20_sig5_zalpha_slope_0.13_gt_lt.root"
    "/sdf/group/hps/users/alspellm/projects/THESIS/ana/isolation_cut_dev/20230724_100MeV/cut_tests/zalpha_gt_lt/signal_55_bkg20_sig5_zalpha_slope_0.15_gt_lt.root"
    "/sdf/group/hps/users/alspellm/projects/THESIS/ana/isolation_cut_dev/20230724_100MeV/cut_tests/zalpha_gt_lt/signal_55_bkg20_sig5_zalpha_slope_0.17_gt_lt.root"
    # Add more input files here
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
#SBATCH --job-name=${Rundir}/${job_name}
#SBATCH --output=${Rundir}/${job_name}.out
#SBATCH --error=${Rundir}/${job_name}.err
#SBATCH --partition=hps

echo "Running on input file ${input_file}"
source /sdf/home/a/alspellm/.bashrc

python3 ${python_script} --infileName ${input_file} --plotsDir ${Rundir}/plots_1d/${filename_without_extension}

EOT
done
