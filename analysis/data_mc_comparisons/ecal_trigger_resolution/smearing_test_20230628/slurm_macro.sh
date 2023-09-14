#!/bin/bash

#SBATCH --partition=hps
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000M
#SBATCH --time=1:00:00

#SBATCH --job-name=fsp
#SBATCH --output=./run/sh/output_%a.log
#SBATCH --error=./run/sh/output_%a.error

OPTIND=1
while getopts ":m:o:s:i:" opt;do
    case ${opt} in
        m) macro=${OPTARG}
            ;;
        o) output_dir=${OPTARG}
            ;;
        s) scratch_dir=${OPTARG}
            ;;
        i) initial=${OPTARG}
            ;;
    esac
done

#job_id=${SLURM_ARRAY_TASK_ID}
job_id=$(( SLURM_ARRAY_TASK_ID -1 + initial ))
echo "job id is: ${job_id}"

source /sdf/home/a/alspellm/.bashrc
echo "output dir: ${output_dir}"
echo "scratch dir: ${scratch_dir}"
output_dir=$(realpath ${output_dir})
macro=$(realpath ${macro})
scratch_dir=$(realpath ${scratch_dir})
echo "output dir: ${output_dir}"
echo "scratch dir: ${scratch_dir}"

mkdir -p ${scratch_dir}/${job_id}
cd ${scratch_dir}/${job_id}

cp ${macro} .

input_file=/sdf/group/hps/users/alspellm/projects/THESIS/ana/ecal_trig_res/20230628/tuples/ecal_trig_res/tritrig-beam_${job_id}_fsp.root
if [ ! -f "$input_file" ]; then
    echo "Input file does not exist: $input_file"
    exit 1
fi
echo "input file: $input_file"

cp ${input_file} .
output_file=./tritrig-beam_${job_id}_fsp_ana.root

root -l -b  <<EOF
gSystem->Load("libevent");
.x fsp_macro.cxx("${input_file}","${output_file}")
EOF

cp ${output_file} ${output_dir}
rm ${scratch_dir}/${job_id}/*
rm -r ${scratch_dir}/${job_id}
