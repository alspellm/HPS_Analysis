#!/bin/bash

OPTIND=1
while getopts s:e: opt; do
    case "${opt}" in
        s) first=${OPTARG};;
        e) end=${OPTARG};;
    esac
done

echo "first: $first"
echo "end: $end"


diff=$((end + 1 - first))
start_arr=1
end_arr=$diff
jobshift=1

if (( first > 1000 )); then
    jobshift=$first
fi
    

array_index="${start_arr}-${end_arr}"
echo "$array_index"
sbatch --array=${array_index} ./slurm_macro.sh -i ${jobshift} -m fsp_macro.cxx -o run/output -s /scratch/alspellm/ana

