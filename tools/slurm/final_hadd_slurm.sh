#!/bin/bash

# Reset OPTIND to its initial value (1)
OPTIND=1
rundir=/scratch/alspellm/hadd
outfilename='final_hadd_default.root'
while getopts r:s:o: flag
do
    case "${flag}" in
        r) rundir=${OPTARG};;
        s) sourcedir=${OPTARG};;
        o) outfilename=${OPTARG};;
    esac
done

export Rundir=$(readlink -f "$rundir")
export Sourcedir=$(readlink -f "$sourcedir")
echo "$Rundir"
echo "$Sourcedir"

# SLURM job submission command
sbatch <<EOT
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=01:00:00  # Adjust the time limit as needed
#SBATCH --mem=4000M
#SBATCH --job-name=final_hadd
#SBATCH --output=${Rundir}/log_final_hadd.out
#SBATCH --error=${Rundir}/log_final_hadd.err
#SBATCH --partition=shared

source /sdf/home/a/alspellm/.bashrc

mkdir -p ${Rundir}
cd $Rundir
hadd ${Rundir}/${outfilename} ${Sourcedir}/*.root 

EOT

