#!/usr/bin/scl enable devtoolset-8 -- /bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=2000M
#SBATCH --array=1
#SBATCH --partition=hps
#SBATCH --job-name=zcutScan

while getopts o:t:s: flag
do
    case "${flag}" in
        o) outfile=${OPTARG};;
        t) lifetime=${OPTARG};;
        s) rundir=${OPTARG};;
    esac

done

export FIRST_ID=-1
export JOB_ID=$(($SLURM_ARRAY_TASK_ID+$FIRST_ID))
source /sdf/home/a/alspellm/.bashrc
runpath=$(readlink -f $rundir)
#export RUNDIR=${runpath}/${JOB_ID}
export RUNDIR=${runpath}/mpifpi3

mkdir -p $RUNDIR
cd $RUNDIR

#python3 /sdf/group/hps/users/alspellm/projects/simps_2019/reach_estimate/reach_estimates/lumiMakeExpSigSimps2019.py -o ${RUNDIR}/${outfile} -t ${lifetime}
python3 /sdf/group/hps/users/alspellm/projects/simps_2019/reach_estimate/reach_estimates/lumiMakeExpSigSimps2019.py -o ${RUNDIR}/${outfile} 

