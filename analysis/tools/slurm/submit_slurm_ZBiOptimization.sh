#!/usr/bin/scl enable devtoolset-8 -- /bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=6000M
#SBATCH --array=1-10
#SBATCH --partition=shared
#SBATCH --job-name=mass_zbi

while getopts o:s:b:z:m:v: flag
do
    case "${flag}" in
        o) outfile=${OPTARG};;
        s) slope=${OPTARG};;
        b) bkg=${OPTARG};;
        z) scan=${OPTARG};;
        m) stepsize=${OPTARG};; 
        v) vd_mass=${OPTARG};;
    esac
done

export JOB_ID=$(($SLURM_ARRAY_TASK_ID))
source /sdf/group/hps/users/alspellm/src/test/hpstr/install/bin/setup.sh
source /sdf/home/a/alspellm/.bashrc

export JOBDIR=/sdf/group/hps/users/alspellm/run/ZBi
export RUNDIR=${JOBDIR}/run/04072023/${vd_mass}/step_size_${stepsize}/${JOB_ID}
mkdir -p $RUNDIR
cd $RUNDIR

#export zalpha_slope=$((0.01*${JOB_ID}))
export zalpha_slope=$(python -c "print 0.005*$JOB_ID")
echo "zalpha slope is ${zalpha_slope}"

export OUTFILE=vdMass_${vd_mass}_zalpha_zbi_zcutscan_optimization_slope_${zalpha_slope}_step_size_${stepsize}.root
echo "Output file is ${OUTFILE}"

hpstr /sdf/group/hps/users/alspellm/src/test/hpstr/processors/config/zbiCutflow_cfg.py -D -o ${OUTFILE} -s ${zalpha_slope} -b ${bkg} -z ${scan} -m ${stepsize} -mass ${vd_mass} > run_log.txt

