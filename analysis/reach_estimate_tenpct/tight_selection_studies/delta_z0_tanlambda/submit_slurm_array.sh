#!/usr/bin/scl enable devtoolset-8 -- /bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=4000M
#SBATCH --array=1-33
#SBATCH --partition=hps
#SBATCH --job-name=zbi_zcuts
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

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
source /sdf/group/hps/users/alspellm/src/hpstr/install/bin/setup.sh
source /sdf/home/a/alspellm/.bashrc

export JOBDIR=/sdf/group/hps/users/alspellm/run/simp_zbi_optimization/delta_z0_tanlambda
export OUTDIR=zbi_opt_cut_eq_35
export signal_mass=$((40 + (JOB_ID-1)*5))
export RUNDIR=${JOBDIR}/${OUTDIR}/mass_${signal_mass}
mkdir -p $RUNDIR
cd $RUNDIR

export OUTFILE=signal_${signal_mass}_logeps2_neg5pt5_run7800_bkg_model_SF_10_tenpct_data_CR.root

echo "Output file is ${OUTFILE}"

hpstr ${JOBDIR}/simp_zbi_config.py -D -o ${OUTFILE} -mass ${signal_mass} > "${OUTFILE%.root}_log.txt"

mv ${RUNDIR}/*.root ${JOBDIR}/${OUTDIR}/
mkdir ${JOBDIR}/${OUTDIR}/logs
mv ${RUNDIR}/*log* ${JOBDIR}/${OUTDIR}/logs/
rm -r ${RUNDIR}
