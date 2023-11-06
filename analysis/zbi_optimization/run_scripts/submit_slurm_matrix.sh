#!/usr/bin/scl enable devtoolset-8 -- /bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=4000M
#SBATCH --array=1-100
#SBATCH --partition=hps
#SBATCH --job-name=z0LR
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

export JOBDIR=/sdf/group/hps/users/alspellm/run/isolation_cut_dev/20230724_100MeV
export OUTDIR=z0tanlambda_left_and_right
export RUNDIR=${JOBDIR}/${OUTDIR}/${JOB_ID}
mkdir -p $RUNDIR
cd $RUNDIR

export param_slope_1=$(python -c "print 0.0+(0.1*(($JOB_ID/10)))")
export param_slope_2=$(python -c "print 0.0+(0.1*(($JOB_ID)%10))")

#export OUTFILE=${vd_mass}_zalpha_ele_${param_slope_1}_pos_${param_slope_2}_bkg20_sig5.root
export OUTFILE=signal_${vd_mass}_z0tanlambda_right_${param_slope_1}_left_${param_slope_2}_bkg20_sig5.root
echo "Output file is ${OUTFILE}"
hpstr /sdf/group/hps/users/alspellm/run/isolation_cut_dev/20230724_100MeV/simp_zbi_config.py -D -o ${OUTFILE} -mass ${vd_mass} -new_vars_params ${param_slope_1} ${param_slope_2} > "${OUTFILE%.root}_log.txt"

mv ${RUNDIR}/*.root ${JOBDIR}/${OUTDIR}/
mkdir ${JOBDIR}/${OUTDIR}/logs
mv ${RUNDIR}/*log* ${JOBDIR}/${OUTDIR}/logs/
rm -r ${RUNDIR}
