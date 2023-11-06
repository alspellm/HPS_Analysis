#!/usr/bin/scl enable devtoolset-8 -- /bin/bash
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=10000M
#SBATCH --array=1-30
#SBATCH --partition=hps
#SBATCH --job-name=zbi

OPTIND=1
# Check if the required arguments are provided
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <rundir> <config> <iterCutsJson> <cutVariables> <newVariables> <newVariableParams> <outputfile> <mass> <nbkg> <scan>"
fi

# Assign command line arguments to variables
rundir="$1"
config="$2"
iterCutsJson="$3"
cutVariables="$4"
newVariables="$5"
newVariableParams="$6"
outputfile="$7"
mass="$8"
nbkg="$9"
scan="${10}"

rundir=$(readlink -f "$rundir")
config=$(readlink -f "$config")
iterCutsJson=$(readlink -f "$iterCutsJson")


# Split the cutVariables, newVariables, and newVariableParams strings into arrays
IFS=' ' read -r -a cutVariablesArray <<< "$cutVariables"
IFS=' ' read -r -a newVariablesArray <<< "$newVariables"
IFS=' ' read -r -a newVariableParamsArray <<< "$newVariableParams"

# Display parameter values
echo "rundir: $rundir"
echo "config: $config"
echo "iterCutsJson: $iterCutsJson"
echo "cutVariables: ${cutVariablesArray[@]}"
echo "mass: $mass"
echo "nbkg: $nbkg"
echo "scan: $scan"

export JOB_ID=$(($SLURM_ARRAY_TASK_ID))
#SBATCH --output=log.${JOB_ID}.out
#SBATCH --error=log.${JOB_ID}.err
source /sdf/group/hps/users/alspellm/src/hpstr/install/bin/setup.sh
source /sdf/home/a/alspellm/.bashrc
#export signal_mass=$((40 + (JOB_ID-1)*5))
export OUTDIR=${rundir}/mass_${mass}
mkdir -p ${OUTDIR}
mkdir -p ${OUTDIR}/logs
cd $OUTDIR

zalphaSlope=$(python -c "print 0.001+(0.001*$JOB_ID)")
newVariableParamsArray=(${zalphaSlope} ${zalphaSlope})
newVariablesArray=("unc_vtx_ele_zalpha" "unc_vtx_pos_zalpha")
outputfile="zalpha_slope_${zalphaSlope}_mass_${mass}.root"
echo "newVariables: ${newVariablesArray[@]}"
echo "newVariableParams: ${newVariableParamsArray[@]}"
echo "outputfile: $outputfile"

# Run the command with the provided arguments
hpstr "$config" -D -cutsJson "$iterCutsJson" -cutVariables "${cutVariablesArray[@]}" -newVariables "${newVariablesArray[@]}" -newVariableParams "${newVariableParamsArray[@]}" -o "$outputfile" -mass "$mass" -b "$nbkg" -z "$scan" > "logs/${outputfile%.root}_log.txt"

