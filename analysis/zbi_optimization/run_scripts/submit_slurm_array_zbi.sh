#!/usr/bin/scl enable devtoolset-8 -- /bin/bash
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=8000M
#SBATCH --array=1-33
#SBATCH --partition=shared
#SBATCH --job-name=z0flat
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# Check if the required arguments are provided
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <rundir> <config> <iterCutsJson> <cutVariables> <newVariables> <newVariableParams> <outputfile> <nbkg> <scan>"
fi

# Assign command line arguments to variables
rundir="$1"
config="$2"
iterCutsJson="$3"
cutVariables="$4"
newVariables="$5"
newVariableParams="$6"
outputfile="$7"
nbkg="$8"
scan="$9"

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
echo "newVariables: ${newVariablesArray[@]}"
echo "newVariableParams: ${newVariableParamsArray[@]}"
echo "outputfile: $outputfile"
echo "nbkg: $nbkg"
echo "scan: $scan"

export JOB_ID=$(($SLURM_ARRAY_TASK_ID))
source /sdf/group/hps/users/alspellm/src/hpstr/install/bin/setup.sh
source /sdf/home/a/alspellm/.bashrc
export mass=$((40 + (JOB_ID-1)*5))
echo "mass: ${mass}"
outputfile=${outputfile}_mass_${mass}.root
export OUTDIR=${rundir}/mass_${mass}
mkdir -p $OUTDIR
mkdir -p ${OUTDIR}/logs
cd $OUTDIR

# Run the command with the provided arguments
hpstr "$config" -D -cutsJson "$iterCutsJson" -cutVariables "${cutVariablesArray[@]}" -newVariables "${newVariablesArray[@]}" -newVariableParams "${newVariableParamsArray[@]}" -o "$outputfile" -mass "$mass" -b "$nbkg" -z "$scan" > "logs/${outputfile%.root}_log.txt"


mv ${OUTDIR}/*log* ${OUTDIR}/logs/

