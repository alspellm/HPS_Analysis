#!/bin/bash

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

# Split the cutVariables, newVariables, and newVariableParams strings into arrays
IFS=' ' read -r -a cutVariablesArray <<< "$cutVariables"
IFS=' ' read -r -a newVariablesArray <<< "$newVariables"
IFS=' ' read -r -a newVariableParamsArray <<< "$newVariableParams"

# Display parameter values
echo "config: $config"
echo "iterCutsJson: $iterCutsJson"
echo "cutVariables: ${cutVariablesArray[@]}"
echo "newVariables: ${newVariablesArray[@]}"
echo "newVariableParams: ${newVariableParamsArray[@]}"
echo "outputfile: $outputfile"
echo "mass: $mass"
echo "nbkg: $nbkg"
echo "scan: $scan"

# Run the command with the provided arguments
hpstr "$config" -D -cutsJson "$iterCutsJson" -cutVariables "${cutVariablesArray[@]}" -newVariables "${newVariablesArray[@]}" -newVariableParams "${newVariableParamsArray[@]}" -o "$outputfile" -mass "$mass" -b "$nbkg" -z "$scan"

