#!/usr/bin/python3
import ROOT as r
import os
import argparse
import subprocess
import re


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directory', type=str)

args = parser.parse_args()

directory = args.directory

total_count = 0
for filename in os.listdir(directory):
    if filename.endswith('.slcio'):
        f = os.path.join(directory, filename)
    else:
        continue
    print("Counting events in ", filename)

    # Your bash command as a string
    command = "lcio_event_counter %s"%(filename)

    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    stdout = result.stdout.decode("utf-8")  # Decode bytes to string
    stderr = result.stderr.decode("utf-8")  # Decode bytes to string

    # Use regex to extract the number at the end of the output
    number_match = re.search(r'\d+', stdout)
    
    output=stdout
    # Use regex to extract the number at the end of the output
    number_match = re.search(r'\d+', output)
             
    extracted_number = int(number_match.group())
    count = extracted_number
    print("Count is ", count)
    total_count = total_count + count

print("Total Event Count of files in %s is %i"%(directory, total_count))

