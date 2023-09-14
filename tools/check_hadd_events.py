#!/usr/bin/python3
import ROOT as r
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directory', type=str)

args = parser.parse_args()

directory = args.directory

total_count = 0
for filename in os.listdir(directory):
    if filename.endswith('.root'):
        f = os.path.join(directory, filename)
    else:
        continue
    print("Counting events in ", filename)
    infile = r.TFile('%s'%(filename),"READ")
    event_h = infile.Get('event_h')
    count = event_h.GetEntries()
    print("Count is ", count)
    infile.Close()
    total_count = total_count + count
print("Total Event Count of files in %s is %i"%(directory, total_count))

