#/bin/env python
import math
import glob
import numpy as np
import ROOT as r
import utilities as utils
import copy
from optparse import OptionParser

Lumi = 10.7 #1/pb

utils.SetStyle()

parser = OptionParser()

parser.add_option("-o", "--outputFile", type="string", dest="outputFile",
            help="Specify the output filename.", metavar="outputFile", default="testOut.root")

(options, args) = parser.parse_args()
outfile = r.TFile("%s"%(options.outputFile),"RECREATE")

mcScale = {}
mcScale['tritrig'] = 1.416e9*Lumi/(50000*9853) #pb2016
mcScale['wab'] = 0.1985e12*Lumi/(100000*9966) #pb2016


ttFile = r.TFile("","READ")
ttTree = ttFile.Get("vtxana_Tight_simpCR/vtxana_Tight_simpSR_tree")
ttTree.SetName("tritrig_Tight_tree")
print("Counting background rate")
Mbin = 10.0
tt_mass_h = r.TH1F("tritrig_SR_mass","tritrig_SR_mass;mass [MeV];events",200,0.0,200.0) 
tt_mass_z_h = r.TH1F("tritrig_SR_mass_zcut","tritrig_SR_mass;mass [MeV];events",200,0.0,200.0) 
for ev in ttTree:
    tt_mass_h.Fill(1000.0*ev.unc_vtx_mass)
ttFile.Close()


wabFile = r.TFile("","READ")
wabTree = wabFile.Get("vtxana_Tight_simpCR/vtxana_Tight_simpSR_tree")
wabTree.SetName("wab_Tight_tree")
print("Counting background rate")
Mbin = 10.0
wab_mass_h = r.TH1F("wab_SR_mass","wab_SR_mass;mass [MeV];events",200,0.0,200.0) 
for ev in wabTree:
    wab_mass_h.Fill(1000.0*ev.unc_vtx_mass)
wabFile.Close()
