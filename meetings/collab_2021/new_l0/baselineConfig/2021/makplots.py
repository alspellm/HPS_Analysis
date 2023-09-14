#!/usr/bin/python3
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import ROOT as r

layerMap = {"L0M01" : "L2T_stereo", "L0M02" : "L1B_stereo", "L0M03" : "L1T_stereo", "L0M04" : "L1B_axial", "L0M05" : "L2T_axial", "L0M06" : "L1T_axial", "L0M08" : "L2B_stereo", "L0M09" : "L2B_axial" }
mg = r.TMultiGraph()
mg.SetTitle("2021_Halfmodule_Production_0V_Bias_Noise")
for infile in sorted(glob.iglob("*anaHD.root")):
    module = infile.split("_")[0]
    tfile = r.TFile("%s"%(infile), "READ")
    tfile.cd()
    noise_g = tfile.Get("noiseSum_g").Clone()
    noise_g.SetTitle("%s"%(layerMap[module]))
    noise_g.SetLineWidth(2)
    mg.Add(noise_g, "lp")
    tfile.Close()

outfile = r.TFile("multigraph.root", "RECREATE")
outfile.cd()
c = r.TCanvas("QA_Noise_Plots_0V_Bias","plot", 1800,1800)
c.cd()
mg.GetXaxis().SetRange(0,512)
mg.GetXaxis().SetTitle("Physical Channel #")
mg.GetYaxis().SetTitle("Noise [ADC Units]")
mg.Draw("A pmc plc")
c.BuildLegend()
mg.Write()
c.SaveAs("halfmodule_production_noise_2021.png")








