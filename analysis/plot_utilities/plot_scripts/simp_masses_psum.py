#!/usr/bin/python3
import ROOT as r
import numpy as np
import glob as glob
import os
import copy

#tstack = r.THStack("tstack","")
colorsMap = { 15:r.kBlue, 14:r.kGreen, 12:r.kOrange, 13:r.kRed, 11:r.kYellow, 8:r.kMagenta, 9:r.kCyan, 10:r.kPink+1, 4:r.kSpring+10, 6:r.kViolet+2, 7:r.kTeal-1, 3:r.kOrange+7, 5:r.kMagenta-3, 2:r.kYellow-3, 1:r.kBlue+2, 0:r.kPink-9}
c = r.TCanvas("c","c",1800,1200)

index = 0
for f in sorted(glob.glob('/sdf/group/hps/users/alspellm/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_*.root')):
    name = os.path.basename(f)
    print(name)
    split = name.split('_')
    mass = int(split[2])

    if mass%10 != 0:
        continue
    #if mass == 20 or mass == 30 or mass == 40 or mass == 200:
    #    continue
    index = index + 1
    if index > 15:
        break

    infile = r.TFile(f,"READ")
    infile.cd()
    psum_h = copy.deepcopy(infile.Get('vtxana_kf_Tight_loose/vtxana_kf_Tight_loose_vtx_Psum_h'))
    psum_h.SetName('%s_mass_VD_mc_reco_Psum'%(mass))
    psum_h.SetTitle('%s_mass_VD;vtx Psum [GeV];Events'%(mass))
    psum_h.Rebin(2)
    psum_h.SetLineWidth(2)
    psum_h.SetLineColor(colorsMap[index])
    psum_h.GetYaxis().SetRangeUser(0.001,10000)

    c.cd()
    if index < 1:
        psum_h.Draw("hist")
    else:
        psum_h.Draw("histsame")
    #tstack.Add(psum_h)

outfile = r.TFile("/sdf/group/hps/users/alspellm/projects/THESIS/ana/reach_estimate_tenpct/components/simp_recon_vtx_psum_2016.root","RECREATE")
outfile.cd()
c.cd()
#tstack.Draw("hist nostack")
legend = c.BuildLegend(0.75,0.75,0.95,0.95)
legend.Draw()
c.Write()

#r.gPad.BuildLegend(0.75,0.75,0.95,0.95,"")

#tstack.Write()



