#!/usr/bin/python3
import ROOT as r
import numpy as np

infilename = 'hadd_tritrig-beam_new_ecal_trig_res_fsp_ana.root'
outfilename = 'cluster_energy_fit_new_ecal_trig_res.root'
infile = r.TFile('%s'%(infilename), "READ")

fid_clustE_h = infile.Get('fiducial_fee_cluster_energy')
fit = r.TF1("fit", "gaus")
fit.SetParameter(1, 2.3)
fit.SetParameter(2, 0.2)
fit.SetRange(2.0,2.5)
fid_clustE_h.Fit("fit","ORQN", "")
fitMean = fit.GetParameter(1)
fitSig = fit.GetParameter(2)
best_chi2 = fit.GetChisquare()/fit.GetNDF()
print("Initial Mean:",fitMean)
print("Initial Sig:",fitSig)

rand = r.TRandom3()
rand.SetSeed(0)
bestMean = fitMean
bestSig = fitSig

best_xmin = bestMean - 2.0*bestSig
best_xmax = bestMean + 2.0*bestSig
for i in range(30):
    testmean = rand.Uniform(2.1, 2.4)
    testsig = rand.Uniform(0.05,0.5)
    print("testmean: ", testmean)
    print("testsig: ", testsig)
    fit.SetParameter(1, testmean)
    fit.SetParameter(2, testsig)
    xmin = testmean - 1.5*testsig
    xmax = testmean + 1.5*testsig
    print("xmin: ", xmin)
    print("xmax: ", xmax)
    fit.SetRange(xmin, xmax)
    fid_clustE_h.Fit("fit","ORQN","")

    testchi2 = fit.GetChisquare()/fit.GetNDF()
    if(testchi2 < best_chi2 and fit.GetParameter(1) > 2.0 and fit.GetParameter(1) < 2.5):
        bestMean = fit.GetParameter(1)
        bestSig = fit.GetParameter(2)
        best_chi2 = testchi2
        best_xmin = bestMean - 1.5*bestSig
        best_xmax = bestMean + 1.5*bestSig
        print("bestMean: ", bestMean)
        print("bestSig: ", bestSig)

print("post loop mean: ", bestMean)
print("post loop sig: ", bestSig)
print("post loop xmin: ", best_xmin)
print("post loop xmax: ", best_xmax)
#xmin = bestMean - 2.0*bestSig
#xmax = bestMean + 2.0*bestSig

fit.FixParameter(1, bestMean)
fit.FixParameter(2, bestSig)
fit.SetRange(best_xmin, best_xmax)
fid_clustE_h.Fit("fit","ORQ","")

text = r.TPaveText(0.50, 0.60, 0.60, 0.70, "NDC")
text.AddText('Param 1: %d'%(fit.GetParameter(1)))
text.AddText('Param 2: %d'%(fit.GetParameter(2)))

outfile = r.TFile('%s'%(outfilename),"RECREATE")
outfile.cd()
c = r.TCanvas("Fiducial_Cluster_Energy","",2500,1440)
c.cd()
fid_clustE_h.Draw()
text.Draw()
c.SaveAs("fiducial_cluster_energy.png")
