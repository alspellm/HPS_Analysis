#/bin/env python
import numpy as np
import ROOT as r
import utilities as utils
import copy
from optparse import OptionParser

Lumi = 10.7 #1/pb

utils.SetStyle()

parser = OptionParser()

parser.add_option("-i", "--inputFile", type="string", dest="inputFile",
    help="Name of file to run on.", metavar="inputFile", default="toys/toys.root")
parser.add_option("-o", "--outputFile", type="string", dest="outputFile",
    help="Specify the output filename.", metavar="outputFile", default="testOut.root")

(options, args) = parser.parse_args()
outfile = r.TFile("%s"%(options.outputFile),"RECREATE")

invMassHistos = {}
inFile = r.TFile("/sdf/group/hps/users/alspellm/projects/simps_2016/reach_estimate/mc/rad/slic/ana/hadd_RADv3_MG5_noXchange_HPS-PhysicsRun2016-Pass2_1_ana.root")
invMassHistos['rad_slic'] = copy.deepcopy(inFile.Get("mcAna/mcAna_mc622Mass_h"))
inFile.Close()
inFile = r.TFile("/sdf/group/hps/users/alspellm/projects/simps_2016/reach_estimate/mc/rad/rad_beam/ana/hadd_rad_beam_ana.root")
invMassHistos['rad'] = copy.deepcopy(inFile.Get("vtxana_radMatchTight_simps/vtxana_radMatchTight_simps_mcMass622_h"))
inFile.Close()

outfile.cd()

#Fix the x axis scaling of rad hist
name = invMassHistos['rad'].GetName()
invMassHistos['rad'].SetName("needs_rescaling")
rad_mev_h = r.TH1F("rescale_rad","rescale_rad",200,0.,200.)
nbins = invMassHistos['rad'].GetXaxis().GetNbins()
for b in range(nbins):
    val = invMassHistos['rad'].GetBinContent(b+1)
    rad_mev_h.SetBinContent(b+1,val)
rad_mev_h.SetName(name)
#rad_mev_h.SetTitle(invMassHistos['rad'].GetTitle())
invMassHistos['rad'] = rad_mev_h

#Scale the histograms
rebin = 1
units = 1 #1/(2 MeV)
#invMassHistos['rad'].Rebin(rebin)
invMassHistos['rad_slic'].Rebin(rebin)
#MC Scale
invMassHistos['rad_slic'].Scale(units*66.36e6*Lumi/(rebin*10000*999))
invMassHistos['rad'].Scale(units*66.36e6*Lumi/(rebin*10000*9967))

'''
before = copy.deepcopy(invMassHistos['rad_slic'])
before.SetName("before")
before.Write()

after = copy.deepcopy(invMassHistos['rad_slic'])
after.SetName("after")
after.GetXaxis().SetRange(1,200)
after.SetBins(200,0.0,200.0)
after.Write()
print("before n bins: ", before.GetNbinsX())
print("after n bins: ", after.GetNbinsX())
'''

#Make N bins the same 
invMassHistos['rad_slic'].GetXaxis().SetRange(1,200)
invMassHistos['rad'].GetXaxis().SetRange(1,200)
invMassHistos['rad_slic'].SetBins(200,0.0,200.0)
invMassHistos['rad'].SetBins(200,0.0,200.0)

print(invMassHistos['rad_slic'].GetNbinsX())
print(invMassHistos['rad'].GetNbinsX())

#Write histos
invMassHistos['rad_slic'].Write()
invMassHistos['rad'].Write()

#canv = r.TCanvas("canv","canv",1600,1200)
#canv.cd()
#invMassHistos['rad'].SetLineColor(r.kBlue)
#invMassHistos['rad_slic'].SetLineColor(r.kRed)
#invMassHistos['rad_slic'].Draw("hist")
#invMassHistos['rad'].Draw("histsame")
#canv.SetLogy(1)
#canv.Write()
#canv.SaveAs(".test.png")

#invMassHistos['rad'].GetXaxis().SetRangeUser(0.05,0.2)
#invMassHistos['rad'].GetYaxis().SetTitle("#frac{dN}{dm} [1 / 2 MeV]")
#invMassHistos['rad'].Draw()
#canv.SetLogy(1)
canv = utils.MakeTotalRadAcc("zeta", ".", [invMassHistos['rad'],invMassHistos['rad_slic']], ['rad', 'rad_slic'],'.png', RatioMin=0.00, RatioMax=0.3, LogY=True)
canv.SaveAs("./totRadAccept.png")
outfile.Write()
