#/bin/env python
import numpy as np
import ROOT as r
import utilities as utils
import copy
from optparse import OptionParser

Lumi = 11.00861649 #1/pb

utils.SetStyle()

parser = OptionParser()

parser.add_option("-i", "--inputFile", type="string", dest="inputFile",
    help="Name of file to run on.", metavar="inputFile", default="toys/toys.root")
parser.add_option("-o", "--outputFile", type="string", dest="outputFile",
    help="Specify the output filename.", metavar="outputFile", default="testOut.root")

(options, args) = parser.parse_args()

invMassHistos = {}
inFile = r.TFile("../coll/rad_anaVtx_Coll.root")
invMassHistos['rad'] = copy.deepcopy(inFile.Get("vtxana_radMatchTight/vtxana_radMatchTight_mcMass622_h"))
inFile.Close()
inFile = r.TFile("../coll/tritrig_anaVtx_Coll.root")
invMassHistos['tritrig'] = copy.deepcopy(inFile.Get("vtxana_Tight/vtxana_Tight_vtx_InvM_h"))
inFile.Close()
inFile = r.TFile("../coll/wab_anaVtx_Coll.root")
invMassHistos['wab'] = copy.deepcopy(inFile.Get("vtxana_Tight/vtxana_Tight_vtx_InvM_h"))
inFile.Close()

#Scale the histograms
rebin = 1
units = 0.5 #1/(2 MeV)
invMassHistos['rad'].Rebin(rebin)
invMassHistos['tritrig'].Rebin(rebin)
invMassHistos['wab'].Rebin(rebin)
invMassHistos['rad'].Scale(units*66.36e6*Lumi/(rebin*10000*9959))
invMassHistos['tritrig'].Scale(units*1.416e9*Lumi/(rebin*50000*9853))
invMassHistos['wab'].Scale(units*0.1985e12*Lumi/(rebin*100000*9966))

triWab_sh = r.THStack("triWab_sh",";m_{e^{+}e^{-}} [GeV];#frac{dN}{dm} [MeV^{-1}]")

triWab_sh.Add(invMassHistos['wab'])
triWab_sh.Add(invMassHistos['tritrig'])

#canv = r.TCanvas("canv","canv",1600,1200)
#canv.cd()

#invMassHistos['rad'].GetXaxis().SetRangeUser(0.05,0.2)
#invMassHistos['rad'].GetYaxis().SetTitle("#frac{dN}{dm} [1 / 2 MeV]")
#invMassHistos['rad'].Draw()
#canv.SetLogy(1)
canv = utils.MakeRadFrac("radFrac", ".", [invMassHistos['rad'],invMassHistos['wab'],invMassHistos['tritrig']], ['rad', 'wab', 'tritrig+wab'],'.png', RatioMin=0.00, RatioMax=0.3, LogY=True)

#canv = utils.MakeStackPlot("radFrac", ".", [invMassHistos['rad'],invMassHistos['tritrig']], ['rad', 'tritrig'],'.png', RatioMin=0.02, RatioMax=0.1, LogY=True)
