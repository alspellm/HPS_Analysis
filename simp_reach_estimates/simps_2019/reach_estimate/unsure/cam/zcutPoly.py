import ROOT as r
from copy import deepcopy
import numpy as np
from optparse import OptionParser

parser = OptionParser()

parser.add_option ("-i", type="string", dest="inFile", default="")
(options,args) = parser.parse_args()

inFilename = options.inFile
inFile = r.TFile(inFilename, "READ")
zcuts_g = inFile.Get("zcuts_g")

canv = r.TCanvas("canv", "canv", 1800, 1000)
canv.cd()

#invm_sigma_h.Fit("pol3", "ES","", frl, fru)
zcuts_g.Fit("pol4", "ES", "", 50.0, 170.0)
zcuts_g.Draw()
canv.SaveAs("zcutPoly.png")
