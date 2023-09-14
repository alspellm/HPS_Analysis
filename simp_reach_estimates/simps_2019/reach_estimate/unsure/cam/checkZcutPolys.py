#/bin/env python

import numpy as np
import ROOT as r
import matplotlib.pyplot as plt

invMasses = [50,60,70,80,90,100,110,120,130,140,150,160,170]
zCut1Vals =[]
zCut2Vals =[]
zCut3Vals =[]
zCut4Vals =[]
outFile = r.TFile("out.root","RECREATE")
outFile.cd()
canvas = r.TCanvas("cc","cc", 1800,800)
canvas.cd()
for mass in invMasses:
    massF = float(mass)
    massFGeV = mass/1000.0
    zCut1 = 12.7252 + 169.564*massFGeV - 5066.71*massFGeV*massFGeV + 39148*pow(massFGeV,3) - 101548*pow(massFGeV,4)
    zCut2 = 17.7702 + 138.166*massFGeV - 5363.29*massFGeV*massFGeV + 44532.4*pow(massFGeV,3) - 120578*pow(massFGeV,4)
    zCut3 = 30.0 + 138.166*massFGeV - 5363.29*massFGeV*massFGeV + 44532.4*pow(massFGeV,3) - 120578*pow(massFGeV,4)
    zCut4 = 21.6343 + 0.612759*massFGeV - 0.014178*massFGeV*massFGeV + 9.58558e-5*pow(massFGeV,3) - 2.12987e-7*pow(massFGeV,4)
    zCut1Vals.append(zCut1)
    zCut2Vals.append(zCut2)
    zCut3Vals.append(zCut3)
    zCut4Vals.append(zCut4)

gr1 = r.TGraph(len(invMasses), np.array([float(x) for x in invMasses]), np.array(zCut1Vals))
gr2 = r.TGraph(len(invMasses), np.array([float(x) for x in invMasses]), np.array(zCut2Vals))
gr3 = r.TGraph(len(invMasses), np.array([float(x) for x in invMasses]), np.array(zCut3Vals))
gr4 = r.TGraph(len(invMasses), np.array([float(x) for x in invMasses]), np.array(zCut4Vals))


gr1.Draw()



gr1.SetLineColor(2)
gr1.Draw()
gr2.Draw("SAME")
gr2.SetLineColor(3)
gr3.Draw("SAME")
gr3.SetLineColor(6)
gr4.Draw("SAME")
gr4.SetLineColor(7)

gr1.GetYaxis().SetRangeUser(0.0, 30.0)

legend = r.TLegend(0.1,0.7,0.48,0.9)
legend.AddEntry(gr1, "12.7252 + 169.564*massFGeV - 5066.71*massFGeV*massFGeV + 39148*pow(massFGeV,3) - 101548*pow(massFGeV,4)","l")
legend.AddEntry(gr2, "17.7702 + 138.166*massFGeV - 5363.29*massFGeV*massFGeV + 44532.4*pow(massFGeV,3) - 120578*pow(massFGeV,4)","l")
legend.AddEntry(gr3, "30.0 + 138.166*massFGeV - 5363.29*massFGeV*massFGeV + 44532.4*pow(massFGeV,3) - 120578*pow(massFGeV,4)","l")
legend.AddEntry(gr4, "21.6343 + 0.612759*massFGeV - 0.014178*massFGeV*massFGeV + 9.58558e-5*pow(massFGeV,3) - 2.12987e-7*pow(massFGeV,4)","l")
legend.Draw()
canvas.Draw()
canvas.Write()


