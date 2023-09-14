#!/usr/bin/python3
import ROOT as r
import numpy as np

def zCut(mass_mev):
    zCut = 69.2555 + -0.916318*float(mass_mev) + 0.00504772*float(mass_mev)*float(mass_mev) + -1.04964e-05*pow(float(mass_mev),3)
    return zCut

outfilename = 'draw_my_zcut.root'
outfile = r.TFile('%s'%(outfilename),"RECREATE")

mass = np.arange(0.02, 0.2, 0.01, dtype=float)
zCut = [zCut(m*1000.0) for m in mass]
print(mass*1000.0)
print(zCut)

zCut = np.array(zCut, dtype=float)

outfile.cd()

gr = r.TGraph(len(mass), mass, zCut)
gr.SetName("zcut")
fit = gr.Fit('pol3','ES')
gr.Draw()
gr.Write()

outfile.Write()


