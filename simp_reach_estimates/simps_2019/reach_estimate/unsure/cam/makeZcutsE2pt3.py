#/bin/env python
import glob
import numpy as np
import ROOT as r
import utilities as utils
import copy
from optparse import OptionParser

def vtxRes(mass):
    res = ( 5.68598 - 7.24924e01*mass + 4.01760e02*pow(mass,2) - 7.38383e02*pow(mass,3) ) #alic
    return res

def massRes(mass):
    res = 2.28198 + 8.83116e-3*mass + 9.26580e-05*mass*mass #alic
    return res

Lumi = 11008.61649e-3
#utils.SetStyle()

parser = OptionParser()

parser.add_option("-i", "--inputFile", type="string", dest="inputFile",
    help="Name of file to run on.", metavar="inputFile", default="toys/toys.root")
parser.add_option("-o", "--outputFile", type="string", dest="outputFile",
    help="Specify the output filename.", metavar="outputFile", default="testOut.root")

(options, args) = parser.parse_args()

r.gROOT.SetBatch(1)

invMasses = [50,60,70,80,90,100,110,120,130,140,150,160,170]

#Calculate the weights
mcScale = {}
mcScale['tritrig'] = 1.416e9*Lumi/(50000*9853) #pb2016
mcScale['wab'] = 0.1985e12*Lumi/(100000*9966) #pb2016


outFile = r.TFile("zcuts.root","RECREATE")

#Get unbinnded MC after selection
ttFile = r.TFile("../coll/tritrig_anaVtx_Coll.root")
ttTree = ttFile.Get("vtxana_Tight/vtxana_Tight_tree")
ttTree.SetName("tritrig_Tight_tree")

#wabFile = r.TFile("../coll/wab_anaVtx_Coll.root")
#wabTree = wabFile.Get("vtxana_Tight/vtxana_Tight_tree")
#wabTree.SetName("wab_Tight_tree")


rand = r.TRandom3()
rand.SetSeed(0)

vtx_z_mass_hh = r.TH2D("vtx_z_mass_hh", "vtx_z_mass_hh", 300, 0, 300, 50, -30.0, 20.0)
vtx_mass_h = r.TH1D("vtx_mass_h", "vtx_mass_h", 300, 0, 300)
for ev in ttTree:
    vtx_mass_h.Fill(1000.0*ev.unc_vtx_mass, mcScale['tritrig'])
    vtx_z_mass_hh.Fill(1000.0*ev.unc_vtx_mass, ev.unc_vtx_z, mcScale['tritrig'])
    pass

zcutTxtFile = open("zcuts.dat","w")
#zcutTxtFile.write("mass/D:zcut/D\n")
zcuts = []
zcutsZres = []
masses = [float(x) for x in invMasses]
for mass in invMasses:
    massF = float(mass)
    print "Running %i MeV with res %f"%(mass, massRes(massF))
    lowMass = massF - 3.0*massRes(massF)/2.0
    highMass = massF + 3.0*massRes(massF)/2.0
    #print "%f\t%f\t%f\t%f"%(lowMass, massF, highMass, massRes(massF))
    zVtx_h = r.TH1D("zVtx%i_h"%mass, "zVtx%i_h"%mass, 150, -50.0, 100.0)
    zCutNorms = []
    norms = []
    #Fill histogram from unbinned data
    for ev in ttTree:
        #print "vtx_mass: %f. lowMass: %f . highMass: %f"%(1000.0*ev.unc_vtx_z, lowMass, highMass)
        if 1000.0*ev.unc_vtx_mass < lowMass: continue
        if 1000.0*ev.unc_vtx_mass > highMass: continue
        zVtx_h.Fill(ev.unc_vtx_z, mcScale['tritrig'])
        pass
    #for ev in wabTree:
    #    if 1000.0*ev.unc_vtx_mass < lowMass: continue
    #    if 1000.0*ev.unc_vtx_mass > highMass: continue
    #    zVtx_h.Fill(ev.unc_vtx_z, mcScale['wab'])
    #    pass
    #fitFunc = r.TF1("fit%i_f"%mass,"[0]*TMath::Exp([1]*x)", -10.0, 90.0)
    #fitFunc = r.TF1("fit%i_f"%mass,"[0]*TMath::Exp([1]*x)+[2]*TMath::Exp([3]*x)", -10.0, 90.0)
    fitFunc = r.TF1("fitfunc","[0]*exp( (((x-[1])/[2])<[3])*(-0.5*(x-[1])^2/[2]^2) + (((x-[1])/[2])>=[3])*(0.5*[3]^2-[3]*(x-[1])/[2]))", -100.0, 100.0)
    gausResult = zVtx_h.Fit("gaus","QS")
    gausParams = gausResult.GetParams()
    gausResult = zVtx_h.Fit("gaus","QS","",gausParams[1]-3.0*gausParams[2],gausParams[1]+3.0*gausParams[2])
    gausParams = gausResult.GetParams()
    tailZ = gausParams[1] + 3.0*gausParams[2]
    bestChi2 = -99.9
    bestParams = [999.9, 999.9, 999.9, 999.9]
    bestFitInit = [999.9, 999.9, 999.9, 999.9]
    #for fitI in range(10):
    #    fitInit = [rand.Uniform(500.0, 1500.0),
    #               rand.Uniform(-2.0,-1.0), 
    #               rand.Uniform(50, 150), 
    #               rand.Uniform(-1.0, -0.5)]
    #    fitFunc.SetParameters(fitInit[0], fitInit[1], fitInit[2], fitInit[3])
    #    #fitResult = zVtx_h.Fit(fitFunc, "LSIM", "", gausParams[1]-2.0*gausParams[2], gausParams[1]+10.0*gausParams[2])
    #    fitResult = zVtx_h.Fit(fitFunc, "LES", "", tailZ, 90.0)
    #    if not fitResult.IsValid(): continue
    #    if fitResult.Chi2() < bestChi2 or bestChi2 < 0.0:
    #        bestChi2 = fitResult.Chi2()
    #        bestParams = fitResult.GetParams()
    #        bestFitInit = fitInit
    #    pass
    fitFunc.SetParameters(gausParams[0], gausParams[1], gausParams[2], 3.0)
    fitResult = zVtx_h.Fit(fitFunc, "QLSIM", "", gausParams[1]-2.0*gausParams[2], gausParams[1]+10.0*gausParams[2])
    #fitFunc.SetParameters(bestFitInit[0], bestFitInit[1], bestFitInit[2], bestFitInit[3])
    #fitResult = zVtx_h.Fit(fitFunc, "LES", "", tailZ, 90.0)
    fullNorm = fitResult.GetParams()[0]
    for normI in range(1, 101):
        nEv = 0.5 #float(normI)*0.01
        zcut = -6.0
        fitFunc.SetParameter(0, fullNorm*normI/100.0)
        testIntegral = fitFunc.Integral(zcut, 90.0)
        while testIntegral > nEv:
            zcut = zcut+0.01
            testIntegral = fitFunc.Integral(zcut, 90.0)
            pass
        norms.append(Lumi*normI/100.0)
        zCutNorms.append(zcut)
        pass
    zcutNorms_g = r.TGraph(len(norms),np.array(norms),np.array(zCutNorms))
    zcutNorms_g.SetName("zcutNorms%i_g"%mass)
    zcutNorms_g.SetTitle(";lumi [1/pb];z_{cut} [mm]")
    outFile.cd()
    zcutNorms_g.Write()
    print "Zcut: %f"%zcut
    zcuts.append(zcut)
    zcutsZres.append(1.5 + 8.0*vtxRes(mass/1000.0))
    zcutTxtFile.write("%f\t%f\n"%(massF, zcut))
    zVtx_h.Write()
    pass
zcutTxtFile.close()

utils.SetStyle()
zcuts_g = r.TGraph(len(masses),np.array(masses),np.array(zcuts))
zcuts_g.SetName("zcuts_g")
zcuts_g.SetTitle(";m_{A'} [MeV];z_{cut} [mm]")

zcutsZres_g = r.TGraph(len(masses),np.array(masses),np.array(zcutsZres))
zcutsZres_g.SetName("zcutsZres_g")
zcutsZres_g.SetTitle(";m_{A'} [MeV];z_{cut} [mm]")

zcutsGeV_g = r.TGraph(len(masses),np.array([x/1000.0 for x in masses]),np.array(zcuts))
zcutsGeV_g.SetName("zcutsGeV_g")
zcutsGeV_g.SetTitle(";m_{A'} [GeV];z_{cut} [mm]")

vtx_mass_h.Write()
vtx_z_mass_hh.Write()
zcutsGeV_g.Write()
zcuts_g.Write()
zcutsZres_g.Write()

outFile.Close()
exit(0)

canv = r.TCanvas("canv", "canv", 1400, 1000)
canv.cd()
massResScaled_ge.SetMinimum(0.0)
massResScaled_ge.SetMaximum(15.0)
massResScaled_ge.Draw("ape")
massRes_ge.Draw("pesame")
fitFunc.Draw("same")
utils.InsertText()
canv.SaveAs("massRes.png")

massRes_ge.Write()
massResScaled_ge.Write()

outFile.Close()
