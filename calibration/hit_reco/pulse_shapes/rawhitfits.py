#!/bin/usr/python
import ROOT as r
import numpy as np


def getAmplitudeIntegralNorm(t1, t2):
    peak_t = 3.0*((t1* (t2**3))**(0.25))
    print("peaktime: ", peak_t)
    A = (t1**2)/((t1-t2)**3) 
    B = (t1-t2)/(t1*t2)
    peakAmp = A * (np.exp(-peak_t / t1) - np.exp(-peak_t/t2) * (1 + peak_t * B + 0.5 * peak_t * peak_t *B*B))
    print("peakAmp: %s peaktime: %s"%(peakAmp, peak_t))
    print("A: %s B: %s "%(A,B))
    print(np.exp(-peak_t / t1))
    print((np.exp(-peak_t/t2)))
    print((1 + peak_t * B + 0.5 * peak_t * peak_t *B*B)) 
    return peakAmp

#########################################################################################
infilename = "/home/alic/HPS/projects/hit_reco/data_events_3.root"
infile = r.TFile("%s"%(infilename), "READ")
print(infile)
infile.cd()


tree = infile.Get("HPS_Event")
nEntries = tree.GetEntries()
r.gSystem.Load("libevent")

pulse = []
noise = []
fitparms = []
#hits[0] = 4320
baseline = 4285 # hits[2]

for i,entry in enumerate(tree):
    if i > 0:
        break
    track = entry.KalmanFullTracks[0]
    hits = track.getSvtHits()
    print(hits)
    rawhit = hits[2].getRawHits()[0]
    print("Channel number: ", rawhit.getLayer()," ", rawhit.getModule()," ", rawhit.getStrip())
    fitT0 = rawhit.getT0(0)
    fitAmp = rawhit.getAmp(0)
    print("fit amp: ", fitAmp)
    print("fit T0: ", fitT0)
    print("NFITS: ", rawhit.getFitN())
    fitparms.append(fitT0)
    fitparms.append(fitAmp)
    adcs = rawhit.getADCs()
    for j,adc in enumerate(adcs):
        if j > 5:
            break
        pulse.append(adc - baseline)
    print(pulse)
    print(fitT0, fitAmp)

##Define fit parameters
t0 = fitparms[0]
t0 = 40.0
fitAmp = fitparms[1]
t1 = 35.0
t2 = 10.0
peakAmp = getAmplitudeIntegralNorm(t1, t2)
#peakAmp = 1.7e-5

#Fit function as found in hps-java (Dropped O(3) term)
fitfuncmod = r.TF1("pulsefit","(TMath::Max(x-[0],0.0)/(x-[0]))*([3])*(([1]^2)/(([1]-[2])^(3))) * ( exp(-(x-[0])/[1]) - ( exp(-(x-[0])/[2]) * ( ((((([1]-[2])/([1]*[2]))*(x-[0]))^(0))) +  ((((([1]-[2])/([1]*[2]))*(x-[0]))^(1))) + ((((([1]-[2])/([1]*[2]))*(x-[0]))^(2))/2) ) ) )")
fitfuncmod.SetParameter(0,t0)
fitfuncmod.FixParameter(1,t1)
fitfuncmod.FixParameter(2,t2)
fitfuncmod.SetParameter(3,fitAmp/peakAmp)
fitfuncmod.SetLineColor(r.kBlack)

#Fit from hps-java
fixedfunc = r.TF1("pulsefit","(TMath::Max(x-[0],0.0)/(x-[0]))*([3])*(([1]^2)/(([1]-[2])^(3))) * ( exp(-(x-[0])/[1]) - ( exp(-(x-[0])/[2]) * ( ((((([1]-[2])/([1]*[2]))*(x-[0]))^(0))) +  ((((([1]-[2])/([1]*[2]))*(x-[0]))^(1))) + ((((([1]-[2])/([1]*[2]))*(x-[0]))^(2))/2) ) ) )")
fixedfunc.SetParameter(0,t0)
fixedfunc.FixParameter(1,t1)
fixedfunc.FixParameter(2,t2)
fixedfunc.FixParameter(3,fitAmp/peakAmp)
fixedfunc.SetLineColor(r.kMagenta)

#T1 = 45
fixedfunc45 = r.TF1("pulsefit","(TMath::Max(x-[0],0.0)/(x-[0]))*([3])*(([1]^2)/(([1]-[2])^(3))) * ( exp(-(x-[0])/[1]) - ( exp(-(x-[0])/[2]) * ( ((((([1]-[2])/([1]*[2]))*(x-[0]))^(0))) +  ((((([1]-[2])/([1]*[2]))*(x-[0]))^(1))) + ((((([1]-[2])/([1]*[2]))*(x-[0]))^(2))/2) ) ) )")
fixedfunc45.SetParameter(0,t0)
fixedfunc45.FixParameter(1,45.0)
fixedfunc45.FixParameter(2,t2)
fixedfunc45.FixParameter(3,fitAmp/peakAmp)
fixedfunc45.SetLineColor(r.kOrange)

#Fit function as defined in Shos Thesis and ref
fitfunc = r.TF1("pulsefit","(TMath::Max(x-[0],0.0)/(x-[0]))*([3])*(([1]^2)/(([1]-[2])^(3))) * ( exp(-(x-[0])/[1]) - ( exp(-(x-[0])/[2]) * ( ((((([1]-[2])/([1]*[2]))*(x-[0]))^(0))) +  ((((([1]-[2])/([1]*[2]))*(x-[0]))^(1))) + ((((([1]-[2])/([1]*[2]))*(x-[0]))^(2))/2) + ((((([1]-[2])/([1]*[2]))*(x-[0]))^(3))/6) ) ) )")
fitfunc.SetParameter(0,t0)
fitfunc.FixParameter(1,t1)
fitfunc.FixParameter(2,t2)
fitfunc.SetParameter(3,fitAmp/peakAmp)
fitfunc.SetLineColor(r.kRed)

#3rd order with T1 = 45
fitfunc45 = r.TF1("pulsefit","(TMath::Max(x-[0],0.0)/(x-[0]))*([3])*(([1]^2)/(([1]-[2])^(3))) * ( exp(-(x-[0])/[1]) - ( exp(-(x-[0])/[2]) * ( ((((([1]-[2])/([1]*[2]))*(x-[0]))^(0))) +  ((((([1]-[2])/([1]*[2]))*(x-[0]))^(1))) + ((((([1]-[2])/([1]*[2]))*(x-[0]))^(2))/2) + ((((([1]-[2])/([1]*[2]))*(x-[0]))^(3))/6) ) ) )")
fitfunc45.SetParameter(0,t0)
fitfunc45.FixParameter(1,45.0)
fitfunc45.FixParameter(2,t2)
fitfunc45.SetParameter(3,fitAmp/peakAmp)
fitfunc45.SetLineColor(r.kGreen)


canvas = r.TCanvas("can","can",1800,800)


#ADC data
ts = np.array(np.arange(0,144,24), dtype=float)
data = np.array(pulse, dtype=float)
pl = r.TGraph(len(data), ts, data)
#pl = r.TGraphErrors(len(data), ts, data, np.array([x*0.0 for x in data], dtype=float),  
pl.SetMarkerStyle(8)
pl.SetMarkerSize(3)
pl.Draw()

#Fit data
pl.Fit(fitfuncmod,"","",0.0,144.0)
pl.Fit(fitfunc,"+","",0.0,144.0)
pl.Fit(fixedfunc,"+","",0.0,144.0)
pl.Fit(fixedfunc45,"+","",0.0,144.0)
pl.Fit(fitfunc45,"+","",0.0,144.0)


canvas.SaveAs("./fit.png")
