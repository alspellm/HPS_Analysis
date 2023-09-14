#!/bin/usr/python
import ROOT as r
import numpy as np
import argparse

def fourPoleFit12():
    func = r.TF1("pulsefit","(TMath::Max(x-[0],0.0)/(x-[0]))*([3])*(([1]^2)/(([1]-[2])^(3))) * ( exp(-(x-[0])/[1]) - ( exp(-(x-[0])/[2]) * ( ((((([1]-[2])/([1]*[2]))*(x-[0]))^(0))) +  ((((([1]-[2])/([1]*[2]))*(x-[0]))^(1))) + ((((([1]-[2])/([1]*[2]))*(x-[0]))^(2))/2) ) ) ) + [4]")
    return func

def fourPole3Tau(t3):
    #A = (t1**2)/((t1-t2)*((t1-t3)**2)) #[10]
    #B = (t2**2)/((t2-t1)((t2-t3)**2)) #[11]
    #C = (1/(t2-t1)) * (((t1**2)/((t1-t3)**2)) - ((t2**2)/((t2-t3)**2)) ) #[12]
    #D = 1/((t3-t1)(t3-t2)) #[13]

    #A = (([1]**2)/(([1]-[2])*(([1]-[5])**2)))
    #B = (([2]**2)/(([2]-[1])*(([2]-[5])**2)))
    #C = ((1/([2]-[1])) * ((([1]**2)/(([1]-[5])**2)) - (([2]**2)/(([2]-[5])**2))))
    #D = (1/(([5]-[1])*([5]-[2]))) 

    #func = r.TF1("3TauFit", "[10]*exp((-x-[0])/[1]) + [11]*exp((-x-[0])/[2]) + ( ([12]+[13]*x) * exp((-x-[0])/[33]))")
    #func = r.TF1("fitfunc3tau", "(TMath::Max(x-[0],0.0)/(x-[0]))*([3])*( (( (([1])**2)/(([1]-[2])*(([1]-[5])**2)) )*exp((-x-[0])/[1])) + (( (([2])**2)/(([2]-[1])*(([2]-[5])**2)) )*exp((-x-[0])/[2]) ) + ( (( (1/([2]-[1])) * (((([1])**2)/(([1]-[5])**2)) - ((([2])**2)/(([2]-[5])**2)) ) ) + (( 1/(([5]-[1])*([5]-[2])) )*x)) * exp((-x-[0])/[5]))) + [4]")

    func = r.TF1("fitfunc3tau", " (TMath::Max(x-[0],0.0/(x-[0]))*[3]) * ( (((([1]**2)/(([1]-[2])*(([1]-[5])**2))))*exp(-(x-[0])/[1])) + (((([2]**2)/(([2]-[1])*(([2]-[5])**2))))*exp(-(x-[0])/[2])) + (( (((1/([2]-[1])) * ((([1]**2)/(([1]-[5])**2)) - (([2]**2)/(([2]-[5])**2))))) + (((1/(([5]-[1])*([5]-[2]))))*(x-[0])))*exp(-(x-[0])/[5])) ) + [4]")

    func.SetParameter(5,t3)

    return func

def setFitAllParams(func, t0, t1, t2, baseline, amp, fixt0 = False, fixt1 = False, fixt2 = False, fixamp = False, fixbaseline=False, color=r.kBlack):
    if fixt0:
        func.FixParameter(0,t0)
    else:
        func.SetParameter(0,t0)

    if fixt1:
        func.FixParameter(1,t1)
    else:
        func.SetParameter(1,t1)

    if fixt2:
        func.FixParameter(2,t2)
    else:
        func.SetParameter(2,t2)

    if fixamp:
        func.FixParameter(3,amp)
    else:
        func.SetParameter(3,amp)

    if fixbaseline:
        func.FixParameter(4, baseline)
    else:
        func.SetParameter(4, baseline)

    func.SetLineColor(color)
    func.SetLineWidth(2)

def getAmplitudeIntegralNorm(t1, t2):
    peak_t = 3.0*((t1* (t2**3))**(0.25))
    A = (t1**2)/((t1-t2)**3) 
    B = (t1-t2)/(t1*t2)
    peakAmp = A * (np.exp(-peak_t / t1) - np.exp(-peak_t/t2) * (1 + peak_t * B + 0.5 * peak_t * peak_t *B*B))
    #print(np.exp(-peak_t / t1))
    #print((np.exp(-peak_t/t2)))
    #print((1 + peak_t * B + 0.5 * peak_t * peak_t *B*B)) 
    return peakAmp


#########################################################################################
parser = argparse.ArgumentParser(description="baseConfig options ")
parser.add_argument("-o", "--outFile", type=str, dest="outFilename", action='store',
                  help="Output filename.", metavar="outFilename", default="out.root")
parser.add_argument('--infile', '-i', type=str, dest="inFilename", action='store',
                  help="Specify the input directory.", metavar="inFilename", default=".")
options = parser.parse_args()

infilename = options.inFilename
outfilename = options.outFilename
infile = r.TFile("%s"%(infilename), "READ")
outfile = r.TFile("%s"%(outfilename),"RECREATE")
infile.cd()


stand_chi2_hh = r.TH2F("Standard_Fit_Chi2_vs_Channel","Standard_Fit_Chi2;channel;Chi2/NDF",512, 0.0, 512.0, 2000,0.0,20000.0) 
fixed_chi2_hh = r.TH2F("Tau1_Tau2_Fixed_Fit_Chi2_vs_Channel","Tau1_Tau2_Fixed_Fit_Chi2;channel;Chi2/NDF",512, 0.0, 512.0, 50000,0.0,500000.0) 
alt_chi2_ch_hh = r.TH2F("Alternative_Fit_Chi2_vs_Channel","Fit_%s_Chi2;channel;Chi2/NDF",512, 0.0, 512.0, 5000,0.0,500000.0) 


stand_tau_hh = r.TH2F("Standard_Fit_Tau1_v_Tau2","Standard_Fit_Tau1_v_Tau2;Tau1;Tau2",1000, 0.0, 100.0, 1000,0.0,100.0) 
stand_tau1_hh = r.TH2F("Standard_Fit_Tau1","Standard_Fit_Tau1;channel;Tau1",512, 0.0, 512.0, 1000,0.0,100.0) 
stand_tau2_hh = r.TH2F("Standard_Fit_Tau2","Standard_Fit_Tau2;channel;Tau2",512, 0.0, 512.0, 1000,0.0,100.0) 


r.gROOT.SetBatch(r.kTRUE)
#r.gStyle.SetPalette(r.kDeepSea)

#Loop over channel 2d histograms
nchannels = 512
for channel in range(0,nchannels):
    hh = None
    canvas = r.TCanvas("ch_%s_pulse_fit"%(channel),"ch_%s_pulse_fit"%(channel),1800,800)
    for key in infile.GetListOfKeys():
        if "channel_%s_"%(channel) in key.GetName():
            hh = infile.Get(key.GetName())
            print(hh.GetName())
    t0 = 0.0
    maxamp = 0.0
    t1 = 45.0
    t2 = 10.0
    baseline = 0
    for tbin in range(1,hh.GetNbinsX()+1):
        h = hh.ProjectionY("projy",int(tbin), int(tbin))
        if tbin < 2:
            baseline = h.GetMean()
        amp = h.GetMean()
        if amp > maxamp:
            maxamp = amp

    outfile.cd()
    fitfunc = fourPoleFit12()
    setFitAllParams(fitfunc,t0,t1,t2,baseline,maxamp, color=r.kGreen)

    fitfuncFixTau = fourPoleFit12()
    setFitAllParams(fitfuncFixTau,t0,35.0,10.0,baseline,maxamp, fixt1=True, fixt2=True,fixbaseline=True, color=r.kRed)
    
    fitfunc3tau = fourPole3Tau(0.05)
    #setFitAllParams(fitfunc3tau,t0,t1,t2,baseline,maxamp, color=r.kMagenta)
    setFitAllParams(fitfunc3tau,26.0,35.0,0.04,baseline,maxamp, color=r.kMagenta)

    #standard fit
    hh.Fit(fitfunc,"q")
    chi2 = fitfunc.GetChisquare()
    ndf = fitfunc.GetNDF()
    chi2ndf = chi2/ndf
    tau1 = fitfunc.GetParameter(1)
    tau2 = fitfunc.GetParameter(2)
    print("chi2: ", chi2, " ndf: ", ndf, " chid/ndf: ", chi2ndf)
    print("tau1: ", tau1, " tau2: ", tau2)
    stand_chi2_hh.Fill(float(channel),chi2ndf,1.0)
    stand_tau_hh.Fill(tau1, tau2, 1.0) 
    stand_tau1_hh.Fill(float(channel), tau1, 1.0)
    stand_tau2_hh.Fill(float(channel), tau2, 1.0)

    hh.Fit(fitfuncFixTau, "+q")
    chi2 = fitfuncFixTau.GetChisquare()
    ndf = fitfuncFixTau.GetNDF()
    chi2ndf = chi2/ndf
    print("Fixed Tau chi2: ", chi2, " ndf: ", ndf, " chid/ndf: ", chi2ndf)
    fixed_chi2_hh.Fill(float(channel),chi2ndf,1.0)
    #hh.Fit(fitfunc3tau, "+")
    hh.Draw("colz")

    canvas.Write()



    #canvas.SaveAs("./ch_%s_pulse_fit.png"%(channel))
stand_chi2_hh.Write()
stand_tau_hh.Write()
stand_tau1_hh.Write()
stand_tau2_hh.Write()
fixed_chi2_hh.Write()


