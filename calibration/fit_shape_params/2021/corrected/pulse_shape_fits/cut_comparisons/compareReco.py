#!/usr/bin/python3
import ROOT as r
import numpy as np
import argparse

parser=argparse.ArgumentParser(description="")
parser.add_argument("-i", type=str, dest="inFilename", help="Input File",default="")
parser.add_argument("-n", type=str, dest="filen", help="file number",default="0")
parser.add_argument("-o", type=str, dest="outFilename", help="output root file",default="pulseFit_ana.root")
options = parser.parse_args()

###############################################################################################
r.gSystem.Load("libevent")

#inFile = r.TFile(options.inFilename,"READ")
outFile = r.TFile(options.outFilename,"RECREATE")
filen = options.filen

histos1d = {}
####### Histograms
chisq_slim_nocg7cut_h = r.TH1F("hit_chisqProb_slimAPV_nocg7cut","hit_chi2_slimAPV_nocg7cut;chisq_prob;entries",100,0,1)
chisq_slim_cg7cut_h = r.TH1F("hit_chisqProb_slimAPV_cg7cut","hit_chi2_slimAPV_cg7cut;chisq_prob;entries",100,0,1)

ampErr_slim_nocg7cut_h = r.TH1F("hit_AmpErr_slimAPV_nocg7cut","hit_AmpErr_slimAPV_nocg7cut;Amp Error (ADC);entries",500,0,1500)
ampErr_slim_cg7cut_h = r.TH1F("hit_AmpErr_slimAPV_cg7cut","hit_AmpErr_slimAPV_cg7cut;Amp Error (ADC);entries",500,0,1500)

t0Err_slim_nocg7cut_h = r.TH1F("hit_t0Err_slimAPV_nocg7cut","hit_t0Err_slimAPV_nocg7cut;t0 Err (ns);entries",150,0,15)
t0Err_slim_cg7cut_h = r.TH1F("hit_t0Err_slimAPV_cg7cut","hit_t0Err_slimAPV_cg7cut;t0 Err (ns);entries",150,0,15)


chisq_thick_nocg7cut_h = r.TH1F("hit_chisqProb_thickAPV_nocg7cut","hit_chi2_thickAPV_nocg7cut;chisq_prob;entries",100,0,1)
chisq_thick_cg7cut_h = r.TH1F("hit_chisqProb_thickAPV_cg7cut","hit_chi2_thickAPV_cg7cut;chisq_prob;entries",100,0,1)

ampErr_thick_nocg7cut_h = r.TH1F("hit_AmpErr_thickAPV_nocg7cut","hit_AmpErr_thickAPV_nocg7cut;Amp Error (ADC);entries",500,0,1500)
ampErr_thick_cg7cut_h = r.TH1F("hit_AmpErr_thickAPV_cg7cut","hit_AmpErr_thickAPV_cg7cut;Amp Error (ADC);entries",500,0,1500)

t0Err_thick_nocg7cut_h = r.TH1F("hit_t0Err_thickAPV_nocg7cut","hit_t0Err_thickAPV_nocg7cut;t0 Err (ns);entries",150,0,15)
t0Err_thick_cg7cut_h = r.TH1F("hit_t0Err_thickAPV_cg7cut","hit_t0Err_thickAPV_cg7cut;t0 Err (ns);entries",150,0,15)

##add histos
histos1d ["hit_chisqProb_slimAPV_nocg7cut"] = chisq_slim_nocg7cut_h
histos1d ["hit_chisqProb_slimAPV_cg7cut"] = chisq_slim_cg7cut_h

histos1d ["hit_AmpErr_slimAPV_nocg7cut"] = ampErr_slim_nocg7cut_h
histos1d ["hit_AmpErr_slimAPV_cg7cut"] = ampErr_slim_cg7cut_h

histos1d ["hit_t0Err_slimAPV_nocg7cut"] = t0Err_slim_nocg7cut_h
histos1d ["hit_t0Err_slimAPV_cg7cut"] = t0Err_slim_cg7cut_h

histos1d ["hit_chisqProb_thickAPV_nocg7cut"] = chisq_thick_nocg7cut_h
histos1d ["hit_chisqProb_thickAPV_cg7cut"] = chisq_thick_cg7cut_h

histos1d ["hit_AmpErr_thickAPV_nocg7cut"] = ampErr_thick_nocg7cut_h
histos1d ["hit_AmpErr_thickAPV_cg7cut"] = ampErr_thick_cg7cut_h

histos1d ["hit_t0Err_thickAPV_nocg7cut"] = t0Err_thick_nocg7cut_h
histos1d ["hit_t0Err_thickAPV_cg7cut"] = t0Err_thick_cg7cut_h

#######

infile_nocg7cut = '/sdf/group/hps/users/alspellm/projects/fit_shape_params/2021/corrected/pulse_shape_fits/cut_comparisons/nocut/reco/hps_14191_data_events_%s_nocut.root'%(filen)
infile_cg7cut = '/sdf/group/hps/users/alspellm/projects/fit_shape_params/2021/corrected/pulse_shape_fits/cut_comparisons/cut_cg7/reco/hps_14191_data_events_%s_cut.root'%(filen)
infiles = {'nocg7cut':infile_nocg7cut,'cg7cut':infile_cg7cut}

for key in infiles:
    infile = r.TFile(infiles[key],"READ")
    infile.cd()
    tree = infile.HPS_Event
    for i,e in enumerate(tree):
        #if i > 100:
        #    break
        if i%1000 == 0:
            print("event: ", i)
        for track in e.KalmanFullTracks:
            trackerhits = track.getSvtHits()
            for trackerhit in trackerhits:
                rawhits = trackerhit.getRawHits()
                for rawhit in rawhits:
                    if rawhit.getFitN() > 1:
                        continue
                    strip = rawhit.getStrip()
                    chisq = rawhit.getChiSq(0)
                    t0Err = rawhit.getT0err(0)
                    ampErr = rawhit.getAmpErr(0)
                    layer = rawhit.getLayer()
                    module = rawhit.getModule()
                    sensor = rawhit.getSensor()

                    if strip%8 != 7:
                        continue

                    if layer < 5:
                        histos1d['hit_chisqProb_slimAPV_%s'%(key)].Fill(chisq)
                        histos1d['hit_AmpErr_slimAPV_%s'%(key)].Fill(ampErr)
                        histos1d['hit_t0Err_slimAPV_%s'%(key)].Fill(t0Err)

                    else:
                        histos1d['hit_chisqProb_thickAPV_%s'%(key)].Fill(chisq)
                        histos1d['hit_AmpErr_thickAPV_%s'%(key)].Fill(ampErr)
                        histos1d['hit_t0Err_thickAPV_%s'%(key)].Fill(t0Err)
    infile.Close()

outFile.cd()
for histo in histos1d.values():
    histo.Write()
outFile.Close()
