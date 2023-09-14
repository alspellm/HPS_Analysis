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
chisq_slim_3125ns_h = r.TH1F("hit_chisqProb_slimAPV_3125ns","hit_chi2_slimAPV_3125ns;chisq_prob;entries",100,0,1)
chisq_slim_3_h = r.TH1F("hit_chisqProb_slimAPV_3ns","hit_chi2_slimAPV_3ns;chisq_prob;entries",100,0,1)

ampErr_slim_3125ns_h = r.TH1F("hit_AmpErr_slimAPV_3125ns","hit_AmpErr_slimAPV_3125ns;Amp Error (ADC);entries",500,0,1500)
ampErr_slim_3_h = r.TH1F("hit_AmpErr_slimAPV_3ns","hit_AmpErr_slimAPV_3ns;Amp Error (ADC);entries",500,0,1500)

t0Err_slim_3125ns_h = r.TH1F("hit_t0Err_slimAPV_3125ns","hit_t0Err_slimAPV_3125ns;t0 Err (ns);entries",150,0,15)
t0Err_slim_3_h = r.TH1F("hit_t0Err_slimAPV_3ns","hit_t0Err_slimAPV_3ns;t0 Err (ns);entries",150,0,15)


chisq_thick_3125ns_h = r.TH1F("hit_chisqProb_thickAPV_3125ns","hit_chi2_thickAPV_3125ns;chisq_prob;entries",100,0,1)
chisq_thick_3_h = r.TH1F("hit_chisqProb_thickAPV_3ns","hit_chi2_thickAPV_3ns;chisq_prob;entries",100,0,1)

ampErr_thick_3125ns_h = r.TH1F("hit_AmpErr_thickAPV_3125ns","hit_AmpErr_thickAPV_3125ns;Amp Error (ADC);entries",500,0,1500)
ampErr_thick_3_h = r.TH1F("hit_AmpErr_thickAPV_3ns","hit_AmpErr_thickAPV_3ns;Amp Error (ADC);entries",500,0,1500)

t0Err_thick_3125ns_h = r.TH1F("hit_t0Err_thickAPV_3125ns","hit_t0Err_thickAPV_3125ns;t0 Err (ns);entries",150,0,15)
t0Err_thick_3_h = r.TH1F("hit_t0Err_thickAPV_3ns","hit_t0Err_thickAPV_3ns;t0 Err (ns);entries",150,0,15)

##add histos
histos1d ["hit_chisqProb_slimAPV_3125ns"] = chisq_slim_3125ns_h
histos1d ["hit_chisqProb_slimAPV_3ns"] = chisq_slim_3_h

histos1d ["hit_AmpErr_slimAPV_3125ns"] = ampErr_slim_3125ns_h
histos1d ["hit_AmpErr_slimAPV_3ns"] = ampErr_slim_3_h

histos1d ["hit_t0Err_slimAPV_3125ns"] = t0Err_slim_3125ns_h
histos1d ["hit_t0Err_slimAPV_3ns"] = t0Err_slim_3_h

histos1d ["hit_chisqProb_thickAPV_3125ns"] = chisq_thick_3125ns_h
histos1d ["hit_chisqProb_thickAPV_3ns"] = chisq_thick_3_h

histos1d ["hit_AmpErr_thickAPV_3125ns"] = ampErr_thick_3125ns_h
histos1d ["hit_AmpErr_thickAPV_3ns"] = ampErr_thick_3_h

histos1d ["hit_t0Err_thickAPV_3125ns"] = t0Err_thick_3125ns_h
histos1d ["hit_t0Err_thickAPV_3ns"] = t0Err_thick_3_h

#######

infile_3125ns = '/sdf/group/hps/users/alspellm/projects/fit_shape_params/2019/corrected/pulse_shape_fits/comparisons/files/3125ns_delay/hps_10030_data_events_%s_3125_withcuts.root'%(filen)
infile_3ns = '/sdf/group/hps/users/alspellm/projects/fit_shape_params/2019/corrected/pulse_shape_fits/comparisons/files/3ns_delay/cut_apvch127/hps_10030_data_events_%s_withcuts.root'%(filen)
infiles = {'3125ns':infile_3125ns,'3ns':infile_3ns}

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

                    if layer == 1:
                        histos1d['hit_chisqProb_slimAPV_%s'%(key)].Fill(chisq)
                        histos1d['hit_AmpErr_slimAPV_%s'%(key)].Fill(ampErr)
                        histos1d['hit_t0Err_slimAPV_%s'%(key)].Fill(t0Err)

                    elif layer == 4:
                        histos1d['hit_chisqProb_thickAPV_%s'%(key)].Fill(chisq)
                        histos1d['hit_AmpErr_thickAPV_%s'%(key)].Fill(ampErr)
                        histos1d['hit_t0Err_thickAPV_%s'%(key)].Fill(t0Err)

                    else:
                        continue
    infile.Close()

outFile.cd()
for histo in histos1d.values():
    histo.Write()
outFile.Close()
