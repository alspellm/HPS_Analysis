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
chisq_nocut_h = r.TH1F("hit_chisqProb_apv127_to_127","hit_chi2_apv127_to_127;chisq_prob;entries",200,0,1)
chisq_withcut_h = r.TH1F("hit_chisqProb_apv127_to_126","hit_chi2_apv127_to_126;chisq_prob;entries",200,0,1)
chisq_withcut2_h = r.TH1F("hit_chisqProb_apv127_to_125","hit_chi2_apv127_to_125;chisq_prob;entries",200,0,1)

ampErr_nocut_h = r.TH1F("hit_AmpErr_apv127_to_127","hit_AmpErr_apv127_to_127;Amp Error (ADC);entries",500,0,1500)
ampErr_withcut_h = r.TH1F("hit_AmpErr_apv127_to_126","hit_AmpErr_apv127_to_126;Amp Error (ADC);entries",500,0,1500)
ampErr_withcut2_h = r.TH1F("hit_AmpErr_apv127_to_125","hit_AmpErr_apv127_to_125;Amp Error (ADC);entries",500,0,1500)

t0Err_nocut_h = r.TH1F("hit_t0Err_apv127_to_127","hit_t0Err_apv127_to_127;t0 Err (ns);entries",150,0,15)
t0Err_withcut_h = r.TH1F("hit_t0Err_apv127_to_126","hit_t0Err_apv127_to_126;t0 Err (ns);entries",150,0,15)
t0Err_withcut2_h = r.TH1F("hit_t0Err_apv127_to_125","hit_t0Err_apv127_to_125;t0 Err (ns);entries",150,0,15)

##add histos
histos1d ["hit_chisqProb_apv127_to_127"] = chisq_nocut_h
histos1d ["hit_chisqProb_apv127_to_126"] = chisq_withcut_h
histos1d ["hit_chisqProb_apv127_to_125"] = chisq_withcut2_h

histos1d ["hit_AmpErr_apv127_to_127"] = ampErr_nocut_h
histos1d ["hit_AmpErr_apv127_to_126"] = ampErr_withcut_h
histos1d ["hit_AmpErr_apv127_to_125"] = ampErr_withcut2_h

histos1d ["hit_t0Err_apv127_to_127"] = t0Err_nocut_h
histos1d ["hit_t0Err_apv127_to_126"] = t0Err_withcut_h
histos1d ["hit_t0Err_apv127_to_125"] = t0Err_withcut2_h
#######

apv127file = '/sdf/group/hps/users/alspellm/projects/fit_shape_params/2019/corrected/pulse_shape_fits/comparisons/files/3ns_delay/cut_apvch127/apv127_to_127/hps_10030_data_events_%s_apv127_to_127.root'%(filen)
apv126file = '/sdf/group/hps/users/alspellm/projects/fit_shape_params/2019/corrected/pulse_shape_fits/comparisons/files/3ns_delay/cut_apvch127/apv127_to_126/hps_10030_data_events_%s_apv127_to_126.root'%(filen)
apv125file = '/sdf/group/hps/users/alspellm/projects/fit_shape_params/2019/corrected/pulse_shape_fits/comparisons/files/3ns_delay/cut_apvch127/apv127_to_125/hps_10030_data_events_%s_apv127_to_125.root'%(filen)
infiles = {'apv127_to_127':apv127file,'apv127_to_126':apv126file, 'apv127_to_125':apv125file}

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
                if len(rawhits) > 1:
                    continue
                for rawhit in rawhits:
                    if rawhit.getFitN() > 1:
                        continue
                    layer = rawhit.getLayer()
                    if layer < 3:
                        continue
                    strip = rawhit.getStrip()
                    if strip%128 != 127:
                        continue
                    #print("Channel: ", strip)
                    chisq = rawhit.getChiSq(0)
                    t0Err = rawhit.getT0err(0)
                    ampErr = rawhit.getAmpErr(0)

                    #print("t0: ",t0Err)
                    #print("amp: ",ampErr)

                    histos1d['hit_chisqProb_%s'%(key)].Fill(chisq)
                    histos1d['hit_AmpErr_%s'%(key)].Fill(ampErr)
                    histos1d['hit_t0Err_%s'%(key)].Fill(t0Err)
    infile.Close()

outFile.cd()
for histo in histos1d.values():
    histo.Write()
outFile.Close()
