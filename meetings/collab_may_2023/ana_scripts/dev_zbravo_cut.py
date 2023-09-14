#!/usr/bin/python3
import ROOT as r
import numpy as np
import root_numpy as rnp
import copy
import pandas as pd


#mass window
mV = 55.0
massRes = 3.0 #MeV for a 55MeV Vd
lowMass = float(mV) - 2.8*massRes/2.0
highMass = float(mV) + 2.8*massRes/2.0
print("Mass Window: ", lowMass, "--", highMass)

#define mass window selection
selection = {'unc_vtx_mass_lt' : highMass/1000., 'unc_vtx_mass_gt' : lowMass/1000.}

'''
#data
data_infile = '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/kf_data/kf_041823/output/hadd_hps_BLPass4_1958_files_recon_4.2_ana_kf.root'
data_treename = 'vtxana_kf_Tight_2016_simp_reach_dev/vtxana_kf_Tight_2016_simp_reach_dev_tree'
data_arr = rnp.root2array(data_infile, data_treename)
df = pd.DataFrame(data_arr)
df = df[(df['unc_vtx_mass'] < selection['unc_vtx_mass_lt']) & (df['unc_vtx_mass'] > selection['unc_vtx_mass_gt']) ]
outfile = r.TFile('simp_data_zbravo.root',"RECREATE")
'''
#mc signal
sig_infile = '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/kf_mc/hadd_mass_55_simp_recon_KF_ana.root'
sig_treename = 'vtxana_kf_Tight_2016_simp_reach_dev/vtxana_kf_Tight_2016_simp_reach_dev_tree'
sig_arr = rnp.root2array(sig_infile, sig_treename)
sig_df = pd.DataFrame(sig_arr)
df = sig_df[(sig_df['unc_vtx_mass'] < selection['unc_vtx_mass_lt']) & (sig_df['unc_vtx_mass'] > selection['unc_vtx_mass_gt']) ]
outfile = r.TFile('simp_signal_zbravo.root',"RECREATE")

plots = {}
#z0 vs recon_z
plots["z0_v_reconz_ele_top_hh"] = r.TH2F("%d_to_%d_MeV_ele_z0_v_reconz_top_hh"%(lowMass,highMass),"%d_%d_MeV_ele_z0_v_reconz_top;recon_z [mm];z0 [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)
plots["z0_v_reconz_pos_top_hh"] = r.TH2F("%d_to_%d_MeV_pos_z0_v_reconz_top_hh"%(lowMass,highMass),"%d_%d_MeV_pos_z0_v_reconz_top;recon_z [mm];z0 [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)
plots["z0_v_reconz_ele_bot_hh"] = r.TH2F("%d_to_%d_MeV_ele_z0_v_reconz_bot_hh"%(lowMass,highMass),"%d_%d_MeV_ele_z0_v_reconz_bot;recon_z [mm];z0 [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)
plots["z0_v_reconz_pos_bot_hh"] = r.TH2F("%d_to_%d_MeV_pos_z0_v_reconz_bot_hh"%(lowMass,highMass),"%d_%d_MeV_pos_z0_v_reconz_bot;recon_z [mm];z0 [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)


plots["zbravo_v_reconz_ele_top_hh"] = r.TH2F("%d_to_%d_MeV_ele_zbravo_v_reconz_top_hh"%(lowMass,highMass),"%d_%d_MeV_ele_zbravo_v_reconz_top;recon_z [mm];zbravo [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)
plots["zbravo_v_reconz_pos_bot_hh"] = r.TH2F("%d_to_%d_MeV_pos_zbravo_v_reconz_bot_hh"%(lowMass,highMass),"%d_%d_MeV_pos_zbravo_v_reconz_bot;recon_z [mm];zbravo [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)
plots["zbravo_v_reconz_ele_bot_hh"] = r.TH2F("%d_to_%d_MeV_ele_zbravo_v_reconz_bot_hh"%(lowMass,highMass),"%d_%d_MeV_ele_zbravo_v_reconz_bot;recon_z [mm];zbravo [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)
plots["zbravo_v_reconz_pos_top_hh"] = r.TH2F("%d_to_%d_MeV_pos_zbravo_v_reconz_top_hh"%(lowMass,highMass),"%d_%d_MeV_pos_zbravo_v_reconz_top;recon_z [mm];zbravo [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)

plots["zbravo_ele_top_pos_bot_hh"] = r.TH2F("%d_to_%d_MeV_zbravo_ele_top_pos_bot_hh"%(lowMass,highMass),"%d_%d_MeV_zbravo_ele_top_pos_bot;pos_bot zbravo;ele_top zbravo"%(lowMass, highMass),100,-4.0,4.0,100,-4.0,4.0)
plots["zbravo_ele_bot_pos_top_hh"] = r.TH2F("%d_to_%d_MeV_zbravo_ele_bot_pos_top_hh"%(lowMass,highMass),"%d_%d_MeV_zbravo_ele_bot_pos_top;pos_top zbravo;ele_bot zbravo"%(lowMass, highMass),100,-4.0,4.0,100,-4.0,4.0)

plots["zbravo_top_ele_h"] = r.TH1F("%d_to_%d_MeV_zbravo_top_ele_h"%(lowMass,highMass),"%d_%d_MeV_zbravo_top_ele_top;zbravo;events"%(lowMass, highMass),100,-4.0,4.0)
plots["zbravo_top_pos_h"] = r.TH1F("%d_to_%d_MeV_zbravo_top_pos_h"%(lowMass,highMass),"%d_%d_MeV_zbravo_top_pos_top;zbravo;events"%(lowMass, highMass),100,-4.0,4.0)
plots["zbravo_bot_ele_h"] = r.TH1F("%d_to_%d_MeV_zbravo_bot_ele_h"%(lowMass,highMass),"%d_%d_MeV_zbravo_bot_ele_bot;zbravo;events"%(lowMass, highMass),100,-4.0,4.0)
plots["zbravo_bot_pos_h"] = r.TH1F("%d_to_%d_MeV_zbravo_bot_pos_h"%(lowMass,highMass),"%d_%d_MeV_zbravo_bot_pos_bot;zbravo;events"%(lowMass, highMass),100,-4.0,4.0)

plots["zbravo_diff_ele_top_hh"] = r.TH2F("%d_to_%d_MeV_zbravo_diff_ele_top_hh"%(lowMass,highMass),"%d_%d_MeV_zbravo_diff_ele_top;recon_z [mm];zbravo [mm]"%(lowMass, highMass),100,-30.0,80.0,200,-8.0,8.0)
plots["zbravo_diff_ele_bot_hh"] = r.TH2F("%d_to_%d_MeV_zbravo_diff_ele_bot_hh"%(lowMass,highMass),"%d_%d_MeV_zbravo_diff_ele_bot;recon_z [mm];zbravo [mm]"%(lowMass, highMass),100,-30.0,80.0,200,-8.0,8.0)

plots["zbravo_sum_ele_top_hh"] = r.TH2F("%d_to_%d_MeV_zbravo_sum_ele_top_hh"%(lowMass,highMass),"%d_%d_MeV_zbravo_sum_ele_top;recon_z [mm];zbravo [mm]"%(lowMass, highMass),100,-30.0,80.0,200,-8.0,8.0)
plots["zbravo_sum_ele_bot_hh"] = r.TH2F("%d_to_%d_MeV_zbravo_sum_ele_bot_hh"%(lowMass,highMass),"%d_%d_MeV_zbravo_sum_ele_bot;recon_z [mm];zbravo [mm]"%(lowMass, highMass),100,-30.0,80.0,200,-8.0,8.0)

#Zbravo defined as difference between z0 and signal Zbravo line (Defined in fit_signal_z0_v_reconz)
Zbravo_top_ele = lambda z0, recon_z : z0 - (-0.039151*recon_z - 0.031282) 
Zbravo_bot_pos = lambda z0, recon_z : z0 - (0.039501*recon_z + 0.004176) 

Zbravo_bot_ele = lambda z0, recon_z : z0 - (0.040086*recon_z + 0.016186) 
Zbravo_top_pos = lambda z0, recon_z : z0 - (-0.037899*recon_z + -0.0094) 

for index, row in df.iterrows():
    #Top ele
    if row['unc_vtx_ele_track_tanLambda'] > 0.0:
        plots["z0_v_reconz_ele_top_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'])

        #1d histo zbravo
        plots["zbravo_top_ele_h"].Fill(Zbravo_top_ele(row['unc_vtx_ele_track_z0'],row['unc_vtx_z']))

        if row['unc_vtx_z'] > 10.0:
            plots["zbravo_v_reconz_ele_top_hh"].Fill(row['unc_vtx_z'], Zbravo_top_ele(row['unc_vtx_ele_track_z0'],row['unc_vtx_z']))
            plots["zbravo_ele_top_pos_bot_hh"].Fill(Zbravo_bot_pos(row['unc_vtx_pos_track_z0'],row['unc_vtx_z']), Zbravo_top_ele(row['unc_vtx_ele_track_z0'],row['unc_vtx_z']))
            #sum and diff
            plots["zbravo_diff_ele_top_hh"].Fill(row['unc_vtx_z'],(Zbravo_top_ele(row['unc_vtx_ele_track_z0'],row['unc_vtx_z']) - (Zbravo_bot_pos(row['unc_vtx_pos_track_z0'],row['unc_vtx_z']))) )
            plots["zbravo_sum_ele_top_hh"].Fill(row['unc_vtx_z'],(Zbravo_top_ele(row['unc_vtx_ele_track_z0'],row['unc_vtx_z']) + (Zbravo_bot_pos(row['unc_vtx_pos_track_z0'],row['unc_vtx_z']))) )


    #Bot ele
    else:
        plots["z0_v_reconz_ele_bot_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'])
        #1d histo zbravo
        plots["zbravo_bot_ele_h"].Fill(Zbravo_bot_ele(row['unc_vtx_ele_track_z0'],row['unc_vtx_z']))

        if row['unc_vtx_z'] > 10.0:
            plots["zbravo_v_reconz_ele_bot_hh"].Fill(row['unc_vtx_z'], Zbravo_bot_ele(row['unc_vtx_ele_track_z0'],row['unc_vtx_z']))
            plots["zbravo_ele_bot_pos_top_hh"].Fill(Zbravo_top_pos(row['unc_vtx_pos_track_z0'],row['unc_vtx_z']), Zbravo_bot_ele(row['unc_vtx_ele_track_z0'],row['unc_vtx_z']))
            #sum and diff
            plots["zbravo_diff_ele_bot_hh"].Fill(row['unc_vtx_z'],(Zbravo_bot_ele(row['unc_vtx_ele_track_z0'],row['unc_vtx_z']) - (Zbravo_top_pos(row['unc_vtx_pos_track_z0'],row['unc_vtx_z']))) )
            plots["zbravo_sum_ele_bot_hh"].Fill(row['unc_vtx_z'],(Zbravo_bot_ele(row['unc_vtx_ele_track_z0'],row['unc_vtx_z']) + (Zbravo_top_pos(row['unc_vtx_pos_track_z0'],row['unc_vtx_z']))) )


    #Top Pos
    if row['unc_vtx_pos_track_tanLambda'] > 0.0:
        plots["z0_v_reconz_pos_top_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'])
        #1d histo zbravo
        plots["zbravo_top_pos_h"].Fill(Zbravo_top_pos(row['unc_vtx_pos_track_z0'],row['unc_vtx_z']))

        if row['unc_vtx_z'] > 10.0:
            plots["zbravo_v_reconz_pos_top_hh"].Fill(row['unc_vtx_z'], Zbravo_top_pos(row['unc_vtx_pos_track_z0'],row['unc_vtx_z']))

    #Bot Pos
    else:
        plots["z0_v_reconz_pos_bot_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'])
        #1d histo zbravo
        plots["zbravo_bot_pos_h"].Fill(Zbravo_bot_pos(row['unc_vtx_pos_track_z0'],row['unc_vtx_z']))

        if row['unc_vtx_z'] > 10.0:
            plots["zbravo_v_reconz_pos_bot_hh"].Fill(row['unc_vtx_z'], Zbravo_bot_pos(row['unc_vtx_pos_track_z0'],row['unc_vtx_z']))

outfile.cd()
for plot in plots.values():
    plot.Write()

