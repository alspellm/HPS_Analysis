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

#data
'''
data_infile = '/sdf/group/hps/users/alspellm/projects/THESIS/impact_param_studies_04182023/data/BLPass4b/ana_04182023/files/hadd_hps_BLPass4b_recon_4.2_simp_ana.root'
data_treename = 'vtxana_kf_Tight_2016_simp_reach_dev/vtxana_kf_Tight_2016_simp_reach_dev_tree'
data_arr = rnp.root2array(data_infile, data_treename)
df = pd.DataFrame(data_arr)
df = df[(df['unc_vtx_mass'] < selection['unc_vtx_mass_lt']) & (df['unc_vtx_mass'] > selection['unc_vtx_mass_gt']) ]
outfile = r.TFile('simp_impact_param_data_04182023.root',"RECREATE")
'''

#mc signal
sig_infile = '/sdf/group/hps/users/alspellm/projects/THESIS/impact_param_studies_04182023/mc/signal/hadd_mass_55_simp_recon_KF_ana.root'
sig_treename = 'vtxana_kf_Tight_2016_simp_reach_dev/vtxana_kf_Tight_2016_simp_reach_dev_tree'
sig_arr = rnp.root2array(sig_infile, sig_treename)
sig_df = pd.DataFrame(sig_arr)
df = sig_df[(sig_df['unc_vtx_mass'] < selection['unc_vtx_mass_lt']) & (sig_df['unc_vtx_mass'] > selection['unc_vtx_mass_gt']) ]
outfile = r.TFile('simp_impact_param_signal_04182023.root',"RECREATE")
plots = {}

#z0 vs recon_z
plots["z0_v_reconz_ele_top_hh"] = r.TH2F("%d_to_%d_MeV_ele_z0_v_reconz_top_hh"%(lowMass,highMass),"%d_%d_MeV_ele_z0_v_reconz_top;recon_z [mm];z0 [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)
plots["z0_v_reconz_pos_top_hh"] = r.TH2F("%d_to_%d_MeV_pos_z0_v_reconz_top_hh"%(lowMass,highMass),"%d_%d_MeV_pos_z0_v_reconz_top;recon_z [mm];z0 [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)
plots["z0_v_reconz_ele_bot_hh"] = r.TH2F("%d_to_%d_MeV_ele_z0_v_reconz_bot_hh"%(lowMass,highMass),"%d_%d_MeV_ele_z0_v_reconz_bot;recon_z [mm];z0 [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)
plots["z0_v_reconz_pos_bot_hh"] = r.TH2F("%d_to_%d_MeV_pos_z0_v_reconz_bot_hh"%(lowMass,highMass),"%d_%d_MeV_pos_z0_v_reconz_bot;recon_z [mm];z0 [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)

#z0 ele vs pos
zcuts = [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0,24.0]
for zcut in zcuts:
    plots["top_ele_z0_v_pos_z0_zcut_%s_hh"%(zcut)] = r.TH2F("%d_to_%d_MeV_top_ele_z0_v_pos_z0__zcut_%s_hh"%(lowMass,highMass,zcut),"%d_%d_MeV_top_ele_z0_v_pos_z0_zcut_%s;e^{+} z0 [mm]; e^{-} z0 [mm]"%(lowMass, highMass,zcut),100,-4.0,4.0,100,-4.0,4.0)
    plots["bot_ele_z0_v_pos_z0_zcut_%s_hh"%(zcut)] = r.TH2F("%d_to_%d_MeV_bot_ele_z0_v_pos_z0__zcut_%s_hh"%(lowMass,highMass,zcut),"%d_%d_MeV_bot_ele_z0_v_pos_z0_zcut_%s;e^{+} z0 [mm]; e^{-} z0 [mm]"%(lowMass, highMass,zcut),100,-4.0,4.0,100,-4.0,4.0)

#z0 vs d0
plots["z0_v_d0_ele_top_hh"] = r.TH2F("%d_to_%d_MeV_ele_z0_v_d0_top_hh"%(lowMass,highMass),"%d_%d_MeV_ele_z0_v_d0_top;d0 [mm];z0 [mm]"%(lowMass, highMass),40,-4.0,4.0,100,-4.0,4.0)
plots["z0_v_d0_ele_bot_hh"] = r.TH2F("%d_to_%d_MeV_ele_z0_v_d0_bot_hh"%(lowMass,highMass),"%d_%d_MeV_ele_z0_v_d0_bot;d0 [mm];z0 [mm]"%(lowMass, highMass),40,-4.0,4.0,100,-4.0,4.0)
plots["z0_v_d0_pos_top_hh"] = r.TH2F("%d_to_%d_MeV_pos_z0_v_d0_top_hh"%(lowMass,highMass),"%d_%d_MeV_pos_z0_v_d0_top;d0 [mm];z0 [mm]"%(lowMass, highMass),40,-4.0,4.0,100,-4.0,4.0)
plots["z0_v_d0_pos_bot_hh"] = r.TH2F("%d_to_%d_MeV_pos_z0_v_d0_bot_hh"%(lowMass,highMass),"%d_%d_MeV_pos_z0_v_d0_bot;d0 [mm];z0 [mm]"%(lowMass, highMass),40,-4.0,4.0,100,-4.0,4.0)

#z0 sum
plots["z0sum_v_reconz_ele_top_hh"] = r.TH2F("%d_to_%d_MeV_z0sum_v_reconz_ele_top_hh"%(lowMass,highMass),"%d_%d_MeV_z0sum_v_reconz_ele_top;recon_z [mm];e^{+}z0+e^{-} z0 [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)
plots["z0sum_v_reconz_ele_bot_hh"] = r.TH2F("%d_to_%d_MeV_z0sum_v_reconz_ele_bot_hh"%(lowMass,highMass),"%d_%d_MeV_z0sum_v_reconz_ele_bot;recon_z [mm];e^{+}z0+e^{-} z0[mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)

#z0 diff
plots["z0dif_v_reconz_ele_top_hh"] = r.TH2F("%d_to_%d_MeV_z0dif_v_reconz_ele_top_hh"%(lowMass,highMass),"%d_%d_MeV_z0dif_v_reconz_ele_top;recon_z [mm];e^{+}z0-e^{-} z0 [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)
plots["z0dif_v_reconz_ele_bot_hh"] = r.TH2F("%d_to_%d_MeV_z0dif_v_reconz_ele_bot_hh"%(lowMass,highMass),"%d_%d_MeV_z0dif_v_reconz_ele_bot;recon_z [mm];e^{+}z0-e^{-} z0 [mm]"%(lowMass, highMass),100,-30.0,80.0,100,-4.0,4.0)

#for index, row in df.iterrows():
        #df[df.columns[0]].count()
for index, row in df.iterrows():
    #Top ele
    if row['unc_vtx_ele_track_tanLambda'] > 0.0:
        plots["z0_v_reconz_ele_top_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'])
        #z0sum
        plots["z0sum_v_reconz_ele_top_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0']+row['unc_vtx_pos_track_z0'])
        #z0dif
        plots["z0dif_v_reconz_ele_top_hh"].Fill(row['unc_vtx_z'],(row['unc_vtx_ele_track_z0']-row['unc_vtx_pos_track_z0']))

        #Recon_Z Cut
        for zcut in zcuts:
            if row['unc_vtx_z'] > zcut: #recon_z slice
                #d0
                #plots["z0_v_d0_ele_top_hh"].Fill(row['unc_vtx_ele_track_d0'],row['unc_vtx_ele_track_z0'])
                #z0 v z0
                plots["top_ele_z0_v_pos_z0_zcut_%s_hh"%(zcut)].Fill(row['unc_vtx_pos_track_z0'],row['unc_vtx_ele_track_z0'])

    #Bot ele
    else:
        plots["z0_v_reconz_ele_bot_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'])
        #z0sum
        plots["z0sum_v_reconz_ele_bot_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0']+row['unc_vtx_pos_track_z0'])
        #z0dif
        plots["z0dif_v_reconz_ele_bot_hh"].Fill(row['unc_vtx_z'],(row['unc_vtx_ele_track_z0']-row['unc_vtx_pos_track_z0']))

        #Recon_Z cut
        for zcut in zcuts:
            if row['unc_vtx_z'] > zcut: #recon_z slice
                #d0
                #plots["z0_v_d0_ele_bot_hh"].Fill(row['unc_vtx_ele_track_d0'],row['unc_vtx_ele_track_z0'])
                #z0 v z0
                plots["bot_ele_z0_v_pos_z0_zcut_%s_hh"%(zcut)].Fill(row['unc_vtx_pos_track_z0'],row['unc_vtx_ele_track_z0'])

    #Top Pos
    if row['unc_vtx_pos_track_tanLambda'] > 0.0:
        plots["z0_v_reconz_pos_top_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'])
        #d0
        if row['unc_vtx_z'] > 8.0: #recon_z slice
            plots["z0_v_d0_pos_top_hh"].Fill(row['unc_vtx_pos_track_d0'],row['unc_vtx_pos_track_z0'])
    #Bot Pos
    else:
        plots["z0_v_reconz_pos_bot_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'])
        #d0
        if row['unc_vtx_z'] > 8.0: #recon_z slice
            plots["z0_v_d0_pos_bot_hh"].Fill(row['unc_vtx_pos_track_d0'],row['unc_vtx_pos_track_z0'])


outfile.cd()
for plot in plots.values():
    plot.Write()

