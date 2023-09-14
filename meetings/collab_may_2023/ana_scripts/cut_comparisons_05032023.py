#!/usr/bin/python3
import ROOT as r
import numpy as np
import root_numpy as rnp
import copy
import pandas as pd

def zbravo_ele(row, z0_correction):
    z0 = row["unc_vtx_ele_track_z0"] + z0_correction
    if row["unc_vtx_ele_track_tanLambda"] > 0.0:
        return z0 - (-0.039151*row["unc_vtx_z"] - 0.031282)
    else:
        return z0 - (0.040086*row["unc_vtx_z"] + 0.016186)

def zbravo_pos(row, z0_correction):
    z0 = row["unc_vtx_pos_track_z0"] + z0_correction
    if row["unc_vtx_pos_track_tanLambda"] > 0.0:
        return z0 - (-0.037899*row["unc_vtx_z"] - 0.0094)
    else:
        return z0 - (0.039501*row["unc_vtx_z"] + 0.004176)

def zbravoAlpha_ele(row, z0_correction,slope):
    zbravo = zbravo_ele(row, z0_correction) 
    if zbravo > 0.0:
        return row["unc_vtx_z"] - ((zbravo)/slope)
    else:
        return row["unc_vtx_z"] - ((zbravo)/-slope)

def zbravoAlpha_pos(row, z0_correction,slope):
    zbravo = zbravo_pos(row, z0_correction) 
    if zbravo > 0.0:
        return row["unc_vtx_z"] - ((zbravo)/slope)
    else:
        return row["unc_vtx_z"] - ((zbravo)/-slope)

def pass_zbravoalpha_gt_cut(row, z0_correction, slope, ele_cut, pos_cut):
    zbravoalpha_ele = zbravoAlpha_ele(row, z0_correction, slope)
    zbravoalpha_pos = zbravoAlpha_pos(row, z0_correction,slope)
    passCut = True
    if zbravoalpha_ele <= ele_cut:
        passCut = False
    if zbravoalpha_pos <= pos_cut:
        passCut = False
    return passCut

def pass_zbravoalpha_lt_cut(row, z0_correction, slope, ele_cut, pos_cut):
    zbravoalpha_ele = zbravoAlpha_ele(row, z0_correction, slope)
    zbravoalpha_pos = zbravoAlpha_pos(row, z0_correction,slope)
    passCut = True
    if zbravoalpha_ele >= ele_cut:
        passCut = False
    if zbravoalpha_pos >= pos_cut:
        passCut = False
    return passCut

def pass_impact_parameter_cut(row, z0_correction): #55 MeV SIMP Signal
    ele_z0 = row["unc_vtx_ele_track_z0"] + z0_correction
    pos_z0 = row["unc_vtx_pos_track_z0"] + z0_correction
    passCut = True
    if ele_z0 > 0.0:
        if ele_z0 < 0.029816*(row["unc_vtx_z"]-3.471875):
            passCut = False
    else:
        if ele_z0 > -0.029530*(row["unc_vtx_z"]-3.471875):
            passCut = False

    if pos_z0 > 0.0:
        if pos_z0 < 0.029816*(row["unc_vtx_z"]-3.471875):
            passCut = False
    else:
        if pos_z0 > -0.029530*(row["unc_vtx_z"]-3.471875):
            passCut = False
    return passCut

def massRes(mass):
    res = 9.73217e-01 + 3.63659e-02*mass + -7.32046e-05*mass*mass #2016 simps alic
    return res
        

#mass window
mV = 55.0
#massRes = 3.0 #MeV for a 55MeV Vd
massRes = massRes(mV)
lowMass = float(mV) - 2.0*massRes/2.0
highMass = float(mV) + 2.0*massRes/2.0
print("Mass Window: ", lowMass, "--", highMass)

#define mass window selection
selection = {'unc_vtx_mass_lt' : highMass/1000., 'unc_vtx_mass_gt' : lowMass/1000.}

#data
data_infile = '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/kf_data/kf_041823/output/hadd_hps_BLPass4_1958_files_recon_4.2_ana_kf.root'
data_treename = 'vtxana_kf_Tight_2016_simp_reach_dev/vtxana_kf_Tight_2016_simp_reach_dev_tree'
data_arr = rnp.root2array(data_infile, data_treename)
df = pd.DataFrame(data_arr)
df = df[(df['unc_vtx_mass'] < selection['unc_vtx_mass_lt']) & (df['unc_vtx_mass'] > selection['unc_vtx_mass_gt']) ]

Lumi = 10.7*0.1 #Data represents ~10% lumi
mcScale = {}
mcScale['tritrig'] = 1.416e9*Lumi/(50000*9853) #pb2016
mcScale['wab'] = 0.1985e12*Lumi/(100000*9966) #pb2016

tritrig_infile = '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/kf_mc/tritrig_beam/tritrig_kf_ana_05122023/final_hadd_tritrigv2-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_KF_ana.root'
tritrig_treename = 'vtxana_kf_Tight_2016_simp_reach_dev/vtxana_kf_Tight_2016_simp_reach_dev_tree'
tritrig_arr = rnp.root2array(tritrig_infile, tritrig_treename)
tdf = pd.DataFrame(tritrig_arr)
tdf = tdf[(tdf['unc_vtx_mass'] < selection['unc_vtx_mass_lt']) & (tdf['unc_vtx_mass'] > selection['unc_vtx_mass_gt']) ]

wab_infile = '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/kf_mc/wab_beam/wab_kf_ana_05022023/output/final_hadd_wabv3-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_KF_ana.root'
wab_treename = 'vtxana_kf_Tight_2016_simp_reach_dev/vtxana_kf_Tight_2016_simp_reach_dev_tree'
wab_arr = rnp.root2array(wab_infile, wab_treename)
wdf = pd.DataFrame(wab_arr)
wdf = wdf[(wdf['unc_vtx_mass'] < selection['unc_vtx_mass_lt']) & (wdf['unc_vtx_mass'] > selection['unc_vtx_mass_gt']) ]

#mc signal
signal_infile = '/sdf/group/hps/users/alspellm/projects/THESIS/cut_dev/ZBi/signal/hadd_mass_55_simp_recon_KF_ana.root'
signal_treename = 'vtxana_kf_radMatchTight_2016_simp_reach_dev/vtxana_kf_radMatchTight_2016_simp_reach_dev_tree'
signal_arr = rnp.root2array(signal_infile, signal_treename)
sigdf = pd.DataFrame(signal_arr)
sigdf = sigdf[(sigdf['unc_vtx_mass'] < selection['unc_vtx_mass_lt']) & (sigdf['unc_vtx_mass'] > selection['unc_vtx_mass_gt'])]


outfile = r.TFile('zbravo_impact_param_comparisons_05072023.root',"RECREATE")
####################################################################

dataPlots2D ={}
dataPlots1D = {}
#data z0 vs recon_z
dataPlots2D["data_z0_v_reconz_hh"] = r.TH2F("data_z0_v_reconz_hh","data_z0_v_reconz;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
dataPlots2D["data_z0_v_reconz_impact_param_cut_hh"] = r.TH2F("data_z0_v_reconz_impact_param_cut_hh","data_z0_v_reconz_impact_param_cut;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
dataPlots2D["data_z0_v_reconz_zbravoalpha_cut_hh"] = r.TH2F("data_z0_v_reconz_zbravoalpha_cut_hh","data_z0_v_reconz_zbravoalpha_cut;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
#dataPlots2D["data_z0_v_reconz_ele_top_hh"] = r.TH2F("data_ele_z0_v_reconz_top_hh","data_ele_z0_v_reconz_top;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
#dataPlots2D["data_z0_v_reconz_ele_bot_hh"] = r.TH2F("data_ele_z0_v_reconz_bot_hh","data_ele_z0_v_reconz_bot;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)

#data unc_vtx_z
dataPlots1D["data_reconz_h"] = r.TH1F("data_reconz_h","data_reconz;recon_z [mm];events",1100,-30.0,80.0)
dataPlots1D["data_reconz_impact_param_cut_h"] = r.TH1F("data_reconz_impact_param_cut_h","data_reconz_impact_param_cut;recon_z [mm];events",1100,-30.0,80.0)
dataPlots1D["data_reconz_zbravoalpha_cut_h"] = r.TH1F("data_reconz_zbravoalpha_cut_h","data_reconz_zbravoalpha_cut;recon_z [mm];events",1100,-30.0,80.0)
#dataPlots1D["data_reconz_ele_top_h"] = r.TH1F("data_reconz_ele_top_h","data_reconz;recon_z [mm];events",1100,-30.0,80.0)
#dataPlots1D["data_reconz_ele_bot_h"] = r.TH1F("data_reconz_ele_bot_h","data_reconz;recon_z [mm];events",1100,-30.0,80.0)
#data zvtx_v_invmass
dataPlots2D["data_reconz_invMass_hh"] = r.TH2F("data_reconz_v_invMass_hh","data_reconz_v_invMass;invMass [GeV];recon_z [mm]",200,0.0,200.0,1100,-30.0,80.0)
dataPlots2D["data_reconz_invMass_impact_param_cut_hh"] = r.TH2F("data_reconz_v_invMass_impact_param_cut_hh","data_reconz_v_invMass_impact_param_cut;invMass [GeV];recon_z [mm]",200,0.0,200.0,1100,-30.0,80.0)
dataPlots2D["data_reconz_invMass_zbravoalpha_cut_hh"] = r.TH2F("data_reconz_v_invMass_zbravoalpha_cut_hh","data_reconz_v_invMass_zbravoalpha_cut;invMass [GeV];recon_z [mm]",200,0.0,200.0,1100,-30.0,80.0)

###################################################################

mcPlots2D ={}
mcPlots1D = {}
#z0 vs recon_z
mcPlots2D["mc_z0_v_reconz_hh"] = r.TH2F("mc_z0_v_reconz_hh","mc_z0_v_reconz;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
mcPlots2D["mc_z0_v_reconz_impact_param_cut_hh"] = r.TH2F("mc_z0_v_reconz_impact_param_cut_hh","mc_z0_v_reconz_impact_param_cut;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
mcPlots2D["mc_z0_v_reconz_zbravoalpha_cut_hh"] = r.TH2F("mc_z0_v_reconz_zbravoalpha_cut_hh","mc_z0_v_reconz_zbravoalpha_cut;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
#mcPlots2D["mc_z0_v_reconz_ele_top_hh"] = r.TH2F("mc_ele_z0_v_reconz_top_hh","mc_ele_z0_v_reconz_top;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
#mcPlots2D["mc_z0_v_reconz_ele_bot_hh"] = r.TH2F("mc_ele_z0_v_reconz_bot_hh","mc_ele_z0_v_reconz_bot;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)

#mc unc_vtx_z
mcPlots1D["mc_reconz_h"] = r.TH1F("mc_reconz_h","mc_reconz;recon_z [mm];events",1100,-30.0,80.0)
mcPlots1D["mc_reconz_impact_param_cut_h"] = r.TH1F("mc_reconz_impact_param_cut_h","mc_reconz_impact_param_cut;recon_z [mm];events",1100,-30.0,80.0)
mcPlots1D["mc_reconz_zbravoalpha_cut_h"] = r.TH1F("mc_reconz_zbravoalpha_cut_h","mc_reconz_zbravoalpha_cut;recon_z [mm];events",1100,-30.0,80.0)
#mcPlots1D["mc_reconz_ele_top_h"] = r.TH1F("mc_reconz_ele_top_h","mc_reconz;recon_z [mm];events",1100,-30.0,80.0)
#mcPlots1D["mc_reconz_ele_bot_h"] = r.TH1F("mc_reconz_ele_bot_h","mc_reconz;recon_z [mm];events",1100,-30.0,80.0)
#mc zvtx_v_invmass
mcPlots2D["mc_reconz_invMass_hh"] = r.TH2F("mc_reconz_v_invMass_hh","mc_reconz_v_invMass;invMass [GeV];recon_z [mm]",200,0.0,200.0,1100,-30.0,80.0)
mcPlots2D["mc_reconz_invMass_impact_param_cut_hh"] = r.TH2F("mc_reconz_v_invMass_impact_param_cut_hh","mc_reconz_v_invMass_impact_param_cut;invMass [GeV];recon_z [mm]",200,0.0,200.0,1100,-30.0,80.0)
mcPlots2D["mc_reconz_invMass_zbravoalpha_cut_hh"] = r.TH2F("mc_reconz_v_invMass_zbravoalpha_cut_hh","mc_reconz_v_invMass_zbravoalpha_cut;invMass [GeV];recon_z [mm]",200,0.0,200.0,1100,-30.0,80.0)

#####################################################################
signalPlots2D ={}
signalPlots1D = {}
#z0 vs recon_z
signalPlots2D["signal_z0_v_reconz_hh"] = r.TH2F("signal_z0_v_reconz_hh","signal_z0_v_reconz;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
signalPlots2D["signal_z0_v_reconz_impact_param_cut_hh"] = r.TH2F("signal_z0_v_reconz_impact_param_cut_hh","signal_z0_v_reconz_impact_param_cut;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
signalPlots2D["signal_z0_v_reconz_zbravoalpha_cut_hh"] = r.TH2F("signal_z0_v_reconz_zbravoalpha_cut_hh","signal_z0_v_reconz_zbravoalpha_cut;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
#signalPlots2D["signal_z0_v_reconz_ele_top_hh"] = r.TH2F("signal_ele_z0_v_reconz_top_hh","signal_ele_z0_v_reconz_top;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)
#signalPlots2D["signal_z0_v_reconz_ele_bot_hh"] = r.TH2F("signal_ele_z0_v_reconz_bot_hh","signal_ele_z0_v_reconz_bot;recon_z [mm];z0 [mm]",1100,-30.0,80.0,800,-4.0,4.0)

#signal unc_vtx_z
signalPlots1D["signal_reconz_h"] = r.TH1F("signal_reconz_h","signal_reconz;recon_z [mm];events",1100,-30.0,80.0)
signalPlots1D["signal_reconz_impact_param_cut_h"] = r.TH1F("signal_reconz_impact_param_cut_h","signal_reconz_impact_param_cut;recon_z [mm];events",1100,-30.0,80.0)
signalPlots1D["signal_reconz_zbravoalpha_cut_h"] = r.TH1F("signal_reconz_zbravoalpha_cut_h","signal_reconz_zbravoalpha_cut;recon_z [mm];events",1100,-30.0,80.0)
#signalPlots1D["signal_reconz_ele_top_h"] = r.TH1F("signal_reconz_ele_top_h","signal_reconz;recon_z [mm];events",1100,-30.0,80.0)
#signalPlots1D["signal_reconz_ele_bot_h"] = r.TH1F("signal_reconz_ele_bot_h","signal_reconz;recon_z [mm];events",1100,-30.0,80.0)
#signal zvtx_v_invmass
signalPlots2D["signal_reconz_invMass_hh"] = r.TH2F("signal_reconz_v_invMass_hh","signal_reconz_v_invMass;invMass [GeV];recon_z [mm]",200,0.0,200.0,1100,-30.0,80.0)
signalPlots2D["signal_reconz_invMass_impact_param_cut_hh"] = r.TH2F("signal_reconz_v_invMass_impact_param_cut_hh","signal_reconz_v_invMass_impact_param_cut;invMass [GeV];recon_z [mm]",200,0.0,200.0,1100,-30.0,80.0)
signalPlots2D["signal_reconz_invMass_zbravoalpha_cut_hh"] = r.TH2F("signal_reconz_v_invMass_zbravoalpha_cut_hh","signal_reconz_v_invMass_zbravoalpha_cut;invMass [GeV];recon_z [mm]",200,0.0,200.0,1100,-30.0,80.0)
######################################################################
z0_correction = 0.1
print("Length of data: ", len(df))
for index, row in df.iterrows():
    if index%100000 == 0:
        print("data event ",index)
    
    #All
    dataPlots2D["data_z0_v_reconz_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'])
    dataPlots2D["data_z0_v_reconz_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'])
    dataPlots1D["data_reconz_h"].Fill(row['unc_vtx_z'])
    dataPlots2D["data_reconz_invMass_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'])

    #impact param cut
    if(pass_impact_parameter_cut(row, z0_correction)):
        dataPlots2D["data_z0_v_reconz_impact_param_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'])
        dataPlots2D["data_z0_v_reconz_impact_param_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'])
        dataPlots1D["data_reconz_impact_param_cut_h"].Fill(row['unc_vtx_z'])
        dataPlots2D["data_reconz_invMass_impact_param_cut_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'])

    #zbravoalpha
    if(pass_zbravoalpha_lt_cut(row, z0_correction, 0.01, 46.0, 44.0) and pass_zbravoalpha_gt_cut(row, z0_correction, 0.01, -18.0, -20.0)):
        dataPlots2D["data_z0_v_reconz_zbravoalpha_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'])
        dataPlots2D["data_z0_v_reconz_zbravoalpha_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'])
        dataPlots1D["data_reconz_zbravoalpha_cut_h"].Fill(row['unc_vtx_z'])
        dataPlots2D["data_reconz_invMass_zbravoalpha_cut_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'])

print("Length of wabs: ", len(wdf))
for index, row in wdf.iterrows():
    if index%100000 == 0:
        print("wab event ",index)
    #All
    mcPlots2D["mc_z0_v_reconz_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'],mcScale['wab'])
    mcPlots2D["mc_z0_v_reconz_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'],mcScale['wab'])
    mcPlots1D["mc_reconz_h"].Fill(row['unc_vtx_z'],mcScale['wab'])
    mcPlots2D["mc_reconz_invMass_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'],mcScale['wab'])

    #impact param cut
    if(pass_impact_parameter_cut(row, z0_correction)):
        mcPlots2D["mc_z0_v_reconz_impact_param_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'],mcScale['wab'])
        mcPlots2D["mc_z0_v_reconz_impact_param_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'],mcScale['wab'])
        mcPlots1D["mc_reconz_impact_param_cut_h"].Fill(row['unc_vtx_z'],mcScale['wab'])
        mcPlots2D["mc_reconz_invMass_impact_param_cut_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'],mcScale['wab'])

    #zbravoalpha cut
    if(pass_zbravoalpha_lt_cut(row, z0_correction, 0.01, 46.0, 44.0) and pass_zbravoalpha_gt_cut(row, z0_correction, 0.01, -18.0, -20.0)):
        mcPlots2D["mc_z0_v_reconz_zbravoalpha_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'],mcScale['wab'])
        mcPlots2D["mc_z0_v_reconz_zbravoalpha_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'],mcScale['wab'])
        mcPlots1D["mc_reconz_zbravoalpha_cut_h"].Fill(row['unc_vtx_z'],mcScale['wab'])
        mcPlots2D["mc_reconz_invMass_zbravoalpha_cut_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'],mcScale['wab'])


print("Length of tritrig: ", len(tdf))
for index, row in tdf.iterrows():
    if index%10000 == 0:
        print("tritrig event ",index)
    #All
    mcPlots2D["mc_z0_v_reconz_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'],mcScale['tritrig'])
    mcPlots2D["mc_z0_v_reconz_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'],mcScale['tritrig'])
    mcPlots1D["mc_reconz_h"].Fill(row['unc_vtx_z'],mcScale['tritrig'])
    mcPlots2D["mc_reconz_invMass_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'],mcScale['tritrig'])

    #impact param cut
    if(pass_impact_parameter_cut(row, z0_correction)):
        mcPlots2D["mc_z0_v_reconz_impact_param_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'],mcScale['tritrig'])
        mcPlots2D["mc_z0_v_reconz_impact_param_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'],mcScale['tritrig'])
        mcPlots1D["mc_reconz_impact_param_cut_h"].Fill(row['unc_vtx_z'],mcScale['tritrig'])
        mcPlots2D["mc_reconz_invMass_impact_param_cut_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'],mcScale['tritrig'])

    #zbravoalpha cut
    if(pass_zbravoalpha_lt_cut(row, z0_correction, 0.01, 46.0, 44.0) and pass_zbravoalpha_gt_cut(row, z0_correction, 0.01, -18.0, -20.0)):
        mcPlots2D["mc_z0_v_reconz_zbravoalpha_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'],mcScale['tritrig'])
        mcPlots2D["mc_z0_v_reconz_zbravoalpha_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'],mcScale['tritrig'])
        mcPlots1D["mc_reconz_zbravoalpha_cut_h"].Fill(row['unc_vtx_z'],mcScale['tritrig'])
        mcPlots2D["mc_reconz_invMass_zbravoalpha_cut_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'],mcScale['tritrig'])

print("Length of signal: ", len(sigdf))
z0_correction = 0.0
for index, row in sigdf.iterrows():
    if index%10000 == 0:
        print("signal event ",index)
    #All
    signalPlots2D["signal_z0_v_reconz_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'])
    signalPlots2D["signal_z0_v_reconz_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'])
    signalPlots1D["signal_reconz_h"].Fill(row['unc_vtx_z'])
    signalPlots2D["signal_reconz_invMass_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'])

    #impact param cut
    if(pass_impact_parameter_cut(row, z0_correction)):
        signalPlots2D["signal_z0_v_reconz_impact_param_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'])
        signalPlots2D["signal_z0_v_reconz_impact_param_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'])
        signalPlots1D["signal_reconz_impact_param_cut_h"].Fill(row['unc_vtx_z'])
        signalPlots2D["signal_reconz_invMass_impact_param_cut_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'])

    #zbravoalpha
    if(pass_zbravoalpha_lt_cut(row, z0_correction, 0.01, 46.0, 44.0) and pass_zbravoalpha_gt_cut(row, z0_correction, 0.01, -18.0, -20.0)):
        signalPlots2D["signal_z0_v_reconz_zbravoalpha_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_ele_track_z0'])
        signalPlots2D["signal_z0_v_reconz_zbravoalpha_cut_hh"].Fill(row['unc_vtx_z'],row['unc_vtx_pos_track_z0'])
        signalPlots1D["signal_reconz_zbravoalpha_cut_h"].Fill(row['unc_vtx_z'])
        signalPlots2D["signal_reconz_invMass_zbravoalpha_cut_hh"].Fill(row['unc_vtx_mass']*1000.0,row['unc_vtx_z'])

outfile.cd()
for plot in dataPlots2D.values():
    plot.Write()
for plot in dataPlots1D.values():
    plot.Write()
for plot in mcPlots2D.values():
    plot.Write()
for plot in mcPlots1D.values():
    plot.Write()

for plot in signalPlots2D.values():
    plot.Write()
for plot in signalPlots1D.values():
    plot.Write()
