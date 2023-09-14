#!/usr/bin/python3
import sys
sys.path.append( '/sdf/group/hps/users/alspellm/projects/THESIS/ana/plot_utils')
import my_plot_utils as utils
import ROOT as r
import numpy as np
import root_numpy as rnp
import pandas as pd
import os
import re
import glob as glob

r.gROOT.SetBatch(1)
colors = utils.getColors()
outdir = '/sdf/group/hps/users/alspellm/projects/THESIS/ana/tight_selection_studies/z0_tanlambda_investigation_20230815/output_files'
outfile = r.TFile('%s/signal_z0tanlambda_investigation_radMatchTight.root'%(outdir),"RECREATE")

directory_path = '/sdf/group/hps/users/alspellm/projects/THESIS/ana/tight_selection_studies/z0_tanlambda_investigation_20230815/simp_sig'

##### OLD AP SAMPLE ########
ap_file = '/sdf/group/hps/users/alspellm/projects/THESIS/ana/tight_selection_studies/z0_tanlambda_investigation_20230815/ap_samp/hadd_mass_75_apsignalv2-beamv6_2500kBunches_displaced_10mm_KF_rereco_target_tracks_24ns_ana.root'
ap_hh = r.TH2F("ap_75_MeV_recon_z_v_z0tanlambda_hh","ap_75_MeV_recon_z_v_z0tanlambda_hh;z0/tan#lambda;recon_z [mm];",150,-400,200,160, -40,120)

tree = 'vtxana_kf_radMatchTight_2016_simp_reach_CR/vtxana_kf_radMatchTight_2016_simp_reach_CR_tree'
#tree = 'vtxana_kf_Tight_2016_simp_reach_CR/vtxana_kf_Tight_2016_simp_reach_CR_tree'
#tree = 'vtxana_kf_Tight_loose/vtxana_kf_Tight_loose_tree'
#tree = 'vtxana_kf_Tight_2016_simp_reach_SR/vtxana_kf_Tight_2016_simp_reach_SR_tree'
#Z0/tanlambda
z0_tanlambda = lambda z0, tanlambda : z0/tanlambda
arr = rnp.root2array(ap_file, tree)
df = pd.DataFrame(arr)
for index, row in df.iterrows():
    ap_hh.Fill(z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda']),row['unc_vtx_z'])
    ap_hh.Fill(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),row['unc_vtx_z'])

# Define a linear fitting function
fit_func = r.TF1("fit_func", "[0] + [1]*x", -200.0, 0.0)
fit_func.SetParameters(1, 1)  # Initial parameter values
fit_result = ap_hh.Fit(fit_func, "QS","", -150.0,0.0)
# Get the fit parameters and their errors
fit_params = fit_result.GetParams()
fit_param_errors = fit_result.GetErrors()
print("AP Fit Parameters:")
print("Intercept:", fit_params[0], "±", fit_param_errors[0])
print("Slope:", fit_params[1], "±", fit_param_errors[1])
ap_hh.Draw("colz")
fit_func.Draw()
outfile.cd()
ap_hh.Write()

############ SIMPS ##############
# Define the pattern to match (e.g., all .txt files)
file_pattern = "hadd_simp*.root"

# Create a list of matching files
matching_files = sorted(glob.glob(os.path.join(directory_path, file_pattern)))

# Loop over each matching file
for i_sig,file_path in enumerate(matching_files):
    print("Matching file:", file_path)

    # Get the basename of the file
    file_basename = os.path.basename(file_path)

    # Define the regular expression pattern
    pattern = r"simp_(\d+)_MeV"
    # Use re.search to find the pattern in the input string
    match = re.search(pattern, file_basename)

    #Set mass based on file name
    mV = 0
    if match:
            mV = int(match.group(1))
    print("mV = ", mV)
    #if mV != 75 or mV != 45:
    #    continue
    print("MAKING PLOTS")
    #mass window
    massRes = 3.0 
    lowMass = float(mV) - 2.8*massRes/2.0
    highMass = float(mV) + 2.8*massRes/2.0
    print("Mass Window: ", lowMass, "--", highMass)

    #DEFINE PLOTS
    plots = {}
    plots["signal_recon_z_v_z0tanlambda_hh"] = r.TH2F("signal_%d_to_%d_MeV_recon_z_v_z0tanlambda_hh"%(lowMass,highMass),"signal_%d_%d_MeV_recon_z_v_z0tanlambda_hh;z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)
    plots["background_recon_z_v_z0tanlambda_hh"] = r.TH2F("background_%d_to_%d_MeV_recon_z_v_z0tanlambda_hh"%(lowMass,highMass),"background_%d_%d_MeV_recon_z_v_z0tanlambda_hh;z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)
    plots["signal_recon_z_v_dz0tanlambda_hh"] = r.TH2F("signal_%d_to_%d_MeV_recon_z_v_dz0tanlambda_hh"%(lowMass,highMass),"signal_%d_%d_MeV_recon_z_v_dz0tanlambda_hh;#Delta z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)
    plots["background_recon_z_v_dz0tanlambda_hh"] = r.TH2F("background_%d_to_%d_MeV_recon_z_v_dz0tanlambda_hh"%(lowMass,highMass),"background_%d_%d_MeV_recon_z_v_dz0tanlambda_hh;#Delta z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)

    plots["signal_recon_z_v_ABSdz0tanlambda_hh"] = r.TH2F("signal_%d_to_%d_MeV_recon_z_v_ABSdz0tanlambda_hh"%(lowMass,highMass),"signal_%d_%d_MeV_recon_z_v_ABSdz0tanlambda_hh;Abs(#Delta z0/tan#lambda);recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)
    plots["background_recon_z_v_ABSdz0tanlambda_hh"] = r.TH2F("background_%d_to_%d_MeV_recon_z_v_ABSdz0tanlambda_hh"%(lowMass,highMass),"background_%d_%d_MeV_recon_z_v_ABSdz0tanlambda_hh;Abs(#Delta z0/tan#lambda);recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)

    #Top/Bot Ele
    plots["signal_recon_z_v_z0tanlambda_topEle_hh"] = r.TH2F("signal_%d_to_%d_MeV_recon_z_v_z0tanlambda_topEle_hh"%(lowMass,highMass),"signal_%d_%d_MeV_recon_z_v_z0tanlambda_topEle_hh;z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)
    plots["signal_recon_z_v_z0tanlambda_botEle_hh"] = r.TH2F("signal_%d_to_%d_MeV_recon_z_v_z0tanlambda_botEle_hh"%(lowMass,highMass),"signal_%d_%d_MeV_recon_z_v_z0tanlambda_botEle_hh;z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)

    plots["background_recon_z_v_z0tanlambda_topEle_hh"] = r.TH2F("background_%d_to_%d_MeV_recon_z_v_z0tanlambda_topEle_hh"%(lowMass,highMass),"background_%d_%d_MeV_recon_z_v_z0tanlambda_topEle_hh;z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)
    plots["background_recon_z_v_z0tanlambda_botEle_hh"] = r.TH2F("background_%d_to_%d_MeV_recon_z_v_z0tanlambda_botEle_hh"%(lowMass,highMass),"background_%d_%d_MeV_recon_z_v_z0tanlambda_botEle_hh;z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)


    plots["signal_recon_z_v_dz0tanlambda_topEle_hh"] = r.TH2F("signal_%d_to_%d_MeV_recon_z_v_dz0tanlambda_topEle_hh"%(lowMass,highMass),"signal_%d_%d_MeV_recon_z_v_dz0tanlambda_topEle_hh;#Delta z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)
    plots["signal_recon_z_v_dz0tanlambda_botEle_hh"] = r.TH2F("signal_%d_to_%d_MeV_recon_z_v_dz0tanlambda_botEle_hh"%(lowMass,highMass),"signal_%d_%d_MeV_recon_z_v_dz0tanlambda_botEle_hh;#Delta z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)

    plots["background_recon_z_v_dz0tanlambda_topEle_hh"] = r.TH2F("background_%d_to_%d_MeV_recon_z_v_dz0tanlambda_topEle_hh"%(lowMass,highMass),"background_%d_%d_MeV_recon_z_v_dz0tanlambda_topEle_hh;#Delta z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)
    plots["background_recon_z_v_dz0tanlambda_botEle_hh"] = r.TH2F("background_%d_to_%d_MeV_recon_z_v_dz0tanlambda_botEle_hh"%(lowMass,highMass),"background_%d_%d_MeV_recon_z_v_dz0tanlambda_botEle_hh;#Delta z0/tan#lambda;recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)

    plots["signal_recon_z_v_ABSdz0tanlambda_topEle_hh"] = r.TH2F("signal_%d_to_%d_MeV_recon_z_v_ABSdz0tanlambda_topEle_hh"%(lowMass,highMass),"signal_%d_%d_MeV_recon_z_v_ABSdz0tanlambda_topEle_hh;Abs(#Delta z0/tan#lambda);recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)
    plots["signal_recon_z_v_ABSdz0tanlambda_botEle_hh"] = r.TH2F("signal_%d_to_%d_MeV_recon_z_v_ABSdz0tanlambda_botEle_hh"%(lowMass,highMass),"signal_%d_%d_MeV_recon_z_v_ABSdz0tanlambda_botEle_hh;Abs(#Delta z0/tan#lambda);recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)
    plots["background_recon_z_v_ABSdz0tanlambda_topEle_hh"] = r.TH2F("background_%d_to_%d_MeV_recon_z_v_ABSdz0tanlambda_topEle_hh"%(lowMass,highMass),"background_%d_%d_MeV_recon_z_v_ABSdz0tanlambda_topEle_hh;Abs(#Delta z0/tan#lambda);recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)
    plots["background_recon_z_v_ABSdz0tanlambda_botEle_hh"] = r.TH2F("background_%d_to_%d_MeV_recon_z_v_ABSdz0tanlambda_botEle_hh"%(lowMass,highMass),"background_%d_%d_MeV_recon_z_v_ABSdz0tanlambda_botEle_hh;Abs(#Delta z0/tan#lambda);recon_z [mm];"%(lowMass, highMass),150,-400,200,160, -40,120)


    #define mass window selection
    selection = {'unc_vtx_mass_lt' : highMass/1000., 'unc_vtx_mass_gt' : lowMass/1000.}
    #tree = 'vtxana_kf_Tight_loose/vtxana_kf_Tight_loose_tree'
    #tree = 'vtxana_kf_Tight_2016_simp_reach_SR/vtxana_kf_Tight_2016_simp_reach_SR_tree'
    tree = 'vtxana_kf_radMatchTight_2016_simp_reach_SR/vtxana_kf_radMatchTight_2016_simp_reach_SR_tree'
    #tree = 'vtxana_kf_Tight_2016_simp_reach_SR/vtxana_kf_Tight_2016_simp_reach_SR_tree'
    #Z0/tanlambda
    z0_tanlambda = lambda z0, tanlambda : z0/tanlambda
    dz0_tanlambda = lambda ele, pos : pos-ele
    abs_dz0_tanlambda = lambda ele, pos : abs(ele-pos)

    ####### Signal #########
    arr = rnp.root2array(file_path, tree)
    df = pd.DataFrame(arr)
    df = df[(df['unc_vtx_mass'] < selection['unc_vtx_mass_lt']) & (df['unc_vtx_mass'] > selection['unc_vtx_mass_gt']) ]

    for index, row in df.iterrows():
        plots["signal_recon_z_v_z0tanlambda_hh"].Fill(z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda']),row['unc_vtx_z'])
        plots["signal_recon_z_v_z0tanlambda_hh"].Fill(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),row['unc_vtx_z'])
        plots["signal_recon_z_v_dz0tanlambda_hh"].Fill(dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])
        plots["signal_recon_z_v_ABSdz0tanlambda_hh"].Fill(abs_dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])
        if row['unc_vtx_ele_track_tanLambda'] > 0.0:
            plots["signal_recon_z_v_z0tanlambda_topEle_hh"].Fill(z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda']),row['unc_vtx_z'])
            plots["signal_recon_z_v_z0tanlambda_topEle_hh"].Fill(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),row['unc_vtx_z'])
            plots["signal_recon_z_v_dz0tanlambda_topEle_hh"].Fill(dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])
            plots["signal_recon_z_v_ABSdz0tanlambda_topEle_hh"].Fill(abs_dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])
        else:
            plots["signal_recon_z_v_z0tanlambda_botEle_hh"].Fill(z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda']),row['unc_vtx_z'])
            plots["signal_recon_z_v_z0tanlambda_botEle_hh"].Fill(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),row['unc_vtx_z'])
            plots["signal_recon_z_v_dz0tanlambda_botEle_hh"].Fill(dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])
            plots["signal_recon_z_v_ABSdz0tanlambda_botEle_hh"].Fill(abs_dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])



    # Define a linear fitting function
    fit_func = r.TF1("fit_func", "[0] + [1]*x", -400.0, 100.0)
    fit_func.SetParameters(1, 1)  # Initial parameter values
    fit_result = plots["signal_recon_z_v_z0tanlambda_hh"].Fit(fit_func, "QS","", -150.0,0.0)
    # Get the fit parameters and their errors
    fit_params = fit_result.GetParams()
    fit_param_errors = fit_result.GetErrors()
    print("Signal Fit Parameters:")
    print("Intercept:", fit_params[0], "±", fit_param_errors[0])
    print("Slope:", fit_params[1], "±", fit_param_errors[1])
    plots["signal_recon_z_v_z0tanlambda_hh"].Draw("colz")
    fit_func.Draw()

    ###### Background ########
    background_path = "/sdf/group/hps/users/alspellm/projects/THESIS/data/2016/hps_7800/20230726_recon/20230726_tuples/20230727_ana/ana_files/hadd_hps_7800_KF_24nsClusterWindow_ana_367_files.root"
    arr = rnp.root2array(background_path, tree)
    df = pd.DataFrame(arr)
    df = df[(df['unc_vtx_mass'] < selection['unc_vtx_mass_lt']) & (df['unc_vtx_mass'] > selection['unc_vtx_mass_gt']) ]

    for index, row in df.iterrows():
        plots["background_recon_z_v_z0tanlambda_hh"].Fill(z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda']),row['unc_vtx_z'])
        plots["background_recon_z_v_z0tanlambda_hh"].Fill(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),row['unc_vtx_z'])
        plots["background_recon_z_v_dz0tanlambda_hh"].Fill(dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])
        plots["background_recon_z_v_ABSdz0tanlambda_hh"].Fill(abs_dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])
        if row['unc_vtx_ele_track_tanLambda'] > 0.0:
            plots["background_recon_z_v_z0tanlambda_topEle_hh"].Fill(z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda']),row['unc_vtx_z'])
            plots["background_recon_z_v_z0tanlambda_topEle_hh"].Fill(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),row['unc_vtx_z'])
            plots["background_recon_z_v_dz0tanlambda_topEle_hh"].Fill(dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])
            plots["background_recon_z_v_ABSdz0tanlambda_topEle_hh"].Fill(abs_dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])
        else:
            plots["background_recon_z_v_z0tanlambda_botEle_hh"].Fill(z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda']),row['unc_vtx_z'])
            plots["background_recon_z_v_z0tanlambda_botEle_hh"].Fill(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),row['unc_vtx_z'])
            plots["background_recon_z_v_dz0tanlambda_botEle_hh"].Fill(dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])
            plots["background_recon_z_v_ABSdz0tanlambda_botEle_hh"].Fill(abs_dz0_tanlambda(z0_tanlambda(row['unc_vtx_pos_track_z0'],row['unc_vtx_pos_track_tanLambda']),z0_tanlambda(row['unc_vtx_ele_track_z0'],row['unc_vtx_ele_track_tanLambda'])),row['unc_vtx_z'])

    # Define a linear fitting function
    fit_func = r.TF1("fit_func", "[0] + [1]*x", -200.0, 0.0)
    fit_func.SetParameters(1, 1)  # Initial parameter values
    fit_result = plots["background_recon_z_v_z0tanlambda_hh"].Fit(fit_func, "QS","", -150.0,0.0)
    fit_params = fit_result.GetParams()
    fit_param_errors = fit_result.GetErrors()
    print("Background Fit Parameters:")
    print("Intercept:", fit_params[0], "±", fit_param_errors[0])
    print("Slope:", fit_params[1], "±", fit_param_errors[1])
    plots["background_recon_z_v_z0tanlambda_hh"].Draw("colz")
    fit_func.Draw("same")

    '''
    fit_func_2 = r.TF1("fit_func_2", "[0] + [1]*x", 0.0, 200.0)
    fit_func_2.SetParameters(1, 1)  # Initial parameter values
    fit_result_2 = plots["background_recon_z_v_z0tanlambda_hh"].Fit(fit_func_2, "QS+","", 0.0,200.0)
    fit_result_2.Draw("same")

    fit_func_3 = r.TF1("fit_func_3", "[0] + [1]*x", -30.0, 30.0)
    fit_func_3.SetParameters(1, 1)  # Initial parameter values
    fit_result_3 = plots["background_recon_z_v_z0tanlambda_hh"].Fit(fit_func_3, "QS+","", -30.0,30.0)
    fit_result_3.Draw("same")
    '''
    outfile.cd()
    for plot in plots.values():
        plot.Write()

