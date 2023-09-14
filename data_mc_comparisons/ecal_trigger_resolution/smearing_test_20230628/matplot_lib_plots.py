import ROOT
import sys
import matplotlib.pyplot as plt
import numpy as np

sys.path.append('/sdf/group/hps/users/alspellm/projects/THESIS/ana/plot_utils/')
import my_plot_utils as utils


# Enable batch mode in PyROOT
ROOT.gROOT.SetBatch(True)

data_path = '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/kf_data/kf_041823/output/hadd_hps_BLPass4_1958_files_recon_4.2_ana_kf.root'
data_dir = 'vtxana_kf_Tight_2016_simp_reach_dev'
lumi_2016 = 10.7 #pb-1
total_lumi = lumi_2016*0.1 #pb-1

old_tritrig_path = '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/kf_mc/tritrig_beam/tritrig_kf_ana_05122023/final_hadd_tritrigv2-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_KF_ana.root'
old_tritrig_dir = 'vtxana_kf_Tight_2016_simp_reach_dev'
old_tritrig_events = 50000*9853 #nevents
old_tritrig_xsec = 1.416e9 #pb

tritrig_path = '/sdf/group/hps/users/alspellm/projects/THESIS/mc/2016/tritrig_beam/20230628_ecal_trig_res/tong_tong/20230628_tuple/20230628_ana/files/hadd_tritrig-beam_ecal_trig_res_ana.root'
tritrig_dir = 'vtxana_kf_Tight_2016_simp_reach_SR'
tritrig_events = 50000*9755 #nevents
tritrig_xsec = 1.416e9 #pb


wab_path = '/sdf/group/hps/users/alspellm/projects/THESIS/mc/2016/wab_beam/ecal_trig_res/ana/hadd_wab-beam_ecal_trig_res_ana.root'
wab_dir = 'vtxana_kf_Tight_2016_simp_reach_SR'
wab_events = 100000*9769 #nevents
wab_xsec = 0.1985e12 #pb

old_wab_path = '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/kf_mc/wab_beam/wab_kf_ana_05022023/output/final_hadd_wabv3-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_KF_ana.root'
old_wab_dir = 'vtxana_kf_Tight_2016_simp_reach_dev'
old_wab_events = 100000*9966 #nevents
old_wab_xsec = 0.1985e12 #pb

#Read histograms
keywords = ["_Psum_","_Esum_","_InvM_h","_ele_p_h","pos_p_h","ele_TanLambda_h","pos_TanLambda_h","ele_Phi_h","pos_Phi_h","vtx_Z_h","_cutflow"]

#Get colors
colors = utils.getColors()

for key in keywords:
    data_h = utils.read_1d_plots_from_root_file(data_path, data_dir, key)[0]
    old_tritrig_h = utils.read_1d_plots_from_root_file(old_tritrig_path, old_tritrig_dir, key)[0]
    old_wab_h = utils.read_1d_plots_from_root_file(old_wab_path, old_wab_dir, key)[0]

    tritrig_h = utils.read_1d_plots_from_root_file(tritrig_path, tritrig_dir, key)[0]
    wab_h = utils.read_1d_plots_from_root_file(wab_path, wab_dir, key)[0]

    #rescale histograms
    old_tritrig_h.Scale(old_tritrig_xsec*total_lumi/(old_tritrig_events))
    old_wab_h.Scale(old_wab_xsec*total_lumi/(old_wab_events))
    tritrig_h.Scale(tritrig_xsec*total_lumi/(tritrig_events))
    wab_h.Scale(wab_xsec*total_lumi/(wab_events))

    #format histograms
    #data
    name = 'data'+data_h.GetName().replace('vtxana_kf','')
    xlabel = data_h.GetXaxis().GetTitle()
    utils.format_TH1(data_h, name=name, title='Data Sample',line_width=2, x_label=xlabel,line_color=colors[0])

    #old tritrig
    name = 'old_tritrig'+old_tritrig_h.GetName().replace('vtxana_kf','')
    xlabel = old_tritrig_h.GetXaxis().GetTitle()
    utils.format_TH1(old_tritrig_h, name=name, title='Tritrig Old EcalTrigRes',line_width=2, x_label=xlabel,line_color=colors[1])

    #old wab
    name = 'old_wab'+old_wab_h.GetName().replace('vtxana_kf','')
    xlabel = old_wab_h.GetXaxis().GetTitle()
    utils.format_TH1(old_wab_h, name=name, title='Wab Old EcalTrigRes',line_width=2, x_label=xlabel,line_color=colors[2])

    #tritrig
    name = 'tritrig'+tritrig_h.GetName().replace('vtxana_kf','')
    xlabel = tritrig_h.GetXaxis().GetTitle()
    utils.format_TH1(tritrig_h, name=name, title='Tritrig New EcalTrigRes',line_width=2, x_label=xlabel,line_color=colors[3])
    #wab
    name = 'wab'+wab_h.GetName().replace('vtxana_kf','')
    xlabel = wab_h.GetXaxis().GetTitle()
    utils.format_TH1(wab_h, name=name, title='Wab New EcalTrigRes',line_width=2, x_label=xlabel,line_color=colors[4])

    #Make tritrig+wab plot
    tritrig_wab_h = tritrig_h.Clone('tritrig_clone')
    tritrig_wab_h.Add(wab_h,1)
    name = 'tritrig+'+wab_h.GetName()
    utils.format_TH1(tritrig_wab_h, name=name, title='Tritrig+Wab New EcalTrigRes',line_width=2, x_label=xlabel,line_color=colors[5])
    print("Tritrig+WAB scaled: ", tritrig_wab_h.GetEntries())

    #Make OLD tritrig+wab plot
    old_tritrig_wab_h = old_tritrig_h.Clone('old_tritrig_clone')
    old_tritrig_wab_h.Add(old_wab_h,1)
    name = 'old_tritrig+'+old_wab_h.GetName()
    utils.format_TH1(old_tritrig_wab_h, name=name, title='Tritrig+Wab Old EcalTrigRes',line_width=2, x_label=xlabel,line_color=colors[6])
    print("old_Tritrig+WAB scaled: ", old_tritrig_wab_h.GetEntries())


    #Plot histograms together
    #utils.SetStyle()
    histograms = [data_h, old_tritrig_h, old_wab_h, old_tritrig_wab_h, tritrig_h, wab_h, tritrig_wab_h]
    utils.plot_TH1s_with_legend(histograms, data_h.GetName().replace('vtxana_kf',''), '.')

    #Plot histogram ratios
    numerators = [data_h, data_h]
    denominators = [tritrig_wab_h, old_tritrig_wab_h]
    names = ["New EcalTrigRes", "Old EcalTrigRes"]
    utils.plot_TH1_ratios_with_legend(histograms, numerators, denominators, names, [colors[0], colors[1]], 'ratios_'+data_h.GetName().replace('vtxana_kf',''), '.')



