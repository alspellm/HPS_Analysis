#!/usr/bin/python3
import ROOT as r
import sys
sys.path.append( '/sdf/group/hps/users/alspellm/projects/THESIS/ana/plot_utils')
import my_plot_utils as utils

selection = 'vtxana_kf_Tight_loose_L1L1'
infile = '/sdf/group/hps/users/alspellm/run/mc/ana/sample_900-9999_hadd_tritrig-beam_ecal_trig_res_ana.root'
outdir = '/sdf/group/hps/users/alspellm/projects/THESIS/ana/tight_selection_studies/z0_tanlambda_investigation_20230815/mc_truth'
outfile = r.TFile("%s/combined_hit_codes.root"%(outdir),"RECREATE")

hc_11 = list(utils.read_2d_plots_from_root_file_dirs(infile,"hc11","z0tanlambda_hh").values())[0]
hc_7 = list(utils.read_2d_plots_from_root_file_dirs(infile,"hc7","z0tanlambda_hh").values())[0]

hc_13 = list(utils.read_2d_plots_from_root_file_dirs(infile,"hc13","z0tanlambda_hh").values())[0]
hc_14 = list(utils.read_2d_plots_from_root_file_dirs(infile,"hc14","z0tanlambda_hh").values())[0]

hc_10 = list(utils.read_2d_plots_from_root_file_dirs(infile,"hc10","z0tanlambda_hh").values())[0]
hc_9 = list(utils.read_2d_plots_from_root_file_dirs(infile,"hc9","z0tanlambda_hh").values())[0]
hc_6 = list(utils.read_2d_plots_from_root_file_dirs(infile,"hc6","z0tanlambda_hh").values())[0]
hc_5 = list(utils.read_2d_plots_from_root_file_dirs(infile,"hc5","z0tanlambda_hh").values())[0]

for i, plot in enumerate(hc_11):
    hc11_hh = plot
    hc7_hh = hc_7[i]
    hc11_hh.Add(hc7_hh)

    outfile.cd()
    name = hc11_hh.GetName().replace('%s'%(selection),'').replace('hc11_1011','')
    hc11_hh.SetName("%s_single_bad_L1_hit"%(name))
    hc11_hh.SetTitle("%s_single_bad_L1_hit"%(name))
    hc11_hh.Write()

for i, plot in enumerate(hc_13):
    hc13_hh = plot
    hc14_hh = hc_14[i]
    hc13_hh.Add(hc14_hh)

    outfile.cd()
    name = hc13_hh.GetName().replace('%s'%(selection),'').replace('hc13_1101','')
    hc13_hh.SetName("%s_single_bad_L2_hit"%(name))
    hc13_hh.SetTitle("%s_single_bad_L2_hit"%(name))

    hc13_hh.Write()

for i, plot in enumerate(hc_10):
    hc10_hh = plot
    hc9_hh = hc_9[i]
    hc6_hh = hc_6[i]
    hc5_hh = hc_5[i]
    hc10_hh.Add(hc9_hh)
    hc10_hh.Add(hc6_hh)
    hc10_hh.Add(hc5_hh)

    outfile.cd()

    name = hc10_hh.GetName().replace('%s'%(selection),'').replace('hc10_1010','')
    hc13_hh.SetName("%s_single_bad_L1_and_L2"%(name))
    hc13_hh.SetTitle("%s_single_bad_L1_and_L2"%(name))
    hc13_hh.Write()
