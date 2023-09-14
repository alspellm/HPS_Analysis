#!/usr/bin/python3
import ROOT as r
import sys
sys.path.append( '/sdf/group/hps/users/alspellm/projects/THESIS/ana/analysis_scripts/plot_utils')
import my_plot_utils as myutils
import hpstr_utilities as utils
import argparse
parser = argparse.ArgumentParser(description=" ")
parser.add_argument('--outDir', '-d', type=str, dest="outDir", action='store',
                            help="Specify the output directory.", metavar="outDir", default=".")

options = parser.parse_args()
outdir = options.outDir 
outfile = r.TFile("%s/data_mc_comparison_plots.root"%(outdir),"RECREATE")

infile_data = '/sdf/group/hps/users/alspellm/projects/THESIS/data/2016/hps_7800/20230726_recon/20230726_tuples/20230727_ana/ana_files/hadd_hps_7800_KF_24nsClusterWindow_ana_367_files.root'
selection = 'vtxana_kf_vtxSelection'

infile_tri = '/sdf/group/hps/users/alspellm/projects/THESIS/mc/2016/tritrig_beam/pass4_2016_mc/20230628_ecal_trig_res/cluster_window_24ns/recon/20230628_tuple/full_hadd_tritrig-beam_ecal_trig_res_target_beamPosCorr_ana.root'
selection = 'vtxana_kf_vtxSelection'

infile_wab = ''
selection = 'vtxana_kf_vtxSelection'

