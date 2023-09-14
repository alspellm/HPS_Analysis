#!/usr/bin/python3
import sys
sys.path.append( '/sdf/group/hps/users/alspellm/projects/THESIS/ana/analysis_scripts/plot_utils')
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

#Define histograms
def histos1D():
    histos = {}
    #momentum
    histos["unc_vtx_track_px"] = r.TH1F('unc_vtx_track_px','unc_vtx_track_px;px [GeV];events',500,-0.5,0.5)
    histos["unc_vtx_track_py"] = r.TH1F('unc_vtx_track_py','unc_vtx_track_py;py [GeV];events',200,-0.1,0.1)
    histos["unc_vtx_track_pz"] = r.TH1F('unc_vtx_track_pz','unc_vtx_track_pz;pz [GeV];events',350,0.0,3.5)
    #histos["unc_vtx_pos_track_px"] = r.TH1F('unc_vtx_pos_track_px','unc_vtx_pos_track_px;px [GeV];events',500,-0.5,0.5)
    #histos["unc_vtx_pos_track_py"] = r.TH1F('unc_vtx_pos_track_py','unc_vtx_pos_track_py;py [GeV];events',200,-0.1,0.1)
    #histos["unc_vtx_pos_track_pz"] = r.TH1F('unc_vtx_pos_track_pz','unc_vtx_pos_track_pz;pz [GeV];events',350,0.0,3.5)

    #covariance
    histos["unc_vtx_cxx"] = r.TH1F('unc_vtx_cxx','unc_vtx_cxx;unc_vtx_cxx;events',200,0.0,0.2)
    histos["unc_vtx_cyy"] = r.TH1F('unc_vtx_cyy','unc_vtx_cyy;unc_vtx_cyy;events',400,0.0,0.04)
    histos["unc_vtx_czz"] = r.TH1F('unc_vtx_czz','unc_vtx_czz;unc_vtx_czz;events',500,0.0,50.0)
    histos["unc_vtx_cyx"] = r.TH1F('unc_vtx_cyx','unc_vtx_cyx;unc_vtx_cyx;events',400,0.0,0.04)
    histos["unc_vtx_czy"] = r.TH1F('unc_vtx_czy','unc_vtx_czy;unc_vtx_czy;events',500,0.0,0.5)
    histos["unc_vtx_czx"] = r.TH1F('unc_vtx_czx','unc_vtx_czx;unc_vtx_czx;events',500,0.0,5.0)

    #vertex momentum
    histos['unc_vtx_px'] = r.TH1F('unc_vtx_px','unc_vtx_px;px [GeV];events',500,-0.5,0.5)
    histos['unc_vtx_py'] = r.TH1F('unc_vtx_py','unc_vtx_py;py [GeV];events',500,-0.5,0.5)
    histos['unc_vtx_pz'] = r.TH1F('unc_vtx_pz','unc_vtx_pz;pz [GeV];events',350,0.0,3.5)

    return histos

def fill1DHisto(histos, name, value, weight):
    histos[name].Fill(value, weight)

def fill2DHisto(histos, name, x, y, weight):
    histos[name].Fill(x, y, weight)


def histos2D():
    histos = {}
    #vertex position
    histos["unc_vtx_x_v_y"] = r.TH2F('unc_vtx_x_v_y','unc_vtx_x_v_y;x[mm];y[mm]',200,-4.0,4.0,200, -4.0, 4.0)
    histos["unc_vtx_x_v_recon_z"] = r.TH2F('unc_vtx_x_v_recon_z','unc_vtx_x_v_recon_z;x[mm];recon_z[mm]',200,-4.0,4.0,200, -20.0, 80.0)
    histos["unc_vtx_y_v_recon_z"] = r.TH2F('unc_vtx_y_v_recon_z','unc_vtx_y_v_recon_z;x[mm];recon_z[mm]',200,-4.0,4.0,200, -20.0, 80.0)

    #covariance
    histos["unc_vtx_cxx_v_recon_z"] = r.TH2F('unc_vtx_cxx_v_recon_z','unc_vtx_cxx_v_recon_z;cxx;recon_z[mm]', 200 , 0.0, 0.2, 200, -20.0, 80.0)
    histos["unc_vtx_cyy_v_recon_z"] = r.TH2F('unc_vtx_cyy_v_recon_z','unc_vtx_cyy_v_recon_z;cyy];recon_z[mm]', 400 , 0.0, 0.04, 200, -20.0, 80.0)
    histos["unc_vtx_czz_v_recon_z"] = r.TH2F('unc_vtx_czz_v_recon_z','unc_vtx_czz_v_recon_z;czz;recon_z[mm]', 500 , 0.0, 50.0, 200, -20.0, 80.0)
    histos["unc_vtx_cyx_v_recon_z"] = r.TH2F('unc_vtx_cyx_v_recon_z','unc_vtx_cyx_v_recon_z;cyx;recon_z[mm]', 400 , 0.0, 0.04, 200, -20.0, 80.0)
    histos["unc_vtx_czy_v_recon_z"] = r.TH2F('unc_vtx_czy_v_recon_z','unc_vtx_czy_v_recon_z;czy;recon_z[mm]', 500 , 0.0, 0.5, 200, -20.0, 80.0)
    histos["unc_vtx_czx_v_recon_z"] = r.TH2F('unc_vtx_czx_v_recon_z','unc_vtx_czx_v_recon_z;czx;recon_z[mm]', 500 , 0.0, 50.0, 200, -20.0, 80.0)
    histos["unc_vtx_cxx_v_cyy"] = r.TH2F('unc_vtx_cxx_v_cyy','unc_vtx_cxx_v_cyy;cxx;cyy', 200 , 0.0, 0.2, 400, 0.0, 0.04)

    #track position vs recon_z
    histos["unc_vtx_track_x_v_recon_z"] = r.TH2F('unc_vtx_track_x_v_recon_z','unc_vtx_track_x_v_recon_z;track x [mm];recon_z[mm]', 1000 , -10.0, 10.0, 200, -20.0, 80.0)
    histos["unc_vtx_track_y_v_recon_z"] = r.TH2F('unc_vtx_track_y_v_recon_z','unc_vtx_track_y_v_recon_z;track y [mm];recon_z[mm]', 1000 , -10.0, 10.0, 200, -20.0, 80.0)
    histos["unc_vtx_track_x_v_track_y"] = r.TH2F('unc_vtx_track_x_v_track_y','unc_vtx_track_x_v_track_y;track x [mm];track y [mm]', 1000 , -10.0, 10.0, 1000, -10.0, 10.0)

    #projected vertex
    histos["unc_vtx_target_projx_v_projy"] = r.TH2F('unc_vtx_target_projx_v_projy','unc_vtx_target_projx_v_projy;vtx projx [mm];vtx projy[mm]', 1000 , -10.0, 10.0, 1000, -10.0, 10.0)
    histos["unc_vtx_target_tracky_v_projy"] = r.TH2F('unc_vtx_target_tracky_v_projy','unc_vtx_target_tracky_v_projy;track y [mm];vtx projy[mm]', 1000 , -10.0, 10.0, 1000, -10.0, 10.0)
    histos["unc_vtx_target_trackx_v_projx"] = r.TH2F('unc_vtx_target_trackx_v_projx','unc_vtx_target_trackx_v_projx;track x [mm];vtx projx[mm]', 1000 , -10.0, 10.0, 1000, -10.0, 10.0)

    return histos

#Defin vertex projection to target
proj_vtx_x = lambda targ_pos, vz, vx, px, pz : vx - (vz-targ_pos)*(px/pz)
proj_vtx_y = lambda targ_pos, vz, vy, py, pz : vy - (vz-targ_pos)*(py/pz)


outfile = r.TFile("v0_proj_ana_run_7800_AtTarget.root","RECREATE")

infilename = '/sdf/group/hps/users/alspellm/projects/THESIS/ana/tight_selection_studies/v0_projection/test_20230906/hps_007800.evio.107_KF_24nsClusterWindow_targetTracks_ana_KF.root'
#infilename = '/sdf/group/hps/users/alspellm/projects/THESIS/ana/tight_selection_studies/v0_projection/test_20230906/hps_007800.evio.107_KF_24nsClusterWindow_AtIP_ana_KF.root'
tree = 'vtxana_kf_Tight_loose/vtxana_kf_Tight_loose_tree'

histos_1d = histos1D()
histos_2d = histos2D()
#Loop over data
arr = rnp.root2array(infilename, tree)
df = pd.DataFrame(arr)
i = 0
for index, row in df.iterrows():
    i = i + 1
    if i%1000 == 0:
        print(i)
    #cov
    cxx = row['unc_vtx_cxx']
    cyy = row['unc_vtx_cyy']
    czz = row['unc_vtx_czz']
    cyx = row['unc_vtx_cyx']
    czy = row['unc_vtx_czy']
    czx = row['unc_vtx_czx']

    #track pos
    ele_track_x = row['unc_vtx_ele_track_x']
    ele_track_y = row['unc_vtx_ele_track_y']
    ele_track_z = row['unc_vtx_ele_track_z']
    pos_track_x = row['unc_vtx_pos_track_x']
    pos_track_y = row['unc_vtx_pos_track_y']
    pos_track_z = row['unc_vtx_pos_track_z']

    #track mom
    ele_track_px = row['unc_vtx_ele_track_px']
    ele_track_py = row['unc_vtx_ele_track_py']
    ele_track_pz = row['unc_vtx_ele_track_pz']
    pos_track_px = row['unc_vtx_pos_track_px']
    pos_track_py = row['unc_vtx_pos_track_py']
    pos_track_pz = row['unc_vtx_pos_track_pz']

    #vertex pos
    vtx_x = row['unc_vtx_x']
    vtx_y = row['unc_vtx_y']
    vtx_z = row['unc_vtx_z']

    #vertex net momentum comps
    vtx_px = ele_track_px + pos_track_px
    vtx_py = ele_track_py + pos_track_py
    vtx_pz = ele_track_pz + pos_track_pz
    target_pos = -4.3
    vtx_projx = proj_vtx_x(target_pos, vtx_z, vtx_x, vtx_px, vtx_pz)
    vtx_projy = proj_vtx_y(target_pos, vtx_z, vtx_y, vtx_py, vtx_pz)

    #Fill 1D histograms
    weight = 1.0
    fill1DHisto(histos_1d, 'unc_vtx_track_px',ele_track_px,weight)
    fill1DHisto(histos_1d, 'unc_vtx_track_py',ele_track_py,weight)
    fill1DHisto(histos_1d, 'unc_vtx_track_pz',ele_track_pz,weight)
    fill1DHisto(histos_1d, 'unc_vtx_track_px',pos_track_px,weight)
    fill1DHisto(histos_1d, 'unc_vtx_track_py',pos_track_py,weight)
    fill1DHisto(histos_1d, 'unc_vtx_track_pz',pos_track_pz,weight)

    fill1DHisto(histos_1d, 'unc_vtx_cxx',cxx,weight)
    fill1DHisto(histos_1d, 'unc_vtx_cyy',cyy,weight)
    fill1DHisto(histos_1d, 'unc_vtx_czz',czz,weight)
    fill1DHisto(histos_1d, 'unc_vtx_cyx',cyx,weight)
    fill1DHisto(histos_1d, 'unc_vtx_czy',czy,weight)
    fill1DHisto(histos_1d, 'unc_vtx_czx',czx,weight)
    fill1DHisto(histos_1d, 'unc_vtx_px',vtx_px,weight)
    fill1DHisto(histos_1d, 'unc_vtx_py',vtx_py,weight)
    fill1DHisto(histos_1d, 'unc_vtx_pz',vtx_pz,weight)

    #Fill 2D histograms
    fill2DHisto(histos_2d, 'unc_vtx_x_v_y', vtx_x, vtx_y, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_x_v_recon_z', vtx_x, vtx_z, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_y_v_recon_z', vtx_y, vtx_z, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_cxx_v_recon_z', cxx, vtx_z, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_cyy_v_recon_z', cyy, vtx_z, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_czz_v_recon_z', czz, vtx_z, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_cyx_v_recon_z', cyx, vtx_z, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_czy_v_recon_z', czy, vtx_z, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_czx_v_recon_z', czx, vtx_z, 1.0)

    fill2DHisto(histos_2d, 'unc_vtx_track_x_v_recon_z', ele_track_x, vtx_z, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_track_x_v_recon_z', pos_track_x, vtx_z, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_track_y_v_recon_z', ele_track_y, vtx_z, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_track_y_v_recon_z', pos_track_y, vtx_z, 1.0)

    fill2DHisto(histos_2d, 'unc_vtx_track_x_v_track_y', ele_track_x, pos_track_x, 1.0)

    fill2DHisto(histos_2d, 'unc_vtx_target_projx_v_projy', vtx_projx, vtx_projy, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_target_tracky_v_projy', ele_track_y, vtx_projy, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_target_tracky_v_projy', pos_track_y, vtx_projy, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_target_trackx_v_projx', ele_track_x, vtx_projx, 1.0)
    fill2DHisto(histos_2d, 'unc_vtx_target_trackx_v_projx', pos_track_x, vtx_projx, 1.0)
    
outfile.cd()
for plot in histos_1d:
    histos_1d[plot].Write()

for plot in histos_2d:
    histos_2d[plot].Write()
    


