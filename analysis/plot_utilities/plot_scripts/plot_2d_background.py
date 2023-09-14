#/bin/env python
import math
import numpy as np
import ROOT as r
import copy
import sys
sys.path.append( '/sdf/group/hps/users/alspellm/projects/THESIS/ana/analysis_scripts/plot_utils')
from optparse import OptionParser
import my_plot_utils as utils

def passDeltaZcut(ele_z0, pos_z0, ele_tl, pos_tl, cut_value):
    val = abs((ele_z0/ele_tl) - (pos_z0/pos_tl))
    if val < cut_value:
        return True
    else:
        return False

def zcut(p0, p1, p2, p3):
    zcut = r.TF1("zcut","[0] + [1]*x + [2]*x*x + [3]*x*x*x",40.0,160.0)
    zcut.FixParameter(0,p0)
    zcut.FixParameter(1,p1)
    zcut.FixParameter(2,p2)
    zcut.FixParameter(3,p3)

    return zcut

bkg_file = '/sdf/group/hps/users/alspellm/projects/THESIS/data/2016/hps_7800/20230726_recon/20230726_tuples/20230727_ana/ana_files/hadd_hps_7800_KF_24nsClusterWindow_ana_367_files.root'

bkg_hh = r.TH2F("Run_7800_Data_in_SR","Run 7800 Data in SR;M_{vtx} [MeV]; Z_{vtx}",190,10.0,200.,160,-20.,60.)

tree_name = 'vtxana_kf_Tight_2016_simp_reach_SR/vtxana_kf_Tight_2016_simp_reach_SR_tree'
infile = r.TFile('%s'%(bkg_file),"READ")
infile.cd()
tree = infile.Get('%s'%(tree_name))
print(tree.GetName())

deltaZcut_lt = 30.0
i = 0
for ev in tree:
    i = i+1
    if i%10000 == 0:
        print('event ',i)
    mass = 1000.0*ev.unc_vtx_mass
    zvtx = ev.unc_vtx_z
    ele_z0 = ev.unc_vtx_ele_track_z0
    pos_z0 = ev.unc_vtx_pos_track_z0
    ele_tl = ev.unc_vtx_ele_track_tanLambda
    pos_tl = ev.unc_vtx_pos_track_tanLambda
    if not passDeltaZcut(ele_z0, pos_z0, ele_tl, pos_tl, deltaZcut_lt):
        continue
    bkg_hh.Fill(mass,zvtx)

#deltaZ = 30.0
zcut_params = [51.6922, -0.443276, 0.000681226, 3.32196e-06]
zcut = zcut(zcut_params[0],zcut_params[1],zcut_params[2],zcut_params[3])

canvas = r.TCanvas('run_7800_bkg_deltaZ','run_7800_bkg_deltaZ',1200,1080)
canvas.cd()
bkg_hh.Draw("colz")
zcut.Draw("same")
text = ['abs(#Delta track z0/tan#lambda) < %s'%(str(deltaZcut_lt))]
utils.InsertText(insertText=text)

canvas.SaveAs('run_7800_bkg_with_deltaZ_lt_%s'%(str(deltaZcut_lt)))



