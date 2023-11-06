#!/usr/bin/python3
import sys
import new_plot_utilities as utils
import ROOT as r
import numpy as np
import glob as glob
import os

style = utils.SetMyStyle(setOptTitle=1)
r.gROOT.SetBatch(1)
colors = utils.getColors()
r.gROOT.SetBatch(1)

#############################################################
for file in sorted(glob.glob('/home/alic/HPS/projects/simps/zbi_opt/zalpha_opt/032423/exp_fit/*.root')):
    filename = os.path.basename(file)
    print(filename)
    slope = filename.split('_')[2]

    signal_zalpha_hh = copy.deepcopy(infile.Get('signal_z0_v_recon_z_alpha_hh'))
    signal_zalpha_hh.SetTitle('signal_z0_v_zalpha_slope_%s'%((slope)))

    tritrig_zalpha_hh = copy.deepcopy(infile.Get('tritrig_z0_v_recon_z_alpha_hh'))
    tritrig_zalpha_hh.SetTitle('tritrig_z0_v_zalpha_slope_%s'%((slope)))

    #signal
    signal_canv = r.TCanvas('signal_z0_v_zalpha_slope_%s'%(slope),'signal_z0_v_zalpha_slope_%s'%(slope), 2560,1440)
    signal_canv.cd()
    signal_canv.SetLogz(1)
    signal_zalpha_hh.Draw("colz")
    signal_zalpha_hh.GetYaxis().SetRangeUser(-4.0,4.0)
    signal_zalpha_hh.GetXaxis().SetRangeUser(-20.0,50.0)
    signal_canv.SaveAs('plots/signal_z0_v_zalpha_slope_%s.png'%(slope))

    #tritrig
    tritrig_canv = r.TCanvas('tritrig_z0_v_zalpha_slope_%s'%(slope),'tritrig_z0_v_zalpha_slope_%s'%(slope), 2560,1440)
    tritrig_canv.cd()
    tritrig_canv.SetLogz(1)
    tritrig_zalpha_hh.Draw("colz")
    tritrig_zalpha_hh.GetYaxis().SetRangeUser(-4.0,4.0)
    tritrig_zalpha_hh.GetXaxis().SetRangeUser(-20.0,50.0)
    tritrig_canv.SaveAs('plots/tritrig_z0_v_zalpha_slope_%s.png'%(slope))




'''
for file in sorted(glob.glob('/home/alic/HPS/projects/simps/zbi_opt/zalpha_opt/032423/exp_fit/*.root')):
    base_filename = os.path.basename(file).replace('.root','')
    print(file)
    slope = base_filename.split('_')[5]
    step_size = base_filename.split('_')[8]
    print(base_filename, slope, step_size)

    infile = r.TFile(file,"READ")
    best_test_cut_ZBi_hh = copy.deepcopy(infile.Get('summary_best_test_cut_ZBi_hh'))
    #best_test_cut_ZBi_hh.RebinY(2)
    best_test_cut_ZBi_hh.SetName('slope_%s_step_%s'%(slope,step_size))
    best_test_cut_ZBi_hh.SetTitle('slope_%s_step_%s;percent signal cut;ZBi'%(slope,step_size))
    best_test_cut_ZBi_hh.SetMarkerStyle(markers[step_size])
    best_test_cut_ZBi_hh.SetMarkerColor(colors[slope])
    plots_hh.append(best_test_cut_ZBi_hh)

canvas = r.TCanvas('zalpha_slope_zbi_comparisons','zalpha_slope_zbi_comparisons',2560,1440)
canvas.cd()

for i,plot in enumerate(plots_hh):
    plot.GetYaxis().SetRangeUser(0.1,6.0)
    if i < 1:
        plot.Draw()
    else:
        plot.Draw("same")
buildLegend(canvas, x1=0.80,y1=0.2,x2=0.95,y2=0.9,textsize=0.02,fill=1)
plots_hh[0].SetTitle('zalpha_slopes_zbi')
canvas.SaveAs('zalpha_slope_comparisons_Zbi.png')
'''
