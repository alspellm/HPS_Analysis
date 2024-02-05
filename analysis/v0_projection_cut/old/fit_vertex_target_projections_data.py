#!/usr/bin/python3
import sys
import new_plot_utilities as utils
import ROOT as r
import numpy as np
import math
import root_numpy as rnp
import pandas as pd
import os
import re
import glob as glob
import json

r.gROOT.SetBatch(1)
colors = utils.getColors()

def runAnalysis(outdir, outfile, run_histos, name, nsigma):
    plots_dir = '%s/%s_plots_nsigma_%s'%(outdir, name, nsigma)
    output_json_file = '%s/data_v0proj_fits_2016_config.json'%(outdir)
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

    outfile.cd()
    outfile.mkdir(name)
    subdir = outfile.GetDirectory(name)
    subdir.cd()

    xpositions = []
    ypositions = []
    xsigmas = []
    ysigmas = []
    run_numbers = []
    rot_angles = []
    for run, histos in runs.items():
        if not outfile.GetDirectory('%s/%s'%(name,str(run))):
            outfile.mkdir('%s/%s'%(name,str(run)))
        rundir = outfile.GetDirectory('%s/%s'%(name,str(run)))
        rundir.cd()
        for key, histo in histos[1].items():
            if key == name:
                histo.Write()
                #Get 1d fit ranges
                projy = histo.ProjectionY('projy_%s'%(histo.GetName()),0, -1, "")
                mean_y, sig_y = gaus1DFit(projy, -0.5, 0.5)
                projx = histo.ProjectionX('projx_%s'%(histo.GetName()),0, -1, "")
                mean_x, sig_x = gaus1DFit(projx, -2.0, 2.0)

                del projy
                del projx

                fitFunc = r.TF2("gaussian", gauss2DFit_Rotated, -5.0,5.0,-1.0, 1.0, 6)
                fitFunc.SetRange(mean_x - (nsigma*sig_x), mean_y - (nsigma*sig_y), mean_x+(nsigma*sig_x), mean_y+(nsigma*sig_y))
                fitFunc.SetParameters(1.0, mean_x, mean_y, sig_x, sig_y, 1.0) 
                histo.Fit(fitFunc, "RS")
                params = fitFunc.GetParameters()
                xpos = params[1]
                ypos = params[2]
                xsigma = params[3]
                ysigma = params[4]
                angle = params[5]
                xrot, yrot = rotat_coords(xpos, ypos, -angle)
                print(xpos, ypos)
                print(xrot, yrot)

                canvas = r.TCanvas("canvas", "Fitted Histogram",2400, 1400)
                histo.GetXaxis().SetRangeUser(-1.5, 1.5)
                histo.GetYaxis().SetRangeUser(-1.5, 1.5)
                histo.Draw("COLZ")
                fitFunc.Draw("SAME")
                canvas.Write()
                canvas.SaveAs('%s/run_%s_%s_2d_gaus_fit_nsigma_%s.png'%(plots_dir,str(run), name, nsigma))
                canvas.Close()

                del fitFunc

                run_numbers.append(run)
                xpositions.append(xpos)
                ypositions.append(ypos)
                xsigmas.append(xsigma)
                ysigmas.append(ysigma)
                #rot_angles.append(math.degrees(angle))
                rot_angles.append(angle*1000)

    subdir.cd()
    #Graph run dependent values
    xmean_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(xpositions, dtype=float))
    utils.formatHisto(xmean_gr, name='%s_x_pos'%(name), title='<x> pos', x_label='Run Number', y_label='Position [mm]', marker_color=colors[0], marker_size=2.0)
    xmean_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    xmean_gr.Write()
    #xsigma
    xsig_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(xsigmas, dtype=float))
    utils.formatHisto(xsig_gr, name='%s_x_sigma'%(name), title='#sigma_{x}', x_label='Run Number', y_label='#sigma_{x} [mm]', marker_color=colors[0], marker_size=2.0, marker_style = 34)
    xsig_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    xsig_gr.Write()
    #y mean
    ymean_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(ypositions, dtype=float))
    utils.formatHisto(ymean_gr, name='%s_y_pos'%(name), title='<y> pos', x_label='Run Number', y_label='Position [mm]', marker_color=colors[1], marker_size=2.0)
    ymean_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    ymean_gr.Write()
    #ysigma
    ysig_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(ysigmas, dtype=float))
    utils.formatHisto(ysig_gr, name='%s_y_sigma'%(name), title='#sigma_{y}', x_label='Run Number', y_label='#sigma_{y} [mm]', marker_color=colors[1], marker_size=2.0, marker_style = 34)
    ysig_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    ysig_gr.Write()

    canvas = r.TCanvas('%s'%(name),'%s'%(name),2400,1440)
    canvas.cd()
    canvas.SetGrid(1,1)
    xmean_gr.Draw("AP")
    ymean_gr.Draw("PSAME")
    xsig_gr.Draw("PSAME")
    ysig_gr.Draw("PSAME")
    legend = canvas.BuildLegend()
    legend.SetFillStyle(0)
    xmean_gr.SetTitle('%s'%(name))
    canvas.SaveAs('%s/%s_nsigma_%s.png'%(plots_dir, name,nsigma))
    #canvas.Write()
    canvas.Close()

    angles_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(rot_angles, dtype=float))
    utils.formatHisto(angles_gr, name='%s_rotation_angles'%(name), title='Rotation Angles', x_label='Run Number', y_label='Rotation Angle [mrad]', marker_color=colors[5], marker_size=2.0)
    angles_gr.Write()

    canvas = r.TCanvas('Rotation_Angles','Rotation_Angles',2400,1440)
    canvas.cd()
    angles_gr.Draw("AP")
    canvas.SaveAs('%s/rotation_angles_nsigma_%s.png'%(plots_dir,nsigma))
    #canvas.Write()
    canvas.Close()


    data = {}
    for i,run in enumerate(run_numbers):
        run_entry = {'target_position':-4.3, 'rotated_mean_x': xpositions[i], 'rotated_mean_y':ypositions[i], 'rotated_sigma_x':xsigmas[i], 
                'rotated_sigma_y':ysigmas[i], 'rotation_angle_mrad': rot_angles[i]}
        data[run] = run_entry
    with open(output_json_file, "w") as json_file:
        json.dump(data, json_file, indent=4)
    
    return run_numbers, xpositions, ypositions, xsigmas, ysigmas

def fill1DHisto(histos, name, value, weight):
    histos[name].Fill(value, weight)

def fill2DHisto(histos, name, x, y, weight):
    histos[name].Fill(x, y, weight)

def gaus2DFit(histo, fitFunc, nsig):
    yproj = histo.ProjectionY('projy_%s'%(histo.GetName()),0, -1, "")
    fity = yproj.Fit("gaus","QRS","",-0.5,0.5)
    ymean = fity.Parameter(1)
    ysig = fity.Parameter(2)
    yproj.Write()
    print('ymean,ysig', ymean, ysig) 
    #fit x gaussian
    xproj = histo.ProjectionX('projx',0, -1, "")
    fitx = xproj.Fit("gaus","QRS","",-2.0,2.0)
    xproj.Write()
    xmean = fitx.Parameter(1)
    xsig = fitx.Parameter(2)
    print('xmean,xsig', xmean, xsig) 

    #Fit 2D gaussian
    xmin = xmean-(nsig*xsig)
    xmax = xmean+(nsig*xsig)
    ymin = ymean-(nsig*ysig)
    ymax = ymean+(nsig*ysig)
    fitFunc.SetParameters(1.0,xmean,xsig,ymean,ysig)
    fitFunc.SetRange(xmin, ymin, xmax, ymax)
    fitresult = histo.Fit(fitFunc,"SR")
    return fitresult

def projectVertex(target_pos, vz, pz, px, py, vx, vy):
    projx = vx - ((vz-target_pos)*(px/pz))
    projy = vy - ((vz-target_pos)*(py/pz))

    return projx, projy

def rotat_coords(x, y, angle):
    x_rotated = x*math.cos(angle) - y*math.sin(angle)
    y_rotated = x*math.sin(angle) + y*math.cos(angle)
    return (x_rotated, y_rotated)

def gauss2DFit_Rotated(x, par):
    amplitude = par[0]
    meanx = par[1]
    meany = par[2]
    sigmax = par[3]
    sigmay = par[4]
    angle = par[5]
    x_rotated = x[0] * math.cos(angle) - x[1] * math.sin(angle)
    y_rotated = x[0] * math.sin(angle) + x[1] * math.cos(angle)

    exponent = -0.5 * ((x_rotated - meanx)**2 / sigmax**2 + (y_rotated - meany)**2 / sigmay**2)
    return amplitude * math.exp(exponent)

def gaus1DFit(histo, xmin, xmax):
    fit = histo.Fit("gaus", "QRS", "", xmin, xmax)
    mean = fit.Parameter(1)
    sig = fit.Parameter(2)
    return mean, sig

####################################################################################
infilename = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/tight_selection_studies/v0_projection/initial_look_20230914/data/hadd_hps_BLPass4_ana_20230919.root'
tree = 'vtxana_kf_Tight_loose/vtxana_kf_Tight_loose_tree'

#Set Style
style = utils.SetMyStyle(setOptStat=1, setOptFit=1)
colors = utils.getColorsHPS()

#Vertex Position histograms
vertex_pos_histos2d = {}

#Loop over data
runs = {}
arr = rnp.root2array(infilename, tree)
df = pd.DataFrame(arr)
i = 0
for index, row in df.iterrows():

    if i%10000 == 0:
        print(i)
    i = i+1
    #if i > 50000:
    #    break
    #Make run dependent histograms
    run = int(row['run_number'])
    if run not in runs:
        histos_1d = {}
        histos_2d = {}
        histos = []
        
        #1d histograms
        histos_1d['unc_vtx_proj_dx'] = r.TH1F('unc_vtx_proj_dx_%s'%(str(run)),'unc_vtx_proj_dx_%s; dx [mm]; Events'%(str(run)),200, -1.0, 1.0)
        histos_1d['unc_vtx_proj_dy'] = r.TH1F('unc_vtx_proj_dy_%s'%(str(run)),'unc_vtx_proj_dy_%s; dy [mm]; Events'%(str(run)),200, -1.0, 1.0)
        histos_1d['unc_vtx_proj_dz'] = r.TH1F('unc_vtx_proj_dz_%s'%(str(run)),'unc_vtx_proj_dz_%s; dz [mm]; Events'%(str(run)),100, -10.0, 10.0)

        histos_2d["unc_vtx_x_v_y"] = r.TH2F('unc_vtx_x_v_y_run_%s'%(str(run)),
                'unc_vtx_x_v_y_run_%s;x[mm];y[mm]'%(str(run)),100,-3.0,3.0,75, -1.5, 1.5)
        histos_2d["unc_vtx_fit_proj"] = r.TH2F('unc_vtx_fit_target_projx_v_projy_run_%s'%(run),'unc_vtx_fit_target_projx_v_projy_%s;vtx projx [mm];vtx projy[mm]'%(run), 100 , -3.0, 3.0, 100, -1.5, 1.5)

        histos = [histos_1d, histos_2d]
        runs[run] = histos

    #vertex pos
    vtx_x = row['unc_vtx_x']
    vtx_y = row['unc_vtx_y']
    vtx_z = row['unc_vtx_z']
    vtx_psum = row['unc_vtx_psum']
    #vertex fit momentum
    vtx_fit_px = row['unc_vtx_px']
    vtx_fit_py = row['unc_vtx_py']
    vtx_fit_pz = row['unc_vtx_pz']

    #Project Vertex using vertex fit momentum
    vtx_fit_projx, vtx_fit_projy = projectVertex(-4.3, vtx_z, vtx_fit_pz, vtx_fit_px, vtx_fit_py, vtx_x, vtx_y)
    histos_2d['unc_vtx_fit_proj'].Fill(vtx_fit_projx, vtx_fit_projy)

    #projection residuals
    dx = vtx_x - vtx_fit_projx
    dy = vtx_y - vtx_fit_projy
    dz = vtx_z - -4.3
    histos_1d['unc_vtx_proj_dx'].Fill(dx)
    histos_1d['unc_vtx_proj_dy'].Fill(dy)
    histos_1d['unc_vtx_proj_dz'].Fill(dz)

    #Fill Histograms
    #angle = 1.3
    #vtx_x_rotated = vtx_x*math.cos(angle) - vtx_y*math.sin(angle)
    #vtx_y_rotated = vtx_x*math.sin(angle) + vtx_y*math.cos(angle)
    #histos_2d['unc_vtx_x_v_y'].Fill(vtx_x_rotated, vtx_y_rotated)
    histos_2d['unc_vtx_x_v_y'].Fill(vtx_x, vtx_y)

nsigma = 1.5
outdir = 'vtxproj_fits_data'
outfile = r.TFile("%s/run_dependent_vtx_proj_fits_nsigma_%s.root"%(outdir,nsigma),"RECREATE")

#Run 2D Gaus Fits
#runAnalysis(outdir,outfile, runs, 'unc_vtx_x_v_y', nsigma)
runAnalysis(outdir,outfile, runs, 'unc_vtx_fit_proj', nsigma)

#Save 1D Histograms
outfile.cd()
dxs = []
dys = []
dzs = []
run_numbers = []
for run, histos in runs.items():
    for key, histo in histos[0].items():
        histo.Write()

