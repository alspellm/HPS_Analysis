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

def projectVertex(target_pos, vz, pz, px, py, vx, vy):
    projx = vx - ((vz-target_pos)*(px/pz))
    projy = vy - ((vz-target_pos)*(py/pz))

    return projx, projy

def runVertex2DFit(vtx_proj_hh, fit_params, mass, outdir, outfile, nsigma=1.5):
    projy = vtx_proj_hh.ProjectionY('projy',0, -1, "")
    print('projy entries: ', projy.GetEntries())
    mean_y, sig_y = gaus1DFit(projy, -0.5, 0.5)
    projx = vtx_proj_hh.ProjectionX('projx',0, -1, "")
    mean_x, sig_x = gaus1DFit(projx, -2.0, 2.0)

    projy.Write()
    projx.Write()

    del projy
    del projx

    outfile.cd()

    fitFunc = r.TF2("gaussian", gauss2DFit_Rotated, -5.0,5.0,-1.0, 1.0, 6)
    fitFunc.SetRange(mean_x - (nsigma*sig_x), mean_y - (nsigma*sig_y), mean_x+(nsigma*sig_x), mean_y+(nsigma*sig_y))
    fitFunc.SetParameters(1.0, mean_x, mean_y, sig_x, sig_y, 1.0) 
    vtx_proj_hh.Fit(fitFunc, "RS")
    params = fitFunc.GetParameters()
    xpos = params[1]
    ypos = params[2]
    xsigma = params[3]
    ysigma = params[4]
    angle = params[5]
    xrot, yrot = rotat_coords(xpos, ypos, -angle)
    print(xpos, ypos)

    canvas = r.TCanvas('simp_%s_mev_vertex_proj'%(mass), "Signal %s MeV Vertex Projection"%(mass),2400, 1400)
    vtx_proj_hh.GetXaxis().SetRangeUser(-1.5, 1.5)
    vtx_proj_hh.GetYaxis().SetRangeUser(-1.5, 1.5)
    vtx_proj_hh.Draw("COLZ")
    fitFunc.Draw("SAME")
    canvas.Write()
    canvas.SaveAs('%s/vertex_projection_rot2dGausFit_signal_%s_mev.png'%(outdir,str(mass)))
    canvas.Close()

    del fitFunc

    fit_params[mass] = [xpos, ypos, xsigma, ysigma, angle]

def graphVertexProjectionFits(outfile, outdir, proj_fits):
    outfile.cd()
    outfile.mkdir('fits_results')
    rundir = outfile.GetDirectory('fits_results')
    rundir.cd()

    mass = []
    xpositions = []
    ypositions = []
    xsigmas = []
    ysigmas = []
    rot_angles = []

    for m,params in proj_fits.items():
       mass.append(int(m)) 
       xpositions.append(params[0])
       ypositions.append(params[1])
       xsigmas.append(params[2])
       ysigmas.append(params[3])
       rot_angles.append(params[4])

    xmean_gr = r.TGraph(len(mass), np.array(mass, dtype=float), np.array(xpositions, dtype=float))
    utils.formatHisto(xmean_gr, name='xpos', title='<x> pos data', x_label=' Mass [MeV]', y_label='Position [mm]', marker_color=colors[0], marker_size=2.0)
    xmean_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    xmean_gr.Write()
    #xsigma
    xsig_gr = r.TGraph(len(mass), np.array(mass, dtype=float), np.array(xsigmas, dtype=float))
    utils.formatHisto(xsig_gr, name='sigma_x', title='#sigma_{x} data', x_label=' Mass [MeV]', y_label='#sigma_{x} [mm]', marker_color=colors[0], marker_size=2.0, marker_style = 34)
    xsig_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    xsig_gr.Write()
    #y mean
    ymean_gr = r.TGraph(len(mass), np.array(mass, dtype=float), np.array(ypositions, dtype=float))
    utils.formatHisto(ymean_gr, name='ypos', title='<y> pos data', x_label=' Mass [MeV]', y_label='Position [mm]', marker_color=colors[1], marker_size=2.0)
    ymean_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    ymean_gr.Write()
    #ysigma
    ysig_gr = r.TGraph(len(mass), np.array(mass, dtype=float), np.array(ysigmas, dtype=float))
    utils.formatHisto(ysig_gr, name='sigma_y', title='#sigma_{y} data', x_label=' Mass [MeV]', y_label='#sigma_{y} [mm]', marker_color=colors[1], marker_size=2.0, marker_style = 34)
    ysig_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    ysig_gr.Write()
    #Fit Rotation Angles
    angles_gr = r.TGraph(len(mass), np.array(mass, dtype=float), np.array(rot_angles, dtype=float))
    utils.formatHisto(angles_gr, name='vtx_proj_fit_angles', title='rot angles data', x_label=' Mass [MeV]', y_label='Rotation Angle', marker_color=colors[5], marker_size=2.0)
    angles_gr.Write()


    plots = [xmean_gr, xsig_gr, ymean_gr, ysig_gr]
    canvas = utils.plotTGraphs(plots, 'vertex_projections_signal_masses',drawOptions='P')
    canvas.cd()
    canvas.SetGrid(1,1)

    legend = utils.makeLegend(canvas, plots)

    canvas.SaveAs('%s/vertex_projections_signal_masses.png'%(outdir))
    canvas.Write()
    canvas.Close()


    canvas = r.TCanvas('simp_vertex_proj_rotations','simp_vertex_proj_rotations',2400,1440)
    canvas.cd()
    angles_gr.Draw("AP")
    canvas.SaveAs('%s/simp_vertex_proj_rotations.png'%(outdir))

    canvas.Write()
    canvas.Close()


####################################################################################
#Set Style
style = utils.SetMyStyle(setOptStat=0, setOptFit=1)
colors = utils.getColorsHPS()
outdir = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/tight_selection_studies/v0_projection/initial_look_20230914/signal/vertex_projection_fits'
output_json_file = '%s/simp_signal_vertex_projections.json'%(outdir)
outfile = r.TFile('signal_vertex_projections.root',"RECREATE")

fit_params = {}
directory = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/tight_selection_studies/v0_projection/initial_look_20230914/signal'
pattern = 'hadd_simp*.root'
file_list = glob.glob(f"{directory}/{pattern}")
tree = 'vtxana_Tight_L1L1_nvtx1/vtxana_Tight_L1L1_nvtx1_tree'

for filepath in file_list:
    filename = os.path.basename(filepath)
    mass = filename.split('_')[3]
    print(filepath)
    print(mass)
    vtx_proj_hh = r.TH2F('vtx_projx_v_projy_simp_%s_mev'%(mass),'vtx_projx_v_projy_%s_mev;vtx projx [mm];vtx projy[mm]'%(mass), 100 , -3.0, 3.0, 100, -1.5, 1.5)

    #Loop over Data
    arr = rnp.root2array(filepath, tree)
    df = pd.DataFrame(arr)
    i = 0
    for index, row, in df.iterrows():
        if i%10000 == 0:
            print(i)
        i = i+1
        #vertex fit momentum
        vtx_x = row['unc_vtx_x']
        vtx_y = row['unc_vtx_y']
        vtx_z = row['unc_vtx_z']
        vtx_fit_px = row['unc_vtx_px']
        vtx_fit_py = row['unc_vtx_py']
        vtx_fit_pz = row['unc_vtx_pz']
        #Project Vertex using vertex fit momentum
        vtx_fit_projx, vtx_fit_projy = projectVertex(-4.3, vtx_z, vtx_fit_pz, vtx_fit_px, vtx_fit_py, vtx_x, vtx_y)
        vtx_proj_hh.Fill(vtx_fit_projx, vtx_fit_projy)

    runVertex2DFit(vtx_proj_hh, fit_params, mass, outdir, outfile, nsigma=1.5)
    
graphVertexProjectionFits(outfile, outdir, fit_params)
print(fit_params)

fit_results = {}
for mass in fit_params.keys():
    entry = {'target_position':-4.3, 'rotated_mean_x': fit_params[mass][0], 'rotated_mean_y':fit_params[mass][1], 'rotated_sigma_x':fit_params[mass][2], 
            'rotated_sigma_y':fit_params[mass][3], 'rotation_angle_mrad': fit_params[mass][4]}
    fit_results[mass] = entry
with open(output_json_file, "w") as json_file:
    json.dump(fit_results, json_file, indent=4)
