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

def runVertex2DFit(vtx_proj_hh, run_fit_params, run, isData, outdir, outfile, nsigma=1.5):
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
    if isData:
        outfile.mkdir('%s_data'%(run))
        rundir = outfile.GetDirectory('%s_data'%(run))
    else:
        outfile.mkdir('%s_mc'%(run))
        rundir = outfile.GetDirectory('%s_mc'%(run))

    rundir.cd()
    print("Run number: ", run)

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

    if isData:
        canvas = r.TCanvas('run_%s_data_vertex_proj'%(run), "Run %s Data Vertex Projection"%(run),2400, 1400)
    else:
        canvas = r.TCanvas('run_%s_mc_vertex_proj'%(run), "Run %s MC Vertex Projection"%(run),2400, 1400)
    vtx_proj_hh.GetXaxis().SetRangeUser(-1.5, 1.5)
    vtx_proj_hh.GetYaxis().SetRangeUser(-1.5, 1.5)
    vtx_proj_hh.Draw("COLZ")
    fitFunc.Draw("SAME")
    canvas.Write()
    if isData:
        canvas.SaveAs('%s/vertex_projection_rot2dGausFit_run_%s_data.png'%(outdir,str(run)))
    else:
        canvas.SaveAs('%s/vertex_projection_rot2dGausFit_run_%s_mc.png'%(outdir,str(run)))
    canvas.Close()

    del fitFunc

    run_fit_params[run] = [xpos, ypos, xsigma, ysigma, angle]

def graphVertexProjectionFitsMCData(outfile, outdir, data_run_proj_fits, mc_run_proj_fits, mc_run):
    outfile.cd()
    outfile.mkdir('fits_results')
    rundir = outfile.GetDirectory('fits_results')
    rundir.cd()

    run_numbers = []
    xpositions = []
    ypositions = []
    xsigmas = []
    ysigmas = []
    rot_angles = []

    for run,params in data_run_proj_fits.items():
       run_numbers.append(int(run)) 
       xpositions.append(params[0])
       ypositions.append(params[1])
       xsigmas.append(params[2])
       ysigmas.append(params[3])
       rot_angles.append(params[4])

    ## DATA
    #Graph run dependent values
    xmean_data_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(xpositions, dtype=float))
    utils.formatHisto(xmean_data_gr, name='xpos_data', title='<x> pos data', x_label='Run Number', y_label='Position [mm]', marker_color=colors[0], marker_size=2.0)
    xmean_data_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    xmean_data_gr.Write()
    #xsigma
    xsig_data_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(xsigmas, dtype=float))
    utils.formatHisto(xsig_gr_data, name='sigma_x_data', title='#sigma_{x} data', x_label='Run Number', y_label='#sigma_{x} [mm]', marker_color=colors[0], marker_size=2.0, marker_style = 34)
    xsig_data_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    xsig_data_gr.Write()
    #y mean
    ymean_data_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(ypositions, dtype=float))
    utils.formatHisto(ymean_data_gr, name='ypos_data', title='<y> pos data', x_label='Run Number', y_label='Position [mm]', marker_color=colors[1], marker_size=2.0)
    ymean_data_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    ymean_data_gr.Write()
    #ysigma
    ysig_data_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(ysigmas, dtype=float))
    utils.formatHisto(ysig_data_gr, name='sigma_y_data', title='#sigma_{y} data', x_label='Run Number', y_label='#sigma_{y} [mm]', marker_color=colors[1], marker_size=2.0, marker_style = 34)
    ysig_data_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    ysig_data_gr.Write()
    #Fit Rotation Angles
    angles_data_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(rot_angles, dtype=float))
    utils.formatHisto(angles_data_gr, name='vtx_proj_fit_angles_data', title='rot angles data', x_label='Run Number', y_label='Rotation Angle', marker_color=colors[5], marker_size=2.0)
    angles_data_gr.Write()

    #MC
    run_numbers = []
    xpositions = []
    ypositions = []
    xsigmas = []
    ysigmas = []
    rot_angles = []

    for run in data_run_proj_fits.keys():
       params = mc_run_proj_fits[mc_run]
       run_numbers.append(int(run)) 
       xpositions.append(params[0])
       ypositions.append(params[1])
       xsigmas.append(params[2])
       ysigmas.append(params[3])
       rot_angles.append(params[4])

    xmean_mc_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(xpositions, dtype=float))
    utils.formatHisto(xmean_mc_gr, name='xpos_mc', title='<x> pos mc', x_label='Run Number', y_label='Position [mm]', marker_color=colors[0], marker_size=2.0, line_width=2)
    xmean_mc_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    xmean_mc_gr.Write()
    #xsigma
    xsig_mc_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(xsigmas, dtype=float))
    utils.formatHisto(xsig_gr_mc, name='sigma_x_mc', title='#sigma_{x} mc', x_label='Run Number', y_label='#sigma_{x} [mm]', marker_color=colors[0], marker_size=2.0, marker_style = 34, line_width=2)
    xsig_mc_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    xsig_mc_gr.Write()
    #y mean
    ymean_mc_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(ypositions, dtype=float))
    utils.formatHisto(ymean_mc_gr, name='ypos_mc', title='<y> pos mc', x_label='Run Number', y_label='Position [mm]', marker_color=colors[1], marker_size=2.0, line_width=2)
    ymean_mc_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    ymean_mc_gr.Write()
    #ysigma
    ysig_mc_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(ysigmas, dtype=float))
    utils.formatHisto(ysig_mc_gr, name='sigma_y_mc', title='#sigma_{y} mc', x_label='Run Number', y_label='#sigma_{y} [mm]', marker_color=colors[1], marker_size=2.0, marker_style = 34, line_width=2)
    ysig_mc_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    ysig_mc_gr.Write()
    #Fit Rotation Angles
    angles_mc_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(rot_angles, dtype=float))
    utils.formatHisto(angles_mc_gr, name='vtx_proj_fit_angles_mc', title='rot angles mc', x_label='Run Number', y_label='Rotation Angle', marker_color=colors[5], marker_size=2.0, line_width=2)
    angles_mc_gr.Write()

    data_plots = [xmean_data_gr, xsig_data_gr, ymean_data_gr, ysig_data_gr]
    canvas = utils.plotTGraphs(data_plots, 'vertex_projections_run_summary',drawOptions='P')
    canvas.cd()
    canvas.SetGrid(1,1)
    #plot MC lines on top
    xmean_mc_gr.Draw("LSAME")
    ymean_mc_gr.Draw("LSAME")
    xsig_mc_gr.Draw("LSAME")
    ysig_mc_gr.Draw("LSAME")

    legend_plots = [xmean_data_gr, xsig_data_gr, ymean_data_gr, ysig_data_gr,
            xmean_mc_gr, xsig_mc_gr, ymean_mc_gr, ysig_mc_gr]

    legend = utils.makeLegend(canvas, legend_plots)

    canvas.SaveAs('%s/vertex_projection_fits_run_summary.png'%(outdir))
    canvas.Write()
    canvas.Close()


    canvas = r.TCanvas('Rotation_Angles','Rotation_Angles',2400,1440)
    canvas.cd()
    angles_data_gr.Draw("AP")
    angles_mc_gr.Draw("LSAME")
    legend = utils.makeLegend(canvas, [angles_data_gr, angles_mc_gr])
    canvas.SaveAs('%s/vertex_projection_fit_rotation_angles.png'%(outdir))

    canvas.Write()
    canvas.Close()


def graphVertexProjectionFits(outfile, outdir, run_numbers, xpositions, ypositions, sigmas_x, sigmas_y):
    outfile.cd()
    outfile.mkdir('fits_results')
    rundir = outfile.GetDirectory('fits_results')
    rundir.cd()
    #Graph run dependent values
    xmean_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(xpositions, dtype=float))
    utils.formatHisto(xmean_gr, name='xpos', title='<x> pos', x_label='Run Number', y_label='Position [mm]', marker_color=colors[0], marker_size=2.0)
    xmean_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    xmean_gr.Write()
    #xsigma
    xsig_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(xsigmas, dtype=float))
    utils.formatHisto(xsig_gr, name='sigma_x', title='#sigma_{x}', x_label='Run Number', y_label='#sigma_{x} [mm]', marker_color=colors[0], marker_size=2.0, marker_style = 34)
    xsig_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    xsig_gr.Write()
    #y mean
    ymean_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(ypositions, dtype=float))
    utils.formatHisto(ymean_gr, name='ypos', title='<y> pos', x_label='Run Number', y_label='Position [mm]', marker_color=colors[1], marker_size=2.0)
    ymean_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    ymean_gr.Write()
    #ysigma
    ysig_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(ysigmas, dtype=float))
    utils.formatHisto(ysig_gr, name='sigma_y', title='#sigma_{y}', x_label='Run Number', y_label='#sigma_{y} [mm]', marker_color=colors[1], marker_size=2.0, marker_style = 34)
    ysig_gr.GetYaxis().SetRangeUser(-0.5, 0.5)
    ysig_gr.Write()

    if isData:
        canvas = r.TCanvas('vertex_projections_run_summary_data',2400,1440)
    else:
        canvas = r.TCanvas('vertex_projections_run_summary_mc',2400,1440)
    canvas.cd()
    canvas.SetGrid(1,1)
    xmean_gr.Draw("AP")
    ymean_gr.Draw("PSAME")
    xsig_gr.Draw("PSAME")
    ysig_gr.Draw("PSAME")
    legend = canvas.BuildLegend()
    legend.SetFillStyle(0)
    xmean_gr.SetTitle('%s'%(name))
    if isData:
        canvas.SaveAs('%s/vertex_projection_fits_run_summary_data.png'%(outdir))
    else:
        canvas.SaveAs('%s/vertex_projections_fits_run_summary_mc.png'%(outdir))
    canvas.Write()
    canvas.Close()

    angles_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(rot_angles, dtype=float))
    utils.formatHisto(angles_gr, name='vertex_projection_fit_rotation_angles', title='Rotation Angles', x_label='Run Number', y_label='Rotation Angle', marker_color=colors[5], marker_size=2.0)
    angles_gr.Write()

    canvas = r.TCanvas('Rotation_Angles','Rotation_Angles',2400,1440)
    canvas.cd()
    angles_gr.Draw("AP")
    if isData:
        canvas.SaveAs('%s/vertex_projection_fit_rotation_angles_data.png'%(outdir))
    else:
        canvas.SaveAs('%s/vertex_projection_fit_rotation_angles_mc.png'%(outdir))

    canvas.Write()
    canvas.Close()


####################################################################################
#Set Style
style = utils.SetMyStyle(setOptStat=0, setOptFit=1)
colors = utils.getColorsHPS()

#MC fits
infile = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/tight_selection_studies/v0_projection/initial_look_20230914/tritrig/hadd_full_tritrig_beam_ana_20231019.root'
root_dir = 'vtxana_Tight_nocuts'
tree = '%s/%s_tree'%(root_dir, root_dir)
nsigma = 1.5
outdir = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/tight_selection_studies/v0_projection/initial_look_20230914/new_output'
outfile = r.TFile("%s/v0_projection_rot2DGaussianFits.root"%(outdir),"RECREATE")
run = '7800'
isData = False
mc_run_fit_params = {}
vtx_proj_hh = r.TH2F('vtx_projx_v_projy_mc_run_%s'%(run),'vtx_projx_v_projy_mc_%s;vtx projx [mm];vtx projy[mm]'%(run), 100 , -3.0, 3.0, 100, -1.5, 1.5)

#Loop over MC
arr = rnp.root2array(infile,'%s'%(tree))
df = pd.DataFrame(arr)
i = 0
for index, row, in df.iterrows():
    if i%10000 == 0:
        print(i)
    i = i+1
    if i > 60000:
        break
    run = int(row['run_number'])
    #vertex fit momentum
    vtx_x = row['unc_vtx_x']
    vtx_y = row['unc_vtx_y']
    vtx_z = row['unc_vtx_z']
    vtx_psum = row['unc_vtx_psum']
    vtx_fit_px = row['unc_vtx_px']
    vtx_fit_py = row['unc_vtx_py']
    vtx_fit_pz = row['unc_vtx_pz']

    #Project Vertex using vertex fit momentum
    vtx_fit_projx, vtx_fit_projy = projectVertex(-4.3, vtx_z, vtx_fit_pz, vtx_fit_px, vtx_fit_py, vtx_x, vtx_y)
    vtx_proj_hh.Fill(vtx_fit_projx, vtx_fit_projy)

runVertex2DFit(vtx_proj_hh, mc_run_fit_params, run, isData, outdir, outfile, nsigma=1.5)
del vtx_fit_proj_hh


#sys.exit()
#Fit data runs
data_run_fit_params = {}
isData = True
directory = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/tight_selection_studies/v0_projection/initial_look_20230914/data'
pattern = 'hadd_hps*.root'

file_list = glob.glob(f"{directory}/{pattern}")

for filepath in file_list:
    del arr
    del df
    filename = os.path.basename(filepath)
    run = str(int(filename.split('_')[2]))
    print(filepath)
    print(run)
    vtx_fit_proj_hh = r.TH2F('vtx_projx_v_projy_data_run_%s'%(run),'vtx_projx_v_projy_data_%s;vtx projx [mm];vtx projy[mm]'%(run), 100 , -3.0, 3.0, 100, -1.5, 1.5)

    #Loop over Data
    arr = rnp.root2array(filepath, '%s/%s_tree'%(root_dir,root_dir))
    df = pd.DataFrame(arr)
    i = 0
    for index, row, in df.iterrows():
        if i%10000 == 0:
            print(i)
        i = i+1
        run = int(row['run_number'])
        if run > 7800:
            break
        #vertex fit momentum
        vtx_x = row['unc_vtx_x']
        vtx_y = row['unc_vtx_y']
        vtx_z = row['unc_vtx_z']
        vtx_psum = row['unc_vtx_psum']
        vtx_fit_px = row['unc_vtx_px']
        vtx_fit_py = row['unc_vtx_py']
        vtx_fit_pz = row['unc_vtx_pz']
        #Project Vertex using vertex fit momentum
        vtx_fit_projx, vtx_fit_projy = projectVertex(-4.3, vtx_z, vtx_fit_pz, vtx_fit_px, vtx_fit_py, vtx_x, vtx_y)
        vtx_fit_proj_hh.Fill(vtx_fit_projx, vtx_fit_projy)

    runVertex2DFit(vtx_proj_hh, data_run_fit_params, run, isData, outdir, outfile, nsigma=1.5)


#graphVertexProjectionFitsMCData(outfile, outdir, data_run_proj_fits, mc_run_proj_fits, 7800)

#data = {}
#for i,run in enumerate(run_numbers):
#    run_entry = {'target_position':-4.3, 'rotated_mean_x': xpositions[i], 'rotated_mean_y':ypositions[i], 'rotated_sigma_x':xsigmas[i], 
#            'rotated_sigma_y':ysigmas[i], 'rotation_angle_mrad': rot_angles[i]}
#    data[run] = run_entry
#with open(output_json_file, "w") as json_file:
#    json.dump(data, json_file, indent=4)
#
#return run_numbers, xpositions, ypositions, xsigmas, ysigmas
