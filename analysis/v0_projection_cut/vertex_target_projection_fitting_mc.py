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

def rotate_coordinates(x, y, angle):
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
    try:
        mean = fit.Parameter(1)
    except:
        return None, None
    sig = fit.Parameter(2)
    return mean, sig

def projectVertex(target_pos, vz, pz, px, py, vx, vy):
    projx = vx - ((vz-target_pos)*(px/pz))
    projy = vy - ((vz-target_pos)*(py/pz))

    return projx, projy

def runVertex2DFit(vtx_proj_hh, run_fit_params, run, outdir, outfile, nsigma=1.5):

    #Get the not-rotated x and y projections to seed the fits
    projy = vtx_proj_hh.ProjectionY('projy',0, -1, "")
    mean_y, sig_y = gaus1DFit(projy, -0.5, 0.5)
    projx = vtx_proj_hh.ProjectionX('projx',0, -1, "")
    mean_x, sig_x = gaus1DFit(projx, -2.0, 2.0)
    projy.Write()
    projx.Write()
    del projy
    del projx

    if mean_y is None or mean_x is None:
        return

    #Make dir for run being fit
    outfile.cd()
    outfile.mkdir('fit_results_run_%s'%(run))
    rundir = outfile.GetDirectory('fit_results_run_%s'%(run))
    rundir.cd()

    #Do rot 2d gaussian fit
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
    xrot, yrot = rotate_coordinates(xpos, ypos, -angle)

    canvas = r.TCanvas('target_proj_vtx_fits_run_%s'%(run), "Run %s Target Vertex Projection Fit"%(run),2400, 1400)
    vtx_proj_hh.GetXaxis().SetRangeUser(-1.5, 1.5)
    vtx_proj_hh.GetYaxis().SetRangeUser(-1.5, 1.5)
    vtx_proj_hh.Draw("COLZ")
    fitFunc.Draw("SAME")
    canvas.Write()
    canvas.Close()

    del fitFunc
    run_fit_params[run] = [xpos, ypos, xsigma, ysigma, angle]

def writeFitResultsJson(fit_results, output_json_file):
    #Save fit results to json file
    json_data = {}
    for key, values in fit_results.items():
        entry = {'target_position':-4.3, 'rotated_mean_x': values[0], 'rotated_mean_y':values[1], 'rotated_sigma_x':values[2], 'rotated_sigma_y':values[3], 'rotation_angle_mrad': values[4]}
        json_data[key] = entry
        with open(output_json_file, "w") as json_file:
            json.dump(json_data, json_file, indent=4)

####################################################################################

def main():
    #Set Style
    style = utils.SetMyStyle(setOptStat=0, setOptFit=1)
    colors = utils.getColorsHPS()

    #Input file
    infile = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/tight_selection_studies/v0_projection/initial_look_20230914/tritrig/hadd_full_tritrig_beam_ana_20231019.root'
    root_dir = 'vtxana_Tight_nocuts'
    tree = '%s/%s_tree'%(root_dir, root_dir)

    #Output file
    outdir = '.'
    outfilename = 'v0proj_rot2DGaussianFits_tritrig'
    outfile = r.TFile("%s/%s.root"%(outdir, outfilename),"RECREATE")
    output_json_file = '%s/%s.json'%(outdir, outfilename)

    #Specify run number
    run = '7800'
    #Hold fit parameter results
    fit_results = {}

    #Init histo that will be fit with rotated Gaussian
    vtx_proj_hh = r.TH2F('vtx_projx_v_projy_mc_run_%s'%(run),'vtx_projx_v_projy_mc_%s;vtx projx [mm];vtx projy[mm]'%(run), 100 , -3.0, 3.0, 100, -1.5, 1.5)

    #Loop over MC
    arr = rnp.root2array(infile,'%s'%(tree))
    df = pd.DataFrame(arr)
    i = 0
    for index, row, in df.iterrows():
        if i%10000 == 0:
            print(i)
        i = i+1
        #Limit 60000 events because of root2array mem limits
        if i > 60000:
            break

        run = int(row['run_number']) #Get run number from file
        vtx_x = row['unc_vtx_x']
        vtx_y = row['unc_vtx_y']
        vtx_z = row['unc_vtx_z']
        vtx_psum = row['unc_vtx_psum']
        vtx_px = row['unc_vtx_px']
        vtx_py = row['unc_vtx_py']
        vtx_pz = row['unc_vtx_pz']

        #Project vertex back to target location
        target_pos = -4.3 #mm
        vtx_fit_projx, vtx_fit_projy = projectVertex(target_pos, vtx_z, vtx_pz, vtx_px, vtx_py, vtx_x, vtx_y)
        #Fill target projected vertex x and y histogram...to be fit with rot 2d gaus
        vtx_proj_hh.Fill(vtx_fit_projx, vtx_fit_projy)

    #Fit the target projected vertex x/y distribution with rotated 2d gaussian.
    runVertex2DFit(vtx_proj_hh, fit_results, run, outdir, outfile, nsigma=1.5)
    writeFitResultsJson(fit_results, output_json_file)

if __name__ == "__main__":
    main()

