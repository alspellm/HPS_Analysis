#!/usr/bin/python3
import vertex_target_projection_fitting_mc as proj
import sys
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

#Input files
fit_results = {}
directory = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/tight_selection_studies/v0_projection/initial_look_20230914/data'
pattern = 'hadd_hps*.root'
file_list = glob.glob(f"{directory}/{pattern}")
root_dir = 'vtxana_Tight_nocuts'
print(file_list)

#Output file
outdir = '.'
outfilename = 'v0proj_rot2DGaussianFits_data'
outfile = r.TFile("%s/%s.root"%(outdir, outfilename),"RECREATE")
output_json_file = '%s/%s.json'%(outdir, outfilename)

arr = None
df = None
for filepath in file_list:
    del arr
    del df
    filename = os.path.basename(filepath)
    run = str(int(filename.split('_')[2]))
    print(filepath)
    print(run)
    vtx_proj_hh = r.TH2F('vtx_projx_v_projy_data_run_%s'%(run),'vtx_projx_v_projy_data_%s;vtx projx [mm];vtx projy[mm]'%(run), 100 , -3.0, 3.0, 100, -1.5, 1.5)

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
        vtx_fit_projx, vtx_fit_projy = proj.projectVertex(-4.3, vtx_z, vtx_fit_pz, vtx_fit_px, vtx_fit_py, vtx_x, vtx_y)
        vtx_proj_hh.Fill(vtx_fit_projx, vtx_fit_projy)

    proj.runVertex2DFit(vtx_proj_hh, fit_results, run, outdir, outfile, nsigma=1.5)

proj.writeFitResultsJson(fit_results, output_json_file)

#Graph results
outfile.cd()
run_numbers = fit_results.keys()
xpositions = [x for x in fit_results.values[0]]

xmean_gr = r.TGraph(len(run_numbers), np.array(run_numbers, dtype=float), np.array(xpositions, dtype=float))
xmean_gr.SetName('fit_mean_x')
xmean_gr.Write()

