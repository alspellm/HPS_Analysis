#!/usr/bin/python3
import sys
sys.path.append( '/sdf/group/hps/users/alspellm/projects/THESIS/ana/analysis_scripts/plot_utils')
import my_plot_utils as utils
import ROOT as r
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

def get_test_cut_values_for_iter(cuts_hh, iteration):
    iter_bin = cuts_hh.GetXaxis().FindBin(iteration)
    projy = cuts_hh.ProjectionY("iteration_%i_cuts"%(iteration),iter_bin,iter_bin,"")
    cuts_values = {}
    for x in range(projy.GetNbinsX()):
        if projy.GetXaxis().GetBinLabel(x+1) == "":
            continue
        cut = projy.GetXaxis().GetBinLabel(x+1)
        val = round(projy.GetBinContent(x+1),2)
        cutid = x+1
        cuts_values[cutid] = (cut,val)
    return cuts_values


################################################################
parser = argparse.ArgumentParser(description="HTML Config")
parser.add_argument('--inputDir', type=str, dest="inputDir", default = ".")
parser.add_argument('--plotsDir', type=str, dest="plotsDir", default = "plots_1d")
options = parser.parse_args()
inputDir = options.inputDir
plotsDir = options.plotsDir

r.gROOT.SetBatch(1)
colors = utils.getColors()

plots_dir = options.plotsDir
os.makedirs(plots_dir, exist_ok=True)
print("Created directory", plots_dir)
logeps2=-5.5

masses = []
zcuts = []
nsigs = []
nbkgs = []
zbis = []
best_cut_values = []
best_cutids = []

cuts = {}
canvases = []
for filenum, infile in enumerate(sorted(os.listdir(options.inputDir))):
    infile = os.path.join(inputDir,infile)
    if infile.endswith('.root'):
        print(infile)
    else:
        continue
    split_name = os.path.basename(infile).split('_')
    mass = split_name[1]
    if int(mass) > 160:
        continue
    print("Signal mass: ", mass)
    masses.append(int(mass))

    #Read and format TH1s
    best_zbi_h = utils.read_plot_from_root_file(infile, 'zbi_processor_best_test_cut_ZBi_h')
    #utils.format_TH1(best_zbi_h, title='%s_MeV_zbi'%(mass))
    best_zcut_h = utils.read_plot_from_root_file(infile, 'zbi_processor_best_test_cut_zcut_h')
    #utils.format_TH1(best_zcut_h, title='%s_MeV_zcut'%(mass))
    best_nsig_h = utils.read_plot_from_root_file(infile, 'zbi_processor_best_test_cut_nsig_h')
    #utils.format_TH1(best_nsig_h, title='%s_MeV_nsig'%(mass))
    best_nbkg_h = utils.read_plot_from_root_file(infile, 'zbi_processor_best_test_cut_nbkg_h')
    #utils.format_TH1(best_nbkg_h, title='%s_MeV_nbkg'%(mass))
    best_cutid_h = utils.read_plot_from_root_file(infile, 'zbi_processor_best_test_cut_id_h')

    #Save 1D Plots
    zbi_c = utils.Make1DPlot('%s_MeV_zbi'%(mass), best_zbi_h, title='%s_MeV_zbi'%(mass), drawOptions="histtext",outdir = plots_dir)
    zcut_c = utils.Make1DPlot('%s_MeV_zcut'%(mass), best_zcut_h,title='%s_MeV_zcut'%(mass),drawOptions="histtext",outdir = plots_dir)
    nsig_c = utils.Make1DPlot('%s_MeV_nsig'%(mass), best_nsig_h,title='%s_MeV_nsig'%(mass), drawOptions="histtext",outdir = plots_dir)
    nbkg_c = utils.Make1DPlot('%s_MeV_nbkg'%(mass), best_nbkg_h,title='%s_MeV_nbkg'%(mass), drawOptions="histtext",outdir = plots_dir)
    canvases.append(zbi_c)
    canvases.extend([zcut_c, nsig_c, nbkg_c])

    best_iteration = 0.0
    best_zbi = -999.9
    use_first_iter = True
    if use_first_iter:
        best_zbi = best_zbi_h.GetBinContent(best_zbi_h.FindBin(0.0))
        best_iteration = 0

    else:
        for nbin in range(1,best_zbi_h.GetXaxis().GetNbins()+1):
            #print(best_zbi_h.GetBinContent(nbin))
            if best_zbi_h.GetBinError(nbin) <= 0:
                continue
            if best_zbi_h.GetBinContent(nbin) > best_zbi:
                best_zbi = best_zbi_h.GetBinContent(nbin)
                best_iteration = (best_zbi_h.GetBinCenter(nbin))
        print("Best ZBi: ", best_zbi)
        print("Best Iteration: ", best_iteration)
        if best_iteration == 0:
            best_iteration = 10
            best_zbi = 0.0

    print("Final best iteration", best_iteration)
    print("Final ZBi: ", best_zbi)

    best_iteration = int(best_iteration)
    best_zbi = round(best_zcut_h.GetBinContent(best_zbi_h.FindBin(best_iteration)),2)
    best_zcut = round(best_zcut_h.GetBinContent(best_zcut_h.FindBin(best_iteration)),2)
    best_nsig = round(best_nsig_h.GetBinContent(best_nsig_h.FindBin(best_iteration)),2)
    best_nbkg = round(best_nbkg_h.GetBinContent(best_nbkg_h.FindBin(best_iteration)),2)
    best_cutid = best_cutid_h.GetBinContent(best_cutid_h.FindBin(best_iteration))
    zcuts.append(best_zcut)
    nsigs.append(best_nsig)
    nbkgs.append(best_nbkg)
    best_cutids.append(best_cutid)
    print("Best Zcut: ", best_zcut)
    print("Best nsig: ", best_nsig)
    print("Best nbkg: ", best_nbkg)

    #Get persistent cuts
    openfile = r.TFile(infile)
    pers_cuts_hh = openfile.Get("zbi_processor_persistent_cuts_hh")
    best_pers_cuts = get_test_cut_values_for_iter(pers_cuts_hh,best_iteration)
    print("Persistent Cuts: ", best_pers_cuts)
    cuts[mass] = best_pers_cuts
    best_cut = best_pers_cuts[best_cutid]
    best_cut_name = best_cut[0]
    best_cut_value = best_cut[1]
    print("Best Cut Value: ", best_cut_value)
    best_cut_values.append(best_cut_value)

    text = ["Cut %s %s"%(best_cut_name, best_cut_value), "Zcut=%s | Nsig=%s | Nbkg=%s"%(str(best_zcut), str(best_nsig), str(best_nbkg))]

    #Read best cut background model
    print("looking for ", "background_zVtx_%s_h"%(best_cut_name))
    print("testCuts_pct_sig_cut_%s.000000"%(str(best_iteration)))
    best_bkg_h = utils.read_plot_from_root_file(infile,"testCutHistos_background_zVtx_%s_h"%(best_cut_name),
            "testCuts_pct_sig_cut_%s.000000"%(str(best_iteration)))
    best_bkg_h.SetName('%s_MeV_bkg_model_%s'%(mass, best_cut_name))
    best_bkg_h.SetTitle('%s_MeV_bkg_model_%s'%(mass, best_cut_name))
    utils.format_TH1(best_bkg_h, line_width=2,line_color = colors[1])
    func_name = best_bkg_h.GetListOfFunctions().At(0).GetName()
    func = best_bkg_h.GetFunction("%s"%(func_name))
    func.SetLineColor(colors[1])

    #Read initial signal distribution
    scaleF = 10
    initbkg_h = utils.read_plot_from_root_file(infile,"background_unc_vtx_z_h",
                        "initial_background")
    initbkg_h.SetName('%s_MeV_initial_background_scaleF_%i'%(mass, scaleF))
    initbkg_h.SetTitle('%s_MeV_initial_background_scaleF_%i'%(mass,scaleF))
    initbkg_h.Scale(scaleF)
    utils.format_TH1(initbkg_h, line_width=2,line_color = colors[2])

    #Read initial signal distribution
    initsig_h = utils.read_plot_from_root_file(infile,"signal_unc_vtx_z_h",
                        "initial_signal")
    initsig_h.SetName('%s_MeV_initial_signal_unscaled'%(mass))
    initsig_h.SetTitle('%s_MeV_initial_signal_unscaled'%(mass))
    utils.format_TH1(initsig_h, line_width=2,line_color = colors[3])

    #Read signal unc_vtx corresponding to this cut
    print("looking for ", "signal_pct_sig_cut_%s.000000"%(str(best_iteration)))
    sig_h = utils.read_plot_from_root_file(infile,"signal_unc_vtx_z_h",
                        "signal_pct_sig_cut_%s.000000"%(str(best_iteration)))
    sig_h.SetName('%s_MeV_signal_unc_vtx_z_unscaled'%(mass))
    sig_h.SetTitle('%s_MeV_signal_unc_vtx_z_unscaled'%(mass))
    utils.format_TH1(sig_h, line_width=2,line_color = colors[4])

    #bkg_model_c = utils.Make1DPlot('%s_MeV_bkg_model_%s'%(mass, best_cut_name), best_bkg_h, 
    bkg_model_c = utils.plot_TH1s_with_legend([initbkg_h, best_bkg_h,initsig_h, sig_h], '%s_MeV_bkg_model_%s'%(mass, best_cut_name),plots_dir, insertText=text, LogY=True)


zcuts = [x for _, x in sorted(zip(masses, zcuts))]
nsigs = [x for _, x in sorted(zip(masses, nsigs))]
nbkgs = [x for _, x in sorted(zip(masses, nbkgs))]
zbis = [x for _, x in sorted(zip(masses, zbis))]
best_cut_values = [x for _, x in sorted(zip(masses, best_cut_values))]
print(masses, best_cut_values)
masses = np.sort(np.array(masses,dtype=float))

#Zcut
zcut_gr = r.TGraph(len(zcuts),np.array(masses, dtype=float),np.array(zcuts,dtype=float))
utils.format_TH1(zcut_gr, name='zcut_positions_0pt5_bkg', title='zcut_positions_0pt5_bkg', x_label='mass [MeV]', 
        y_label='zcut [mm]', line_width=2, line_color = 2, marker_style = 8)
fitFunc = r.TF1("fitFunc", "[0] + [1]*x +[2]*x*x + [3]*x*x*x")
fitFunc.SetParameters(1.0,1.0,1.0,1.0)
zcut_gr.Fit(fitFunc,"SIM","",40.0,160.0)
print(fitFunc.GetParameter(0))
params = [0.0] * 4
for i in range(len(params)):
    params[i] = fitFunc.GetParameter(i)
param_text = [f"Fit Parameters:\n[0] = {params[0]}\n [1] = {params[1]}\n [2] = {params[2]}\n [3] = {params[3]}"]
zcut_gr_c = utils.Make1DPlot('zcut_positions_0pt5_bkg', zcut_gr, outdir = plots_dir,insertText=param_text)
print("ZCUT VALUES FOR EACH MASS POINT")
print("masses: ", masses)
print("zcuts: ", zcuts)

nsig_gr = r.TGraph(len(nsigs),np.array(masses, dtype=float),np.array(nsigs,dtype=float))
utils.format_TH1(zcut_gr, name='nsig', title='nsig', x_label='mass [MeV]', 
        y_label='zcut [mm]', line_width=2, line_color = 4, marker_style = 8)
nsig_gr_c = utils.Make1DPlot('nsig_positions_0pt5_bkg', nsig_gr, outdir = plots_dir)

nbkg_gr = r.TGraph(len(nbkgs),np.array(masses, dtype=float),np.array(nbkgs,dtype=float))
utils.format_TH1(zcut_gr, name='nbkg', title='nbkg', x_label='mass [MeV]', 
        y_label='zcut [mm]', line_width=2, line_color = 6, marker_style = 8)
nbkg_gr_c = utils.Make1DPlot('nbkg_positions_0pt5_bkg', nbkg_gr, outdir = plots_dir)

cuts_by_name = {}
for mass, cut_ids in cuts.items():
    for cutid, (name, value) in cut_ids.items():
        cuts_by_name.setdefault(name,[]).append((mass,value))
print('cuts_by_name:', cuts_by_name)
for cut, values in cuts_by_name.items():
    xvals = []
    yvals = []
    for pair in values:
        print('mass:', pair[0])
        print('val: ', pair[1])
        xvals.append(float(pair[0]))
        yvals.append(pair[1])

    yvals = [x for _, x in sorted(zip(xvals, yvals))]
    xvals = np.sort(xvals)
    print('xvals:', xvals)
    print('yvals: ', yvals)
    gr = r.TGraph(len(xvals),np.array(xvals, dtype=float),np.array(yvals, dtype=float))
    utils.format_TH1(gr, name='%s'%(cut), title='%s'%(cut), x_label='mass [MeV]', 
            y_label='cut value', line_width=2, line_color = 6, marker_style = 8)
    gr_c = utils.Make1DPlot('%s_vs_mass'%(cut), gr, outdir = plots_dir)

#cut_gr = r.TGraph(len(best_cut_values),np.array(masses, dtype=float),np.array(best_cut_values,dtype=float))
#utils.format_TH1(zcut_gr, name='cut values', title='cut values', x_label='mass [MeV]', 
#        y_label='cut value', line_width=2, line_color = 6, marker_style = 8)
#nbkg_gr_c = utils.Make1DPlot('nbkg_positions_0pt5_bkg', nbkg_gr, outdir = plots_dir)

#matplotlib plots
fig, ax1 = plt.subplots()
ax1.set_xlabel('mass [MeV]')
ax1.set_ylabel('zcut [mm]', color='blue')
ax1.plot(masses, zcuts, color='blue', label='zcut')
ax1.tick_params(axis='y', labelcolor='blue')

ax2 = ax1.twinx()
ax2.set_ylabel('Nsig & Nbkg', color='red')
ax2.plot(masses, nsigs, color='red', label='Nsig')
ax2.plot(masses, nbkgs, color='green', label='Nbkg')
ax2.tick_params(axis='y', labelcolor='red')

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
lines = lines1 + lines2
labels = labels1 + labels2
ax1.legend(lines, labels, loc='upper left')

plt.title('Invariant Mass Zcut Positions')
plt.savefig('%s/invariant_mass_zcuts.png'%(plots_dir),dpi=300)

        

