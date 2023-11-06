#!/usr/bin/python3
import sys
import new_plot_utilities as utils
import ROOT as r
import os
import argparse
import numpy as np

################################################################
def old_get_test_cut_values_for_iter(cuts_hh, iteration):
    iter_bin = cuts_hh.GetXaxis().FindBin(iteration)
    projy = cuts_hh.ProjectionY("iteration_%i_cuts"%(iteration),iter_bin,iter_bin,"")
    cuts_values = {}
    for x in range(projy.GetNbinsX()):
        if projy.GetXaxis().GetBinLabel(x+1) == "":
            continue
        cut = projy.GetXaxis().GetBinLabel(x+1)
        val = round(projy.GetBinContent(x+1),2)
        cuts_values[cut] = val
    return cuts_values

def get_test_cut_values_for_iter(cuts_hh, iteration):
    iter_bin = cuts_hh.GetXaxis().FindBin(iteration)
    cuts_values = {}
    for y in range(cuts_hh.GetNbinsY()):
        if cuts_hh.GetYaxis().GetBinLabel(y+1) == "":
            continue
        cut = cuts_hh.GetYaxis().GetBinLabel(y+1)
        val = round(cuts_hh.GetBinContent(iter_bin,y+1),2)
        cuts_values[cut] = val
    return cuts_values


################################################################
parser = argparse.ArgumentParser(description="HTML Config")
parser.add_argument('--infileName', type=str, dest="infileName", default = "")
parser.add_argument('--plotsDir', type=str, dest="plotsDir", default = "plots_1d")
options = parser.parse_args()

style = utils.SetMyStyle(setOptTitle=1, setOptStat=1)
r.gROOT.SetBatch(1)
colors = utils.getColors()

infileName = options.infileName
print("Input File: ", infileName)
plots_dir = options.plotsDir
os.makedirs(plots_dir, exist_ok=True)
print("Created directory", plots_dir)

plots = utils.read_1d_plots_from_root_file_dirs(infileName,"testCuts","background_zVtx")

#Get the test cuts
pers_cuts_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_persistent_cuts_hh")
#Zbis
zbis_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_best_test_cut_ZBi_hh")

#Test Cuts
testcuts_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_test_cuts_values_hh")
#Test Cut ZBis
testcuts_zbi_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_test_cuts_ZBi_hh")
#Test Cut Zcut
testcuts_zcut_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_test_cuts_zcut_hh")
#Test Cut Nsig
testcuts_nsig_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_test_cuts_nsig_hh")
#Test Cut nbkg
testcuts_nbkg_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_test_cuts_nbkg_hh")

# Open the ROOT file
root_file = r.TFile.Open(infileName,"READ")

iterations = []
nsigs = []
nbkgs = []
zcuts = []
zbis = []
for key, values in plots.items():
    iteration = int(key.split('_')[-1].split('.')[0])
    iterations.append(iteration)
    test_cuts = get_test_cut_values_for_iter(testcuts_hh, iteration)
    test_cuts_zbi = get_test_cut_values_for_iter(testcuts_zbi_hh, iteration)
    test_cuts_zcut = get_test_cut_values_for_iter(testcuts_zcut_hh, iteration)
    test_cuts_nsig = get_test_cut_values_for_iter(testcuts_nsig_hh, iteration)
    test_cuts_nbkg = get_test_cut_values_for_iter(testcuts_nbkg_hh, iteration)
    text = []
    best_zbi = -9999.9
    best_cut = ''
    for key, value in test_cuts.items():
        zbi = str(test_cuts_zbi[key])
        zcut = str(test_cuts_zcut[key])
        nsig = str(test_cuts_nsig[key])
        nbkg = str(test_cuts_nbkg[key])
        line = "Test Cut %s:%s | %s [zbi] | zcut %s [mm]. Nsig=%s | Nbkg=%s"%(key, value, zbi, zcut, nsig, nbkg)
        text.append(line)
        if float(zbi) > best_zbi:
            best_cut = key
    nsigs.append(test_cuts_nsig[best_cut])
    nbkgs.append(test_cuts_nbkg[best_cut])
    zcuts.append(test_cuts_zcut[best_cut])
    zbis.append(test_cuts_zbi[best_cut])

    histos = []
    for n, value in enumerate(values):
        h = value
        name = h.GetName()
        var_name = name.replace('testCutHistos_background_zVtx','').replace('_h','')
        h.SetName(var_name)
        h.SetTitle('iter_%s%s'%(iteration,var_name))
        utils.formatHisto(h,line_width = 3, line_color = colors[n])
        func_name = h.GetListOfFunctions().At(0).GetName()
        func = h.GetFunction("%s"%(func_name))
        func.SetLineColor(colors[n])
        func.SetLineWidth(3)
        func.Draw()
        histos.append(h)

    canvas_name = 'iteration_%s_test_cuts'%(str(iteration)) 
    canvas = utils.plotTH1s(histos, canvas_name,LogY=True, xmin=-30.0, xmax=80.0,ymin=0.1)
    legend = utils.makeLegend(canvas, histos)
    textbox = utils.makeLatexBox(canvas, text, text_x=0.2, text_size=0.025)
    canvas.SaveAs('%s/%s.png'%(plots_dir,canvas_name))

nsigs_g = utils.makeTGraph('Nsig', np.array(iterations,dtype=float), np.array(nsigs,dtype=float),title='Nsig',xlabel='Iteration', ylabel='Nsig', lineColor=colors[0], lineWidth=3)
canvas = utils.plotTGraphs([nsigs_g],'Nsig')
canvas.SaveAs('%s/%s.png'%(plots_dir,'a0_Nsig'))

nbkgs_g = utils.makeTGraph('Nbkg', np.array(iterations,dtype=float), np.array(nbkgs,dtype=float),title='Nbkg',xlabel='Iteration',ylabel='Nbkg',lineColor=colors[1], lineWidth=3)
canvas = utils.plotTGraphs([nbkgs_g],'Nbkg')
canvas.SaveAs('%s/%s.png'%(plots_dir,'a1_Nbkg'))

zcuts_g = utils.makeTGraph('Zcut', np.array(iterations,dtype=float), np.array(zcuts,dtype=float),title='Zcut',xlabel='Iteration',lineColor=colors[2],lineWidth=3)
canvas = utils.plotTGraphs([zcuts_g],'Zcut')
canvas.SaveAs('%s/%s.png'%(plots_dir,'a2_Zcut'))

zbis_g = utils.makeTGraph('ZBi', np.array(iterations,dtype=float), np.array(zbis,dtype=float),title='ZBi',xlabel='Iteration', ylabel='ZBi', lineColor=colors[3], lineWidth=3)
canvas = utils.plotTGraphs([zbis_g],'ZBi')
canvas.SaveAs('%s/%s.png'%(plots_dir,'a3_ZBi'))

ymin = min(min(nsigs), min(nbkgs), min(zcuts), min(zbis))
ymax = max(max(nsigs),max(nbkgs), max(zcuts), max(zbis))
nsigs_g.GetYaxis().SetRangeUser(0.1,ymax)
#nsigs_g.GetYaxis().SetNdivisions(500)
graphs = [nsigs_g, nbkgs_g, zcuts_g, zbis_g]
canvas = utils.plotTGraphs(graphs, 'Iteration Results', LogY=True)
canvas.SetGrid()
canvas.SetGridx()
legend = utils.makeLegend(canvas, graphs)
canvas.SaveAs('%s/%s.png'%(plots_dir,'a4_iteration_results'))




