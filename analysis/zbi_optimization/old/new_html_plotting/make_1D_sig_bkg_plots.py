#!/usr/bin/python3
import sys
import new_plot_utilities as utils
import ROOT as r
import os
import argparse

################################################################
def old_get_test_cut_values_for_iter(cuts_hh, iteration):
    iter_bin = cuts_hh.GetXaxis().FindBin(iteration)
    projy = cuts_hh.ProjectionY("iteration_%i_cuts"%(iteration),iter_bin,iter_bin,"")
    cuts_values = {}
    for x in range(projy.GetNbinsX()):
        if projy.GetXaxis().GetBinLabel(x) == "":
            continue
        cut = projy.GetXaxis().GetBinLabel(x)
        print("cut:",cut)
        print("value;",value)
        print("bin:",x)
        val = round(projy.GetBinContent(x),2)
        cuts_values[cut] = val
    return cuts_values

def get_test_cut_values_for_iter(cuts_hh, iteration):
    iter_bin = cuts_hh.GetXaxis().FindBin(iteration)
    cuts_values = {}
    for y in range(cuts_hh.GetNbinsY()):
        if cuts_hh.GetYaxis().GetBinLabel(y+1) == "":
            continue
        cut = cuts_hh.GetYaxis().GetBinLabel(y+1)
        val = round(cuts_hh.GetBinContent(iter_bin,y+1))
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

signal_1d_plots = utils.read_1d_plots_from_root_file_dirs(infileName,"signal","_h")
background_1d_plots = utils.read_1d_plots_from_root_file_dirs(infileName,"background","_h")
#print(signal_1d_plots.keys())

#Get persistent cuts
pers_cuts_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_persistent_cuts_hh")

#Get best cut ZBi from previous iteration
zbis_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_best_test_cut_ZBi_hh")

#Test Cut Values from previous iteration
testcuts_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_test_cuts_values_hh")

#Test Cut ZBis from previous iteration
testcuts_zbi_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_test_cuts_ZBi_hh")

#Test Cut Zcut
testcuts_zcut_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_test_cuts_zcut_hh")

#Test Cut Nsig
testcuts_nsig_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_test_cuts_nsig_hh")

#Test Cut nbkg
testcuts_nbkg_hh = utils.read_plot_from_root_file(infileName, "zbi_processor_test_cuts_nbkg_hh")


#sys.exit()
for key, values in signal_1d_plots.items():
    if "initial" in key:
        iteration = "initial"
    else:
        iteration = int(key.split('_')[-1].split('.')[0])

    text = []
    if iteration != "initial":
        pers_cuts = get_test_cut_values_for_iter(pers_cuts_hh, iteration)
        for cut_key, value in pers_cuts.items():
           line = "Active:%s %s"%(cut_key, str(value))
           text.append(line) 
        if iteration > 0:
            #Read test cut values from previous iteration, which correspond to the sig/bkg plots for this iteration
            #test_cuts = get_test_cut_values_for_iter(testcuts_hh, iteration-1)
            test_cuts_zbi = get_test_cut_values_for_iter(testcuts_zbi_hh, iteration-1)
            test_cuts_zcut = get_test_cut_values_for_iter(testcuts_zcut_hh, iteration-1)
            test_cuts_nsig = get_test_cut_values_for_iter(testcuts_nsig_hh, iteration-1)
            test_cuts_nbkg = get_test_cut_values_for_iter(testcuts_nbkg_hh, iteration-1)
            max_zbi = -999.9
            max_zbi_cut = ""
            for zbi_key, value in test_cuts_zbi.items():
                if value > max_zbi:
                    max_zbi = value
                    max_zbi_cut = zbi_key
            zbi = float(test_cuts_zbi[max_zbi_cut])
            zcut = float(test_cuts_zcut[max_zbi_cut])
            nsig = float(test_cuts_nsig[max_zbi_cut])
            nbkg = float(test_cuts_nbkg[max_zbi_cut])
            #text.append("Iteration has peak zbi %s at zcut %s [mm]. Nsig=%s | Nbkg=%s"%(zbi, zcut, nsig, nbkg))
            text.append("ZBi:%f | Nsig:%f | Nbkg:%f"%(zbi, nsig, nbkg))
            text.append("Zcut: %f [mm]"%(zcut))

    for n,value in enumerate(values):
        signal_h = value
        if signal_h.GetEntries() < 1:
            continue
        signal_name = signal_h.GetName()
        var_name = signal_name.replace('signal_','')

        bkg_key = key.replace('signal','background')
        bkg_h = background_1d_plots[bkg_key][n]
        bkg_name = bkg_h.GetName()

        #format histos
        utils.formatHisto(signal_h,title='%s_%s'%(iteration,signal_h.GetTitle()),line_width = 3, line_color = colors[0])
        utils.formatHisto(bkg_h,title='%s_%s'%(iteration, bkg_h.GetTitle()),line_width = 3, line_color = colors[1])

        canvas_name = 'iteration_%s_%s'%(str(iteration), var_name) 
        print(canvas_name)
        canvas = utils.plotTH1s([signal_h, bkg_h], canvas_name, LogY=True) 
        legend = utils.makeLegend(canvas, [signal_h, bkg_h])
        latext = utils.makeLatexBox(canvas, text)
        canvas.SaveAs('%s/%s.png'%(plots_dir, canvas_name))

