#!/usr/bin/python3
import ROOT as r
import numpy as np
from array import array
import sys
sys.path.append( '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/plot_utilities')
import plot_utilities as utils
import json
from tabulate import tabulate

#Integrate Histograms
def integrateHisto(histo, cutname, cutval, cut_type):

    total_integral = histo.Integral(0, histo.GetNbinsX()+1)
    firstbin = -99999.9
    lastbin = 99999.9

    if cut_type == 'abs':
        firstbin = histo.FindBin(-cutval)
        lastbin = histo.FindBin(cutval)
        if firstbin == 0:
            firstbin = 1
        if lastbin > histo.GetNbinsX():
            lastbin = histo.GetNbinsX()

    elif cut_type == 'lt':
        firstbin = 0
        lastbin = int(histo.FindBin(cutval))
        if lastbin > histo.GetNbinsX():
            lastbin = histo.GetNbinsX()

    elif cut_type == 'gt':
        firstbin = int(histo.FindBin(cutval))
        lastbin = int(histo.GetNbinsX()+1)
        if firstbin == 0:
            firstbin = 1 


    print('firstbinval:',histo.GetBinCenter(firstbin))
    print('lastbinval:',histo.GetBinCenter(lastbin))
    print('firstbin:',firstbin)
    print('lastbin:',lastbin)
    integral = round(float(histo.Integral(firstbin, lastbin)),2)
    nentries = float(histo.GetEntries())
    eff = round(integral/total_integral,3)
    print('integral',integral, 'entries', nentries, 'eff', eff)

    return integral, eff
    

r.gROOT.SetBatch(1)
style = utils.SetMyStyle()
colors = utils.getColorsHPS()

#r.gROOT.SetStyle("myStyle")
#r.gROOT.ForceStyle()

#data
infile = '../data/hadd_hps_BLPass4_10pct_preselection_nminus1_20230914.root'
#tritrig
tritrig = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/tritrig_beam/hadd_tritrig_beam_preselection_20230915.root'
#wab
wab = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/wab_beam/hadd_wab_beam_preselection_20230915.root'

#signal 60MeV
sig_60 = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/signal/hadd_signal_60_MeV_beam_preselection_20230915.root'
sig_100 = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/signal/hadd_signal_60_MeV_beam_preselection_20230915.root'

#MC Scaling
Lumi = 0.107 #1/pb
mcScale = {}
mcScale['tritrig'] = 1.416e9*Lumi/(50000*9456) #pb2016
mcScale['wab'] = 0.1985e12*Lumi/(100000*9688) #pb2016


#selection
preselection = '../configs/vertexSelection_2016_simp_preselection.json'

#Define cut variables for each cut in json
cuts = {'eleTrkTime_lt': 'ele_time',
        'posTrkTime_lt': 'pos_time',
        'eleposCluTimeDiff_lt': 'ele_pos_clusTimeDiff',
        'eleTrkCluTimeDiff_lt': 'ele_track_clus_dt',
        'posTrkCluTimeDiff_lt': 'pos_track_clus_dt',
        'eleTrkChi2Ndf_lt': 'ele_chi2ndf',
        'posTrkChi2Ndf_lt': 'pos_chi2ndf',
        'eleMom_lt': 'ele_p',
        'eleMom_gt': 'ele_p',
        'posMom_gt': 'pos_p',
        'eleN2Dhits_gt': 'ele_track_n2dhits',
        'posN2Dhits_gt': 'pos_track_n2dhits',
        'maxVtxMom_lt': 'vtx_Psum',
        'chi2unc_lt':'vtx_chi2'
        }

latex_cuts = {'eleTrkTime_lt': '$|e^{-} Track_{t}| < ',
        'posTrkTime_lt': '$|e^{+} Track_{t}| <',
        'eleposCluTimeDiff_lt': '$\Delta_{t}(cluster_{e^{-}},cluster_{e^{+}} < ',
        'eleTrkCluTimeDiff_lt': '$e^{-}\Delta_{t}(track,cluster) < ',
        'posTrkCluTimeDiff_lt': '$e^{+}\Delta_{t}(track,cluster) < ',
        'eleTrkChi2Ndf_lt': '$e^{-} Track \Chi^2/n.d.f. < ',
        'posTrkChi2Ndf_lt': '$e^{+} Track \Chi^2/n.d.f. < ',
        'eleMom_lt': '$p_{e^-} < ',
        'eleMom_gt': '$p_{e^-} > ',
        'posMom_gt': '$p_{e^+} > ',
        'eleN2Dhits_gt': '$N_{2d hits} on e^{-}_{Track} >= ',
        'posN2Dhits_gt': '$N_{2d hits} on e^{+}_{Track} >= ',
        'maxVtxMom_lt': '$p_{e^{-}+e^{+}} < ',
        'chi2unc_lt':'$Vtx_{\Chi^2} < '
        }

#Store values for latex table
table_entries = {}
row_labels = []

f = open(preselection)
sel = json.load(f)
i = 0
for cut, info in sel.items():
    cutval = float(info['cut'])
    print(cut)
    print('cutval: ', cutval)
    directory = 'vtxana_kf_%s_selection'%(cut)

    cut_type = 'null'
    if 'Time' in cut:
        cut_type = 'abs'
    elif '_lt' in cut:
        cut_type = 'lt'
    elif '_gt' in cut:
        cut_type = 'gt'
 
    #skip pair1
    if cut == 'Pair1_eq':
        continue

    i = i +1
    #if i > 2:
    #    continue

    #Format row labels for latex
    latex_label = latex_cuts[cut]
    #Append cut value
    latex_label = latex_label+'%s$'%(str(cutval))
    row_labels.append(latex_label)

    #data
    data_plot = utils.read_1d_plots_from_root_file(infile, directory, '%s_h'%(cuts[cut]))[0]
    utils.formatHisto(data_plot, line_color=colors[0],name='Data', title='1% data')

    #tritrig
    tritrig_plot = utils.read_1d_plots_from_root_file(tritrig, directory, '%s_h'%(cuts[cut]))[0]
    utils.formatHisto(tritrig_plot, line_color=colors[1],name='Tritrig-Beam', title='Tritrig+Beam')

    #wab
    wab_plot = utils.read_1d_plots_from_root_file(wab, directory, '%s_h'%(cuts[cut]))[0]
    utils.formatHisto(wab_plot, line_color=colors[2],name='WAB-Beam', title='WAB+Beam')

    #Signal
    sig60_plot = utils.read_1d_plots_from_root_file(sig_60, directory, '%s_h'%(cuts[cut]))[0]
    utils.formatHisto(sig60_plot, line_color=colors[3],name='Signal-60 MeV', title='Signal_60MeV')
    sig100_plot = utils.read_1d_plots_from_root_file(sig_100, directory, '%s_h'%(cuts[cut]))[0]
    utils.formatHisto(sig100_plot, line_color=colors[4],name='Signal-100 MeV', title='Signal_100MeV')

    #Scale tritrig+wab+beam
    #tritrig_plot.Scale(mcScale['tritrig'])
    #wab_plot.Scale(mcScale['wab'])
    mc_bkg_plot = tritrig_plot.Clone()
    mc_bkg_plot.Scale(mcScale['tritrig'])
    wab_clone_plot = wab_plot.Clone()
    wab_clone_plot.Scale(mcScale['wab'])
    mc_bkg_plot.Add(wab_clone_plot)
    utils.formatHisto(mc_bkg_plot, line_color=colors[5],name='Tritrig-WAB-Beam', title='tritri+wab+beam')

    #collect plots
    plots = [data_plot, tritrig_plot, wab_plot, mc_bkg_plot, sig60_plot, sig100_plot]

    #Integrate histogram
    text = []
    for plot in plots:
        name = plot.GetName()
        integral, eff = integrateHisto(plot, cut, cutval, cut_type)
        text.append('%s integral: %s | eff: %s'%(plot.GetTitle(),integral, eff))
        if name in table_entries:
            table_entries[name][0].append(integral)
            table_entries[name][1].append(eff)
        else:
            table_entries[name] = [ [integral], [eff] ]


    #Scale Histograms
    for plot in plots:
        plot.Scale(1.0/float(plot.Integral(-1, plot.GetNbinsX()+1)))

    #Plots
    logY = True
    freezeX = True
    no_logy = ['N2Dhits','TrkCluTimeDiff']
    #if any(substring in cut for substring in no_logy):
    #    logY = False
    no_freezeX = ['TrkCluTimeDiff','eleposCluTimeDiff_lt']
    if any(substring in cut for substring in no_freezeX):
        freezeX = False

    #c = utils.plot_TH1s_with_legend(plots, cut, insertText=text, text_x=0.6, text_y=0.9,text_size=0.02, 
            #freezeXaxis=freezeX, LogY=logY, save=False)
    c = utils.plot_TH1s_with_legend(plots, cut, freezeXaxis=freezeX, LogY=logY, save=False)
    c.Range( 0., -10., 1., 10. )
    line = r.TLine(cutval, c.GetUymin(), cutval, c.GetUymax())
    line.SetLineColor(2)  # Set the line color (optional)
    line.SetLineWidth(3)
    line.Draw("same")
    if cut_type == 'abs':
        line2 = r.TLine(-cutval, c.GetUymin(), -cutval, c.GetUymax())
        line2.SetLineColor(2)  # Set the line color (optional)
        line2.SetLineWidth(3)
        line2.Draw("same")

    c.SaveAs('../plots/%s.png'%(cut))
    c.Close()

#Make latex table
table = []
headers = []
for key in table_entries.keys():
    headers.extend([f"{key} Integral", f"{key} Eff"])
table.append(headers)

for i, row_label in enumerate(row_labels):
    row = [row_label]
    for key, values in table_entries.items():
        for value in values:
            row.append(value[i])
    table.append(row)

# Generate the LaTeX table
latex_table = tabulate(table, headers="firstrow", tablefmt="latex_raw")
print(latex_table)

print('\n Latex Table Format:\n')

latex_table = tabulate(table, headers="firstrow", tablefmt="latex")
print(latex_table)

#lTable = utils.makeLatexTable(table)







