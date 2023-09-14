#!/usr/bin/python3

import ROOT as r
import numpy as np
import glob as glob
import os
import copy as copy

def getStatsYPositions(nplots):
    start = 1.0
    end = 0.1
    height = (start-end)/nplots
    if height > 0.1:
        height = 0.1
    #ypositions = list(range(start,end,height))
    ypositions = np.arange(end,start,height)
    ypositions=np.flip(ypositions)
    print(ypositions)
    return ypositions,height

def setStatsBox(plot,ypos,height,linecolor=r.kBlack,color=r.kWhite,alpha=0.8):
    r.gPad.Update()
    stats = plot.GetListOfFunctions().FindObject("stats")
    stats.SetFillColorAlpha(color,alpha)
    stats.SetLineColor(linecolor)
    r.gStyle.SetStatY(ypos)
    r.gStyle.SetStatH(height)
    r.gStyle.SetStatW(0.1)
    r.gStyle.SetStatX(0.9)
    r.gPad.Update()

def buildLegend(canvas,x1=0.60, y1=0.45,x2=0.80,y2=0.7,textsize=0.02,separation=0.1,fill=0,border=0):
    #canvas.cd()
    legend = canvas.BuildLegend(x1,y1,x2,y2)
    legend.Draw()
    legend.SetTextSize(textsize)
    #legend.SetEntrySeparation(separation)
    if(fill==0):
        legend.SetFillStyle(fill)
    else:
        legend.SetFillColorAlpha(r.kGray,0.5)
    legend.SetBorderSize(border)

def readGraphsFromRootDir(infile, directory, keyword="none"):
    graphs = {}
    rootdir = infile.Get(directory)
    rootdir.cd()
    for key in rootdir.GetListOfKeys():
        if "TGraph" not in key.GetClassName():
            continue
        if keyword != "none":
            if keyword not in key.GetName():
                continue
        graph = rootdir.Get(key.GetName())
        name = graph.GetName().replace(" ","_").replace(".","")
        graph.SetName(name)
        graphs[name] = graph
    return graphs

def read1DPlotsFromRootDir(infile, directory, keyword="none"):
    plots = {}
    rootdir = infile.Get(directory)
    rootdir.cd()
    for key in rootdir.GetListOfKeys():
        if "TH1" not in key.GetClassName():
            continue
        if keyword != "none":
            if keyword not in key.GetName():
                continue
        plot = rootdir.Get(key.GetName())
        name = plot.GetName().replace(" ","_").replace(".","")
        plot.SetName(name)
        plots[name] = plot

    return plots

def format1DPlot(plot, name, title = "",linecolor = 1, linewidth = 2, linestyle = 1):
    plot.SetLineWidth(linewidth)
    plot.SetLineColor(linecolor)
    plot.SetLineStyle(linestyle)
    plot.SetName(name)
    if title != "":
        plot.SetTitle(title)

def overlay1DGraphs(graphs,canvas_name, graphs_dir, logY=False, drawLegend = True, legx1=0.60, legx2=0.8,legy1=0.45,legy2=0.7,drawStats = True, setXmin='',setXmax='',setYmin='',setYmax='',draw_options=''):

    ymax = -99999
    ymin = 0.00001

    if setYmin != '':
        ymin = setYmin
    if setYmax != '':
        ymax = setYmax

    c = r.TCanvas("%s"%(canvas_name),"%s"%(canvas_name),2560,1440)
    c.cd()

    #set yaxis
    for i,graph in enumerate(graphs):
        graph.GetYaxis().SetRangeUser(ymin,ymax*1.1)
        if setXmin != '' and setXmax != '':
            graph.GetXaxis().SetRangeUser(setXmin,setXmax)
        if i < 1:
            graph.Draw("%s"%(draw_options))
        else:
            graph.Draw("%ssame"%(draw_options))

    if(drawLegend):
        buildLegend(c,x1=legx1,y1=legy1,x2=legx2,y2=legy2)

    if logY:
        c.SetLogy()

    graphs[0].SetTitle(canvas_name)
    c.SaveAs("./%s/%s.png"%(graphs_dir,canvas_name))


def overlay1DPlots(plots,canvas_name, plots_dir, logY=False, drawLegend = True, drawStats = True, setXmin='',setXmax='',setYmin='',setYmax='',draw_options=''):
    #statsbox formatting
    statspositions,height = getStatsYPositions(len(plots))

    ymax = -99999
    ymin = 0.00001

    if setYmin != '':
        ymin = setYmin
    if setYmax != '':
        ymax = setYmax
    else:
        for i,plot in enumerate(plots):
            ymaxi = plot.GetMaximum()
            if ymaxi > ymax:
                ymax = ymaxi

    c = r.TCanvas("%s"%(canvas_name),"%s"%(canvas_name),2560,1440)
    c.cd()

    #set yaxis
    for i,plot in enumerate(plots):
        plot.GetYaxis().SetRangeUser(ymin,ymax*1.1)
        if setXmin != '' and setXmax != '':
            plot.GetXaxis().SetRangeUser(setXmin,setXmax)
        if i < 1:
            plot.Draw("%s"%(draw_options))
            if(drawStats):
                setStatsBox(plot,statspositions[i],height, linecolor=colors[i])
        else:
            plot.Draw("%ssames"%(draw_options))
            if(drawStats):
                setStatsBox(plot,statspositions[i],height,linecolor=colors[i])

    if(drawLegend):
        buildLegend(c)

    if logY:
        c.SetLogy()
    plots[0].SetTitle(canvas_name)
    c.SaveAs("./%s/%s.png"%(plots_dir,canvas_name))

def getColors():
    colors = [r.kBlue, r.kGreen, r.kOrange, r.kRed, r.kMagenta, r.kCyan, r.kPink+1, r.kBlack, r.kViolet+2, r.kTeal-1, r.kOrange+7, r.kMagenta-3, r.kBlue+2 ,r.kPink-9]
    return colors

colors = getColors()
#markers = {'0.005':2, '0.01':3,'0.015':4,'0.02':5,'0.025':24,'0.03':25,'0.035':26,'0.04':27,'0.045':28,'0.05':30}
#colors = {'0.005':1,'0.01':4,'0.02':6}

#colors = {'0.005':r.kBlue, '0.01':r.kGreen,'0.015':r.kOrange,'0.02':r.kRed,'0.025':r.kMagenta,'0.03':r.kBlack,'0.035':r.kViolet+2,'0.04':r.kTeal-1,'0.045':r.kBlue+2,'0.05':r.kPink-9}
colors = {'0.012':r.kBlue, '0.014':r.kGreen,'0.016':r.kOrange,'0.018':r.kRed,'0.022':r.kMagenta,'0.024':r.kBlack,'0.026':r.kViolet+2,'0.028':r.kTeal-1,'0.02':r.kBlue+2,'0.03':r.kPink-9}
markers = {'0.005':34,'0.01':34,'0.02':4}

#############################################################
r.gROOT.SetBatch(1)
#Get best zbi for each iteration
cut_map = {2:'ele_track_p_lt',3:'ele_track_p_gt',5:'pos_track_p_lt',6:'pos_track_p_gt',7:'vtx_psum_gt',8:'vtx_psum_lt',9:'ele_clust_E_lt',10:'pos_clust_E_lt',11:'ele_track_zalpha_lt',12:'pos_track_zalpha_lt'}

plots_hh = []
for file in sorted(glob.glob('/sdf/group/hps/users/alspellm/run/cut_dev/zalpha_slope_scan_05032023/zalpha_slope_scan_05082023_results/*.root')):
    base_filename = os.path.basename(file).replace('.root','')
    print(file)
    slope = base_filename.split('_')[7]
    step_size = base_filename.split('_')[10]
    print(base_filename, slope, step_size)

    infile = r.TFile(file,"READ")
    best_test_cut_ZBi_hh = copy.deepcopy(infile.Get('zbi_processor_best_test_cut_ZBi_hh'))
    #best_test_cut_ZBi_hh.RebinY(2)
    best_test_cut_ZBi_hh.SetName('slope_%s_step_%s'%(slope,step_size))
    best_test_cut_ZBi_hh.SetTitle('slope_%s_step_%s;percent signal cut;ZBi'%(slope,step_size))
    best_test_cut_ZBi_hh.SetMarkerStyle(markers[step_size])
    best_test_cut_ZBi_hh.SetMarkerColor(colors[slope])
    plots_hh.append(best_test_cut_ZBi_hh)

canvas = r.TCanvas('zbravoalpha_slope_zbi_comparisons','zalphabravo_slope_zbi_comparisons',2560,1440)
canvas.cd()

for i,plot in enumerate(plots_hh):
    plot.GetYaxis().SetRangeUser(0.1,6.0)
    if i < 1:
        plot.Draw()
    else:
        plot.Draw("same")
buildLegend(canvas, x1=0.80,y1=0.2,x2=0.95,y2=0.9,textsize=0.02,fill=1)
plots_hh[0].SetTitle('zbravoalpha_slopes_zbi')
canvas.SaveAs('zalpha_slope_comparisons_Zbi_05082023.png')

