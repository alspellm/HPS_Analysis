#!/usr/bin/python3

import ROOT as r
import numpy as np

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
    legend.SetFillStyle(fill)
    #legend.SetFillColorAlpha(r.kGray,0.7)
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

#############################################################

infilename = '/home/alic/HPS/projects/simps/zbi_opt/zalpha_opt/032423/exp_fit/zalpha_zbi_zcutscan_optimization_slope_0.025_step_size_0.01.root'

infile_1 = r.TFile('%s'%(infilename),"READ")

#Get best zbi for each iteration
cut_map = {3:'ele_track_p_lt',4:'ele_track_p_gt',5:'pos_track_p_lt',6:'pos_track_p_gt',7:'vtx_psum_gt',8:'vtx_psum_lt',9:'ele_clust_E_lt',10:'pos_clust_E_lt',11:'ele_track_zalpha_lt',12:'pos_track_zalpha_lt'}
best_zbi_hh = infile_1.Get("summary_best_test_cut_ZBi_hh")

r.gROOT.SetBatch(1)
#Loop over all iterations, make plots per iteration
step_size = float(infilename.replace('.root','').split('step_size_')[-1])
for iteration in range(1,99):
    pct_sig = 100*step_size*iteration
    print(pct_sig)

    #Get best zbi for iteration
    proj_bin = best_zbi_hh.GetXaxis().FindBin(pct_sig)
    print(proj_bin)
    best_cut_id = best_zbi_hh.ProjectionY("proj",proj_bin,proj_bin).GetMaximum()
    best_cut = cut_map[best_cut_id]

    #Plot background models
    bkg_plots = read1DPlotsFromRootDir(infile_1,"testCuts_pct_sig_cut_%i.000000"%(iteration),keyword="tritrig_zVtx")
    plots = []
    for i,plot in enumerate(bkg_plots.values()):
        format1DPlot(plot,plot.GetName().replace('cutHistos_tritrig_zVtx_',''),plot.GetName().replace('cutHistos_tritrig_zVtx_',''),linecolor=colors[i])
        #If this cut was the best ZBi for the iteration, set special marker
        if best_cut in plot.GetName():
            plot.SetMarkerStyle(34)
            plot.SetMarkerSize(3)
            plot.SetMarkerColor(colors[i])
        plot.GetFunction(plot.GetListOfFunctions()[0].GetName()).SetLineColor(colors[i])
        plots.append(plot)
    overlay1DPlots(plots,'tritrig_zVtx_pct_sig_cut_%i'%(pct_sig),'test_plots',setXmin=-5.0,setXmax=80.0,setYmin=1.0,logY=True)

    #Read graphs
    graphs = []
    zbi_graphs = readGraphsFromRootDir(infile_1, "testCuts_pct_sig_cut_%i.000000"%(pct_sig),keyword="zbi")
    #c = r.TCanvas('c_%i'%(i),'c_%i'%(i),2560,1440)
    #c.cd()
    for i,graph in enumerate(zbi_graphs.values()):
        print(graph.GetName())
        #graph.Draw()
        format1DPlot(graph, graph.GetName().replace('unc_vtx','').replace('_g','').replace('zcut_vs_',''),graph.GetName().replace('unc_vtx','').replace('_g','').replace('zcut_vs_',''),linecolor=colors[i], linestyle=3)
        #If this cut was the best ZBi for the iteration, set special marker
        if best_cut in graph.GetName():
            plot.SetMarkerStyle(34)
            plot.SetMarkerSize(3)
            plot.SetMarkerColor(colors[i])
        graphs.append(graph)


    #c.SaveAs('test.png')
    nsig_graphs = readGraphsFromRootDir(infile_1, "testCuts_pct_sig_cut_%i.000000"%(pct_sig),keyword="nsig")
    for i,graph in enumerate(nsig_graphs.values()):
        format1DPlot(graph, graph.GetName().replace('unc_vtx','').replace('_g','').replace('zcut_vs_',''),graph.GetName().replace('unc_vtx','').replace('_g','').replace('zcut_vs_',''),linecolor=colors[i], linestyle=1)
        #If this cut was the best ZBi for the iteration, set special marker
        if best_cut in graph.GetName():
            plot.SetMarkerStyle(34)
            plot.SetMarkerSize(3)
            plot.SetMarkerColor(colors[i])
        graphs.append(graph)

    nsig_graphs = readGraphsFromRootDir(infile_1, "testCuts_pct_sig_cut_%i.000000"%(pct_sig),keyword="nbkg")
    for i,graph in enumerate(nsig_graphs.values()):
        format1DPlot(graph, graph.GetName().replace('unc_vtx','').replace('_g','').replace('zcut_vs_',''),graph.GetName().replace('unc_vtx','').replace('_g','').replace('zcut_vs_',''),linecolor=colors[i], linestyle=2)
        #If this cut was the best ZBi for the iteration, set special marker
        if best_cut in graph.GetName():
            plot.SetMarkerStyle(34)
            plot.SetMarkerSize(3)
            plot.SetMarkerColor(colors[i])
        graphs.append(graph)

    overlay1DGraphs(graphs,'pct_sig_cut_%i_graphs'%(pct_sig),'test_plots',legx1=0.2,legx2=0.3,legy1=0.1,legy2=0.8,setXmin=0.0,setXmax=80.0,setYmin = 0.1, setYmax = 200.0, drawStats=False,logY=True)
