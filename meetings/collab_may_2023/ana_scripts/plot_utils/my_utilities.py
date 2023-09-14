#!/usr/bin/python3
import ROOT as r
import numpy as np
import copy
from difflib import SequenceMatcher

my_colors = [r.kGreen, r.kBlue,  r.kRed, r.kOrange, r.kYellow, r.kMagenta, r.kCyan, r.kPink+1, r.kSpring+10, r.kViolet+2, r.kTeal-1, r.kOrange+7, r.kMagenta-3, r.kYellow-3, r.kBlue+2 ,r.kPink-9]
def getColors():
    my_colors = [r.kGreen, r.kBlue,  r.kRed, r.kOrange, r.kYellow, r.kMagenta, r.kCyan, r.kPink+1, r.kSpring+10, r.kViolet+2, r.kTeal-1, r.kOrange+7, r.kMagenta-3, r.kYellow-3, r.kBlue+2 ,r.kPink-9]
    return my_colors

bottomStatsPos = 0
nplots = 0

def getRootSubDirs(infile, keyword = ""):
    subdirs = []
    for key in infile.GetListOfKeys():
        if keyword in key.GetName():
            subdirs.append(key.GetName())
    return subdirs

def read2DPlotsFromRootDir(infile, directory, keywords=[]):
    infile.cd()
    plots = {}
    rootdir = infile.Get(directory)
    rootdir.cd()
    for key in rootdir.GetListOfKeys():
        if "TH2" not in key.GetClassName():
            continue
        found = True
        if len(keywords) > 0:
            found = False
            for keyword in keywords:
                if keyword not in key.GetName():
                    continue
                else:
                    found = True
        if found != True:
            continue
        plot = copy.deepcopy(rootdir.Get(key.GetName()))
        if plot.GetEntries() == 0:
            continue
        name = plot.GetName().replace(" ","_").replace(".","")
        plot.SetName(name)
        plots[name] = plot

    return plots

def read1DPlotsFromRootDir(infile, directory, keywords=[],scale=1.0):
    infile.cd()
    print(infile.GetName())
    print("Looking for directory ", directory)
    plots = {}
    rootdir = infile.Get(directory)
    print(rootdir)
    rootdir.cd()
    for key in rootdir.GetListOfKeys():
        if "TH1" not in key.GetClassName():
            continue
        found = True
        if len(keywords) > 0:
            found = False
            for keyword in keywords:
                if keyword not in key.GetName():
                    continue
                else:
                    found = True
        if found != True:
            continue
        plot = copy.deepcopy(rootdir.Get(key.GetName()))
        if plot.GetEntries() == 0:
            continue
        name = plot.GetName().replace(" ","_").replace(".","")
        plot.SetName(name)
        if scale != 1.0:
            plot.Scale(scale)
        plots[name] = plot

    return plots

def format1DPlot(plot, name, title = "",linecolor = 1, linewidth = 2, linestyle = 1, markerstyle=0, markersize=1):
    if markerstyle > 0:
        plot.SetMarkerStyle(markerstyle)
        plot.SetMarkerSize(markersize)
    plot.SetLineWidth(linewidth)
    plot.SetLineColor(linecolor)
    plot.SetLineStyle(linestyle)
    plot.SetName(name)
    if title != "":
        plot.SetTitle(title)

def overlay1DPlots(plots,canvas_name, plots_dir,colors=my_colors,drawStyles=[],statsbox=False):
    #statsbox formatting
    global bottomStatsPos
    global nplots
    nplots = len(plots)

    statspositions,height = getStatsYPositions(len(plots))
    bottomStatsPos = statspositions[nplots-1] - height 

    ymax = -99999
    xmin = 999999
    xmax = -999999
    for i,plot in enumerate(plots):
        print("PLOT ", i)
        ymaxi = plot.GetMaximum()
        if ymaxi > ymax:
            ymax = ymaxi
        ixmin = plot.GetBinLowEdge(plot.FindFirstBinAbove(0))
        ixmax = plot.GetBinLowEdge(plot.FindLastBinAbove(0))
        if ixmin < xmin:
            xmin = ixmin
        if ixmax > xmax:
            xmax = ixmax

    c = r.TCanvas("%s"%(canvas_name),"%s"%(canvas_name),1800,1000)
    c.cd()

    #set yaxis
    for i,plot in enumerate(plots):
        plot.GetYaxis().SetRangeUser(0.00001,ymax*1.1)
        plot.GetXaxis().SetRangeUser(xmin*0.9,xmax*1.1)
        if len(drawStyles) > 0:
            if i < 1:
                plot.Draw("%s"%(drawStyles[i]))
                setStatsBox(plot,statspositions[i],height, linecolor=colors[i],statsbox=statsbox)
            else:
                plot.Draw("sames%s"%(drawStyles[i]))
                setStatsBox(plot,statspositions[i],height,linecolor=colors[i],statsbox=statsbox)
        else:
            if i < 1:
                plot.Draw("hist")
                setStatsBox(plot,statspositions[i],height, linecolor=colors[i],statsbox=statsbox)
            else:
                plot.Draw("histsames")
                setStatsBox(plot,statspositions[i],height,linecolor=colors[i],statsbox=statsbox)

    buildLegend(c)
    plots[0].SetTitle(canvas_name)
    #plots[0].SetTitle(canvas_name)
    c.SaveAs("%s/%s.png"%(plots_dir,canvas_name))

def getStatsYPositions(nplots):
    start = 1.0
    end = 0.1
    height = (start-end)/nplots
    if height > 0.1:
        height = 0.1
    #ypositions = list(range(start,end,height))
    ypositions = np.arange(end,start,height)
    ypositions=np.flip(ypositions)
    return ypositions,height

def setStatsBox(plot,ypos,height,linecolor=r.kBlack,color=r.kWhite,alpha=0.8,statsbox=False):
    r.gPad.Update()
    stats = plot.GetListOfFunctions().FindObject("stats")
    plot.SetStats(statsbox)
    if statsbox:
        stats.SetFillColorAlpha(color,alpha)
        stats.SetLineColor(linecolor)
        r.gStyle.SetStatY(ypos)
        r.gStyle.SetStatH(height)
        r.gStyle.SetStatW(0.1)
        r.gStyle.SetStatX(0.9)
        r.gPad.Update()

def buildLegend(canvas,x1=0.50, y1=0.7,x2=0.80,y2=0.9,textsize=0.025,separation=0.05,fill=0,border=0,statsbox=False):
    canvas.cd()
    if statsbox:
        height = 0.025*nplots
        if bottomStatsPos > 0:
            y2 = bottomStatsPos
            y1 = y2 - 1.2*height 
    legend = canvas.BuildLegend(x1,y1,x2,y2)
    legend.Draw()
    legend.SetTextSize(textsize)
    #legend.SetEntrySeparation(separation)
    legend.SetFillStyle(0)
    #legend.SetFillColorAlpha(r.kGray,0.0)
    legend.SetBorderSize(border)

def buildSimpleLegend(canvas,x1=0.7,y1=0.7,x2=0.9,y2=0.9,textsize=0.020,separation=0.02,fill=0,border=0):
    legend=canvas.BuildLegend(x1,y1,x2,y2)
    legend.Draw()
    legend.SetFillColorAlpha(r.kGray,0.1)
    legend.SetBorderSize(border)

def makeRatioPlot(plots,ratioPairs,ratioNames,canvas_name,plots_dir,colors=my_colors,drawStyles=[],LogX=False,LogY=False,RatioMin=0.01, RatioMax=2.0,statsbox=False):
    #statsbox formatting
    global bottomStatsPos
    global nplots
    nplots = len(plots)

    statspositions,height = getStatsYPositions(len(plots))
    bottomStatsPos = statspositions[nplots-1] - height 

    ymax = -99999
    for i,plot in enumerate(plots):
        ymaxi = plot.GetMaximum()
        if ymaxi > ymax:
            ymax = ymaxi

    c = r.TCanvas("%s"%(canvas_name),"%s"%(canvas_name),2200,1400)
    c.SetMargin(0,0,0,0)
    top = r.TPad("top","top",0,0.42,1,1)
    bot = r.TPad("bot","bot",0,0,1,0.38)
    #top = r.TPad("top","top",0,0.4,1,1)
    #bot = r.TPad("bot","bot",0,0,1,0.4)

    if LogX:
        top.SetLogx(1)
        bot.SetLogx(1)
    if LogY:
        top.SetLogy(1)
        bot.SetLogy(1)

    top.Draw()
    top.SetBottomMargin(0.0)
    #top.SetTopMargin(r.gStyle.GetPadTopMargin()*topScale)
    bot.Draw()
    bot.SetTopMargin(0)
    bot.SetBottomMargin(0.1)
    top.cd()
    plotsProperties=[]

    #set yaxis
    for i,plot in enumerate(plots):
        plot.GetYaxis().SetRangeUser(0.00001,ymax*1.1)
        if len(drawStyles) > 0:
            if i < 1:
                plot.Draw("%s"%(drawStyles[i]))
                setStatsBox(plot,statspositions[i],height, linecolor=colors[i],statsbox=statsbox)
            else:
                plot.Draw("sames%s"%(drawStyles[i]))
                setStatsBox(plot,statspositions[i],height,linecolor=colors[i],statsbox=statsbox)
        else:
            if i < 1:
                plot.Draw("hist")
                setStatsBox(plot,statspositions[i],height, linecolor=colors[i],statsbox=statsbox)
            else:
                plot.Draw("histsames")
                setStatsBox(plot,statspositions[i],height,linecolor=colors[i],statsbox=statsbox)

    buildLegend(top)

    #-------------Ratio---------------#
    bot.cd()
    ratioPlots = []
    for i,pair in enumerate(ratioPairs):
        numerator_idx, denom_idx = pair
        print("Numerator: ", plots[numerator_idx].GetName())
        print("Denominator: ", plots[denom_idx].GetName())
        numerator = plots[numerator_idx].Clone("numerator_%s"%(plots[numerator_idx].GetName()))
        numerator.SetStats(0)
        numerator.GetYaxis().SetRangeUser(RatioMin,RatioMax)
        numerator.SetNdivisions(508)
        numerator.GetYaxis().SetDecimals(True)
        numerator.Draw("axis")
        numerator.SetTitle(ratioNames[i])
        
        denominator = plots[denom_idx].Clone("denom_%s"%(plots[denom_idx].GetName()))
        denominator.SetStats(0)
        numerator.Divide(denominator)
        ratioPlots.append(numerator)

    for i,plot in enumerate(ratioPlots):
        if i < 1:
            plot.Draw("ep")
        else:
            plot.Draw("ep same")
            print("drawing sames")
        #numerator.SetTitle('%s_divided_by_%s'%(numstring,denstring))
    buildSimpleLegend(bot)

    plots[0].SetTitle(canvas_name)
    c.SaveAs("%s/%s.png"%(plots_dir,canvas_name))
