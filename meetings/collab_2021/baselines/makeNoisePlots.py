import ROOT as r
import argparse
import numpy as np

def getColors():
    colors = [r.kBlue, r.kGreen, r.kOrange, r.kRed, r.kYellow, r.kMagenta, r.kCyan, r.kPink+1, r.kSpring+10, r.kViolet+2, r.kTeal-1, r.kOrange+7, r.kMagenta-3, r.kYellow-3, r.kBlue+2 ,r.kPink-9]
    return colors

def getRootSubDirs(infile, keyword = ""):
    subdirs = []
    for key in infile.GetListOfKeys():
        if "TDirectoryFile" in key.GetClassName(): 
            if keyword in key.GetName():
                subdirs.append(key.GetName())
    return subdirs

def readHistogram(infile, histoname, subDir=None):
    #print("Looking for %s/%s"%(subDir,histoname))
    infile.cd()
    histo = None
    keys = None
    if subDir is not None:
        subDir = infile.Get("%s"%(subDir))
    else:
        subDir = infile

    for key in subDir.GetListOfKeys():
        if histoname in key.GetName():
            histo = subDir.Get("%s"%(key.GetName()))

    #if histo is not None:
        #print("Found %s"%(histo.GetName()))
    #else:
    #    print("FAILED TO FIND HISTOGRAM %s/%s IN FILE"%(subDir,key.GetName()))

    return histo

def formatHistogram(h,linecolor=1,linewidth=1,name=None,title=None,xlabel=None, ylabel=None,fillcolor=None):
    h.SetLineColor(linecolor)
    h.SetLineWidth(linewidth)
    if name is not None:
        h.SetName("%s"%(name))
    if title is not None:
        h.SetTitle("%s"%(title))
    if xlabel is not None:
        h.GetXaxis().SetTitle("%s"%(xlabel))
    if ylabel is not None:
        h.GetYaxis().SetTitle("%s"%(ylabel))
    if fillcolor is not None:
        h.SetFillColor(fillcolor)

def buildLegend(canvas,x1=0.7, y1=0.24,x2=0.82,y2=0.88,textsize=0.025,separation=1.0,fill=0,border=0):
    #canvas.cd()
    legend = canvas.BuildLegend(x1,y1,x2,y2)
    legend.Draw()
    legend.SetTextSize(textsize)
    legend.SetEntrySeparation(separation)
    legend.SetFillStyle(fill)
    legend.SetBorderSize(border)


##################################################################################################
parser = argparse.ArgumentParser(description="baseConfig options ")
parser.add_argument('--infile','-i', type=str, dest="inFilename", metavar='infiles',
                    help="Input files, specify on or more.")
parser.add_argument("-o", "--outFile", type=str, dest="outFilename", action='store',
                  help="Output filename.", metavar="outFilename", default="out.root")
parser.add_argument('--outDir', '-d', type=str, dest="outDir", action='store',
                  help="Specify the output directory.", metavar="outDir", default=".")
parser.add_argument('--run', '-r', type=str, dest="run", action='store',
                  help="Specify run number.", metavar="run", default="")
options = parser.parse_args()


r14191 = "/home/alic/HPS/projects/collab_2021/baselines/hps_14191_offline_baseline_fits_ana.root"
r14552 = "/home/alic/HPS/projects/collab_2021/baselines/hps_14552_offline_baseline_fits_ana.root"
r14596 = "/home/alic/HPS/projects/collab_2021/baselines/hps_14596_offline_baseline_fits_ana.root"
r14691 = "/home/alic/HPS/projects/collab_2021/baselines/hps_14691_offline_baseline_fits_ana.root"
r14710 = "/home/alic/HPS/projects/collab_2021/baselines/hps_14710_offline_baseline_fits_ana.root"
runs = {"14191":r14191,"14552":r14552,"14596":r14691,"14710":r14710}

outfile = r.TFile("%s"%(options.outFilename), "RECREATE")
outdir = options.outDir
run = options.run
colors = getColors()

keywords = ["L0", "L1"]

allplots = []
runplots = {}
for run in runs:
    plots = []
    infile = r.TFile("%s"%(runs[run]),"READ")
    for key in keywords:
        layers = getRootSubDirs(infile, keyword = key)
        for layer in layers:
            findname = "baseline"
            plot = readHistogram(infile, findname, subDir=layer)
            basename = plot.GetName()
            plot.SetName("%s_run_%s"%(basename,run))
            plots.append(plot)
            allplots.append(plot)
    runplots[run] = plots

sensPlots = {}
skip = []
for plot in allplots:
    plots = []
    split = plot.GetName().split("_")
    sensor = split[0] + split[1]
    run = split[-1]
    if sensor in skip:
        continue
    plots.append(plot)
    for plot2 in allplots:
        if plot == plot2:
            continue
        split2 = plot2.GetName().split("_")
        sensor2 = split2[0] + split2[1]
        run2 = split2[-1]
        if sensor == sensor2:
            skip.append(sensor2)
            plots.append(plot2)
    sensPlots[sensor] = plots 

outfile.cd()
for sensor in sensPlots:
    #print(sensor)
    plots = sensPlots[sensor]
    graphs = []
    for plot in plots:
        #print(plot.GetName())
        nbins = plot.GetXaxis().GetNbins()
        noises = []
        channels = []
        for b in range(nbins):
            noise = plot.GetErrorY(b)
            noises.append(noise)
            channels.append(b)
        gr = r.TGraph(nbins,np.array(channels, dtype=float),np.array(noises, dtype=float))
        gr.SetName("%s"%(plot.GetName()))
        gr.SetTitle("run_%s"%(plot.GetName().split("_")[-1]))
        graphs.append(gr)

    #Set title
    r.gStyle.SetOptTitle(0);
    t = r.TPaveText(0.0, 0.9, 0.3, 1.0, "brNDC")
    t.Draw()

    c = r.TCanvas("%s_Noise_wBeam"%(sensor),"%s_Noise_wBeam"%(sensor), 1800, 1000)
    c.cd()
    for i,graph in enumerate(graphs):
        print(graph.GetName())
        if i < 1:
            graph.Draw("ap")
        else:
            graph.Draw("samep")
        graph.SetMarkerStyle(8)
        graph.SetMarkerSize(1)
        graph.SetMarkerColor(colors[i])
    buildLegend(c, textsize=0.025, y1=0.75,y2=0.85)
    graphs[0].SetTitle("%s_Noise_wBeam"%(sensor))
    c.SetTitle("%s_Noise_wBeam;Channel,Noise (ADC)"%(sensor))
    c.Draw()
    c.Write()
        #graphs.append(gr)
outfile.Write()

            



#Loop over baseline plots for each run, and plot noise for each sensor on same plot
#for run in runplots:
#    plots = runplots[run]
#    print(run)
#    for plot in plots:
#       sensor = plot.GetName().split("_")[0] + plot.GetName().split("_")[1]
#       print(plot.GetName())
#       nbins = plot.GetXaxis().GetNbins()
#       noises = []
#       channels = []
#       for b in range(nbins):
#           noise = plot.GetErrorY(b)
#           noises.append(noise)
#           channels.append(b)
#       gr = r.TGraph(nbins,np.array(channels, dtype=float),np.array(noises, dtype=float))






