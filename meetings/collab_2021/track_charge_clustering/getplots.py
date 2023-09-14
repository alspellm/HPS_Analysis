#!usr/bin/python3
import ROOT as r

inFiles = ["hps_014191_40_trkClus.root", "hps_014191_40_trkClus10hits.root", "hps_014191_40_trkClus12hits.root", "hps_014552_10_trkClus.root", "hps_014552_10_trkClus10hits.root", "hps_014552_10_trkClus12hits.root", "hps_014596_40_trkClus.root", "hps_014596_40_trkClus10hits.root", "hps_014596_40_trkClus12hits.root", "hps_014691_150_trkClus.root", "hps_014691_150_trkClus10hits.root", "hps_014691_150_trkClus12hits.root"]  
rundir = "/home/alic/HPS/projects/hit_reco/track_charge_clustering"

#layers = ["L0", "L1", "L3","L6"]
layers = ["L0", "L1", "L5"]


for f in inFiles:
    s = f.split("_")
    run = s[0]+"_"+ s[1]
    hits = s[3].replace(".root","")
    plots = []
    projxplots = []
    projyplots = []

    infile = r.TFile("%s/%s"%(rundir,f),"READ")
    anaDir = infile.Get("anaClusOnTrk")
    for histo in anaDir.GetListOfKeys():
        if histo.GetClassName().find("TH2") == -1:
            continue
        histoname = histo.GetName()
        for layer in layers:
            if histoname.find(layer) != -1 and histoname.find("charge_vs_stripPos") != -1:
                plot = anaDir.Get("%s"%(histoname))
                histoname = histoname.replace("anaClusOnTrk_","") + "_" + run + "_" + hits 
                plot.SetName("%s"%(histoname))
                plot.SetTitle("%s"%(histoname))

                ymax = plot.GetMean(2) + 2* plot.GetStdDev(2)
                plot.GetYaxis().SetRangeUser(0.0, ymax)
                fb =0
                lb = 640
                if(layer == "L0" or layer == "L1"):
                    for half in ("left","right"):
                        plotclone = plot.Clone()
                        if half == "left":
                            fb = 0
                            lb = 100
                        elif half == "right":
                            fb = 400
                            lb = 512

                        plotclone.GetXaxis().SetRange(fb,lb)
                        plotclone.SetName("%s_%s"%(half, histoname))
                        plotclone.SetTitle("%s_%s"%(half,histoname))
                        print(half)
                        print(plotclone.GetName())
                        plots.append(plotclone)

                        projx = plotclone.ProjectionX('%s_projX'%(plotclone.GetName()),0,plotclone.GetNbinsY(), "e")
                        projx.SetLineWidth(2)
                        projx.SetLineColor(3)
                        projxplots.append(projx)

                        projy = plotclone.ProjectionY('%s_projY'%(plotclone.GetName()),fb,lb, "e")
                        projy.SetLineWidth(2)
                        projy.SetLineColor(4)
                        projyplots.append(projy)

                else:
                    if plot.GetMean(1) > 300:
                        fb = 300
                        lb = 640
                    else:
                        fb = 0
                        lb = 300

                    plot.GetXaxis().SetRange(fb, lb)
                    plots.append(plot)

                    projx = plot.ProjectionX('%s_projX'%(plot.GetName()),0,plot.GetNbinsY(), "e")
                    projx.SetLineWidth(2)
                    projx.SetLineColor(3)
                    projxplots.append(projx)

                    projy = plot.ProjectionY('%s_projY'%(plot.GetName()),fb,lb, "e")
                    projy.SetLineWidth(2)
                    projy.SetLineColor(4)
                    projyplots.append(projy)
            else:
                continue

    for p in plots:
        dir1 = hits
        split = p.GetName().split("_") 
        layer = ""
        if "left" in p.GetName() or "right" in p.GetName():
            layer = split[1][0:2]
        else:
            layer = split[0][0:2]
        c = r.TCanvas("%s"%(p.GetName()),"%s"%(p.GetName()) , 1200, 800)
        c.cd()
        p.Draw("colz")
        c.Draw()

        c.SaveAs("./images/%s/%s/%s.png"%(dir1,layer, p.GetName() ))
        c.Close()

    for x in projxplots:
        dir1 = hits
        split = x.GetName().split("_") 
        layer = ""
        if "left" in x.GetName() or "right" in x.GetName():
            layer = split[1][0:2]
        else:
            layer = split[0][0:2]

        c = r.TCanvas("projx_%s"%(x.GetName()),"projx_%s"%(x.GetName()) , 1200, 800)
        c.cd()
        x.Draw("LP")
        c.Draw()

        c.SaveAs("./images/xprojections/%s/%s/%s.png"%(layer,run, x.GetName() ))
        c.Close()

    for x in projyplots:
        dir1 = hits
        split = x.GetName().split("_") 
        layer = ""
        if "left" in x.GetName() or "right" in x.GetName():
            layer = split[1][0:2]
        else:
            layer = split[0][0:2]

        c = r.TCanvas("projy_%s"%(x.GetName()),"projy_%s"%(x.GetName()) , 1200, 800)
        c.cd()
        x.Draw("LP")
        c.Draw()

        c.SaveAs("./images/yprojections/%s/%s/%s.png"%(layer,run, x.GetName() ))
        c.Close()

    infile.Close()





        


