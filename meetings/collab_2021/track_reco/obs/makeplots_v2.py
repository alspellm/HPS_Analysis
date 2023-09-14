import ROOT as r

colors = [r.kRed, r.kOrange-3, r.kYellow+2, r.kSpring+9, r.kGreen+3,r.kTeal+1,r.kBlue+2, r.kAzure+10,r.kBlue+1,r.kViolet-6,r.kMagenta,r.kPink+2, r.kGreen-1, r.kOrange-1, r.kViolet+10, r.kCyan+3,]
#colors  = [r.kBlue+2, r.kCyan+2, r.kRed+2,r.kOrange+10,r.kYellow+2,r.kGreen-1,r.kAzure-2,r.kGreen-8,r.kOrange+3,r.kYellow+2,r.kRed+2,r.kMagenta, r.kCyan+3, r.kViolet-6, r.kTeal+1, r.kPink+2]

infilename = "./output/hps_14552_track_ana.root"
outfilename = "./output/hps_14552_HitCodeStacks.root"
infile = r.TFile("%s"%(infilename),"READ")
outfile = r.TFile("%s"%(outfilename), "RECREATE")

hcodes = []
for code in infile.GetListOfKeys():
    if "hc" in code.GetName():
        hcodes.append(code.GetName())
print(hcodes)

charges = ["Pos", "Ele"]
halves = ["top", "bot"]
params = ["TanLambda", "Z0", "d0", "Phi", "p_h"]

histonames = []
for param in params:
    for charge in charges:
        for half in halves:
            name = half+charge+"_"+param
            histonames.append(name)

for hist in histonames:
    stack = []
    if "TanLambda" in hist:
        if "bot" in hist:
            continue
        else:
            hist = hist.replace("top",'')
    print("Find ",hist)
    for code in hcodes:
        intcode = int(code.split("_")[0].replace("hc",''))
        print(intcode)
        infile.cd()
        cdir = infile.Get("%s"%(code))
        for key in cdir.GetListOfKeys():
            if hist in key.GetName():
                print("Found histogram ", key.GetName())
                histo = cdir.Get("%s"%(key.GetName()))
                #histo.SetLineColor(colors[intcode])
                histo.SetLineColor(r.kBlack)
                histo.SetLineWidth(1)
                histo.SetFillColor(colors[intcode])
                histo.SetName("%s"%(code))
                histo.SetTitle("%s"%(code))
                stack.append(histo)

    #Sort histos by stats
    sortedStack = []
    stats = []
    for histo in stack:
        stats.append(histo.GetEntries())
    sortedStats = sorted(stats)
    for stat in sortedStats:
        idx = stats.index(stat)
        sortedStack.append(stack[idx])

    outfile.cd()
    print("Building stack for ",hist)
    hstack = r.THStack("hs_%s"%(hist), "%s_Stacked_HitCodes"%(hist))
    for histo in sortedStack:
        hstack.Add(histo)
    c = r.TCanvas("c_%s"%(hist), "%s"%(hist), 1800, 900)
    c.cd()
    r.gStyle.SetPalette(r.kOcean)
    hstack.Draw("HIST")
    hstack.GetXaxis().SetTitle("%s"%(hist.split("_")[1]))
    #hstack.Write()
    legend = c.BuildLegend()
    legend.Draw()
    legend.SetTextSize(0.018)
    legend.SetEntrySeparation(1.0)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    c.Write()
outfile.Write()
