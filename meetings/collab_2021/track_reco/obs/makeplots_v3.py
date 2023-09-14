import ROOT as r

#colors = [r.kRed, r.kOrange-3, r.kYellow+2, r.kSpring+9, r.kGreen+3,r.kTeal+1,r.kBlue+2, r.kAzure+10,r.kBlue+1,r.kViolet-6,r.kMagenta,r.kPink+2, r.kGreen-1, r.kOrange-1, r.kViolet+10, r.kCyan+3,]
#colors  = [r.kBlue+2, r.kCyan+2, r.kRed+2,r.kOrange+10,r.kYellow+2,r.kGreen-1,r.kAzure-2,r.kGreen-8,r.kOrange+3,r.kYellow+2,r.kRed+2,r.kMagenta, r.kCyan+3, r.kViolet-6, r.kTeal+1, r.kPink+2]
colors = [r.kBlue, r.kGreen, r.kOrange, r.kRed, r.kYellow, r.kMagenta, r.kCyan, r.kPink+1, r.kSpring+10, r.kViolet+2, r.kTeal-1, r.kOrange+7, r.kMagenta-3, r.kYellow-3, r.kBlue+2 ,r.kPink-9]
colorsMap = { 15:r.kBlue, 14:r.kGreen, 12:r.kOrange, 13:r.kRed, 11:r.kYellow, 8:r.kMagenta, 9:r.kCyan, 10:r.kPink+1, 4:r.kSpring+10, 6:r.kViolet+2, 7:r.kTeal-1, 3:r.kOrange+7, 5:r.kMagenta-3, 2:r.kYellow-3, 1:r.kBlue+2, 0:r.kPink-9}
colors.reverse()
run = "run_14191"
infilename = "./output/hps_14191_ana_track.root"
outfilename = "./output/hps_14191_HitCodeStacks.root"
infile = r.TFile("%s"%(infilename),"READ")
outfile = r.TFile("%s"%(outfilename), "RECREATE")
outdir = "./images"

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

histostacks = []

for hist in histonames:
    stack = []
    print("Find ",hist)
    for code in hcodes:
        intcode = int(code.split("_")[0].replace("hc",''))
        infile.cd()
        cdir = infile.Get("%s"%(code))
        for key in cdir.GetListOfKeys():
            if hist in key.GetName():
                print("Found histogram ", key.GetName())
                histo = cdir.Get("%s"%(key.GetName()))
                #histo.SetLineColor(colors[intcode])
                histo.SetLineColor(r.kBlack)
                histo.SetLineWidth(1)
                #histo.SetFillColor(colors[intcode])
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
    hstack = r.THStack("%s_%s_hitcodes"%(hist,run), "%s_%s_hitcodes"%(hist, run))
    for i,histo in enumerate(sortedStack):
        intcode = int(histo.GetName().split("_")[0].replace("hc",''))
        histo.SetFillColor(colorsMap[intcode])
        hstack.Add(histo)
    histostacks.append(hstack)

tl_Ele = []
tl_Pos = []
other = []
for stack in histostacks:
    if "Ele_TanLambda" in stack.GetName():
        tl_Ele.append(stack)
    elif "Pos_TanLambda" in stack.GetName():
        tl_Pos.append(stack)
    else:
        other.append(stack)


for hstack in other:
    hist = hstack.GetName()
    print("HIST!!!! ", hist)
    c = r.TCanvas("c_%s"%(hist), "%s"%(hist), 1800, 900)
    c.cd()
    hstack.Draw("HIST")
    hstack.GetXaxis().SetTitle("%s"%(hist.split("_")[1]))
    xmin = hstack.GetXaxis().GetXmin()
    xmax = hstack.GetXaxis().GetXmax()
    ymax = 1.2*hstack.GetYaxis().GetXmax()
    if "Phi" in hist:
        xmin = -0.18
        xmax = 0.18
    elif "Z0" in hist:
        xmin = -2.0
        xmax = 2.0
    hstack.GetXaxis().SetRangeUser(xmin,xmax)
    hstack.GetYaxis().SetRangeUser(0.0,ymax)
    legend = c.BuildLegend(0.7,0.24,0.82,0.88)
    legend.Draw()
    legend.SetTextSize(0.018)
    legend.SetEntrySeparation(1.0)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    c.Write()
    c.SaveAs("%s/%s.png"%(outdir,hstack.GetName()))
    c.Close()

c = r.TCanvas("%s_%s"%("Ele_TanLambda", run), "%s_%s"%("Ele_TanLambda",run), 1800, 900)
c.cd()
ymax = 0
for i,hstack in enumerate(tl_Ele):
    if i < 1:
        hstack.Draw("HIST")
        legend = c.BuildLegend(0.7,0.24,0.82,0.88)
        legend.Draw()
        legend.SetTextSize(0.018)
        legend.SetEntrySeparation(1.0)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        y = 1.2*hstack.GetYaxis().GetXmax()
        if y > ymax:
            ymax = y
    else:
        hstack.Draw("HISTsame")
        y = 1.2*hstack.GetYaxis().GetXmax()
        if y > ymax:
            ymax = y
    hstack.GetXaxis().SetTitle("%s"%(hstack.GetName().split("_")[1]))
    #xmin = hstack.GetXaxis().GetXmin()
    #xmax = hstack.GetXaxis().GetXmax()
    ymax = 1.2*hstack.GetYaxis().GetXmax()
    hstack.GetXaxis().SetRangeUser(-0.1,0.1)
    hstack.GetYaxis().SetRangeUser(0.0,ymax)
c.Write()
c.SaveAs("%s/%s.png"%(outdir,c.GetName()))
c.Close()

c = r.TCanvas("%s_%s"%("Pos_TanLambda", run), "%s_%s"%("Pos_TanLambda", run), 1800, 900)
ymax = 0
for i,hstack in enumerate(tl_Pos):
    if i < 1:
        hstack.Draw("HIST")
        legend = c.BuildLegend(0.7,0.24,0.82,0.88)
        legend.Draw()
        legend.SetTextSize(0.018)
        legend.SetEntrySeparation(1.0)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        y = 1.2*hstack.GetYaxis().GetXmax()
        if y > ymax:
            ymax = y
    else:
        hstack.Draw("HISTsame")
        y = 1.2*hstack.GetYaxis().GetXmax()
        if y > ymax:
            ymax = y
    hstack.GetXaxis().SetTitle("%s"%(hstack.GetName().split("_")[1]))
    xmin = hstack.GetXaxis().GetXmin()
    xmax = hstack.GetXaxis().GetXmax()
    ymax = 1.2*hstack.GetYaxis().GetXmax()
    #hstack.GetXaxis().SetRangeUser(xmin,xmax)
    hstack.GetXaxis().SetRangeUser(-0.1,0.1)
    hstack.GetYaxis().SetRangeUser(0.0,ymax)
c.Write()
c.SaveAs("%s/%s.png"%(outdir,c.GetName()))
c.Close()
outfile.Write()
