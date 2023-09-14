#!/usr/bin/python3
import ROOT as r
from optparse import OptionParser

def buildLegend(canvas,x1=0.7, y1=0.24,x2=0.94,y2=0.88,textsize=0.04,separation=1.0,fill=0,border=0):
    #canvas.cd()
    legend = canvas.BuildLegend(x1,y1,x2,y2)
    legend.Draw()
    legend.SetTextSize(textsize)
    legend.SetEntrySeparation(separation)
    legend.SetFillStyle(fill)
    legend.SetBorderSize(border)

parser = OptionParser()
parser.add_option("-i", "--inputFile", type="string", dest="inputFile",
                help="Name of file to run on.", metavar="inputFile", default="expSigRate.root")
parser.add_option("-o", "--outputFile", type="string", dest="outputFile",
                help="Name of file to run on.", metavar="outputFile", default="excContours.root")
(options, args) = parser.parse_args()
infile = r.TFile(options.inputFile,"READ")
outFile = r.TFile(options.outputFile,"RECREATE")
#outFile = r.TFile("contourPlots.root","RECREATE")

colorsMap = { 15:r.kBlue, 14:r.kGreen, 12:r.kOrange, 13:r.kRed, 11:r.kYellow, 8:r.kMagenta, 9:r.kCyan, 10:r.kPink+1, 4:r.kSpring+10, 6:r.kViolet+2, 7:r.kTeal-1, 3:r.kOrange+7, 5:r.kMagenta-3, 2:r.kYellow-3, 1:r.kBlue+2, 0:r.kPink-9}

histos = []
infile.cd()
for key in infile.GetListOfKeys():
    if "excContour_low_Lumi" in key.GetName():
        histo = infile.Get("%s"%(key.GetName()))
        print(histo.GetName())
        histos.append(histo)



outFile.cd()
histos.reverse()
c = r.TCanvas("contour_plots","exc_contours_4pt55_4Pi; m_{A'} (MeV); log_{10}(\epsilon^2)",1800,900)
c.cd()
legend = r.TLegend(0.8,0.24,0.94,0.40);
legend.SetBorderSize(0)
for i,histo in enumerate(histos):
    col = colorsMap[i]
    histo.SetLineWidth(2)
    histo.SetLineColor(col)
    split = histo.GetName().split("_") 
    lumi = split[3]
    name="Lumi: " + str(lumi) + " (1/pb)"
    print(histo.GetName())
    if i < 1:
        histo.Draw()
    else:
        histo.Draw("same")
    legend.AddEntry(histo,"%s"%(name),"l")

#buildLegend(c)
c.SetLogy()
legend.Draw()
c.Write()
infile.Close()
outFile.Close()

