#!/usr/bin/python3
import ROOT as r
import copy
from array import array
from copy import deepcopy
import glob
import os
import sys
from tabulate import tabulate

def addHisto1D(name, nbins, xmin, xmax, title='', xlabel='', ylabel=''):
    h = r.TH1F('%s'%(name),'%s;%s;%s'%(title, xlabel, ylabel),nbins, xmin, xmax)
    return h

def addHisto2D(name, nbinsx, xmin, xmax, nbinsy, ymin, ymax, title='', xlabel='', ylabel=''):
    h = r.TH2F('%s'%(name),'%s;%s;%s'%(title, xlabel, ylabel),nbinsx, xmin, xmax, nbinsy, ymin, ymax)
    return h

def drawTLine(canvas, x_value, line_color=2):
    canvas.cd()
    line = r.TLine(x_value, canvas.GetUymin(), x_value, canvas.GetUymax())
    line.SetLineColor(line_color)  # Set the line color (optional)
    line.Draw()

def applyGeneralHPSConfig():
    #General configuration
    bottomFraction = 0.4
    bottomScale = 1./bottomFraction
    topScale = 1./(1. - bottomFraction)
    r.TProfile.Approximate(True)

def SetMyStyle(tsize=0.035, tzsize=0.035, font=42, setOptTitle=0, setOptStat=0, setOptFit=0):
    print("SETTING MY STYLE")

    colors = getColorsHPS()
    r.gROOT.SetBatch(1)

    myStyle = r.TStyle("myStyle", "my style")

    # Set your custom attributes here
    myStyle.SetOptTitle(setOptTitle)
    myStyle.SetOptStat(setOptStat)
    myStyle.SetOptFit(setOptFit)
    myStyle.SetTitleFont(font)
    myStyle.SetTitleSize(tsize)
    #myStyle.SetTitleX(0.5)
    #myStyle.SetTitleY(0.98)

    #Set legend text size
    myStyle.SetLegendTextSize(0.02)

    # Set the title text color to black
    myStyle.SetTitleTextColor(r.kBlack)

    # use plain black on white colors
    icol = 0
    myStyle.SetFrameBorderMode(icol)
    myStyle.SetCanvasBorderMode(icol)
    myStyle.SetPadBorderMode(icol)
    myStyle.SetPadColor(icol)
    myStyle.SetCanvasColor(icol)
    myStyle.SetStatColor(icol)

    # set the paper & margin sizes
    myStyle.SetPaperSize(20, 26)
    myStyle.SetPadTopMargin(0.10)
    myStyle.SetPadRightMargin(0.05)
    myStyle.SetPadBottomMargin(0.10)
    myStyle.SetPadLeftMargin(0.10)

    myStyle.SetTextSize(tsize) 
    myStyle.SetLabelFont(font, "x")
    myStyle.SetTitleFont(font, "x")
    myStyle.SetLabelFont(font, "y")
    myStyle.SetTitleFont(font, "y")
    myStyle.SetLabelFont(font, "z")
    myStyle.SetTitleFont(font, "z")

    myStyle.SetLabelSize(tsize, "x")
    myStyle.SetTitleSize(tsize, "x")
    myStyle.SetLabelSize(tsize, "y")
    myStyle.SetTitleSize(tsize, "y")
    myStyle.SetLabelSize(tzsize, "z")
    myStyle.SetTitleSize(tzsize, "z")

    myStyle.SetTitleOffset(1.0, "y")
    myStyle.SetTitleOffset(1.15, "x")

    #use bold lines and markers
    myStyle.SetMarkerSize(1.0)
    myStyle.SetMarkerStyle(8)
    myStyle.SetMarkerColor(1)
    myStyle.SetLineColor(1)
    myStyle.SetHistLineWidth(3)
    #myStyle.SetLineStyleString(2, "[12 12]")  # postscript dashes

    # put tick marks on top and RHS of plots
    #myStyle.SetPadTickX(1)
    #myStyle.SetPadTickY(1)

    r.gROOT.SetStyle("myStyle")
    r.gROOT.ForceStyle()

    NRGBs = 5
    NCont = 255

    stops = array("d", [0.00, 0.34, 0.61, 0.84, 1.00])
    red = array("d", [0.00, 0.00, 0.87, 1.00, 0.51])
    green = array("d", [0.00, 0.81, 1.00, 0.20, 0.00])
    blue = array("d", [0.51, 1.00, 0.12, 0.00, 0.00])
    r.TColor.CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont)

    return myStyle

def InsertText(insertText=[], text_x=0.2, text_y=0.8, line_spacing=0.03, text_size=0.025, Hps=True):

    drawText = insertText

    latex = r.TLatex()
    latex.SetTextFont(42)
    latex.SetTextSize(text_size)
    latex.SetTextAlign(12)
    latex.SetTextColor(r.kBlack)

    if (Hps):
        latex.DrawLatexNDC(text_x, text_y,'#bf{#it{HPS}} Internal')
        text_y = text_y - line_spacing

    for line in drawText:
        latex.DrawLatexNDC(text_x, text_y,line)
        text_y = text_y - line_spacing
    return latex

def read_2d_plots_from_root_file_dirs(file_path, root_dir_key="",keyword=""):
    # List to store the matching 1D plots
    dir_plots = {}

    # Open the ROOT file
    root_file = r.TFile.Open(file_path)

    #Loop over all directories matching key
    for dir_key in root_file.GetListOfKeys():
        dir_obj = dir_key.ReadObj()

        # Check if the object is TDir
        if isinstance(dir_obj, r.TDirectory):
            if root_dir_key not in dir_obj.GetName():
                continue

            print("Navigating to directory ", dir_obj.GetName())
            # Get the directory within the ROOT file
            root_dir = root_file.GetDirectory(dir_obj.GetName())

            plots = []

            # Loop over all objects in the directory
            for key in root_dir.GetListOfKeys():
                obj = key.ReadObj()

                # Check if the object is a 1D histogram
                if isinstance(obj, r.TH2) and obj.GetDimension() == 2:
                    # Check if the object name contains the keyword (if provided)
                    if keyword and keyword not in obj.GetName():
                        continue
                    print("Copying plot", obj.GetName())
                    # Create a deepcopy of the plot and append to the list
                    plots.append(copy.deepcopy(obj))

            dir_plots[dir_obj.GetName()] = plots

    # Close the ROOT file
    root_file.Close()

    return dir_plots

def plot_2d_plots_side_by_side(hist1, hist2, canvas_name, save_directory, insertText=[], text_x=0.4, text_y=0.3, 
        line_spacing=0.03, text_size=0.02):
    # Set transparent fill style for the statistics box
    r.gStyle.SetOptStat(1)  # Show statistics box
    r.gStyle.SetStatStyle(0)  # Set fill style to transparent
    canvas = r.TCanvas("%s"%(canvas_name), "%s"%(canvas_name), 1800, 900 )
    pad1 = r.TPad("pad1","Pad 1", 0.01, 0.01, 0.49, 0.99)
    pad2 = r.TPad("pad2","Pad 2", 0.51, 0.01, 0.99, 0.99)
    # Set margins and draw pads
    #pad1.SetBottomMargin(0)  # No margin at the bottom for pad1
    #pad2.SetTopMargin(0)     # No margin at the top for pad2
    pad1.SetRightMargin(0.15)
    pad1.SetLeftMargin(0.05)
    pad2.SetRightMargin(0.15)
    pad2.SetLeftMargin(0.05)

    pad1.Draw()
    pad2.Draw()

    pad1.cd()
    hist1.Draw("colz")

    pad2.cd()
    hist2.Draw("colz")

    canvas.Update()

    InsertText(insertText)
    canvas.Update()

    file_name = save_directory + "/" + canvas_name + ".png"
    canvas.SaveAs(file_name)

    # Clean up
    #canvas.Close()

def read_1d_plots_from_root_file_dirs(file_path, root_dir_key="", keyword=""):
    # List to store the matching 1D plots
    dir_plots = {}

    # Open the ROOT file
    root_file = r.TFile.Open(file_path)

    #Loop over all directories matching key
    for dir_key in root_file.GetListOfKeys():
        dir_obj = dir_key.ReadObj()

        # Check if the object is TDir
        if isinstance(dir_obj, r.TDirectory):
            if root_dir_key not in dir_obj.GetName():
                continue

            print("Navigating to directory ", dir_obj.GetName())
            # Get the directory within the ROOT file
            root_dir = root_file.GetDirectory(dir_obj.GetName())

            plots = []

            # Loop over all objects in the directory
            for key in root_dir.GetListOfKeys():
                obj = key.ReadObj()

                # Check if the object is a 1D histogram
                if isinstance(obj, r.TH1) and obj.GetDimension() == 1:
                    # Check if the object name contains the keyword (if provided)
                    if keyword and keyword not in obj.GetName():
                        continue
                    print("Copying plot", obj.GetName())
                    # Create a deepcopy of the plot and append to the list
                    plots.append(copy.deepcopy(obj))

            dir_plots[dir_obj.GetName()] = plots

    # Close the ROOT file
    root_file.Close()

    return dir_plots

def read_plot_from_root_file(file_path, name, root_dir=""):
    # Open the ROOT file
    root_file = r.TFile(file_path)

    # Get the directory within the ROOT file
    root_dir = root_file.GetDirectory(root_dir)

    # Check if the directory exists
    if not root_dir:
        print(f"Failed to find directory: {root_dir}")
        root_file.Close()
        return None

    plot = copy.deepcopy(root_dir.Get("%s"%(name)))

    root_file.Close()

    return plot

def read_1d_plots_from_root_file(file_path, root_dir="", keyword=""):
    # Open the ROOT file
    root_file = r.TFile.Open(file_path)

    # Check if the file is open
    if not root_file.IsOpen():
        print(f"Failed to open ROOT file: {file_path}")
        return []

    # Get the directory within the ROOT file
    root_dir = root_file.GetDirectory(root_dir)

    # Check if the directory exists
    if not root_dir:
        print(f"Failed to find directory: {root_dir}")
        root_file.Close()
        return []

    # List to store the matching 1D plots
    plots = []

    # Loop over all objects in the directory
    for key in root_dir.GetListOfKeys():
        obj = key.ReadObj()

        # Check if the object is a 1D histogram
        if isinstance(obj, r.TH1) and obj.GetDimension() == 1:
            # Check if the object name contains the keyword (if provided)
            if keyword and keyword not in obj.GetName():
                continue

            # Create a deepcopy of the plot and append to the list
            plots.append(copy.deepcopy(obj))

    # Close the ROOT file
    root_file.Close()

    return plots

def formatHisto(histogram, name=None, title=None, x_label=None, y_label=None,
        line_width=None, line_color=None, marker_style=None,
        marker_size=None, marker_color=None, line_style=None):
    if name is not None:
        histogram.SetName(name)
    if title is not None:
        histogram.SetTitle(title)
    if x_label is not None:
        histogram.GetXaxis().SetTitle(x_label)
    if y_label is not None:
        histogram.GetYaxis().SetTitle(y_label)
    if line_width is not None:
        histogram.SetLineWidth(line_width)
    if line_color is not None:
        histogram.SetLineColor(line_color)
    if marker_style is not None:
        histogram.SetMarkerStyle(marker_style)
    if marker_size is not None:
        histogram.SetMarkerSize(marker_size)
    if marker_color is not None:
        histogram.SetMarkerColor(marker_color)
    if line_style is not None:
        histogram.SetLineStyle(line_style)

def getMarkersHPS():
    markers = [r.kFullCircle, r.kFullTriangleUp, r.kFullSquare, r.kOpenSquare, r.kOpenTriangleUp, r.kOpenCircle, r.kFullCircle, r.kOpenSquare, r.kFullSquare, r.kOpenTriangleUp, r.kOpenCircle, r.kFullCircle, r.kOpenSquare, r.kFullSquare, r.kOpenTriangleUp, r.kOpenCircle, r.kFullCircle, r.kOpenSquare, r.kFullSquare, r.kOpenTriangleUp, r.kOpenCircle, r.kFullCircle, r.kOpenSquare, r.kFullSquare, r.kOpenTriangleUp, r.kOpenCircle, r.kFullCircle, r.kOpenSquare, r.kFullSquare, r.kOpenTriangleUp]
    return markers

def getColorsHPS():
    colors = [r.kBlue+2, r.kCyan+2, r.kRed+2, r.kOrange+10, r.kYellow+2, r.kGreen-1, r.kAzure-2, r.kGreen-8, r.kOrange+3, r.kYellow+2, r.kRed+2, r.kBlue+2, r.kGreen-8, r.kOrange+3, r.kYellow+2, r.kRed+2, r.kBlue+2, r.kGreen-8, r.kOrange+3, r.kYellow+2, r.kRed+2, r.kBlue+2, r.kGreen-8, r.kOrange+3, r.kYellow+2, r.kRed+2, r.kBlue+2, r.kGreen-8, r.kOrange+3]
    return colors

def getColors():
    # Array to store the colors
    colors = []

    # List of colors to choose from
    available_colors = [
            r.kBlack,
            r.kBlue,
            r.kRed,
            r.kGreen,
            r.kMagenta,
            r.kOrange,
            r.kTeal,
            r.kSpring,
            r.kGray,
            ]

    for color in available_colors:
        colors.append(color)

    print("Colors: ", colors)
    return colors

def format_multiStats(n):
    # Define the position and size of each statistics box
    box_x = 0.7  # X-coordinate of the top-right corner of the first box
    box_y = 0.7  # Y-coordinate of the top-right corner of the first box
    box_width = 0.2  # Width of each box
    box_height = 0.15  # Height of each box
    box_margin = 0.05  # Margin between each box

    # Calculate the positions for each statistics box
    box_positions = []
    for i in range(n):
        x = box_x
        y = box_y - (box_height + box_margin) * i
        box_positions.append((x, y))

    return box_positions 

def plot_TH1s_with_legend(histograms, canvas_name, drawOptions='hist', save_directory = '.',setStats=False,freezeXaxis=True,legx1=0.85,legy1=0.7,legx2=0.9,legy2=0.9, clear_legend=True, LogX=False, LogY=False,LogY_min=0.5, insertText=[],text_x=0.6, text_y=0.6, text_size = 0.03, line_spacing=0.03, save=False):
    # Create a canvas
    canvas = r.TCanvas(canvas_name, canvas_name, 2560, 1440)

    # Find the maximum x and y values among all histograms
    min_x = min(h.GetBinLowEdge(h.FindFirstBinAbove(0.0)) for h in histograms)
    max_x, max_y = max(h.GetBinLowEdge(h.FindLastBinAbove(0.0)) + h.GetBinWidth(0) for h in histograms), max(h.GetMaximum() for h in histograms)

    min_y = 1e10
    for h in histograms:
        local_miny = 1e10
        min_y_l = h.GetBinContent(h.FindFirstBinAbove(0.0))
        min_y_u = h.GetBinContent(h.FindLastBinAbove(0.0))
        if min_y_l < min_y_u:
            local_miny = min_y_l
        else:
            local_miny = min_y_u
        if local_miny < min_y:
            min_y = local_miny

    # Create a legend corresponding to each histogram
    legend = buildLegend(legx1, legy1, legx2, legy2, clear_legend)

    for i,histogram in enumerate(histograms):
        if(freezeXaxis == False):
            # Adjust the axis ranges for all histograms
            histogram.SetAxisRange(0.9*min_x, 1.1*max_x, "X")
            # Set the same maximum and minimum for both axes
            histogram.GetXaxis().SetRangeUser(0.9*min_x,1.1*max_x)
        histogram.SetAxisRange(min_y, 1.1 * max_y, "Y")
        histogram.GetYaxis().SetRangeUser(min_y, 1.1 * max_y)
        #if LogY and max_y > 0 and min_y <= 0:
        #    histogram.SetAxisRange(0.5, 1.1 * max_y, "Y")
        #    histogram.GetYaxis().SetRangeUser(LogY_min, 1.1 * max_y)

        # Plot the histogram on the canvas
        if i < 1:
            histogram.Draw('%s'%(drawOptions))
        else:
            histogram.Draw("%sSAME"%(drawOptions))

        if setStats == False:
            histogram.SetStats(0)
        else:
            stat_box = histogram.GetListOfFunctions().FindObject("stats")
            x, y = format_multiStats(len(histograms))
            stat_box.SetX1NDC(x)
            stat_box.SetY1NDC(y)
            stat_box.SetX2NDC(x + box_width)
            stat_box.SetY2NDC(y - box_height)
            stat_box.SetTextSize(0.02)  # Adjust the text size of the statistics box

        #Add legend entries
        legend.AddEntry(histogram, histogram.GetTitle(), "l")

    # Draw the legend on the canvas
    legend.Draw()

    if LogX:
        canvas.SetLogx(1)
    if LogY:
        canvas.SetLogy(1)

    # Save the canvas as a PNG file
    histograms[0].SetTitle(canvas_name)

    #Insert Text
    InsertText(insertText, text_x=text_x, text_y=text_y,text_size=text_size)
    canvas.Update()

    if save:
        file_name = save_directory + "/" + canvas_name + ".png"
        canvas.SaveAs(file_name)

    #return canvas, legend
    return deepcopy(canvas)

def plot_TH1_ratios_with_legend(histograms, numerators, denominators, ratioNames, ratioColors, canvas_name, save_directory,ratioMin=0.01, ratioMax=2.0, setStats=False,legx1=0.7,legy1=0.7,legx2=0.9,legy2=0.9, clear_legend=True, LogX=False, LogY=False, save=False):
    # Create a canvas
    canvas = r.TCanvas(canvas_name, canvas_name, 2560, 1440)
    canvas.SetMargin(0,0,0,0)
    top = r.TPad("top","top",0,0.42,1,1)
    bot = r.TPad("bot","bot",0,0,1,0.38)

    if LogX:
        top.SetLogx(1)
        bot.SetLogx(1)
    if LogY:
        top.SetLogy(1)
        bot.SetLogy(1)

    top.Draw()
    top.SetBottomMargin(0.0)
    bot.Draw()
    bot.SetTopMargin(0)
    bot.SetBottomMargin(0.1)
    top.cd()

    # Find the maximum x and y values among all histograms
    max_x, max_y = max(h.GetBinLowEdge(h.FindLastBinAbove(0)) for h in histograms), max(h.GetMaximum() for h in histograms)
    # Find the minimumg x and y values among all histograms
    min_x, min_y = min(h.GetBinLowEdge(h.FindFirstBinAbove(0)) for h in histograms), min(h.GetMinimum() for h in histograms)

    #Build Legend
    legend = buildLegend(legx1, legy1, legx2, legy2, clear_legend)

    for histogram in histograms:

        # Adjust the axis ranges for all histograms
        histogram.SetAxisRange(0.9*min_x, 1.1*max_x, "X")
        histogram.SetAxisRange(min_y, 1.1 * max_y, "Y")
        # Set the same maximum and minimum for both axes
        histogram.GetXaxis().SetRangeUser(0.9*min_x,1.1*max_x)
        histogram.GetYaxis().SetRangeUser(min_y, 1.1 * max_y)

        # Plot the histogram on the canvas
        histogram.Draw("hist SAME")

        if setStats == False:
            histogram.SetStats(0)
        else:
            stat_box = histogram.GetListOfFunctions().FindObject("stats")
            x, y = format_multiStats(len(histograms))
            stat_box.SetX1NDC(x)
            stat_box.SetY1NDC(y)
            stat_box.SetX2NDC(x + box_width)
            stat_box.SetY2NDC(y - box_height)
            stat_box.SetTextSize(0.02)  # Adjust the text size of the statistics box

        legend.AddEntry(histogram, histogram.GetTitle(), "l")

    # Draw the legend on the canvas
    legend.Draw()

    histograms[0].SetTitle(canvas_name)

    #-------------Ratio---------------#

    ratio_plots = []
    bot.cd()
    for i, plot in enumerate(numerators):
        numerator = plot.Clone("numerator_%s"%plot.GetName())
        denominator = denominators[i].Clone("denominator_%s"%denominators[i].GetName())
        numerator.SetStats(0)
        numerator.GetYaxis().SetRangeUser(ratioMin,ratioMax)
        numerator.SetNdivisions(508)
        numerator.GetYaxis().SetDecimals(True)
        numerator.Draw("axis")
        numerator.SetTitle(ratioNames[i])

        denominator.SetStats(0)
        numerator.Divide(denominator)
        numerator.SetLineColor(ratioColors[i])
        ratio_plots.append(numerator)

    for i, plot in enumerate(ratio_plots):
        if i < 1:
            plot.Draw("ep")
        else:
            plot.Draw("ep same")

    bot_legend = bot.BuildLegend(legx1, legy1, legx2, legy2)
    bot_legend.Draw()
    if clear_legend:
        bot_legend.SetFillStyle(0)
        bot_legend.SetFillColor(0)
        bot_legend.SetLineColor(0)
        bot_legend.SetBorderSize(0)

    # Save the canvas as a PNG file
    if save:
        file_name = save_directory + "/" + canvas_name + ".png"
        canvas.SaveAs(file_name)

    #return canvas, legend, bot_legend
    return deepcopy(canvas)


def Make2DPlot(canvas_name, histo, title="", xtitle="", ytitle="", ztitle="", insertText=[], zmin="", zmax="", outdir='.', save=False):
    oFext = ".png"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if title:
        histo.SetTitle(title)

    can = r.TCanvas('%s'%(canvas_name),'%s'%(canvas_name),2500,1440)
    #can.SetTitle('%s'%(canvas_name))
    #can.SetRightMargin(0.2)

    #histolist[ih].GetZaxis().SetRangeUser(zmin,zmax)
    if xtitle:
        histo.GetXaxis().SetTitle(xtitle)
    #histo.GetXaxis().SetTitleSize(
    #        histo.GetXaxis().GetTitleSize()*0.7)
    #histo.GetXaxis().SetLabelSize(
    #        histo.GetXaxis().GetLabelSize()*0.75)
    #histo.GetXaxis().SetTitleOffset(
    #        histo.GetXaxis().GetTitleOffset()*0.8)

    #histo.GetYaxis().SetTitleSize(
    #        histo.GetYaxis().GetTitleSize()*0.7)
    #histo.GetYaxis().SetLabelSize(
    #        histo.GetYaxis().GetLabelSize()*0.75)
    #histo.GetYaxis().SetTitleOffset(
    #        histo.GetYaxis().GetTitleOffset()*1.7)
    if ytitle:
        histo.GetYaxis().SetTitle(ytitle)

    histo.Draw("colz")

    InsertText(insertText)

    if save:
        can.SaveAs(outdir+"/"+name+oFext)
    return deepcopy(can)

def Make1DPlot(canvas_name, histo, title="",xtitle="", ytitle="", ztitle="", drawOptions = "", insertText=[], outdir='.', LogY=False,save=True):
    oFext = ".png"
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if title:
        histo.SetTitle(title)

    can = r.TCanvas('%s'%(canvas_name),'%s'%(canvas_name),2500,1440)
    #can.SetName('%s'%(canvas_name))
    #can.SetTitle('%s'%(canvas_name))
    can.SetRightMargin(0.2)

    #xmax = histo.GetBinLowEdge(histo.FindLastBinAbove(0))
    #ymax = histo.GetMaximum()
    #xmin = histo.GetBinLowEdge(histo.FindFirstBinAbove(0))
    #ymin = histo.GetMinimum()
    #if setymin:
    #    ymin = set_ymin

    if xtitle:
        histo.GetXaxis().SetTitle(xtitle)
    #histo.GetXaxis().SetTitleSize(
    #        histo.GetXaxis().GetTitleSize()*0.7)
    #histo.GetXaxis().SetLabelSize(
    #        histo.GetXaxis().GetLabelSize()*0.75)
    #histo.GetXaxis().SetTitleOffset(
    #        histo.GetXaxis().GetTitleOffset()*0.8)

    if ytitle:
        histo.GetYaxis().SetTitle(ytitle)
    #histo.GetYaxis().SetTitleSize(
    #        histo.GetYaxis().GetTitleSize()*0.7)
    #histo.GetYaxis().SetLabelSize(
    #        histo.GetYaxis().GetLabelSize()*0.75)
    #histo.GetYaxis().SetTitleOffset(
    #        histo.GetYaxis().GetTitleOffset()*1.7)

    histo.Draw("%s"%(drawOptions))

    if LogY:
        canvas.SetLogy(1)

    InsertText(insertText)

    if save:
        can.SaveAs(outdir+"/"+canvas_name+oFext)
    return deepcopy(can)

def saveCanvasesToPDF(canvases=[], pdf_name='my_canvases.pdf'):
    outfile = r.TFile('%s'%(pdf_name),"RECREATE")
    for canvas in canvases:
        canvas.Print('%s'%(pdf_name))
    outfile.Close()

def buildLegend(x1=0.7,y1=0.7,x2 =0.9,y2=0.9,clear_legend=True):
    legend = r.TLegend(x1,y1,x2,y2)
    # Set the legend to transparent (clear) if the option is specified
    if clear_legend:
        legend.SetFillStyle(0)
        legend.SetFillColor(0)
        legend.SetLineColor(0)
        legend.SetBorderSize(0)
    return legend

def draw_centered_title(canvas, title_text):
    title = r.TLatex()
    #title.SetTextSize(text_size) 

    # Calculate the position for centering the title at the top
    canvas.Update()
    title_x = -0.1 + canvas.GetLeftMargin() + (1 - canvas.GetRightMargin() - canvas.GetLeftMargin()) / 2.0
    title_y = 1 - canvas.GetTopMargin() + 0.01 

    # Draw the title
    title.DrawLatexNDC(title_x, title_y, title_text)

def makeLatexTable(data = []):
    table = tabulate(data, headers="firstrow", tablefmt="latex")
    return table



