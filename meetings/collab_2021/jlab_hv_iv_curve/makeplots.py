import numpy as np
import matplotlib.pyplot as plt
import ROOT as r
import sys
np.set_printoptions(threshold=sys.maxsize)

colors = [r.kRed, r.kOrange-3, r.kYellow+2, r.kSpring+9, r.kGreen+3,r.kTeal+1,r.kCyan+3,r.kAzure+10,r.kBlue+1,r.kViolet-6,r.kMagenta,r.kPink+2]
lines = [1, 3, 9, 10]
styles = []
for line in lines:
    for color in colors:
        style = [line,color]
        styles.append(style)

outfile = r.TFile("./iv_curve_out.root","RECREATE")
voltMap = {}
currentMap = {}
date = ""
with open('svtBiasMainIV_curve_20211105.txt') as f:
    header = f.readline().split()
    #print(header)
    ncol = len(header)
    data = []
    for i in range(ncol):
        data.append([])
    lines = f.readlines()[1:]
    for x in lines:
        date = x.split()[0]
        for i,c in enumerate(x.split()[1:]):
            data[i].append(c)

    for i, c in enumerate(header):
        data[i].insert(0,c)

    for row in data[1:]:
        key = row[0]
        vals = np.array((row[1:]), dtype=float)
        if "i_rd" in key:
            currentMap[key] = vals

        elif "v_sens" in key:
            voltMap[key] = vals

outfile.cd()
graphs = []
mg = r.TMultiGraph()
styleiter = 0
for vkey in voltMap:
    print(vkey)
    for ikey in currentMap:
        if vkey.replace("v_sens",'') in ikey:
            voltages = []
            currents = []
            voltArray = np.array(np.array((voltMap[vkey]), dtype = int), dtype = float)
            currentArray = (currentMap[ikey])

            repeat = [-9999]
            maxindex = 99999
            for i,v in enumerate(voltArray):
                print("V: ", v, " | I: ", currentArray[i])
                if v in repeat:
                    continue
                if v < max(repeat):
                    #maxindex = max(voltArray.where(max(repeat)))
                    maxindex = np.where(voltArray == max(repeat)))
                    continue
                repeat.append(v)

                #indices = [index for index, element in enumerate(voltArray) if element == v]
                indices = [index for index, element in enumerate(voltArray) if element == v and currentArray[index] > 0 and index < maxindex]
                #Smooth out spikes
                index = max(indices)
                voltage = voltArray[index]
                current = currentArray[index]
                print("At Max Index: ", voltArray[index], " : ", currentArray[index])
                print("Before Max Index: ", voltArray[index-1], " : ", currentArray[index-1])
                print("After Max Index: ", voltArray[index+1], " : ", currentArray[index+1])
                if(index > 1 and index < len(currentArray)-1):
                    #if currentArray[index-1] < 5*current:
                    if current > 2*currentArray[index-1]:
                        #voltage = voltArray[index-1]
                        current = currentArray[index-1]
                        print("Change A: ", current)

                    elif current > 2*currentArray[index+1]:
                        #voltage = voltArray[index+1]
                        current = currentArray[index+1]
                        print("Change B: ", current)

                print("Selected: ", voltage, " : ", current)
                voltages.append(voltage)
                currents.append(current)
                #voltages.append(voltArray[index])
                #currents.append(currentArray[index])
            for j,volt in enumerate(voltages):
                print("CORRECTED V: ", volt, " | I: ", currents[j])
                    
            gr = r.TGraph(len(voltages), np.array(voltages, dtype = float), np.array(currents, dtype = float))
            gr.SetLineWidth(2)
            gr.SetLineStyle(int(styles[styleiter][0]))
            gr.SetLineColor((styles[styleiter][1]))
            gr.GetYaxis().SetRangeUser(0.0,0.000035)
            styleiter = styleiter + 1
            gr.SetName("%s"%(vkey.replace("v_sens",'')))
            gr.SetTitle("%s;HV Bias (V);Readout Current (I)"%(vkey.replace("v_sens",'')))
            graphs.append(gr)
            mg.Add(gr, 'lp')

c1 = r.TCanvas("iv_curve_custom_colz","iv_curve_custom_colz", 1800,800)
c1.cd()
for i,graph in enumerate(graphs):
    if i < 1:
        graph.Draw()
    else:
        graph.Draw("same")
    graph.Write()
c1.SetTitle("SVT_MAIN_HV_IV_Curve")
c1.Draw()
legend = c1.BuildLegend(0.13,0.24,0.25,0.88)
legend.SetEntrySeparation(0.8)
legend.SetTextSize(0.015)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.Draw()
c1.Write()
c1.Close()

c = r.TCanvas("title","title", 1800,800)
c.cd()
mg.GetXaxis().SetRange(0,350)
mg.GetXaxis().SetTitle("HV Bias Voltage (V)")
mg.GetYaxis().SetTitle("Readout Current (A)")
mg.Draw("A pmc plc")
legend = c.BuildLegend(0.13,0.24,0.25,0.88)
legend.SetEntrySeparation(0.8)
legend.SetTextSize(0.015)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
c.Draw()
c.Write()
c.SaveAs("./test.png")


