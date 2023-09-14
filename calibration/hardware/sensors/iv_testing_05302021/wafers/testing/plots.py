#!/usr/bin/python3
import pandas as pd
import openpyxl
import glob
import ROOT as r
import numpy as np
import os

import utilities as utils



def parseInputName():

    basename = os.path.basename(filename)


myMarkers = [1, r.kFullCircle, r.kFullTriangleUp, r.kFullSquare, r.kFullStar, r.kFullDiamond, r.kFullCross]

directory = '/home/alic/HPS/projects/hardware/sensors/iv_testing_05302021/wafers/testing'

dotw = ["Mon", "Tue", "Wed", "Thur", "Fri", "Sat", "Sun"]

graphs = []
legends = []
modifierLineStyles = []
sensorColorStyles = []
sensorMarkerStyles = []
wafers = ["W04","W07","W08","W12","W14"]
nwafers = len(wafers)
plotnames = [None] * nwafers

for i in range(nwafers):
    graphs.append([])
    legends.append([])
    modifierLineStyles.append([])
    sensorColorStyles.append([])
    sensorMarkerStyles.append([])

colors = []
modifiers = ["baked", "UVjig", "UVshort", "CNM"]
check = []

for filename in sorted(glob.iglob("%s/*"%(directory))):
    print("file: ",filename)
    basename = os.path.basename(filename).split("_",5)
    print("basename: ", basename)
    sensorid = basename[0]+"_"+basename[1]
    print("sensor: ", sensorid)

    #Get Wafer and Sensor number
    wafer = sensorid.split("_",2)[0]
    sensor = sensorid.split("_",2)[1]
    mod = basename[2]
    print("wafer_sensor: ", wafer+sensor)
    print("mod: ", mod)

    #Set Plot line styles, markers, and colors based on Wafer, Sensor, Modifier
    linestyle = 1
    linecolor = int(sensor.split("S",2)[1])
    modindex = -1
    for modifier in modifiers:
        if mod in modifier:
            modindex = modifiers.index(mod)
            break
    if modindex < 0:
        markerstyle = myMarkers[0] 
    else:
        markerstyle = myMarkers[modindex]

    if modindex < 0:
        name = wafer+"_"+sensor
    else:
        name = wafer+"_"+sensor+"_"+modifiers[modindex]

    check.append(name)
    count = 0
    for e in check:
        if e == name:
            count = count + 1
    if count > 1:
        name = name+"_"+str(count)
    print("plot name: ", name)
    legend = name
    continue
    
    #Read in data from xlsx file
    df = pd.read_excel(filename, header=None,engine='openpyxl')
    n = len(df.index)
    x = np.double(df.iloc[:, 0].to_numpy())
    y = df.iloc[:,1].to_numpy()
    gr = r.TGraph(n, x, y)
    gr.SetTitle("Wafer%s_iv_curves"%(wafer))

    #Append graphs by wafer number
    graphs[int(wafer)-1].append(gr)
    legends[int(wafer)-1].append(legend)
    modifierLineStyles[int(wafer)-1].append(linestyle)
    sensorColorStyles[int(wafer)-1].append(linecolor)
    sensorMarkerStyles[int(wafer)-1].append(markerstyle)

    if wafer+"_"+"iv_curves" not in plotnames:
        plotnames[int(wafer)-1] = ("wafer%s_iv_curves"%(wafer))

for i, (sensorLegends, sensorGraphs, linestyles, sensorColors, markerStyles) in enumerate(zip(legends,graphs, modifierLineStyles, sensorColorStyles, sensorMarkerStyles)):
    if plotnames[i] is None:
        continue
    
    utils.Make1DGraphs(plotnames[i], "./", sensorGraphs,sensorColors,markerStyles, sensorLegends, ".png",xtitle="Voltage (V)",ytitle="Leakage Current (A)", ymax=0.000025, ymin = 0.0 , legLocation=[0.7,0.7], xlimits=[0,200], lineStyles=linestyles, LogY=True)


