#!/usr/bin/python3
import pandas as pd
import openpyxl
import glob
import ROOT as r
import numpy as np
import os

import utilities as utils

directory = '/home/alic/HPS/projects/hardware/sensors/iv_testing_05302021/wafers_1-6/testing'

dotw = ["Mon", "Tue", "Wed", "Thur", "Fri", "Sat", "Sun"]

graphs = []
legends = []
modifierLineStyles = []
sensorColorStyles = []
sensorMarkerStyles = []
nwafers = 6
plotnames = [None] * nwafers

for i in range(nwafers):
    graphs.append([])
    legends.append([])
    modifierLineStyles.append([])
    sensorColorStyles.append([])
    sensorMarkerStyles.append([])

colors = []
modifiers = ["CNM", "baked", "baked2", "baked-stored", "UVjig", "UVshort"]

for filename in sorted(glob.iglob("%s/*.xls*"%(directory))):
    print(filename)
    infile = os.path.basename(filename).split("_",3)
    basename = os.path.basename(filename)
    sensorid = "unknown"
    #Split file name by day of the week, which is standard output of Probe station
    for d in dotw:
        if d in basename:
            sensorid = basename.split(d,1)[0]
    if "CNM" in basename:
        sensorid = basename.split(".")[0]
    print(sensorid)
            
    if sensorid is "unknown":
        print("ERROR! INPUT FILE NAME FORMAT INCORRECT. NOT COMPATIBLE WITH PLOTTER")
        break
    #Get Wafer and Sensor number
    split_sensor = sensorid.split("_",3)
    print(split_sensor)
    wafer = split_sensor[0].replace('W','')
    sensor = split_sensor[1].replace('S','')
    mod = split_sensor[2].replace('_','')
    number = 0
    print(wafer, sensor, mod)
    if len(split_sensor) > 3:
        if split_sensor[3] not in '':
            number = int(split_sensor[3].replace('_',''))

    #Set Plot line styles, markers, and colors based on Wafer, Sensor, Modifier
    linestyle = 1
    markerstyle = utils.myMarkers[number]
    linecolor = utils.mycolors[int(sensor)]

    print("is modifier ", mod, "in list: ", modifiers)
    if mod in modifiers:
        linestyle = utils.myLineStyles[1 + modifiers.index(mod)]
        if number == 0:
            legend = "W%s_S%s_%s"%(wafer, sensor, mod)
        else:
            legend = "W%s_S%s_%s_%i"%(wafer, sensor, mod, number)
    else:
        legend = "W%s_S%s_reception"%(wafer, sensor)
    
    #Read in data from xlsx file
    x = [0]
    y = [0]
    n = 0
    if os.path.splitext(basename)[1] == '.xlsx':
        df = pd.read_excel(filename, header=None,engine='openpyxl')
        n = len(df.index)
        x = np.double(df.iloc[:, 0].to_numpy())
        y = df.iloc[:,1].to_numpy()
    elif os.path.splitext(basename)[1] == '.xls':
        df = pd.read_excel(filename)
        n = len(df.index)
        x = np.double(-1 * df.iloc[:, 1].to_numpy())
        y = -1*(df.iloc[:,0].to_numpy())

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
    
    utils.Make1DGraphs(plotnames[i], "./", sensorGraphs,sensorColors,markerStyles, sensorLegends, ".png",xtitle="Voltage (V)",ytitle="Leakage Current (A)", ymax=0.00015, ymin = 0.0 , legLocation=[0.78,0.65, 0.899,0.2], xlimits=[0,220], lineStyles=linestyles, LogY=True)


