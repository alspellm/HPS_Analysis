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
modifiers = ["baked", "UVjig", "UVshort", "CNM"]

for filename in sorted(glob.iglob("%s/*.xlsx"%(directory))):
    print(filename)
    infile = os.path.basename(filename).split("_",3)
    basename = os.path.basename(filename)
    sensorid = "unknown"
    #Split file name by day of the week, which is standard output of Probe station
    for d in dotw:
        if d in basename:
            sensorid = basename.split(d,1)[0]
    if sensorid is "unknown":
        print("ERROR! INPUT FILE NAME FORMAT INCORRECT. NOT COMPATIBLE WITH PLOTTER")
        break

    #Get Wafer and Sensor number
    split_sensor = sensorid.split("_",3)
    wafer = split_sensor[0].replace('W','')
    sensor = split_sensor[1].replace('S','')
    mod = split_sensor[2].replace('_','')
    number = 0
    if len(split_sensor) > 3:
        if split_sensor[3] not in '':
            number = int(split_sensor[3].replace('_',''))

    #Set Plot line styles, markers, and colors based on Wafer, Sensor, Modifier
    linestyle = 1
    markerstyle = utils.markers[number]
    linecolor = int(sensor)

    if mod in modifiers:
        linestyle = 2 + modifiers.index(mod)
        if number == 0:
            legend = "W%s_S%s_%s"%(wafer, sensor, mod)
        else:
            legend = "W%s_S%s_%s_%i"%(wafer, sensor, mod, number)
    else:
        legend = "W%s_S%s_factory"%(wafer, sensor)
    
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


