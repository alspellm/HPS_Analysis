#!/usr/bin/python3
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import pandas as pd

plots = {}

colors = ["black", "r", "g", "b", "m", "gray", "orange", "saddlebrown", "darkslategray", "aqua", "purple", "indigo", "springgreen"]
markers = ['.', 'x', '<', '>', '*', 'D']
cases = ["CNM", "initial", "UVjig", "UVshort", "baked"]

layerMap = {"L0M01" : "L2T_stereo", "L0M02" : "L1B_stereo", "L0M03" : "L1T_stereo", "L0M04" : "L1B_axial", "L0M05" : "L2T_axial", "L0M06" : "L1T_axial", "L0M08" : "L2B_stereo", "L0M09" : "L2B_axial" }


case_styles = {}
for i,case in enumerate(cases):
    case_styles[case] = markers[i]


for infile in sorted(glob.iglob("*.xls*")):
    name = os.path.splitext(infile)[0]
    ext = os.path.splitext(infile)[1]
    module = name.split("_")[1]
    wafer = name.split("_")[0]
    case = name.split("_")[2]

    x,y,n = [0], [0], 0
    if ext == ".xlsx":
        df = pd.read_excel(infile, header=None, engine='openpyxl')
        n = len(df.index)
        x = np.double(df.iloc[:, 0].to_numpy())
        y = df.iloc[:,1].to_numpy()
    elif ext == ".xls":
        df = pd.read_excel(infile)
        n = len(df.index)
        x = np.double(-1 * df.iloc[:, 1].to_numpy())
        y = -1*(df.iloc[:,0].to_numpy())

    plots[module+" "+wafer+" "+case] = np.transpose(np.array(list(zip(x,y))))

for key in plots:
    module = int((key.split(" ")[0][3:5]))
    string = layerMap[key.split(" ")[0]]
    wafer = key.split(" ")[1]
    case = key.split(" ")[2]
    col = colors[module]
    style = case_styles[case]
    plt.plot(plots[key][0], plots[key][1], '%s-'%(style), color=col, markersize = 12, label='%s'%(string+" "+ wafer+ " "+case))
    plt.legend(prop={'size' : 10})

plt.title("Installed L0 Wafer IV Curves", fontsize=20)
plt.xlabel('Volts (V)', fontsize=20)
plt.ylabel('Leakage Current (A)', fontsize=20)

plt.xticks([x for x in range(0,250,10)])
plt.yticks([y for y in range(0,25,1)])
plt.yscale('log')

plt.grid()
plt.show()


