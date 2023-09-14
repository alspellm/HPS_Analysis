#!/usr/bin/python3
import numpy as np
import glob
import os
import matplotlib.pyplot as plt

layerMap = {"L0M01" : "L2T_stereo", "L0M02" : "L1B_stereo", "L0M03" : "L1T_stereo", "L0M04" : "L1B_axial", "L0M05" : "L2T_axial", "L0M06" : "L1T_axial", "L0M08" : "L2B_stereo", "L0M09" : "L2B_axial" }
plots = {}
for infile in sorted(glob.iglob("*.dat")):
    module = infile.split("_")[0]
    print(module)

    data = []
    with open(infile) as f:
        lines = (line for line in f if not line.startswith('V'))
        data = np.transpose(np.loadtxt(lines))
        plots[module] = data

for module in plots:
    plt.plot(plots[module][0], plots[module][1], 'o-', label='%s'%(layerMap[module]))
    plt.legend(prop={'size' : 20})

plt.title("2021 Installed Slim Edge HalfModule IV Curves", fontsize=20)
plt.xlabel('Volts (V)', fontsize=20)
plt.ylabel('Leakage Current (uA)', fontsize=20)

plt.xticks([x for x in range(0,310,10)])
plt.yticks([y for y in range(0,25,1)])
#plt.yscale('log')

plt.grid()
plt.show()


