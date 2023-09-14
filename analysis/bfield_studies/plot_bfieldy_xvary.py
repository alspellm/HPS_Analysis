#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import csv

#Read in svt sensor locations global (x,y,z)
sensor_file = '/sdf/group/hps/users/bravo/run/geoDump/HPS-PhysicsRun2016-Pass2.txt'
'''
#column_names = ['sensor', 'x', 'y', 'z']
#sens = pd.read_csv(sensor_file, delimiter=' ', index_col=False, names=column_names)
#print(sens.head())
sens = pd.read_csv(sensor_file, sep='^', header=None)
for index,row in sens.iterrows():
    row = row.str.replace('[','').str.replace(']','')
    row = row.str.strip()

    row = pd.DataFrame(row)
    print(row.head())
    '''
sens_pos_z = []
with open(sensor_file, 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        row = ' '.join(row)
        row = row.replace('[','').replace(']','')
        row = row.strip()
        row = row.split()
        sens_pos_z.append(float(row[3]))
print(sens_pos_z)

file_path = '209acm2_5kg_corrected_unfolded_scaled_1.04545_v4.dat'
column_names = ['x','y','z','bx','by','bz']
data = pd.read_csv(file_path, delimiter=' ', skiprows=9, index_col=False, names=column_names)

cmap = plt.cm.get_cmap('Set1')
plt.figure()

zoffset = 457.2 #mm
xmasks = [-200.0,-100.0,-75.0,-50.0,-25.0,0.0,25.0,50.0,75.0,100.0,200.0]
for e,xmask in enumerate(xmasks):
    mask = (data['x'] == xmask) & (data['y'] == 0.0)
    mask_data = data[mask]

    z_arr = []
    bfy_arr = []

    for index, row in mask_data.iterrows():
        z = row['z']
        if z < -10.0-zoffset or z > 1500.0-zoffset:
            continue
        bfy = 1000.0*row['by']
        z_arr.append(z+zoffset)
        bfy_arr.append(bfy)

    color = cmap(e % cmap.N)
    plt.plot(np.array(z_arr), np.array(bfy_arr),linewidth = 0.5, label='xpos %s'%(xmask))

for i,line in enumerate(sens_pos_z):
    plt.axvline(line, color ='black',linewidth=0.5, linestyle='--')

ymin,ymax = plt.ylim()
line_height = (ymin - ymax)/2
plt.axvline(-4.3, color='red',linewidth=1.0, label='target', ymin=ymin+line_height/(ymax-ymin), ymax=ymax-line_height /(ymax-ymin))

plt.axvline(0.0, color='black',linewidth=1.0, label='origin', ymin=ymin+line_height/(ymax-ymin), ymax=ymax-line_height /(ymax-ymin))

plt.axvline(1390.0, color='green',linewidth=1.0, label='ecal', ymin=ymin+line_height/(ymax-ymin), ymax=ymax-line_height /(ymax-ymin))

plt.xlabel('Global Z [mm]')
plt.ylabel('BField Y [T]')
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.locator_params(axis='y', nbins=20)
plt.title('2016 SVT BField Y versus Z')
plt.legend()

plt.savefig('2016_svt_bfield.png',dpi=300)

