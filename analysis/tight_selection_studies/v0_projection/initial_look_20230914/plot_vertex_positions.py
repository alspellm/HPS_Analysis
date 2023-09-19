import ROOT as r
import numpy as np
import plot_utilities as utils
import root_numpy as rnp
import pandas as pd
import os
import re
import glob as glob
import json

db_conds = '/sdf/group/hps/users/alspellm/src/hpstr/analysis/data/beamspot_positions_2016.json'

conds = json.load(open(db_conds))

runs = []
xpositions = []
ypositions = []
for run in conds:
    print(run)
    runs.append(int(run))
    x = conds[run]['beamspot_x']
    y = conds[run]['beamspot_y']
    xpositions.append(x)
    ypositions.append(y)

style = utils.SetMyStyle()
colors = utils.getColorsHPS()

x_gr = r.TGraph(len(runs), np.array(runs,dtype=float), np.array(xpositions,dtype=float))
utils.formatHisto(x_gr, name='x_pos_db_corr',title='<x> pos db corr',x_label='Run Number', y_label='Database Position Correction [mm]',marker_color=colors[0], marker_size=2.0)
x_gr.GetYaxis().SetRangeUser(-0.3,0.3)

y_gr = r.TGraph(len(runs), np.array(runs,dtype=float), np.array(ypositions,dtype=float))
utils.formatHisto(y_gr, name='y_pos_db_corr',title='<y> pos db corr',x_label='Run Number', y_label='Database Position Correction [mm]',marker_color=colors[1], marker_size=2.0)
y_gr.GetYaxis().SetRangeUser(-0.3,0.3)

canvas = r.TCanvas('c','c',2400,1440)
canvas.cd()
canvas.SetGrid(1,1)
x_gr.Draw("AP")
y_gr.Draw("PSAME")
legend = canvas.BuildLegend()
legend.SetFillStyle(0)

#canvas.Update()
canvas.SaveAs('test_plot.png')
canvas.Close()

del style



