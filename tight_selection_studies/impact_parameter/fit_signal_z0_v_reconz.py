#!/usr/bin/python3
import ROOT as r
import numpy as np
import copy

def getHistoFromFile(inFile, name = ''):
    inFile.cd()
    for key in inFile.GetListOfKeys():
        plotName = key.GetName()
        print(plotName)
        if (plotName.find(name) != -1):
            print("FOUND ", plotName)
            return inFile.Get(plotName)

def getHistosFromFile(inFile, histoType = "TH2D", name = ""):
    inFile.cd()
    histos = []
    for key in inFile.GetListOfKeys():
        plotName = key.GetName()
        plotType = key.GetClassName()
        if plotType.find(histoType) != -1 and (plotName.find(name) != -1):
            histos.append(inFile.Get(plotName))
    return histos

def linearFitZ0vReconZ(outfile, z0_v_reconz_hh, gr_name="graph"):
    means = []
    reconz = []
    entries = []
    for i in range(1, z0_v_reconz_hh.GetNbinsX()+1, 3):
        projy_h = z0_v_reconz_hh.ProjectionY("proj_y_%i"%(i), i, i+2,"")
        if(projy_h.GetEntries() < 40):
            continue
        entries.append(projy_h.GetEntries())
        fit = projy_h.Fit("gaus","QES","",projy_h.GetBinCenter(projy_h.FindFirstBinAbove(0.0)),projy_h.GetBinCenter(projy_h.FindLastBinAbove(0.0)))
        means.append(fit.Parameter(1))

        z_avg = 0.0
        for j in range(1,3+1):
            z_avg = z_avg + z0_v_reconz_hh.GetXaxis().GetBinCenter(i)
        z_avg = z_avg/3.0
        reconz.append(z_avg)

    #Convert to np arrays
    means = np.array(means, dtype=float)
    reconz = np.array(reconz,dtype=float)
    errors = [1/(x) for x in entries]
    errors = np.array(errors,dtype=float)

    outfile.cd()

    #Plot Tgraph
    gr = r.TGraphErrors(len(means), reconz, means,np.zeros(len(means)),errors )
    gr.SetName("%s"%(gr_name))

    fit_func = r.TF1("linear_fit","[0]*x + [1]", -40.0,60.0)
    fit = gr.Fit("linear_fit","QES","",-40.0,60.0)
    fit.Draw()
    gr.Write()

    print("%s fit results: "%(gr_name))
    print("slope = %f"%(fit.Parameter(0)))
    print("y-inter = %f"%(fit.Parameter(1)))

###################################
outfile = r.TFile('z0_v_zrecon_signal_fit.root',"RECREATE")

#signal in
signal_infile = r.TFile('/sdf/group/hps/users/alspellm/projects/THESIS/impact_param_studies_04182023/test/simp_impact_param_signal_04182023.root',"READ")
ele_z0_v_reconz_top_hh = getHistoFromFile(signal_infile, name ='ele_z0_v_reconz_top')
pos_z0_v_reconz_bot_hh = getHistoFromFile(signal_infile, name ='pos_z0_v_reconz_bot')
ele_z0_v_reconz_bot_hh = getHistoFromFile(signal_infile, name ='ele_z0_v_reconz_bot')
pos_z0_v_reconz_top_hh = getHistoFromFile(signal_infile, name ='pos_z0_v_reconz_top')

#Top ele_z0_v_reconz_top
linearFitZ0vReconZ(outfile, ele_z0_v_reconz_top_hh, gr_name="signal_ele_z0_v_reconz_top_linear_fit")
#y = -0.039151*x + -0.031282  #Fit Results 04/19/2023

#Top pos_z0_v_reconz_bot
linearFitZ0vReconZ(outfile, pos_z0_v_reconz_bot_hh, gr_name="signal_pos_z0_v_reconz_bot_linear_fit")
#y = 0.039501*x + 0.004176 #Results 04/19/2023

#Top ele_z0_v_reconz_bot
linearFitZ0vReconZ(outfile, ele_z0_v_reconz_bot_hh, gr_name="signal_ele_z0_v_reconz_bot_linear_fit")
# y = 0.040086*x + 0.016186

#Top pos_z0_v_reconz_top
linearFitZ0vReconZ(outfile, pos_z0_v_reconz_top_hh, gr_name="signal_pos_z0_v_reconz_top_linear_fit")
# y = -0.037899 + -0.009400

