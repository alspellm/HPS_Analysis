#!/usr/bin/python3
import ROOT as r
import numpy as np
import root_numpy as rnp
import copy
import pandas as pd


def massRes(mass):
    res = 9.73217e-01 + 3.63659e-02*mass + -7.32046e-05*mass*mass #2016 simps alic
    return res

#mass window
mass = 55.0
#massRes = 3.0 #MeV for a 55MeV Vd
massRes = massRes(mass)
lowMass = float(mass) - 2.0*massRes/2.0
highMass = float(mass) + 2.0*massRes/2.0
print("Mass Window: ", lowMass, "--", highMass)

#define mass window selection
selection = {'unc_vtx_mass_lt' : highMass/1000., 'unc_vtx_mass_gt' : lowMass/1000. , 'unc_vtx_psum_gt' : 1.0}

#data
signal_infile = '/sdf/group/hps/users/alspellm/projects/THESIS/cut_dev/ZBi/signal/hadd_mass_55_simp_recon_KF_ana.root'
signal_treename = 'vtxana_kf_radMatchTight_2016_simp_reach_dev/vtxana_kf_radMatchTight_2016_simp_reach_dev_tree'
signal_arr = rnp.root2array(signal_infile, signal_treename)
df = pd.DataFrame(signal_arr)
df = df[(df['unc_vtx_mass'] < selection['unc_vtx_mass_lt']) & (df['unc_vtx_mass'] > selection['unc_vtx_mass_gt']) & (df['unc_vtx_psum'] > selection['unc_vtx_psum_gt']) ]

Lumi = 10.7*0.1 #Data represents ~10% lumi
mcScale = {}
mcScale['tritrig'] = 1.416e9*Lumi/(50000*9853) #pb2016
mcScale['wab'] = 0.1985e12*Lumi/(100000*9966) #pb2016

outfile = r.TFile('%s_impactparameter_cut.root'%(mass),"RECREATE")

###################################################################

signalPlots2D ={}
signalPlots1D = {}
#z0 vs recon_z
signalPlots2D["signal_z0_v_reconz_hh"] = r.TH2F("signal_z0_v_reconz_hh","signal_z0_v_reconz;recon_z [mm];z0 [mm]",220,-30.0,80.0,800,-4.0,4.0)
signalPlots1D["signal_impact_parameter_up_h"] = r.TH1F("signal_impact_parameter_up_h","signal_impact_parameter_up_h;recon_z[mm];events",220, -30.0, 80.0)
signalPlots1D["signal_impact_parameter_down_h"] = r.TH1F("signal_impact_parameter_down_h","signal_impact_parameter_down_h;recon_z[mm];events",220, -30.0, 80.0)
#####################################################################
alpha = 0.15
#Fill impact parameter distribution
for index, row in df.iterrows():
    signalPlots2D["signal_z0_v_reconz_hh"].Fill(row['unc_vtx_z'], row['unc_vtx_ele_track_z0']) 
    signalPlots2D["signal_z0_v_reconz_hh"].Fill(row['unc_vtx_z'], row['unc_vtx_pos_track_z0']) 

hh = signalPlots2D["signal_z0_v_reconz_hh"]
for i in range(hh.GetNbinsX()):
    projy = hh.ProjectionY("projy_bin_%i"%(i+1),i+1,i+1)
    if projy.GetEntries() < 1:
        continue

    #z0 > 0.0
    start_bin = projy.FindBin(0.0)
    end_bin = projy.FindLastBinAbove(0.0)
    refIntegral = projy.Integral(start_bin,end_bin)
    cutz0_bin = start_bin
    testIntegral = refIntegral
    while(testIntegral > (1.0 - alpha)*refIntegral and cutz0_bin < end_bin-1):
        cutz0_bin = cutz0_bin + 1
        testIntegral = projy.Integral(cutz0_bin, end_bin)

    cutz0_up = projy.GetXaxis().GetBinLowEdge(cutz0_bin)
    for j in range(int(refIntegral)):
        signalPlots1D["signal_impact_parameter_up_h"].Fill(hh.GetXaxis().GetBinCenter(i+1), cutz0_up/refIntegral)

    #z0 < 0.0
    end_bin = projy.FindFirstBinAbove(0.0)
    start_bin = projy.FindBin(0.0)
    refIntegral = projy.Integral(end_bin, start_bin)
    cutz0_bin = start_bin
    testIntegral = refIntegral
    while(testIntegral > (1.0 - alpha)*refIntegral and cutz0_bin > end_bin+1):
        cutz0_bin = cutz0_bin - 1
        testIntegral = projy.Integral(end_bin, cutz0_bin)

    cutz0_down = projy.GetXaxis().GetBinUpEdge(cutz0_bin)
    for j in range(int(refIntegral)):
        signalPlots1D["signal_impact_parameter_down_h"].Fill(hh.GetXaxis().GetBinCenter(i+1), cutz0_down/refIntegral)

outfile.cd()
#Fit upper cut
fitFunc_up = r.TF1("linear_fit","[0]*(x-[1])",5.0,70.0)
fitResult_up = signalPlots1D["signal_impact_parameter_up_h"].Fit("linear_fit","QS","",5.0,70.0)
fitFunc_up.Draw()
slope_up = fitResult_up.GetParams()[0]
xinter_up = fitResult_up.GetParams()[1]


#Fit lower cut
fitFunc_down = r.TF1("linear_fit","[0]*(x-[1])",5.0,70.0)
fitResult_down = signalPlots1D["signal_impact_parameter_down_h"].Fit("linear_fit","QS","",5.0,70.0)
fitFunc_down.Draw()
slope_down = fitResult_down.GetParams()[0]
xinter_down = fitResult_down.GetParams()[1]

xintercept = .5*(xinter_up+xinter_down)


print("%s MeV Signal Impact Parameter Up Cut: Slope = %f, X-intercept = %f"%(mass, slope_up, xinter_up))
print("%s MeV Signal Impact Parameter down Cut: Slope = %f, X-intercept = %f"%(mass, slope_down, xinter_down))
print("Average X-intercept = %f"%(xintercept))

for plot in signalPlots2D.values():
    plot.Write()
for plot in signalPlots1D.values():
    plot.Write()

'''
55.0 MeV Signal Impact Parameter Up Cut: Slope = 0.027313, X-intercept = 3.833276
55.0 MeV Signal Impact Parameter down Cut: Slope = -0.027450, X-intercept = 4.469586
Average X-intercept = 4.151431
'''
'''
55.0 MeV Signal Impact Parameter Up Cut: Slope = 0.027572, X-intercept = 3.705004
55.0 MeV Signal Impact Parameter down Cut: Slope = -0.027837, X-intercept = 4.498615
Average X-intercept = 4.101810
'''
'''
5.0 MeV Signal Impact Parameter Up Cut: Slope = 0.029816, X-intercept = 3.346113
55.0 MeV Signal Impact Parameter down Cut: Slope = -0.029530, X-intercept = 3.597638
Average X-intercept = 3.471875
'''
