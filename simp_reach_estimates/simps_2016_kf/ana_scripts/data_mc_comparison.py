#!/usr/bin/python3
import my_utilities as utils
import ROOT as r

r.gROOT.SetBatch(1)
colors = utils.getColors()

mc_gbl_name = '/home/alic/HPS/projects/simp_analysis/mc/tritrig_beam/ana/CR/hadd_tritrigv2-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_ana_GBL.root'
mc_kf_name = '/home/alic/HPS/projects/simp_analysis/mc/tritrig_beam/ana/CR/hadd_tritrigv2-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_ana_KF.root'

data_gbl_name = '/home/alic/HPS/projects/simp_analysis/data/7800/hadd_hps_007800_ana_GBL.root'
data_kf_name = '/home/alic/HPS/projects/simp_analysis/data/7800/hadd_hps_007800_ana_KF.root'

mc_gbl = r.TFile(mc_gbl_name,"READ")
mc_kf = r.TFile(mc_kf_name,"READ")

data_gbl = r.TFile(data_gbl_name,"READ")
data_kf = r.TFile(data_kf_name,"READ")

#Read histograms
keywords = ["_Psum_","_InvM_h"]
mc_gbl_plots1d = utils.read1DPlotsFromRootDir(mc_gbl,"vtxana_gbl_Tight_noPsum",keywords) 
mc_kf_plots1d = utils.read1DPlotsFromRootDir(mc_kf,"vtxana_kf_Tight_noPsum",keywords) 

data_gbl_plots1d = utils.read1DPlotsFromRootDir(data_gbl,"vtxana_gbl_Tight_noPsum",keywords) 
data_kf_plots1d = utils.read1DPlotsFromRootDir(data_kf,"vtxana_kf_Tight_noPsum",keywords) 

Lumi = (167.17/388)*387*0.001 #1/pb
#units = 1/11.66
units = 1/12

for name in mc_gbl_plots1d:
    print(name)
    cname = 'mc_data_comp_'+name.replace('vtxana_gbl_Tight_','')
    plot1 = mc_gbl_plots1d[name]
    plot1.Scale(1.416e9*Lumi*units/(50000*4163))
    #plot1.Scale(0.01)
    name1 = 'tritrig_mc_'+plot1.GetName().replace('vtxana_','')
    utils.format1DPlot(plot1, name1,title=name1,linecolor=colors[0])

    plot2 = data_gbl_plots1d[name]
    name2 = 'data_'+plot2.GetName().replace('vtxana_','')
    plot2.SetName(name2)
    utils.format1DPlot(plot2, name2,title=name2,linecolor=colors[1], linestyle=1)

    plot3 = mc_kf_plots1d[name.replace('gbl','kf')]
    plot3.Scale(1.416e9*Lumi*units/(50000*4163))
    name3 = 'tritrig_mc_'+plot3.GetName().replace('vtxana_','')
    utils.format1DPlot(plot3, name3,title=name3,linecolor=colors[2])

    plot4 = data_kf_plots1d[name.replace('gbl','kf')]
    name4 = 'data_'+plot4.GetName().replace('vtxana_','')
    plot4.SetName(name4)
    utils.format1DPlot(plot4, name4,title=name4,linecolor=colors[3], linestyle=1)

    plots = [plot1,plot2,plot3,plot4]
    utils.overlay1DPlots(plots,cname,'/home/alic/HPS/projects/simp_analysis/scripts/plots/images')

#Scaling
#tritrig_beam Scale 1.416e9*Lumi/(50000*n)
#data_scale -> total lumi for 7800 = 167.17 (1/nb) for 388 files. I have 387 files.

#invMassHistos['rad'].Scale(units*66.36e6*Lumi/(rebin*10000*9967))
#invMassHistos['tritrig'].Scale(units*1.416e9*Lumi/(rebin*50000*9853))
#invMassHistos['wab'].Scale(units*0.1985e12*Lumi/(rebin*100000*9966))
