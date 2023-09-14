#!/usr/bin/python3
import sys
sys.path.append( '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/ana_scripts/plot_utils')
import my_utilities as utils
import ROOT as r

r.gROOT.SetBatch(1)
colors = utils.getColors()

lumi_2016 = 10.7 #pb-1

tritrig_name = '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/kf_mc/tritrig_beam/tritrig_kf_ana_05122023/final_hadd_tritrigv2-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_KF_ana.root'
tritrig_events = 50000*9853 #nevents
tritrig_xsec = 1.416e9 #pb

wab_name = '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/kf_mc/wab_beam/wab_kf_ana_05022023/output/final_hadd_wabv3-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_KF_ana.root'
wab_events = 100000*9966 #nevents
wab_xsec = 0.1985e12 #pb

data_name = '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/kf_data/kf_041823/output/hadd_hps_BLPass4_1958_files_recon_4.2_ana_kf.root'
total_lumi = lumi_2016*0.1 #pb-1

#Read histograms
keywords = ["_Psum_","_InvM_h","_ele_p_h","pos_p_h","ele_TanLambda_h","pos_TanLambda_h","ele_Phi_h","pos_Phi_h","vtx_Z_h","_cutflow"]
#keywords = [""]
#vtxselection = 'Tight_2016_simp_reach_dev_CR'
vtxselection = 'Tight_2016_simp_reach_dev'
#vtxselection = 'vtxSelection'

#Load root files
data = r.TFile(data_name,"READ")
tritrig = r.TFile(tritrig_name,"READ")
wab = r.TFile(wab_name,"READ")

kfselection = 'vtxana_kf_%s'%(vtxselection)
tritrig_h = utils.read1DPlotsFromRootDir(tritrig,kfselection,keywords)
wab_h = utils.read1DPlotsFromRootDir(wab,kfselection,keywords)
data_h = utils.read1DPlotsFromRootDir(data,kfselection,keywords)

for i,pname in enumerate(data_h):
    canv_name = 'KF'+pname.replace("vtxana_kf",'')

    #histos
    data = data_h[pname]
    print("Data Entries: ", data.GetEntries())
    name = 'data'+data.GetName().replace("vtxana_kf",'')
    utils.format1DPlot(data,name,title=name,linewidth=1,linecolor=1,markerstyle=8,markersize=0.5)

    tritrig = tritrig_h[pname]
    print("Tritrig Entries: ", tritrig.GetEntries())
    print("Tritrig Scaled Entries: ", tritrig.GetEntries()*(tritrig_xsec*total_lumi/(tritrig_events)))
    tritrig.Scale(tritrig_xsec*total_lumi/(tritrig_events))
    name = 'tritrig'+tritrig.GetName().replace("vtxana_kf",'')
    utils.format1DPlot(tritrig,name,title=name,linecolor=colors[1])

    wab = wab_h[pname]
    print("WAB Entries: ", wab.GetEntries())
    print("WAB Entries: ", wab.GetEntries()*wab_xsec*total_lumi/(wab_events))
    wab.Scale(wab_xsec*total_lumi/(wab_events))
    name = 'wab'+wab.GetName().replace("vtxana_kf",'')
    utils.format1DPlot(wab,name,title=name,linecolor=colors[2])

    tritrig_wab = tritrig.Clone('tritrig_clone')
    tritrig_wab.Add(wab,1)
    name = 'tritrig+WAB'+data.GetName().replace("vtxana_kf",'').replace("data",'')
    utils.format1DPlot(tritrig_wab,name,title=name,linecolor=colors[5])
    print("Tritrig+WAB scaled: ", tritrig_wab.GetEntries())

    plots = [data,tritrig,wab,tritrig_wab]
    #utils.overlay1DPlots(plots,canv_name,'/home/alic/HPS/projects/simp_analysis/plots/sample_0_comparisons/081622/vertexSelection/%s'%(kfselection),colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])
    utils.overlay1DPlots(plots,canv_name,'/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/ana_scripts/plots/',colors=[1,colors[1],colors[2],colors[5]],drawStyles=["Ehist","hist","hist","hist"])


    #ratio plots
    ratioPairs = [(0,3)]
    ratioNames = ["Data:Tritrig+Wab+Beam"]
    #utils.makeRatioPlot(plots,ratioPairs,ratioNames,"ratios_"+canv_name, '/home/alic/HPS/projects/simp_analysis/plots/sample_0_comparisons/081622/vertexSelection/%s'%(kfselection),colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])
    utils.makeRatioPlot(plots,ratioPairs,ratioNames,"ratios_"+canv_name, '/sdf/group/hps/users/alspellm/projects/collaboration_meetings/may_2023/ana_scripts/plots/',colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])


#outfile.Close()

#hist_units_dic
#tritrig_beam Scale 1.416e9*Lumi/(50000*n)
#data_scale -> total lumi for 7800 = 167.17 (1/nb) for 388 files. I have 387 files.

#invMassHistos['rad'].Scale(units*66.36e6*Lumi/(rebin*10000*9967))
#invMassHistos['tritrig'].Scale(units*1.416e9*Lumi/(rebin*50000*9853))
#invMassHistos['wab'].Scale(units*0.1985e12*Lumi/(rebin*100000*9966))
