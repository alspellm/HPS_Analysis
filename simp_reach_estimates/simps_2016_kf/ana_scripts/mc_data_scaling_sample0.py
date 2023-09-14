#!/usr/bin/python3
import sys
sys.path.append( '/sdf/group/hps/users/alspellm/projects/THESIS/simp_reach_estimates/simps_2016_kf/ana_scripts/plot_utils')
import my_utilities as utils
import ROOT as r

r.gROOT.SetBatch(1)
colors = utils.getColors()

runKF = True
runGBL = False

#tritrig_kf_name = '/home/alic/HPS/projects/simp_analysis/mc/tritrig_beam/ana/hadd_tritrigv2-beamv6_rerecon_ana_KF.root'
#tritrig_gbl_name = '/home/alic/HPS/projects/simp_analysis/mc/tritrig_beam/ana/hadd_tritrigv2-beamv6_rerecon_ana_GBL.root'

#Old SeedTracker
#tritrig_gbl_name = '/home/alic/HPS/projects/simp_analysis/mc/seedtracker/ana/hadd_tritrigv2-beamv6_seedtracker_ana_GBL.root'


#wab_kf_name = '/home/alic/HPS/projects/simp_analysis/mc/wab_beam/ana/hadd_wabv3-beamv6_ana_KF.root'
#wab_gbl_name = '/home/alic/HPS/projects/simp_analysis/mc/wab_beam/ana/hadd_wabv3-beamv6_ana_GBL.root'

#Old SeedTracker
#wab_gbl_name = '/home/alic/HPS/projects/simp_analysis/mc/seedtracker/ana/hadd_wabv3-beamv6_seedtracker_ana_GBL.root'


#data_kf_name = '/home/alic/HPS/projects/simp_analysis/data/sample0/hps_sample0_data_rerecon_ana_KF.root'
#data_gbl_name = '/home/alic/HPS/projects/simp_analysis/data/sample0/hps_sample0_data_rerecon_ana_GBL.root'

#Old SeedTracker
#data_gbl_name = '/home/alic/HPS/projects/simp_analysis/data/seedtracker/ana/hadd_data_seedtracker_ana_GBL.root'

total_lumi = 16.67 + 18.33 + 18.86 + 18.54 + 18.90 #nb-1
total_lumi = total_lumi * 0.001 #pb-1

tritrig_kf_name = '/sdf/group/hps/users/alspellm/projects/THESIS/simp_reach_estimates/simps_2016_kf/mc/tritrig_beam/ana/final_hadd_tritrigv2-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_976_KF_CR.root'
tritrig_events = 50000*9853 #nevents
tritrig_xsec = 1.416e9 #pb

wab_kf_name = '/sdf/group/hps/users/alspellm/projects/THESIS/simp_reach_estimates/simps_2016_kf/mc/wab_beam/ana/final_hadd_wabv3-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_KF_ana_CR.root'
wab_events = 100000*9966 #nevents
wab_xsec = 0.1985e12 #pb

data_kf_name = '/sdf/group/hps/users/alspellm/projects/THESIS/simp_reach_estimates/simps_2016_kf/data/hadd_sample0_KF_ana_CR.root'
data_lumi = 0.0913 #pb-1
lumi_2016 = 10.7 #pb-1

#Read histograms
keywords = ["_Psum_","_InvM_h","_ele_p_h","pos_p_h","ele_TanLambda_h","pos_TanLambda_h","ele_Phi_h","pos_Phi_h","vtx_Z_h","_cutflow"]
#keywords = [""]
#kfselection = 'vtxana_kf_vertexSelection_loose'

vtxselection = 'vtxSelection'
if runKF:

    #Load root files
    data_kf = r.TFile(data_kf_name,"READ")
    tritrig_kf = r.TFile(tritrig_kf_name,"READ")
    wab_kf = r.TFile(wab_kf_name,"READ")

    kfselection = 'vtxana_kf_%s'%(vtxselection)
    tritrig_kf_h = utils.read1DPlotsFromRootDir(tritrig_kf,kfselection,keywords)
    wab_kf_h = utils.read1DPlotsFromRootDir(wab_kf,kfselection,keywords)
    data_kf_h = utils.read1DPlotsFromRootDir(data_kf,kfselection,keywords)

    for i,pname in enumerate(data_kf_h):
        canv_name = 'KF'+pname.replace("vtxana_kf",'')

        #histos
        data = data_kf_h[pname]
        name = 'Data_KF'+data.GetName().replace("vtxana_kf",'')
        utils.format1DPlot(data,name,title=name,linewidth=1,linecolor=1,markerstyle=8,markersize=0.5)

        tritrig = tritrig_kf_h[pname]
        tritrig.Scale(tritrig_xsec*total_lumi/(tritrig_events))
        name = 'tritrig_KF'+tritrig.GetName().replace("vtxana_kf",'')
        utils.format1DPlot(tritrig,name,title=name,linecolor=colors[1])

        wab = wab_kf_h[pname]
        wab.Scale(wab_xsec*total_lumi/(wab_events))
        name = 'WAB_KF'+wab.GetName().replace("vtxana_kf",'')
        utils.format1DPlot(wab,name,title=name,linecolor=colors[2])

        tritrig_wab = tritrig.Clone('tritrig_clone')
        tritrig_wab.Add(wab,1)
        name = 'tritrig+WAB'+data.GetName().replace("vtxana_kf",'').replace("Data",'')
        utils.format1DPlot(tritrig_wab,name,title=name,linecolor=colors[5])

        plots = [data,tritrig,wab,tritrig_wab]
        #utils.overlay1DPlots(plots,canv_name,'/home/alic/HPS/projects/simp_analysis/plots/sample_0_comparisons/081622/vertexSelection/%s'%(kfselection),colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])
        utils.overlay1DPlots(plots,canv_name,'/sdf/group/hps/users/alspellm/projects/THESIS/simp_reach_estimates/simps_2016_kf/ana_scripts/plots/',colors=[1,colors[1],colors[2],colors[5]],drawStyles=["Ehist","hist","hist","hist"])


        #ratio plots
        ratioPairs = [(0,3)]
        ratioNames = ["Data:Tritrig+Wab+Beam"]
        #utils.makeRatioPlot(plots,ratioPairs,ratioNames,"ratios_"+canv_name, '/home/alic/HPS/projects/simp_analysis/plots/sample_0_comparisons/081622/vertexSelection/%s'%(kfselection),colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])
        utils.makeRatioPlot(plots,ratioPairs,ratioNames,"ratios_"+canv_name, '/sdf/group/hps/users/alspellm/projects/THESIS/simp_reach_estimates/simps_2016_kf/ana_scripts/plots',colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])

if runGBL:

    #Load root files
    data_gbl = r.TFile(data_gbl_name,"READ")
    tritrig_gbl = r.TFile(tritrig_gbl_name,"READ")
    wab_gbl = r.TFile(wab_gbl_name,"READ")

    #gblselection = 'vtxana_gbl_vertexSelection'
    gblselection = 'vtxana_gbl_%s'%(vtxselection)
    tritrig_gbl_h = utils.read1DPlotsFromRootDir(tritrig_gbl,gblselection,keywords)
    wab_gbl_h = utils.read1DPlotsFromRootDir(wab_gbl,gblselection,keywords)
    data_gbl_h = utils.read1DPlotsFromRootDir(data_gbl,gblselection,keywords)
    keywords = [""]

    for i,pname in enumerate(data_gbl_h):
        canv_name = 'gbl'+pname.replace("vtxana_gbl",'')

        #histos
        data = data_gbl_h[pname]
        #name = 'Data_gbl'+data.GetName().replace("vtxana_gbl",'')
        name = 'SeedTracker_Data_gbl'+data.GetName().replace("vtxana_gbl",'')
        utils.format1DPlot(data,name,title=name,linewidth=1,linecolor=1,markerstyle=8,markersize=0.5)

        tritrig = tritrig_gbl_h[pname]
        tritrig.Scale(tritrig_xsec*total_lumi/(tritrig_events))
        #name = 'tritrig_gbl'+tritrig.GetName().replace("vtxana_gbl",'')
        name = 'SeedTracker_tritrig_gbl'+tritrig.GetName().replace("vtxana_gbl",'')
        utils.format1DPlot(tritrig,name,title=name,linecolor=colors[1])

        wab = wab_gbl_h[pname]
        wab.Scale(wab_xsec*total_lumi/(wab_events))
        #name = 'WAB_gbl'+wab.GetName().replace("vtxana_gbl",'')
        name = 'SeedTracker_WAB_gbl'+wab.GetName().replace("vtxana_gbl",'')
        utils.format1DPlot(wab,name,title=name,linecolor=colors[2])

        tritrig_wab = tritrig.Clone('tritrig_clone')
        tritrig_wab.Add(wab,1)
        #name = 'tritrig+WAB'+data.GetName().replace("vtxana_gbl",'').replace("Data",'')
        name = 'SeedTracker_tritrig+WAB'+data.GetName().replace("vtxana_gbl",'').replace("Data",'')
        utils.format1DPlot(tritrig_wab,name,title=name,linecolor=colors[5])

        plots = [data,tritrig,wab,tritrig_wab]
        #utils.overlay1DPlots(plots,canv_name,'/home/alic/HPS/projects/simp_analysis/plots/sample_0_comparisons/081622/vertexSelection/%s'%(gblselection),colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])
        #utils.overlay1DPlots(plots,canv_name,'/home/alic/HPS/projects/simp_analysis/plots/sample_0_comparisons/081622/vertexSelection/',colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])
        utils.overlay1DPlots(plots,canv_name,'/home/alic/HPS/projects/simp_analysis/plots/sample_0_comparisons/081622/vertexSelection_SeedTracker/',colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])


        #ratio plots
        ratioPairs = [(0,3),(2,1)]
        ratioNames = ["Data:Tritrig+Wab+Beam","WAB+Beam:Tritrig+Beam"]
        #utils.makeRatioPlot(plots,ratioPairs,ratioNames,"ratios_"+canv_name, '/home/alic/HPS/projects/simp_analysis/plots/sample_0_comparisons/081622/vertexSelection/%s'%(gblselection),colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])
        #utils.makeRatioPlot(plots,ratioPairs,ratioNames,"ratios_"+canv_name, '/home/alic/HPS/projects/simp_analysis/plots/sample_0_comparisons/081622/vertexSelection/',colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])
        utils.makeRatioPlot(plots,ratioPairs,ratioNames,"ratios_"+canv_name, '/home/alic/HPS/projects/simp_analysis/plots/sample_0_comparisons/081622/vertexSelection_SeedTracker/',colors=[1,colors[1],colors[2],colors[5]],drawStyles=["EP","hist","hist","hist"])



#outfile.Close()

#hist_units_dic
#tritrig_beam Scale 1.416e9*Lumi/(50000*n)
#data_scale -> total lumi for 7800 = 167.17 (1/nb) for 388 files. I have 387 files.

#invMassHistos['rad'].Scale(units*66.36e6*Lumi/(rebin*10000*9967))
#invMassHistos['tritrig'].Scale(units*1.416e9*Lumi/(rebin*50000*9853))
#invMassHistos['wab'].Scale(units*0.1985e12*Lumi/(rebin*100000*9966))
