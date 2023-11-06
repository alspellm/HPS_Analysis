#!/usr/bin/python3
import ROOT as r
import numpy as np
from array import array
import sys
import plot_utilities as utils
import json
from tabulate import tabulate

class AnalysisTools:
    @staticmethod
    def getCutType(cut):
        if 'Time' in cut:
            cut_type = 'abs'
        elif '_lt' in cut:
            cut_type = 'lt'
        elif '_gt' in cut:
            cut_type = 'gt'
        else:
            return ''
        return cut_type
    @staticmethod
    def integrateHisto(histo, cutval=None, cut_type=''):
        total_integral = histo.Integral(0, histo.GetNbinsX()+1)
        firstbin = -99999.9
        lastbin = 99999.9

        if cut_type == 'abs':
            firstbin = histo.FindBin(-cutval)
            lastbin = histo.FindBin(cutval)
            if firstbin == 0:
                firstbin = 1
            if lastbin > histo.GetNbinsX():
                lastbin = histo.GetNbinsX()

        elif cut_type == 'lt':
            firstbin = 0
            lastbin = int(histo.FindBin(cutval))
            if lastbin > histo.GetNbinsX():
                lastbin = histo.GetNbinsX()

        elif cut_type == 'gt':
            firstbin = int(histo.FindBin(cutval))
            lastbin = int(histo.GetNbinsX()+1)
            if firstbin == 0:
                firstbin = 1 

        elif cut_type == '':
            firstbin = 0
            lastbin = int(histo.GetNbinsX()+1)

        integral = round(float(histo.Integral(firstbin, lastbin)),2)
        eff = round(integral/total_integral,3)
        print('integral',integral, 'entries', nentries, 'eff', eff)

        return integral, eff


class PreselectionNminus1:

    def __init__(self, data, tritrig, wab, sigA, sigB, preselection, cutvariables, latex_cutvariables, plots_dir, mcScale):
        self.data = data
        self.tritrig = tritrig
        self.wab = wab
        self.sigA = sigA
        self.sigB = sigB
        self.preselection = preselection
        self.cutvariables = cutvariables
        self.latex_cutvariables = latex_cutvariables
        self.plots_dir = plots_dir
        self.mcScale = mcScale

    def runAnalysis(self):
        tools = AnalysisTools()

        #Store values for latex table
        table_entries = {}
        row_labels = []

        f = open(self.preselection)
        sel = json.load(f)
        i = 0
        for cut, info in sel.items():

            #skip pair1
            if cut == 'Pair1_eq':
                continue

            cutvar = self.cutvariables[cut]
            cutval = float(info['cut'])
            print(cut)
            print('cutval: ', cutval)
            directory = 'vtxana_kf_%s_selection'%(cut)

            cut_type = tools.getCutType(cut)

            #Format row labels for latex
            latex_label = self.latex_cutvariables[cut]
            #Append cut value
            latex_label = '$'+latex_label+'%s$'%(str(cutval))
            row_labels.append(latex_label)

            #data
            data_plot = utils.read_1d_plots_from_root_file(self.data, directory, '%s_h'%(cutvar))[0]
            utils.formatHisto(data_plot, line_color=colors[0],name='Data', title='1% data')

            #tritrig
            tritrig_plot = utils.read_1d_plots_from_root_file(self.tritrig, directory, '%s_h'%(cutvar))[0]
            utils.formatHisto(tritrig_plot, line_color=colors[1],name='Tritrig-Beam', title='Tritrig+Beam')

            #wab
            wab_plot = utils.read_1d_plots_from_root_file(self.wab, directory, '%s_h'%(cutvar))[0]
            utils.formatHisto(wab_plot, line_color=colors[2],name='WAB-Beam', title='WAB+Beam')

            #Signal
            sig60_plot = utils.read_1d_plots_from_root_file(self.sigA, directory, '%s_h'%(cutvar))[0]
            utils.formatHisto(sig60_plot, line_color=colors[3],name='Signal-60 MeV', title='Signal_60MeV')
            sig100_plot = utils.read_1d_plots_from_root_file(self.sigB, directory, '%s_h'%(cutvar))[0]
            utils.formatHisto(sig100_plot, line_color=colors[4],name='Signal-100 MeV', title='Signal_100MeV')

            #Scale tritrig+wab+beam
            mc_bkg_plot = tritrig_plot.Clone()
            mc_bkg_plot.Scale(self.mcScale['tritrig'])
            wab_clone_plot = wab_plot.Clone()
            wab_clone_plot.Scale(self.mcScale['wab'])
            mc_bkg_plot.Add(wab_clone_plot)
            utils.formatHisto(mc_bkg_plot, line_color=colors[5],name='Tritrig-WAB-Beam', title='tritri+wab+beam')

            #collect plots
            plots = [data_plot, tritrig_plot, wab_plot, mc_bkg_plot, sig60_plot, sig100_plot]

            #Integrate histogram
            for plot in plots:
                name = plot.GetName()
                integral, eff = tools.integrateHisto(plot, cutval, cut_type)
                if name in table_entries:
                    table_entries[name][0].append(integral)
                    table_entries[name][1].append(eff)
                else:
                    table_entries[name] = [ [integral], [eff] ]

            #Scale Histograms
            for plot in plots:
                plot.Scale(1.0/float(plot.Integral(-1, plot.GetNbinsX()+1)))

            #Plots
            logY = True
            freezeX = True
            no_logy = ['N2Dhits','TrkCluTimeDiff']
            #if any(substring in cut for substring in no_logy):
            #    logY = False
            no_freezeX = ['TrkCluTimeDiff','eleposCluTimeDiff_lt']
            if any(substring in cut for substring in no_freezeX):
                freezeX = False

            c = utils.plot_TH1s_with_legend(plots, cut, freezeXaxis=freezeX, LogY=logY, save=False)
            c.Range( 0., -10., 1., 10. )
            line = r.TLine(cutval, c.GetUymin(), cutval, c.GetUymax())
            line.SetLineColor(2)  # Set the line color (optional)
            line.SetLineWidth(3)
            line.Draw("same")
            if cut_type == 'abs':
                line2 = r.TLine(-cutval, c.GetUymin(), -cutval, c.GetUymax())
                line2.SetLineColor(2)  # Set the line color (optional)
                line2.SetLineWidth(3)
                line2.Draw("same")

            c.SaveAs('%s/%s.png'%(self.plots_dir,cut))
            c.Close()

        #Make latex table
        table = []
        headers = []
        for key in table_entries.keys():
            headers.extend([f"{key} Integral", f"{key} Eff"])
        table.append(headers)

        for i, row_label in enumerate(row_labels):
            row = [row_label]
            for key, values in table_entries.items():
                for value in values:
                    row.append(value[i])
            table.append(row)

        # Generate the LaTeX table
        latex_table = tabulate(table, headers="firstrow", tablefmt="latex_raw")
        print(latex_table)
        print('\n Latex Table Format:\n')
        stdoutOrigin=sys.stdout
        sys.stdout = open('latex_table.txt','w')
        latex_table = tabulate(table, headers="firstrow", tablefmt="latex")
        print(latex_table)
        sys.stdout.close()
        sys.stdout=stdoutOrigin

class PreselectionCutflow:

    def __init__(self, data, tritrig, wab, sigA, sigB, preselection, cutvariables, latex_cutvariables, plots_dir, mcScale):
        self.data = data
        self.tritrig = tritrig
        self.wab = wab
        self.sigA = sigA
        self.sigB = sigB
        self.preselection = preselection
        self.cutvariables = cutvariables
        self.latex_cutvariables = latex_cutvariables
        self.plots_dir = plots_dir
        self.mcScale = mcScale

    def runAnalysis(self):
        tools = AnalysisTools()
        row_labels = []

        #Store cut plots for each plot type
        invMass = {'data':{}, 'tritrig':{}, 'wab':{}, 'mc_bkg':{}, 'signal_60':{}, 'signal_100':{}}
        reconZ = {'data':{}, 'tritrig':{}, 'wab':{}, 'mc_bkg':{}, 'signal_60':{}, 'signal_100':{}}
        psum = {'data':{}, 'tritrig':{}, 'wab':{}, 'mc_bkg':{}, 'signal_60':{}, 'signal_100':{}}

        f = open(self.preselection)
        load_sel = json.load(f)
        #insert preprocessing with no cuts
        sel = {
            'Preprocessing': {
                'cut': 1,
                'id': -1,
                'info': 'Preprocessing'
            }
        }
        for key, value in load_sel.items():
            sel[key] = value
        
        for i, (cut, info) in enumerate(sel.items()):
            latex_cut = self.latex_cutvariables[cut] 
            cutvar = self.cutvariables[cut]
            cutval = float(info['cut'])
            #Format row labels for latex
            latex_label = self.latex_cutvariables[cut]
            #Append cut value
            latex_label = '$'+latex_label+'%s$'%(str(cutval))
            row_labels.append(latex_label)

            print('cut:', cut)
            print('cutval: ', cutval)
            print('latex_label',latex_label)

            if cut == 'Preprocessing':
                directory = 'vtxana_kf_vtxSelection'
            else:
                directory = 'vtxana_kf_%s_selection_inclusive'%(cut)

        #Mass
            ##data
            data_mass = utils.read_1d_plots_from_root_file(self.data, directory, 'vtx_InvM_h')[0]
            utils.formatHisto(data_mass, line_color=colors[i],name='Data_InvM_%s'%(cut), title='\$%s %s\$'%(latex_cut, cutval))
            ##tritrig
            tritrig_mass = utils.read_1d_plots_from_root_file(self.tritrig, directory, 'vtx_InvM_h')[0]
            utils.formatHisto(tritrig_mass, line_color=colors[i],name='Tritrig_Beam_InvM_%s'%(cut), title='\$%s %s\$'%(latex_cut, cutval))
            ##wab
            wab_mass = utils.read_1d_plots_from_root_file(self.wab, directory, 'vtx_InvM_h')[0]
            utils.formatHisto(wab_mass, line_color=colors[i],name='WAB_Beam_InvM_%s'%(cut), title='\$%s %s\$'%(latex_cut, cutval))
            ##Signal
            sig60_mass = utils.read_1d_plots_from_root_file(self.sigA, directory, 'vtx_InvM_h')[0]
            utils.formatHisto(sig60_mass, line_color=colors[i],name='Signal_60_MeV_%s'%(cut), title='\$%s %s\$'%(latex_cut, cutval))
            sig100_mass = utils.read_1d_plots_from_root_file(self.sigB, directory, 'vtx_InvM_h')[0]
            utils.formatHisto(sig100_mass, line_color=colors[i],name='Signal_100_MeV_%s'%(cut), title='\$%s %s\$'%(latex_cut, cutval))
            ##Scale tritrig+wab+beam
            mc_bkg_mass = tritrig_mass.Clone()
            mc_bkg_mass.Scale(self.mcScale['tritrig'])
            wab_mass_clone = wab_mass.Clone()
            wab_mass_clone.Scale(self.mcScale['wab'])
            mc_bkg_mass.Add(wab_mass_clone)
            utils.formatHisto(mc_bkg_mass, line_color=colors[i],name='Tri_WAB_Beam_InvM_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))

            #collect plots
            invMass['data'][cut] = data_mass
            invMass['tritrig'][cut] = tritrig_mass
            invMass['wab'][cut] = wab_mass
            invMass['signal_60'][cut] = sig60_mass
            invMass['signal_100'][cut] = sig100_mass
            invMass['mc_bkg'][cut] = mc_bkg_mass

        #Recon Z
            ##data
            data_z = utils.read_1d_plots_from_root_file(self.data, directory, 'vtx_Z_h')[0]
            utils.formatHisto(data_z, line_color=colors[i],name='Data_Z_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))
            ##tritrig
            tritrig_z = utils.read_1d_plots_from_root_file(self.tritrig, directory, 'vtx_Z_h')[0]
            utils.formatHisto(tritrig_z, line_color=colors[i],name='Tritrig_Beam_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))
            ##wab
            wab_z = utils.read_1d_plots_from_root_file(self.wab, directory, 'vtx_Z_h')[0]
            utils.formatHisto(wab_z, line_color=colors[i],name='WAB_Beam_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))
            ##Signal
            sig60_z = utils.read_1d_plots_from_root_file(self.sigA, directory, 'vtx_Z_h')[0]
            utils.formatHisto(sig60_z, line_color=colors[i],name='Signal_60_MeV_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))
            sig100_z = utils.read_1d_plots_from_root_file(self.sigB, directory, 'vtx_Z_h')[0]
            utils.formatHisto(sig100_z, line_color=colors[i],name='Signal_100_MeV_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))
            ##Scale tritrig+wab+beam
            mc_bkg_z = tritrig_z.Clone()
            mc_bkg_z.Scale(self.mcScale['tritrig'])
            wab_z_clone = wab_z.Clone()
            wab_z_clone.Scale(self.mcScale['wab'])
            mc_bkg_z.Add(wab_z_clone)
            utils.formatHisto(mc_bkg_z, line_color=colors[i],name='Tri_WAB_Beam_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))

            #collect plots
            reconZ['data'][cut] = data_z
            reconZ['tritrig'][cut] = tritrig_z
            reconZ['wab'][cut] = wab_z
            reconZ['signal_60'][cut] = sig60_z
            reconZ['signal_100'][cut] = sig100_z
            reconZ['mc_bkg'][cut] = mc_bkg_z


        #Psum
            ##data
            data_Psum = utils.read_1d_plots_from_root_file(self.data, directory, 'vtx_Psum_h')[0]
            utils.formatHisto(data_Psum, line_color=colors[i],name='Data_Psum_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))
            ##tritrig
            tritrig_Psum = utils.read_1d_plots_from_root_file(self.tritrig, directory, 'vtx_Psum_h')[0]
            utils.formatHisto(tritrig_Psum, line_color=colors[i],name='Tritrig_Beam_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))
            ##wab
            wab_Psum = utils.read_1d_plots_from_root_file(self.wab, directory, 'vtx_Psum_h')[0]
            utils.formatHisto(wab_Psum, line_color=colors[i],name='WAB_Beam_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))
            ##Signal
            sig60_Psum = utils.read_1d_plots_from_root_file(self.sigA, directory, 'vtx_Psum_h')[0]
            utils.formatHisto(sig60_Psum, line_color=colors[i],name='Signal_60_MeV_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))
            sig100_Psum = utils.read_1d_plots_from_root_file(self.sigB, directory, 'vtx_Psum_h')[0]
            utils.formatHisto(sig100_Psum, line_color=colors[i],name='Signal_100_MeV_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))
            ##Scale tritrig+wab+beam
            mc_bkg_Psum = tritrig_Psum.Clone()
            mc_bkg_Psum.Scale(self.mcScale['tritrig'])
            wab_Psum_clone = wab_Psum.Clone()
            wab_Psum_clone.Scale(self.mcScale['wab'])
            mc_bkg_Psum.Add(wab_Psum_clone)
            utils.formatHisto(mc_bkg_Psum, line_color=colors[i],name='Tri_WAB_Beam_%s'%(cut), title='\$%s %s\$'%(latex_cut,cutval))

            #collect plots
            psum['data'][cut] = data_Psum
            psum['tritrig'][cut] = tritrig_Psum
            psum['wab'][cut] = wab_Psum
            psum['signal_60'][cut] = sig60_Psum
            psum['signal_100'][cut] = sig100_Psum
            psum['mc_bkg'][cut] = mc_bkg_Psum

        #Store values for latex table
        table_entries = {}
        
        for entry in invMass:
            plots = []
            cuts = invMass[entry]
            initial_integral = 0.0
            for cut, plot in cuts.items():
                plots.append(plot)

                #Integrate histogram
                name = plot.GetTitle()
                integral = plot.Integral(0,plot.GetNbinsX()+1)
                if cut == 'Preprocessing':
                    initial_integral = integral
                    name = 'Preprocessing'
                eff = integral/initial_integral
                print("Cut: ", cut)
                print("Eff: ", eff)

                if entry in table_entries:
                    table_entries[entry][0].append(integral)
                    table_entries[entry][1].append(eff)
                else:
                    table_entries[entry] = [ [integral], [eff] ]

            c = utils.plot_TH1s_with_legend(plots, entry, legx1=0.7, legx2=0.9, legy1=0.5, legy2=0.9, line_spacing=0.05, save=False)
            c.SaveAs('%s/%s.png'%(plots_dir,'%s_invM_presele_cutflow'%(entry)))

        #Make latex table
        table = []
        headers = []
        print(table_entries)
        for key in table_entries.keys():
            headers.extend([f"{key} Integral", f"{key} Eff"])
        table.append(headers)

        for i, row_label in enumerate(row_labels):
            row = [row_label]
            for key, values in table_entries.items():
                for value in values:
                    row.append(value[i])
            table.append(row)

        # Generate the LaTeX table
        latex_table = tabulate(table, headers="firstrow", tablefmt="latex_raw")
        print(latex_table)
        print('\n Latex Table Format:\n')
        stdoutOrigin=sys.stdout
        sys.stdout = open('preselection_cutflow.txt','w')
        latex_table = tabulate(table, headers="firstrow", tablefmt="latex")
        print(latex_table)
        sys.stdout.close()
        sys.stdout=stdoutOrigin

        for entry in reconZ:
            plots = []
            cuts = reconZ[entry]
            for cut, plot in cuts.items():
                plots.append(plot)

            c = utils.plot_TH1s_with_legend(plots, entry, legx1=0.7, legx2=0.9, legy1=0.5, legy2=0.9, line_spacing=0.05, LogY=True, save=False)
            c.SaveAs('%s/%s.png'%(plots_dir,'%s_reconZ_presele_cutflow'%(entry)))

        for entry in psum:
            plots = []
            cuts = psum[entry]
            for cut, plot in cuts.items():
                plots.append(plot)

            c = utils.plot_TH1s_with_legend(plots, entry, legx1=0.7, legx2=0.9, legy1=0.5, legy2=0.9, line_spacing=0.05, save=False)
            c.SaveAs('%s/%s.png'%(plots_dir,'%s_psum_presele_cutflow'%(entry)))


if __name__ == "__main__":

    r.gROOT.SetBatch(1)
    style = utils.SetMyStyle()
    colors = utils.getColorsHPS()

    #selection
    preselection = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/configs/vertexSelection_2016_simp_preselection.json'

    #Define cut variables for each cut in json
    cuts = {'Preprocessing': '',
            'Pair1_eq': '',
            'eleTrkTime_lt': 'ele_time',
            'posTrkTime_lt': 'pos_time',
            'eleposCluTimeDiff_lt': 'ele_pos_clusTimeDiff',
            'eleTrkCluTimeDiff_lt': 'ele_track_clus_dt',
            'posTrkCluTimeDiff_lt': 'pos_track_clus_dt',
            'eleTrkChi2Ndf_lt': 'ele_chi2ndf',
            'posTrkChi2Ndf_lt': 'pos_chi2ndf',
            'eleMom_lt': 'ele_p',
            'eleMom_gt': 'ele_p',
            'posMom_gt': 'pos_p',
            'eleN2Dhits_gt': 'ele_track_n2dhits',
            'posN2Dhits_gt': 'pos_track_n2dhits',
            'maxVtxMom_lt': 'vtx_Psum',
            'chi2unc_lt':'vtx_chi2'
            }

    latex_cuts = {'Preprocessing': '',
            'Pair1_eq': 'Pairs1 Trigger ',
            'eleTrkTime_lt': '|e^{-} Track_{t}| < ',
            'posTrkTime_lt': '|e^{+} Track_{t}| <',
            'eleposCluTimeDiff_lt': '\Delta_{t}(cluster_{e^{-}},cluster_{e^{+}} < ',
            'eleTrkCluTimeDiff_lt': 'e^{-}\Delta_{t}(track,cluster) < ',
            'posTrkCluTimeDiff_lt': 'e^{+}\Delta_{t}(track,cluster) < ',
            'eleTrkChi2Ndf_lt': 'e^{-} Track \chi^2/n.d.f. < ',
            'posTrkChi2Ndf_lt': 'e^{+} Track \chi^2/n.d.f. < ',
            'eleMom_lt': 'p_{e^-} < ',
            'eleMom_gt': 'p_{e^-} > ',
            'posMom_gt': 'p_{e^+} > ',
            'eleN2Dhits_gt': 'N_{2d hits} e^{-}_{Track} > ',
            'posN2Dhits_gt': 'N_{2d hits} e^{+}_{Track} > ',
            'maxVtxMom_lt': 'p_{e^{-}+e^{+}} < ',
            'chi2unc_lt':'Vtx_{\chi^2} < '
            }

    #data
    data = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/data/hadd_hps_BLPass4_10pct_preselection_nminus1_20230914.root'
    #tritrig
    tritrig = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/tritrig_beam/hadd_tritrig_beam_preselection_20230915.root'
    #wab
    wab = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/wab_beam/hadd_wab_beam_preselection_20230915.root'

    #signal 60MeV
    sig_60 = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/signal/hadd_signal_60_MeV_beam_preselection_20230915.root'
    sig_100 = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/signal/hadd_signal_60_MeV_beam_preselection_20230915.root'

    plots_dir = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/plots'

    #MC Scaling
    Lumi = 0.107 #1/pb
    mcScale = {}
    mcScale['tritrig'] = 1.416e9*Lumi/(50000*9456) #pb2016
    mcScale['wab'] = 0.1985e12*Lumi/(100000*9688) #pb2016

    #preselectionAna = PreselectionNminus1(data, tritrig, wab, sig_60, sig_100, preselection, cuts, latex_cuts, plots_dir, mcScale)
    #preselectionAna.runAnalysis()


#Cutflow
#data
    data = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/data/hadd_hps_BLPass4b_preselection_cutflow_ana_20230922.root'
    #tritrig
    tritrig = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/tritrig_beam/hadd_tritrig-beam_preselection_cutflow_ana_20230920.root'
    #wab
    wab = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/wab_beam/hadd_wab_beam_preselection_cutflow_ana_20230920.root'

    #signal 60MeV
    sig_60 = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/signal/hadd_simp_60_MeV_preselection_cutflow.root'
    sig_100 = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/signal/hadd_simp_100_MeV_preselection_cutflow.root'

    plots_dir = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/plots/cutflow'
    cutflowAna = PreselectionCutflow(data, tritrig, wab, sig_60, sig_100, preselection, cuts, latex_cuts, plots_dir, mcScale)
    cutflowAna.runAnalysis()







