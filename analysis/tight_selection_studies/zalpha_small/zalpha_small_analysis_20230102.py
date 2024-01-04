import os
import uproot
import ROOT as r
import math
import awkward as ak
import hist
from hist import Hist
import sys
import argparse
import hist
from hist import loc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
from SimpEquations_validated import SimpEquations as simpeqs
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import mpl_plot_utilities as mplutils
#get_ipython().run_line_magic('load_ext', 'autoreload')
#get_ipython().run_line_magic('autoreload', '2')

def calculateZBi(n_on, n_off, tau):
    P_Bi = r.TMath.BetaIncomplete(1./(1.+tau),n_on,n_off+1) #why plus 1?
    Z_Bi = 2.**(0.5)*r.TMath.ErfInverse(1-2*P_Bi)
    return Z_Bi
    
def radiativeFraction(mass_mev):
    radF = -1.04206e-01 + 9.92547e-03*mass_mev + -1.99437e-04*pow(mass_mev,2) + 1.83534e-06*pow(mass_mev,3) + -7.93138e-9*pow(mass_mev,4) + 1.30456e-11*pow(mass_mev,5) #alic 2016 simps kf 11/15/22
    return radF

def radiativeAcceptance(mass_mev):
    acc = ( -7.35934e-01 + 9.75402e-02*mass_mev + -5.22599e-03*pow(mass_mev,2) + 1.47226e-04*pow(mass_mev,3) + -2.41435e-06*pow(mass_mev,4) + 2.45015e-08*pow(mass_mev,5) + -1.56938e-10*pow(mass_mev,6) + 6.19494e-13*pow(mass_mev,7) + -1.37780e-15*pow(mass_mev,8) + 1.32155e-18*pow(mass_mev,9) ) #alic 2016 simps kf 11/15/22 
    return acc

#Calculated in 'makeMassRes.py'
def massRes(mass):
    res = 1.06314 + 3.45955e-02*mass + -6.62113e-05*pow(mass,2) # 2016 simps kf 11/15/22
    return res

#Calculate efficiency vertex
def calculateTotalApProduction(mass_ap, eps, radFrac, radAcc, dNdm_CR):
    apProduction = (3.*(137./2.)*3.14159)* (mass_ap * eps*eps * radFrac * dNdm_CR)/radAcc
    return apProduction

def calculateExpectedSignal(mass, mass_ap, mass_pid, mass_lepton, fpid, alpha_dark, eps, selEffZ_h, target_pos, apProduction, signal_meanEnergyGeV, rho=True):
    ctau = simpeqs.getCtau(mass_ap, mass_pid, mass, eps, alpha_dark, fpid, mass_lepton, rho)
    gcTau = ctau*simpeqs.gamma(mass/1000.0, signal_meanEnergyGeV)
    effVtx = 0.0
    for z in range(len(selEffZ_h.axes[0].centers)):
        zz = selEffZ_h.axes[0].centers[z]
        if zz < target_pos:
            continue
        effVtx += (math.exp((target_pos - zz)/gcTau)/gcTau)*selEffZ_h.values()[z]
    br_VPi = 0.0
    if rho:
        br_VPi = simpeqs.br_Vrho_pi(mass_ap, mass_pid, mass, alpha_dark, fpid)
    else:
        br_VPi = simpeqs.br_Vphi_pi(mass_ap, mass_pid, mass, alpha_dark, fpid)
    br_V_ee = 1.0
    expSignal = apProduction * effVtx * br_VPi * br_V_ee
    return expSignal

def getTBranchArrays(events, branches=[],mass_low=0.0, mass_high=999.9):
    arrays = events.arrays(branches,f"(unc_vtx_mass*1000.0 > {mass_low}) & (unc_vtx_mass*1000.0 < {mass_high})")
    return arrays

def applyCorrection(array, variable, correction):
    array[variable] = array[variable] + correction
    
def v0ProjSigCut(array, cut_value):
    condition = array["unc_vtx_proj_sig"] < cut_value
    array = array[condition]
    return array
    
def deltaZCut(array, par0, par1, mass_mev):
    cut = par0 + par1*mass_mev
    condition = array["unc_vtx_deltaZ"] < cut
    array = array[condition]
    return array

def zalphaCut(array, cut_value):
    condition = array["unc_vtx_zalpha_max"] < cut_value
    array = array[condition]
    return array

def flatZ0Cut(array, cut_value):
    condition = (array["unc_vtx_min_z0"])> cut_value
    array = array[condition]
    return array

def zalphaTransformation(slope, recon_z, z0):
    condition = z0 > 0.0
    zalpha = ak.where(condition, recon_z - (z0/slope), recon_z + (z0/slope))
    return zalpha

def defineHistos1d(category, color='blue'):
    histos = {}
    histos["unc_vtx_proj_sig"] = mplutils.defHist1d(f"unc_vtx_proj_sig_{category}",100,0,10,xlabel="V0 Projection Significance", label=category, color=color)
    histos["unc_vtx_deltaZ"] = mplutils.defHist1d(f"unc_vtx_deltaZ_{category}",50,0,10,xlabel="deltaZ", label=category, color=color)
    histos["unc_vtx_z"] = mplutils.defHist1d(f"unc_vtx_z_{category}",120,-20,100,xlabel="recon z [mm]", logY=True, label=category, color=color)
    histos["unc_vtx_min_z0"] = mplutils.defHist1d(f"unc_vtx_min_z0_{category}",1000,0,5,xlabel="min z0 [mm]", logY=False, label=category, color=color)
    histos["unc_vtx_zalpha_max"] = mplutils.defHist1d(f"unc_vtx_zalpha_max_{category}",1000,-10000,100,xlabel="zalpha max", logY=False, label=category, color=color)
    histos["vd_true_vtx_z"] = mplutils.defHist1d(f"vd_true_vtx_z_{category}",200,-50.3,149.7,xlabel="truth vtx z [mm]", logY=True, label=category, color=color)
    histos["vd_true_vtx_energy"] = mplutils.defHist1d(f"vd_true_vtx_energy_{category}",250,0.0,2.5,xlabel="truth energy [GeV]", logY=False, label=category, color=color)
    return histos

def defineHistos2d(category):
    histos = {}
    histos["recon_z_vs_track_z0"] = mplutils.defHist2d(f'recon_z_vs_track_z0_{category}', 160,-20,60,500,-4,4, xlabel='recon z [mm]',ylabel='track z0 [mm]')
    return histos

def fillVarHistos1d(histos_dict, array, reset=False):
    for key, histo in histos_dict.items():
        if reset == True:
            histo.reset()
        if key in array.fields:
            histo.fill(array[key])
            
def fillHistos2d(histos_dict, array, reset=False):
    for key, histo in histos_dict.items():
        if reset == True:
            histo.reset()
    #track z0
    histos_dict["recon_z_vs_track_z0"].fill(array["unc_vtx_z"],array["unc_vtx_ele_track_z0"])
    histos_dict["recon_z_vs_track_z0"].fill(array["unc_vtx_z"],array["unc_vtx_pos_track_z0"])
    
def saveHistsToROOT(outfile, histos_dict, subdir=''):
    for key, histo in histos_dict.items():
        mplutils.writeHistToROOT(outfile, histo, subdir=subdir)


def findCutValueLT(initial_histo, cut_fraction):
    initial_integral = initial_histo[::sum].value
    test_int = 0.0
    start_bin = initial_histo.shape[0]-1
    while test_int < (cut_fraction)*initial_integral:
        test_int = initial_histo[start_bin:hist.overflow:sum].value
        if test_int >= (cut_fraction)*initial_integral:
            break
        else:
            start_bin = start_bin - 1
        if start_bin == 0:
            break
    return initial_histo.axes[0].centers[start_bin]


def findCutValueGT(initial_histo, cut_fraction):
    initial_integral = initial_histo[::sum].value
    test_int = 0.0
    start_bin = 1
    while test_int < (cut_fraction)*initial_integral:
        test_int = initial_histo[:start_bin:sum].value
        if test_int >= (cut_fraction)*initial_integral:
            break
        else:
            start_bin = start_bin + 1
        if start_bin >= initial_histo.shape[0]-1:
            start_bin = initial_histo.shape[0]-1
            break
    return initial_histo.axes[0].centers[start_bin]

def getMaxAbsZ0(ele_z0, pos_z0):
    condition = abs(ele_z0) > abs(pos_z0)
    max_z0 = ak.where(condition, abs(ele_z0), abs(pos_z0))
    return max_z0
    
def getMinAbsZ0(ele_z0, pos_z0):
    condition = abs(ele_z0) < abs(pos_z0)
    min_z0 = ak.where(condition, abs(ele_z0), abs(pos_z0))
    return min_z0

################################# RUN ######################################

parser = argparse.ArgumentParser(description='Process some inputs.')
    
parser.add_argument('--mass', type=int, default=55, help='Specify the mass as an integer (default: 55)')
parser.add_argument('--logeps2', type=float, default=-5.5, help='Specify logeps2 as a float (default: -5.5)')

args = parser.parse_args()

print(f"Mass: {args.mass}")
print(f"logeps2: {args.logeps2}")

#Specify mixing
logeps2 = args.logeps2
eps2 = pow(10,logeps2)
eps = np.sqrt(eps2)

#SIMP Parameters
alpha_dark = 0.01
mass_ratio_ap_to_vd = 1.66
mass_ratio_ap_to_pid = 3.0
ratio_mpi_to_fpi = 4.0*3.14159
mass_lepton = 0.511

#Configure output files
outdir = 'run_20230103'
outfilename = 'zalpha_slope_study'
if not os.path.exists(outdir):
    os.mkdir(outdir)
    
## **Configure High-Z Cuts** ##
deltaZ_par0 = 18.5972
deltaZ_par1 = 0.159555
v0projSig_cutvalue = 2.0

#Loop over mass range
for mass in range(args.mass, args.mass+1):
    plots_dir = f'{outdir}/{int(mass)}_plots'
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)
    mass = float(mass)
    mass_nsigma = 3.5
    mass_low = float(mass) - mass_nsigma*massRes(float(mass))/2.0
    mass_high = float(mass) + mass_nsigma*massRes(float(mass))/2.0
    print(f'MASS WINDOW: {mass_low} - {mass_high} MeV')
    mass_ap = mass*mass_ratio_ap_to_vd
    mass_low_ap = mass_ratio_ap_to_vd*mass_low
    mass_high_ap = mass_ratio_ap_to_vd*mass_high
    mass_pid = mass_ap/mass_ratio_ap_to_pid
    fpid = mass_pid/ratio_mpi_to_fpi

    #RadFrac and RadAcc eval at Ap mass
    radFrac = radiativeFraction(mass_ap) 
    radAcc = radiativeAcceptance(mass_ap) 
    
    # Instantiating pdf document
    #PDF = PdfPages(f'{outdir}/{pdf_name}.pdf')
    
    #init outfile
    outfile = uproot.recreate(f'./{outdir}/mass_{mass}_{outfilename}.root')
    
    #background/data
    print("Reading Data")
    bkg_branches = ["unc_vtx_mass","unc_vtx_z","unc_vtx_proj_sig","unc_vtx_deltaZ","unc_vtx_ele_track_z0","unc_vtx_pos_track_z0","unc_vtx_ele_track_tanLambda"]
    bkg_filename = '/sdf/group/hps/users/alspellm/projects/THESIS/data/2016/BLPass4c_20231006/ana_20231019/full_hadd_blpass4c_ana.root'
    bkg_subdir = 'vtxana_Tight_2016_simp_reach_SR'
    bkg_treename = 'vtxana_Tight_2016_simp_reach_SR_tree'
    bkg_subdirCR = 'vtxana_Tight_2016_simp_reach_CR'
    bkg_treenameCR = 'vtxana_Tight_2016_simp_reach_CR_tree'
    bkg_branchesCR = ["unc_vtx_mass"]
    bkg_sf_SR = 10.0
    bkg_sf_CR = 10.0
    print("Finished Reading Data")

    #Signal
    print("Reading Signal MC")
    signal_branches = ["unc_vtx_mass","unc_vtx_z","unc_vtx_proj_sig","unc_vtx_deltaZ","unc_vtx_ele_track_z0","unc_vtx_pos_track_z0","unc_vtx_ele_track_tanLambda","unc_vtx_pos_track_tanLambda","vd_true_vtx_z","vd_true_vtx_energy"]
    signal_subdir = 'vtxana_radMatchTight_2016_simp_reach_SR'
    signal_treename = 'vtxana_radMatchTight_2016_simp_reach_SR_tree'
    signal_filename = '/sdf/group/hps/users/alspellm/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/hps-java_v5pt2pt1/pass4/recon_20231009/ana_20231020/hadd_simp_signal_%s_MeV_beam_ana.root'%(int(mass))
    signal_sf = 1.0
    print("Finished Reading Signal MC")
    
    #Truth Signal
    print("Reading Signal Truth MC")
    signal_truth_filename = '/sdf/group/hps/users/alspellm/projects/THESIS/mc/2016/simps/slic/20230713_slic/20230724_slic_ana/ana_files/hadd_simp_%s_MeV_rot_slic_mcana.root'%(int(mass))
    truth_mc_selection = 'vtxana_mc_radMatchTight_2016_simp_reach_SR'
    truth_mc_tree='vtxana_mc_radMatchTight_2016_simp_reach_SR_tree'
    pdgid = 625
    print("Finished Reading Signal Truth MC")

    #Signal truth recon z
    sbudir = "vtxana_radMatchTight_2016_simp_reach_SR"
    signal_truth_file = uproot.open(signal_truth_filename)
    truth_z_h = signal_truth_file["mcAna/mcAna_mc625Z_h"].to_hist()
    signal_truth_file.close
    signal_truthE_h = mplutils.defHist1d('signal_truthE',250,0.0,2.5,title='truth energy',label='truth energy',xlabel='truth energy [GeV]',ylabel='MC Events',logY=False, color='blue')
    simZ_h = mplutils.defHist1d('signal_simZ',200,-50.3,149.7,title='signal_simZ',label='signal simZ',xlabel='true zvtx [mm]',ylabel='MC Events',logY=False, color='blue')
    simZ_h.fill(truth_z_h.axes[0].centers, weight=truth_z_h.counts())

    ## **LOAD DATA** ##
    #Read background for mass window
    mass_selection_SR = f"(unc_vtx_mass*1000.0 > {mass_low}) & (unc_vtx_mass*1000.0 < {mass_high})"
    bkg_arrays = mplutils.readTBranchAwk(bkg_filename, bkg_subdir, bkg_treename, bkg_branches, mass_selection_SR)

    #Get background rate in control region, eval at Ap mass, NOT Vd mass
    mass_selection_CR = f"(unc_vtx_mass*1000.0 > {mass_low_ap}) & (unc_vtx_mass*1000.0 < {mass_high_ap})"
    bkg_arraysCR = mplutils.readTBranchAwk(bkg_filename, bkg_subdirCR, bkg_treenameCR, bkg_branchesCR, mass_selection_CR)

    dNdm = bkg_sf_CR*len(bkg_arraysCR["unc_vtx_mass"])/(mass_high_ap-mass_low_ap)
    print("Background Rate dNdm:",dNdm)

    ## **LOAD SIGNAL** ##
    #Read signal for mass window
    signal_arrays = mplutils.readTBranchAwk(signal_filename, signal_subdir, signal_treename, signal_branches, mass_selection_SR)

    #Apply z0 corrections
    applyCorrection(signal_arrays,'unc_vtx_ele_track_z0', -0.058)
    applyCorrection(signal_arrays,'unc_vtx_pos_track_z0', -0.098)


    ## **INITIALIZE HISTOGRAMS** ##
    print("Initializing Histograms")
    #Define Initial background and signal histograms
    bkg_histos1d = defineHistos1d('bkg') 
    signal_histos1d = defineHistos1d('signal', color='red')
    
    #Fill initial histograms
    fillVarHistos1d(bkg_histos1d, bkg_arrays)
    fillVarHistos1d(signal_histos1d, signal_arrays)
    
    #Fill 2d histograms
    signal_histos2d = defineHistos2d('signal')
    bkg_histos2d = defineHistos2d('bkg')
    fillHistos2d(bkg_histos2d, bkg_arrays)
    fillHistos2d(signal_histos2d, signal_arrays)

    #Write and plot 
    saveHistsToROOT(outfile, bkg_histos1d, subdir='initial')
    saveHistsToROOT(outfile, bkg_histos2d, subdir='initial')
    saveHistsToROOT(outfile, signal_histos1d, subdir='initial')
    saveHistsToROOT(outfile, signal_histos2d, subdir='initial')
    
    vmin = 1
    vmax = 12000
    for key, histo in signal_histos2d.items():
        fig, ax = mplutils.plotHist2d(histo, text='initial', textpos=[0.7,0.8], text_ax=True, logZ=True)
        plt.savefig(f'{plots_dir}/{histo.metadata.get("title")}_initial.png')
        
    for key, histo in bkg_histos2d.items():
        fig = mplutils.plotHist2d(histo, text='initial', textpos=[0.7,0.8], text_ax=True, vmin=vmin, vmax=vmax, logZ=True)
        plt.savefig(f'{plots_dir}/{histo.metadata.get("title")}_initial.png')
        
    plt.close('all')
    
    #Apply v0proj and deltaZ cuts
    print("Applying high-z cuts")
    bkg_arrays = v0ProjSigCut(bkg_arrays, v0projSig_cutvalue)
    bkg_arrays = deltaZCut(bkg_arrays, deltaZ_par0, deltaZ_par1, mass)
    signal_arrays = v0ProjSigCut(signal_arrays, v0projSig_cutvalue)
    signal_arrays = deltaZCut(signal_arrays, deltaZ_par0, deltaZ_par1, mass) 
    
    n_slopes = 21
    figures = {}
    #for n,zalpha_slope in enumerate(range(1,n_slopes+1)):
    #for n, zalpha_slope in enumerate([0.0001,0.0002, 0.0004, 0.0006,0.0008,0.001,0.0015,0.002,0.004,0.020,0.03,0.04]):
    #for n, zalpha_slope in enumerate([0.0001, 0.005, 0.01]):
    for n in range(0,n_slopes):
        if n == 0:
            zalpha_slope = 0.0001
        else:
            zalpha_slope = n/500.0
        plt.close('all')
        zalpha_subdir = f'slope_{zalpha_slope}'
        zalpha_plots_dir = f'{plots_dir}/{zalpha_subdir}'
        if not os.path.exists(zalpha_plots_dir):
            os.mkdir(zalpha_plots_dir)
        print("Running Zalpha Slope:",zalpha_slope)

        #Copy arrays
        iter_signal_arrays = signal_arrays
        iter_bkg_arrays = bkg_arrays

        #Add signal zalpha max
        max_zalpha_arr_sig = np.maximum(zalphaTransformation(zalpha_slope,iter_signal_arrays["unc_vtx_z"],iter_signal_arrays["unc_vtx_ele_track_z0"]),
                                    zalphaTransformation(zalpha_slope,iter_signal_arrays["unc_vtx_z"],iter_signal_arrays["unc_vtx_pos_track_z0"]))
        iter_signal_arrays["unc_vtx_zalpha_max"] = max_zalpha_arr_sig

        #Add bkg zalpha max
        max_zalpha_arr_bkg = np.maximum(zalphaTransformation(zalpha_slope,iter_bkg_arrays["unc_vtx_z"],iter_bkg_arrays["unc_vtx_ele_track_z0"]),
                                    zalphaTransformation(zalpha_slope,iter_bkg_arrays["unc_vtx_z"],iter_bkg_arrays["unc_vtx_pos_track_z0"]))
        iter_bkg_arrays["unc_vtx_zalpha_max"] = max_zalpha_arr_bkg

        minzalpha = ak.min(ak.concatenate([max_zalpha_arr_sig, max_zalpha_arr_bkg]))-10.0
        #maxzalpha = ak.max(ak.concatenate([max_zalpha_arr_sig, max_zalpha_arr_bkg]))+10.0
        maxzalpha = 100
        #nbins = int((maxzalpha-minzalpha)/10.0)
        nbins = 1000

        signal_histos1d["unc_vtx_zalpha_max"] = mplutils.defHist1d(f"unc_vtx_zalpha_max_signal",nbins,minzalpha,maxzalpha,xlabel="zalpha max", logY=False, label='signal')
        bkg_histos1d["unc_vtx_zalpha_max"] = mplutils.defHist1d(f"unc_vtx_zalpha_max_bkg",nbins,minzalpha,maxzalpha,xlabel="zalpha max", logY=False, label='bkg')

        #Fill histograms before iterating
        fillVarHistos1d(bkg_histos1d, bkg_arrays,reset=True)
        fillVarHistos1d(signal_histos1d, signal_arrays, reset=True)

        #Fill 2d histograms
        fillHistos2d(bkg_histos2d, bkg_arrays, reset=True)
        fillHistos2d(signal_histos2d, signal_arrays, reset=True)

        #Write and plot 
        mod = f'{zalpha_subdir}/pre_iter'
        saveHistsToROOT(outfile, bkg_histos1d, mod)
        saveHistsToROOT(outfile, bkg_histos2d, mod)
        saveHistsToROOT(outfile, signal_histos1d, mod)
        saveHistsToROOT(outfile, signal_histos2d, mod)

        for key, histo in signal_histos2d.items():
            fig, ax = mplutils.plotHist2d(histo, text=f'iter-1 \n slope {zalpha_slope}', textpos=[0.7,0.8], text_ax=True, vmin=1, vmax=20000, logZ=True)
            #fig = plotHisto2D(histo, text='iter-1', text_x=0.8, text_y=0.8, text_ax=True)
            plt.savefig(f'{zalpha_plots_dir}/{histo.metadata.get("title")}_iter-1.png')

        for key, histo in bkg_histos2d.items():
            fig, ax = mplutils.plotHist2d(histo, text=f'iter-1 \n slope {zalpha_slope}', textpos=[0.7,0.8], text_ax=True, vmin=1, vmax=20000, logZ=True)
            #fig = plotHisto2D(histo, text='iter-1', text_x=0.8, text_y=0.8, text_ax=True)
            plt.savefig(f'{zalpha_plots_dir}/{histo.metadata.get("title")}_iter-1.png')

        plt.close('all')

        #Graphs to track performance
        nsigs = []
        nbkgs = []
        zbis = []
        cut_values = []
        effs = []

        max_iter = 95
        step_size = 0.01
        initial_histo = outfile[f'{zalpha_subdir}/pre_iter/unc_vtx_zalpha_max_signal'].to_hist()
        for iteration in range(0,max_iter):
            if iteration%10 == 0:
                print("Iteration: ", iteration)
            #print('Iteration',iteration)
            iter_subdir=f'{zalpha_subdir}/iter_{iteration}'
            cut_fraction = iteration*step_size
            cut_value = findCutValueLT(initial_histo, cut_fraction)
            cut_values.append(cut_value)

            #Apply cut to signal and background
            if iteration > 0:
                iter_signal_arrays = zalphaCut(iter_signal_arrays,cut_value)
                iter_bkg_arrays = zalphaCut(iter_bkg_arrays,cut_value)

            #Signal efficiency starting from no cut
            eff = len(iter_signal_arrays['unc_vtx_zalpha_max'])/initial_histo[::sum].value
            effs.append(eff)

            #Fill histograms for this iteration
            #mod = f'iter_{iteration}_slope_{zalpha_slope}'
            mod = f'slope_{zalpha_slope}/iter_{iteration}'
            fillVarHistos1d(bkg_histos1d, iter_bkg_arrays,reset=True)
            fillVarHistos1d(signal_histos1d, iter_signal_arrays, reset=True)
            #Fill 2d histograms
            fillHistos2d(bkg_histos2d, iter_bkg_arrays, reset=True)
            fillHistos2d(signal_histos2d, iter_signal_arrays, reset=True)

            saveHistsToROOT(outfile, bkg_histos1d, mod)
            saveHistsToROOT(outfile, bkg_histos2d, mod)
            saveHistsToROOT(outfile, signal_histos1d, mod)
            saveHistsToROOT(outfile, signal_histos2d, mod)

            for key, histo in signal_histos2d.items():
                fig, ax = mplutils.plotHist2d(histo, text=f'iteration_{iteration} \n slope {zalpha_slope}', textpos=[0.7,0.8], text_ax=True, vmin=1, vmax = 20000,  logZ=True)
                #fig = plotHisto2D(histo, text=f'iteration_{iteration}', text_x=0.8, text_y=0.8, text_ax=True)
                plt.savefig(f'{zalpha_plots_dir}/{histo.metadata.get("title")}_iter{iteration}.png')

            for key, histo in bkg_histos2d.items():
                fig, ax = mplutils.plotHist2d(histo, text=f'iteration_{iteration} \n slope {zalpha_slope}', textpos=[0.7,0.8], text_ax=True, vmin=1, vmax=20000, logZ=True)
                #fig = plotHisto2D(histo, text=f'iteration_{iteration}', text_x=0.8, text_y=0.8, text_ax=True)
                plt.savefig(f'{zalpha_plots_dir}/{histo.metadata.get("title")}_iter{iteration}.png')

            plt.close('all')

            #Count remaining background beyond zposition
            target_pos = -4.3 #mm
            nbkg = bkg_sf_SR*bkg_histos1d["unc_vtx_z"][loc(target_pos):hist.overflow:sum] + 0.5 #add half background event to keep nbkg > 0 for calculation
            nbkgs.append(nbkg)
            #print("Nbkg:",nbkg)

            #Get signal selection efficiency F(z)
            weights=signal_histos1d["vd_true_vtx_z"].counts()/simZ_h.counts()
            nan_mask = np.isnan(weights)
            weights = ak.where(nan_mask, 0, weights)
            selEffZ_h = mplutils.defHist1d(f"sel_eff_signal",200,-50.3,149.7,xlabel="truth vtx z [mm]", logY=True, label='signal', color='green')
            selEffZ_h.fill(signal_histos1d["vd_true_vtx_z"].axes[0].centers, 
                           weight=weights)
            mplutils.writeHistToROOT(outfile,selEffZ_h, subdir=iter_subdir)

            #Calculate Expected Signal
            signal_meanEnergyGeV = 1.4

            apProduction = calculateTotalApProduction(mass_ap, eps, radFrac, radAcc, dNdm)
            nsig_rho = calculateExpectedSignal(mass, mass_ap, mass_pid, mass_lepton, fpid, alpha_dark, eps, selEffZ_h, target_pos, apProduction, signal_meanEnergyGeV, True)
            nsig_phi = calculateExpectedSignal(mass, mass_ap, mass_pid, mass_lepton, fpid, alpha_dark, eps, selEffZ_h, target_pos, apProduction, signal_meanEnergyGeV, False)
            nsig_total = signal_sf*(nsig_rho+nsig_phi)
            nsigs.append(nsig_total)
            zbi = calculateZBi(nsig_total+nbkg, nbkg, 1.0)
            zbis.append(zbi)

            plt.close('all')

        xvalues = [i+1 for i in range(len(zbis))]
        fig, ax = plt.subplots(figsize=(15,10))
        #fig, ax = plt.subplots()
        ax.plot(xvalues, nsigs, label='Nsig',color='blue')
        ax.plot(xvalues, nbkgs, label='Nbkg',color='red')
        ax.set_yscale('log')
        ax.set_ylabel('Nsig and Nbkg', fontsize = 20)
        ax.set_xlabel("Iteration Number", fontsize = 20)
        ax.tick_params(axis='y', labelsize = 20)
        ax.tick_params(axis='x', labelsize = 20)
        ax.set_ylim(bottom=0.0, top=1000)
        ax.yaxis.grid(True, linestyle='--', which='major', color='gray', alpha=0.7)

        ax2 = ax.twinx()
        ax2.plot(xvalues, zbis, label='ZBi',color='green')
        ax2.set_ylabel('ZBi', color='green', fontsize = 20)
        ax2.tick_params(axis='y',labelcolor='green', labelsize = 20)
        ax2.set_ylim(bottom=0.1, top=8)
        ax2.text(1,6,f'zalpha_slope_{zalpha_slope}')

        fig.legend(loc='center')
        figures[f'zalpha_slope_{zalpha_slope}'] = fig
        fig.savefig(f'{zalpha_plots_dir}/mass_{mass}_zalpha_slope_{zalpha_slope}_results.png')
        plt.close('all')
        
        #Save data as TGraphs
        graphs = {}
        graphs['nsig'] = r.TGraph(len(xvalues), np.array(xvalues,dtype=float), np.array(nsigs,dtype=float))
        graphs['nbkg'] = r.TGraph(len(xvalues), np.array(xvalues,dtype=float), np.array(nbkgs,dtype=float))
        graphs['zbis'] = r.TGraph(len(xvalues), np.array(xvalues,dtype=float), np.array(zbis,dtype=float))
        graphs['cut_values'] = r.TGraph(len(xvalues), np.array(xvalues,dtype=float), np.array(cut_values,dtype=float))
        graphs['eff'] = r.TGraph(len(xvalues), np.array(xvalues,dtype=float), np.array(effs,dtype=float))
        for key, graph in graphs.items():
            outfile[f'{zalpha_subdir}/{key}_g'] = graph
            
    with PdfPages(f'{outdir}/mass_{mass}_zalpha_slope_plots.pdf') as pdf:
        for fig in figures.values():
            fig.axes[0].set_ylim(bottom=1,top=500000)
            pdf.savefig(fig)
            plt.close('all')

