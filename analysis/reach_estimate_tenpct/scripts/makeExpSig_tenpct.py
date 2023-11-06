#/bin/env python
import sys
import plot_utilities as utils
import math
import glob
import numpy as np
import ROOT as r
import copy
from optparse import OptionParser
from SimpEquations import SimpEquations

#Calculated in 'makeRadFrac.py'
def radiativeFraction(mass):
    radF = -1.04206e-01 + 9.92547e-03*mass + -1.99437e-04*pow(mass,2) + 1.83534e-06*pow(mass,3) + -7.93138e-9*pow(mass,4) + 1.30456e-11*pow(mass,5) #alic 2016 simps kf 11/15/22
    return radF

#Calculated in 'makeTotRadAcc.py'
def radiativeAcceptance(mass):
    acc = ( -7.35934e-01 + 9.75402e-02*mass + -5.22599e-03*pow(mass,2) + 1.47226e-04*pow(mass,3) + -2.41435e-06*pow(mass,4) + 2.45015e-08*pow(mass,5) + -1.56938e-10*pow(mass,6) + 6.19494e-13*pow(mass,7) + -1.37780e-15*pow(mass,8) + 1.32155e-18*pow(mass,9) ) #alic 2016 simps kf 11/15/22 
    return acc

#Calculated by running 'vtxhProcess.py' in hpstr, then using 'makeVtxResolution.py'
def vtxRes(mass):
    mass = mass/1000.0 #cnv MeV to GeV 
    res = 8.09149 + -1.54072e02*mass + 1.25624e03*pow(mass,2) + -3.43499e03*pow(mass,3) # 2016 simps kf 11/15/22
    return res

#Calculated in 'makeMassRes.py'
def massRes(mass):
    res = 1.06314 + 3.45955e-02*mass + -6.62113e-05*pow(mass,2) # 2016 simps kf 11/15/22
    return res

def countDiffBackgroundData(m_Ap, infile_data, tree_name, massRes):
    dNdm = 0.0
    datafile = r.TFile("%s"%(infile_data),"READ")
    tree = datafile.Get("%s/%s_tree"%(tree_name,tree_name))
    tree.SetName("data_%s_tree"%(tree_name))
    print("Counting background rate")
    Mbin = 30.0
    for ev in tree:
        if 1000.0*ev.unc_vtx_mass > m_Ap + (Mbin/2): continue
        if 1000.0*ev.unc_vtx_mass < m_Ap - (Mbin/2): continue
        dNdm += 1.0
        pass

    datafile.Close()

    dNdm = dNdm/Mbin
    print("Background Rate: %f"%dNdm)

    return dNdm


def countDiffBackgroundMC(m_Ap, infile_tritrig, infile_wab, tree_name, tritrig_mcScale, wab_mcScale, massRes):
    
    dNdm = 0.0
    #tritrig
    ttFile = r.TFile("%s"%(infile_tritrig),"READ")
    ttTree = ttFile.Get("%s/%s_tree"%(tree_name,tree_name))
    ttTree.SetName("tritrig_%s_tree"%(tree_name))
    print("Counting background rate")
    Mbin = 30.0
    for ev in ttTree:
        if 1000.0*ev.unc_vtx_mass > m_Ap + (Mbin/2): continue
        if 1000.0*ev.unc_vtx_mass < m_Ap - (Mbin/2): continue
        dNdm += tritrig_mcScale
        pass
    ttFile.Close()

    #WAB
    wabFile = r.TFile("%s"%(infile_wab))
    wabTree = wabFile.Get("%s/%s_tree"%(tree_name, tree_name))
    wabTree.SetName("wab_Tight_tree")
    for ev in wabTree:
        if 1000.0*ev.unc_vtx_mass > m_Ap + (Mbin/2): continue
        if 1000.0*ev.unc_vtx_mass < m_Ap - (Mbin/2): continue
        dNdm += wab_mcScale
        pass
    wabFile.Close()

    dNdm = dNdm/Mbin
    print("Background Rate: %f"%dNdm)
    return dNdm

#This cut is a good first high-z background cut, helps get more realistic zcut
def pass_deltaZ(ele_z0, ele_tanlambda, pos_z0, pos_tanlambda, cut_value):
    value = abs((ele_z0/ele_tanlambda) - (pos_z0/pos_tanlambda) )
    if value < cut_value:
        return True
    else:
        return False

#Suppress plots
r.gROOT.SetBatch(1)

#2016 Lumi from Golden runs
Lumi = 10.7 #1/pb
mcScale = {}
mcScale['tritrig'] = 1.416e9*Lumi/(50000*9853) #pb2016
mcScale['wab'] = 0.1985e12*Lumi/(100000*9966) #pb2016

#utils.SetStyle()

parser = OptionParser()

parser.add_option("-o", "--outputFile", type="string", dest="outputFile",
        help="Specify the output filename.", metavar="outputFile", default="expSigRate.root")
parser.add_option("-x", "--exclusionContourFile", type="string", dest="exContFile",
        help="Specify the output exlusion contour filename.", metavar="exContFile", default="exContour.txt")
parser.add_option("-p", "--darkPionDecayConstRatio", type=float, dest="darkPionDecayConstRatio",
        help="Specify the Dark Pion Decay Constant (3 or 4pi).", metavar="darkPionDecayConstRatio", default=12.556)
parser.add_option("-z", "--zcutFile", type="string", dest="zcutFile",
        help="Name of file containing zcut values.", metavar="zcutFile", default="/sdf/group/hps/users/alspellm/projects/THESIS/simp_reach_estimates/simps_2016/reach_estimate/makeComponents/components/zcuts.dat")

(options, args) = parser.parse_args()

outFile = r.TFile(options.outputFile,"RECREATE")

zCutVals = []
dNdms = []
#MC Generated Dark Vector Masses
#invMasses = [40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200]
invMasses = [40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160]
#invMasses = [40, 50, 60, 70, 80, 200]

#Simp Equations Config
year = 2016
alpha_dark = 0.01
mass_ratio_Ap_to_Vd = 1.66
mass_ratio_Ap_to_Pid = 3.0
ratio_mPi_to_fPi = float(options.darkPionDecayConstRatio)
lepton_mass = 0.511
target_pos = -4.3
simpeqs = SimpEquations(year, alpha_dark, mass_ratio_Ap_to_Vd, mass_ratio_Ap_to_Pid, ratio_mPi_to_fPi, lepton_mass)

# SIMP MC is generated using fixed mass ratio of Dark Vector to A'
ap_invMasses = [round(x*(mass_ratio_Ap_to_Vd),1) for x in invMasses]

#Read zcuts from file if specified
#zcuts = {}
#zcutFile = open(options.zcutFile,"r")
#for line in zcutFile:
#    lineList = line.split()
#    zcuts[float(lineList[0])] = float(lineList[1])
#    pass
#print(zcuts)

#Initialize Histograms
nMasses = len(ap_invMasses)
lowM = float(ap_invMasses[0] - 5.0)
highM = float(ap_invMasses[-1] + 5.0)
logEps2_min = -1400
logEps2_max = -100
epsBins = 650
logEps2_min = -1400
logEps2_max = -100
epsBins = 620
eps_min = -10.005
eps_max = -3.905
#apProd_hh = r.TH2D("apProd_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, 620, -10.005, -3.905)
#Nsig_hh = r.TH2D("Nsig_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, 620, -10.005, -3.905)
#Nsig_vd_hh = r.TH2D("Nsig_vd_hh", ";m_{V_{D}'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, 620, -10.005, -3.905)
#effVtx_rho_hh = r.TH2D("effVtx_rho_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, 620, -10.005, -3.905)
#effVtx_phi_hh = r.TH2D("effVtx_phi_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, 620, -10.005, -3.905)
apProd_hh = r.TH2D("apProd_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM,epsBins, eps_min, eps_max)
Nsig_hh = r.TH2D("Nsig_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, epsBins, eps_min, eps_max)
effVtx_rho_hh = r.TH2D("effVtx_rho_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, epsBins, eps_min, eps_max)

effVtx_phi_hh = r.TH2D("effVtx_phi_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, epsBins, eps_min, eps_max)
gcTau_rho_hh = r.TH2D("gcTau_rho_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, epsBins, eps_min, eps_max)
gcTau_phi_hh = r.TH2D("gcTau_phi_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, epsBins, eps_min, eps_max)

#Vd mass plots
nMasses = len(invMasses)
lowM_vd = float(invMasses[0] - 5.0)
highM_vd = float(invMasses[-1] + 5.0)
Nsig_vd_hh = r.TH2D("Nsig_m_{VD}_hh", ";m_{V_{D}'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM_vd, highM_vd, epsBins, eps_min, eps_max)
Nsig_rho_hh = r.TH2D("Nsig_rho_m_{VD}_hh", ";m_{m_{#rhoD}'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM_vd, highM_vd, epsBins, eps_min, eps_max)
Nsig_phi_hh = r.TH2D("Nsig_phi_m_{VD}_hh", ";m_{#phi_{D}'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM_vd, highM_vd, epsBins, eps_min, eps_max)
effVtx_rho_vd_hh = r.TH2D("effVtx_rho_m_{VD}_hh", ";m_{VD} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM_vd, highM_vd, epsBins, eps_min, eps_max)
effVtx_phi_vd_hh = r.TH2D("effVtx_phi_m_{VD}_hh", ";m_{VD} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM_vd, highM_vd, epsBins, eps_min, eps_max)


#Exclusion Contours for Nsig mean
upExContourMass = []
upExContourEps2 = []
downExContourMass = []
downExContourEps2 = []

#Exclusion Contours for NSig Upper
upExContourMass_up = []
upExContourEps2_up = []
downExContourMass_up = []
downExContourEps2_up = []

#Exclusion Contours for NSig Lower
upExContourMass_low = []
upExContourEps2_low = []
downExContourMass_low = []
downExContourEps2_low = []

#Looping over A' masses, NOT VECTOR MASSES
for m_Vd in invMasses:

    m_Ap = round(m_Vd * mass_ratio_Ap_to_Vd,1) 

    print("Running A' mass = %i MeV"%m_Ap)

    #low stats causes failure when totalRadiativeAcceptance = 0. Skip those masses
    radAcc = radiativeAcceptance(m_Ap)
    if radAcc < 0:
        radAcc = 0.0
        continue
    
    radFrac = radiativeFraction(m_Ap)

    #First grab the pretrigger vtx z distribution
    vdSimFilename = "/sdf/group/hps/users/alspellm/projects/THESIS/mc/2016/simps/slic/20230713_slic/20230724_slic_ana/ana_files/hadd_simp_%i_MeV_rot_slic_mcana.root"%(m_Vd)
    vdSimFile = r.TFile(vdSimFilename)
    vdSimZ_hcp = copy.deepcopy(vdSimFile.Get("mcAna/mcAna_mc625Z_h") )
    vdSimFile.Close()
    vdSimZ_hcp.SetName("vdSimZ%i_hcp"%m_Vd)
    vdSimZ_h = r.TH1F("vdSimZ%i_h"%m_Vd, ";true z_{vtx} [mm];MC Events", 200, -50.3, 149.7)
    for i in range(201):
        vdSimZ_h.SetBinContent(i, vdSimZ_hcp.GetBinContent(i))
        pass
    outFile.cd()
    vdSimZ_h.Write()

    #Next count the differential background rate in 1 MeV bin
    infile_data = '/sdf/group/hps/users/alspellm/projects/THESIS/data/2016/BLPass4_rereco_20230823/ana_20230823/output/final_hadd_data_BLPass4_ReRecon_20230828.root'
    infile_tritrig = '/sdf/group/hps/users/alspellm/projects/THESIS/mc/2016/tritrig_beam/pass4_2016_mc/rerecon_kf_v5_1/simps_2016_kf/final_hadd_tritrigv2-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_976_KF_CR.root'
    infile_wab = '/sdf/group/hps/users/alspellm/projects/THESIS/mc/2016/wab_beam/pass4_2016_mc/rerecon_kf_v5_1/simps_2016_kf/final_hadd_wabv3-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_KF_ana_CR.root'
    #dNdm = countDiffBackgroundMC(m_Ap, infile_tritrig, infile_wab,"vtxana_kf_vertexSelection_Tight_CR", mcScale['tritrig'], mcScale['wab'],
    dNdm = 10.0*countDiffBackgroundData(m_Ap, infile_data,"vtxana_kf_Tight_2016_simp_reach_CR",massRes(float(m_Vd)))
    dNdms.append(dNdm)

    #Next get flat tuple from anaVtx and fill eff_vtx numerator
    lowMass = float(m_Vd) - 2.8*massRes(float(m_Vd))/2.0
    highMass = float(m_Vd) + 2.8*massRes(float(m_Vd))/2.0

    #### Zcut ####
    #zCut = 69.2555 + -0.916318*m_vdF + 0.00504772*m_vdF*m_vdF + -1.04964e-05*pow(m_vdF,3) #simp kf 11/15/22
    #zCut = 72.6454 + -1.08843*m_vdF + 0.00676883*m_vdF*m_vdF + -1.52583e-05*pow(m_vdF,3) #simp gbl old
    #zCut = 9.71425 + -0.140865*m_vdF + 0.000441817*math.pow(m_vdF,2) + -4.73974e-07*math.pow(m_vdF,3) #GBL
    #zCut = 69.2555 + -0.916318*float(m_Vd) + 0.00504772*float(m_Vd)*float(m_Vd) + -1.04964e-05*pow(float(m_Vd),3) #simp kf 11/15/22
    zCut = 51.6922 + -0.443276*(float(m_Vd)) + 0.000681226*math.pow(float(m_Vd),2) + 3.32196e-06*math.pow(float(m_Vd),3) #simp kf with abs(delta z0/tanlambda = 30.0, nbkg = 0.5) 
    #zCut = zcuts[m_Vd] 
    zCutVals.append(zCut)
    vdSelZ_h = r.TH1F("vdSelZ%i_h"%m_Vd, ";true z_{vtx} [mm];MC Events", 200, -50.3, 149.7)
    vdSelNoZ_h = r.TH1F("vdSelNoZ%i_h"%m_Vd, ";true z_{vtx} [mm];MC Events", 200, -50.3, 149.7)
    vdFilename = "/sdf/group/hps/users/alspellm/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_%i_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"%(m_Vd)
    vd_energy_h = r.TH1F("vd_%s_truth_energy"%(m_Vd),"vd_%s_truth_energy;Truth Energy [GeV];Events"%(m_Vd),250,0.0,2.5)
    vdFile = r.TFile("%s"%(vdFilename),"READ")
    tree_name = 'vtxana_kf_radMatchTight_2016_simp_reach_SR'
    vdTree = vdFile.Get("%s/%s_tree"%(tree_name, tree_name))
    print("Counting Signal")
    for ev in vdTree:
        if 1000.0*ev.unc_vtx_mass > highMass: continue
        if 1000.0*ev.unc_vtx_mass < lowMass: continue
        vd_energy_h.Fill(ev.vd_true_vtx_energy)
       
        #Apply high-z cut used to set the zcut values 08/31/2023
        if not pass_deltaZ(ev.unc_vtx_ele_track_z0, ev.unc_vtx_ele_track_tanLambda, 
                ev.unc_vtx_pos_track_z0, ev.unc_vtx_pos_track_tanLambda, 30.0):
            continue

        if ev.vd_true_vtx_z > 135.0: continue
        vdSelNoZ_h.Fill(ev.vd_true_vtx_z)
        if ev.unc_vtx_z < zCut: continue
        vdSelZ_h.Fill(ev.vd_true_vtx_z)
        pass
    vdFile.Close()

    #Get mean VD Energy
    E_Vd = vd_energy_h.GetMean()
    print("Dark Vector Mass %i Energy: %f"%(m_Vd, E_Vd))


    #Make the efficiencies
    vdEffVtxZ_gae = r.TGraphAsymmErrors(vdSelZ_h, vdSimZ_h, "shortest")
    vdEffVtxZ_gae.SetName("vdEffVtxZ%i_gae"%m_Vd)
    vdEffVtxZ_e = r.TEfficiency(vdSelZ_h, vdSimZ_h)
    vdEffVtxZ_e.SetName("vdEffVtxZ%i_e"%m_Vd)
    effCalc_h = vdEffVtxZ_e

    outFile.cd()
    vdSelNoZ_h.Write()
    vdSelZ_h.Write()
    vdEffVtxZ_gae.Write()
    vdEffVtxZ_e.Write()

    prevRate = 0.0
    prevRate_low = 0.0
    prevRate_up = 0.0

    excThr = 2.3
    print("Calculate expected signal rate")

    epsilons = []
    effVtxs = []
    gctaus = []
    apProduced = []
    NSigs = []
    checks = []

    length = 0
    for logEps2 in range(logEps2_min, logEps2_max):
        length = length+1
        Ap_prod = 0.0
        Nsig_rho = 0.0
        Nsig_phi = 0.0
        Nsig = 0.0

        logEps2 = (round(logEps2/100.0,2))
        eps2 = pow(10, logEps2)
        eps = float(np.sqrt(eps2))
        epsilons.append(logEps2)
        
        results_rho = simpeqs.expectedSignalCalculation(float(m_Vd), eps, True, False, E_Vd, 
                radFrac, radAcc, dNdm, effCalc_h, target_pos, zCut)

        results_phi = simpeqs.expectedSignalCalculation(float(m_Vd), eps, False, True, E_Vd, 
                radFrac, radAcc, dNdm, effCalc_h, target_pos, zCut)

        #Total A' production at target
        tot_apProd = results_rho["tot_apProd"]
        #Rho meson
        Nsig_rho = results_rho["signal"]
        effVtx_rho = results_rho["effVtx"]
        gcTau_rho = results_rho["gcTau"]
        br_VPi_rho = results_rho["br_VPi"]
        '''
        print("Nsig_rho: ", Nsig_rho)
        print("effVtx_rho: ", effVtx_rho)
        print("gcTau_rho: ", gcTau_rho)
        print("br_VPi_rho: ", br_VPi_rho)
        '''
        #Phi meson
        Nsig_phi = results_phi["signal"]
        effVtx_phi = results_phi["effVtx"]
        gcTau_phi = results_phi["gcTau"]
        br_VPi_phi = results_phi["br_VPi"]
        '''
        print("Nsig_phi: ", Nsig_phi)
        print("effVtx_phi: ", effVtx_phi)
        print("gcTau_phi: ", gcTau_phi)
        print("br_VPi_phi: ", br_VPi_phi)
        '''
        Nsig = Nsig_rho + Nsig_phi

        NSigs.append(Nsig)

        #Fill histos for signal counts, etc.
        apProd_hh.Fill(m_Ap, logEps2, tot_apProd)
        Nsig_hh.Fill(m_Ap, logEps2, Nsig)
        effVtx_rho_hh.Fill(m_Ap, logEps2, effVtx_rho)
        effVtx_phi_hh.Fill(m_Ap, logEps2, effVtx_phi)
        gcTau_rho_hh.Fill(m_Ap, logEps2, gcTau_rho)
        gcTau_phi_hh.Fill(m_Ap, logEps2, gcTau_phi)

        #Fill VD mass histograms
        Nsig_vd_hh.Fill(float(m_Vd), logEps2, Nsig)
        Nsig_rho_hh.Fill(float(m_Vd), logEps2, Nsig_rho)
        Nsig_phi_hh.Fill(float(m_Vd), logEps2, Nsig_phi)
        effVtx_rho_vd_hh.Fill(float(m_Vd), logEps2, effVtx_rho)
        effVtx_phi_vd_hh.Fill(float(m_Vd), logEps2, effVtx_phi)

        #Build contour by checking for signal to cross above or below detection threshold (2.3 Events)
        #We build three sets of contours, for NSig upper, NSig mean, and NSig lower
        if prevRate < excThr and Nsig > excThr:
            downExContourMass.append(m_Ap)
            downExContourEps2.append(logEps2)
        if prevRate > excThr and Nsig < excThr:
            upExContourMass.append(m_Ap)
            upExContourEps2.append(logEps2)
        prevRate = Nsig
        pass

    #debug
    '''
    effVtxEps_g = r.TGraph(len(epsilons),np.array(epsilons),np.array(effVtxs))
    effVtxEps_g.SetName("effVtx_eps_%i"%(m_Vd))
    effVtxEps_g.Draw()
    effVtxEps_g.Write()

    gamctauEps_g = r.TGraph(len(epsilons),np.array(epsilons),np.array(gctaus))
    gamctauEps_g.SetName("gctau_eps_%i"%(m_Vd))
    gamctauEps_g.Draw()
    gamctauEps_g.Write()

    totApEps_g = r.TGraph(len(epsilons),np.array(epsilons),np.array(apProduced))
    totApEps_g.SetName("produced As_%i"%(m_Ap))
    totApEps_g.Draw()
    totApEps_g.Write()

    NSigEps_g = r.TGraph(len(epsilons),np.array(epsilons),np.array(NSigs))
    NSigEps_g.SetName("Nsig_%i"%(m_Ap))
    NSigEps_g.Draw()
    NSigEps_g.Write()
    '''
    pass
#Contour Plots
upExContourMass.reverse()
upExContourEps2.reverse()
exContourMass = upExContourMass + downExContourMass
exContourEps2 = upExContourEps2 + downExContourEps2
exContourEps = [math.sqrt(pow(10,x)) for x in exContourEps2]
if(len(exContourEps) > 0):
    contOutFile = open("%s_mean.txt"%(options.exContFile),"w")
    for i in range(len(exContourMass)):
        contOutFile.write("%f\t%E\n"%(exContourMass[i], exContourEps[i]))
        pass
    contOutFile.close()
    exContourEps.append(exContourEps[0])
    exContourMass.append(exContourMass[0])
    excContour_g = r.TGraph(len(exContourMass), np.array(exContourMass), np.array(exContourEps))
    excContour_g.SetName("excContour_Lumi_%s_g"%str(Lumi))
    excContour_g.Write()

zCuts_g = r.TGraph(len(invMasses), np.array([float(x) for x in invMasses]), np.array(zCutVals))
zCuts_g.SetName("zCuts_g")
zCuts_g.Write()

dNdm_g = r.TGraph(len(ap_invMasses), np.array([float(x) for x in ap_invMasses]), np.array(dNdms))
dNdm_g.SetName("dNdm_g")
dNdm_g.Write()


#Format Plots
write_text = ["2016 MC Reach Estimate","Tritrig+WAB+Beam Scaled Lumi %s"%(str(Lumi)),"aD=%s"%(str(alpha_dark))]
c_apProd = utils.Make2DPlot(apProd_hh.GetName(), apProd_hh, apProd_hh.GetXaxis().GetTitle(), apProd_hh.GetYaxis().GetTitle(), insertText=write_text)
c_Nsig = utils.Make2DPlot(Nsig_hh.GetName(), Nsig_hh, Nsig_hh.GetXaxis().GetTitle(), Nsig_hh.GetYaxis().GetTitle(), insertText=write_text)
c_gcTau_rho = utils.Make2DPlot(gcTau_rho_hh.GetName(), gcTau_rho_hh, gcTau_rho_hh.GetXaxis().GetTitle(), gcTau_rho_hh.GetYaxis().GetTitle(), insertText=write_text)
c_gcTau_phi = utils.Make2DPlot(gcTau_phi_hh.GetName(), gcTau_phi_hh, gcTau_phi_hh.GetXaxis().GetTitle(), gcTau_phi_hh.GetYaxis().GetTitle(), insertText=write_text)
c_effVtx_rho = utils.Make2DPlot(effVtx_rho_hh.GetName(), effVtx_rho_hh, effVtx_rho_hh.GetXaxis().GetTitle(), effVtx_rho_hh.GetYaxis().GetTitle(), insertText=write_text)
c_effVtx_phi = utils.Make2DPlot(effVtx_phi_hh.GetName(), effVtx_phi_hh, effVtx_phi_hh.GetXaxis().GetTitle(), effVtx_phi_hh.GetYaxis().GetTitle(), insertText=write_text)
c_Nsig_vd = utils.Make2DPlot(Nsig_vd_hh.GetName(), Nsig_vd_hh, Nsig_vd_hh.GetXaxis().GetTitle(), Nsig_vd_hh.GetYaxis().GetTitle(), insertText=write_text)
c_Nsig_rho = utils.Make2DPlot(Nsig_rho_hh.GetName(), Nsig_rho_hh, Nsig_rho_hh.GetXaxis().GetTitle(), Nsig_rho_hh.GetYaxis().GetTitle(), insertText=write_text)
c_Nsig_phi = utils.Make2DPlot(Nsig_phi_hh.GetName(), Nsig_phi_hh, Nsig_phi_hh.GetXaxis().GetTitle(), Nsig_phi_hh.GetYaxis().GetTitle(), insertText=write_text)
c_effVtx_rho_vd = utils.Make2DPlot(effVtx_rho_vd_hh.GetName(), effVtx_rho_vd_hh, effVtx_rho_vd_hh.GetXaxis().GetTitle(), effVtx_rho_vd_hh.GetYaxis().GetTitle(), insertText=write_text)
c_effVtx_phi_vd = utils.Make2DPlot(effVtx_phi_vd_hh.GetName(), effVtx_phi_vd_hh, effVtx_phi_vd_hh.GetXaxis().GetTitle(), effVtx_phi_vd_hh.GetYaxis().GetTitle(), insertText=write_text)

c_apProd.Write()
c_Nsig.Write()
c_gcTau_rho.Write()
c_gcTau_phi.Write()
c_effVtx_rho.Write()
c_effVtx_phi.Write()
c_Nsig_vd.Write()
c_Nsig_rho.Write()
c_Nsig_phi.Write()
c_effVtx_rho_vd.Write()
c_effVtx_phi_vd.Write()
#apProd_hh.Write()
#Nsig_hh.Write()
#gcTau_rho_hh.Write()
#gcTau_phi_hh.Write()
#effVtx_rho_hh.Write()
#effVtx_phi_hh.Write()
#Nsig_vd_hh.Write()
#Nsig_rho_hh.Write()
#Nsig_phi_hh.Write()
#effVtx_rho_vd_hh.Write()
#effVtx_phi_vd_hh.Write()
outFile.Close()
