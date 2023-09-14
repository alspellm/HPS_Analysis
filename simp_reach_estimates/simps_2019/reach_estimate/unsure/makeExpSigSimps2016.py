#/bin/env python
import math
import glob
import numpy as np
import ROOT as r
import utilities as utils
import copy
from optparse import OptionParser

def radFrac(mass):
    radF = ( -1.92497e-01 + 1.47144e-02*mass + -2.91966e-04*pow(mass,2) + 2.65603e-06*pow(mass,3) + -1.12471e-8*pow(mass,4) + 1.74765e-11*pow(mass,5) + 2.235718e-15*pow(mass,6)) #alic 2016 simps 
    return radF

def totRadAcc(mass):
    zeta = ( -7.93151e-01 + 1.04324e-01*mass + -5.55225e-03*pow(mass,2) + 1.55480e-04*pow(mass,3) + -2.53281e-06*pow(mass,4) + 2.54558e-08*pow(mass,5) + -1.60877e-10*pow(mass,6) + 6.24627e-13*pow(mass,7) + -1.36375e-15*pow(mass,8) + 1.28312e-18*pow(mass,9) ) #alic 2016 simps 
    return zeta

def vtxRes(mass):
    mass = mass/1000.0 #cnv MeV to GeV 
    res = ( 5.99358 - 8.22402e01*mass + 4.91751e02*pow(mass,2) + -9.98972e02*pow(mass,3)) # 2016 simps
    return res

def massRes(mass):
    res = 9.73217e-01 + 3.63659e-02*mass + -7.32046e-05*mass*mass #2016 simps alic
    return res

def rate_2pi(m_Ap,m_pi,m_V,alpha_D):
    coeff = (2*alpha_D/3) * m_Ap
    pow1 = (1-(4*(m_pi)**2/(m_Ap**2)))**(3/2.)
    pow2 = ((m_V**2)/((m_Ap**2)-(m_V**2)))**2
    return coeff * pow1 * pow2

def rate_Vpi(m_Ap,m_pi,m_V,alpha_D,f_pi,rho,phi):
    x = m_pi/m_Ap
    y = m_V/m_Ap
    pi = 3.14159
    coeff = alpha_D*Tv(rho,phi)/(192*(pi**4))
    return coeff * (m_Ap/m_pi)**2 * (m_V/m_pi)**2 * (m_pi/f_pi)**4 * m_Ap*(Beta(x,y))**(3./2.)

def br_Vpi(m_Ap,m_pi,m_V,alpha_D,f_pi,rho,phi):
    rate = rate_Vpi(m_Ap,m_pi,m_V,alpha_D,f_pi,rho,phi) + rate_2pi(m_Ap,m_pi,m_V,alpha_D)
    if(2*m_V < m_Ap): rate = rate_Vpi(m_Ap,m_pi,m_V,alpha_D,f_pi,rho,phi) + rate_2pi(m_Ap,m_pi,m_V,alpha_D) + rate_2V(m_Ap,m_V,alpha_D)
    return rate_Vpi(m_Ap,m_pi,m_V,alpha_D,f_pi,rho,phi)/rate

def br_2V(m_Ap,m_pi,m_V,alpha_D,f_pi,rho,phi):
    if(2*m_V >= m_Ap): return 0.
    rate = rate_Vpi(m_Ap,m_pi,m_V1,alpha_D,f_pi,rho,phi) + rate_2pi(m_Ap,m_pi,m_V1,alpha_D) + rate_2V(m_Ap,m_V1,alpha_D)
    return rate_2V(m_Ap,m_V1,alpha_D)/rate

def Tv(rho,phi):
    if rho:
        return 3/4.
    elif phi:
        return 3/2.
    else:
        return 18

def Beta(x,y):
    return (1+y**2-x**2-2*y)*(1+y**2-x**2+2*y)

def rate_2V(m_Ap,m_V,alpha_D):
    r = m_V/m_Ap
    return alpha_D/6 * m_Ap * f(r)

def f(r):
    num = 1 + 16*r**2 - 68*r**4 - 48*r**6
    den = (1-r**2) ** 2
    return num/den * (1-4*r**2)**0.5

def rate_2l(m_Ap,m_pi,m_V,eps,alpha_D,f_pi,m_l,rho):
    alpha = 1/137.
    pi = 3.14159
    coeff = 16*pi*alpha_D*alpha*eps**2*f_pi**2/(3*m_V**2)
    term1 = (m_V**2/(m_Ap**2 - m_V**2))**2
    term2 = (1-(4*m_l**2/m_V**2))**0.5
    term3 = 1+(2*m_l**2/m_V**2)
    const = 1
    if rho:
        const = 2
    return coeff * term1 * term2 * term3 * m_V * const

def getCtau(m_Ap,m_pi,m_V,eps,alpha_D,f_pi,m_l,rho):
    c = 3.00e10 #cm/s
    hbar = 6.58e-22 #MeV*sec
    rate = rate_2l(m_Ap,m_pi,m_V,eps,alpha_D,f_pi,m_l,rho)#MeV
    tau = hbar/rate
    ctau = c*tau
    return ctau

def Vdistribution(z,targZ,gammact):
    return np.exp(targZ/gammact-1/gammact*z)/gammact

def gamma(m_V, E_V):
    gamma = E_V/m_V
    return gamma

#2016 Lumi from Golden runs
Lumi = 10.7 #1/pb
mcScale = {}
mcScale['tritrig'] = 1.416e9*Lumi/(50000*9853) #pb2016
mcScale['wab'] = 0.1985e12*Lumi/(100000*9966) #pb2016

utils.SetStyle()

parser = OptionParser()

parser.add_option("-i", "--inputFile", type="string", dest="inputFile",
        help="Name of file to run on.", metavar="inputFile", default="toys/toys.root")
parser.add_option("-o", "--outputFile", type="string", dest="outputFile",
        help="Specify the output filename.", metavar="outputFile", default="expSigRate.root")
parser.add_option("-z", "--zcutFile", type="string", dest="zcutFile",
        help="Name of file containing zcut values.", metavar="zcutFile", default="/data/src/hpstr/run/reach/zcut/zcuts.dat")

(options, args) = parser.parse_args()

zCutVals = []
dNdms = []
invMasses = [40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200]
ap_invMasses = [round(x*(3./1.8),1) for x in invMasses]

zcuts = {}
zcutFile = open(options.zcutFile,"r")
for line in zcutFile:
    lineList = line.split()
    zcuts[float(lineList[0])] = float(lineList[1])
    pass
print(zcuts)


outFile = r.TFile(options.outputFile,"RECREATE")

nMasses = len(invMasses)
lowM = float(ap_invMasses[0] - 5.0)
highM = float(ap_invMasses[-1] + 5.0)
apProd_hh = r.TH2D("apProd_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, 620, -10.005, -3.905)
Nsig_hh = r.TH2D("Nsig_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, 620, -10.005, -3.905)
effVtx_hh = r.TH2D("effVtx_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, 620, -10.005, -3.905)
gcTau_hh = r.TH2D("gcTau_hh", ";m_{A'} [MeV];log_{10}(#epsilon^{2})", nMasses, lowM, highM, 620, -10.005, -3.905)

upExContourMass = []
upExContourEps2 = []
downExContourMass = []
downExContourEps2 = []

#Looping over A' masses, NOT VECTOR MASSES
for m_Ap in ap_invMasses:

    print("Running A' mass = %i MeV"%m_Ap)

    #SIMP Params in MeV units
    m_ApGeV = m_Ap/1000.0
    m_pi = m_Ap/3.0
    m_vdI = int(round(m_Ap*(1.8/3.0),0))
    m_vdF = float(m_vdI)
    m_vdFGeV = m_vdF/1000.0
    alpha_D = 0.001
    #f_pi = m_pi/3.
    f_pi = m_pi/(4*math.pi)
    m_l = 0.511

    #First grab the pretrigger vtx z distribution
    vdSimFilename = "/sdf/group/hps/users/alspellm/projects/simps_2016/reach_estimate/mc/simps/slic/simp_masses/%i/hadd_mass_%i_simp_mcAna.root"%(m_vdI, m_vdI)
    vdSimFile = r.TFile(vdSimFilename)
    vdSimZ_hcp = copy.deepcopy(vdSimFile.Get("mcAna/mcAna_mc625Z_h") )
    vdSimFile.Close()
    vdSimZ_hcp.SetName("vdSimZ%i_hcp"%m_vdI)
    vdSimZ_h = r.TH1F("vdSimZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 200, -50.3, 149.7)
    for i in range(201):
        vdSimZ_h.SetBinContent(i, vdSimZ_hcp.GetBinContent(i))
        pass
    outFile.cd()
    vdSimZ_h.Write()

    #Next count the differential background rate in 1 MeV bin
    dNdm = 0.0
    #tritrig
    ttFile = r.TFile("/sdf/group/hps/users/alspellm/projects/simps_2016/reach_estimate/mc/tritrig_beam/ana/hadd_tritrigv2-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_ana.root")
    ttTree = ttFile.Get("vtxana_Tight_simpCR/vtxana_Tight_simpCR_tree")
    ttTree.SetName("tritrig_Tight_tree")
    print("Counting background rate")
    Mbin = 10.0
    for ev in ttTree:
        if 1000.0*ev.unc_vtx_mass > m_Ap + (Mbin/2): continue
        if 1000.0*ev.unc_vtx_mass < m_Ap - (Mbin/2): continue
        dNdm += mcScale['tritrig']
        pass
    ttFile.Close()

    #WAB
    wabFile = r.TFile("/sdf/group/hps/users/alspellm/projects/simps_2016/reach_estimate/mc/wab_beam/ana/hadd_wabv3-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_ana.root")
    wabTree = wabFile.Get("vtxana_Tight_simpCR/vtxana_Tight_simpCR_tree")
    wabTree.SetName("wab_Tight_tree")
    for ev in wabTree:
        if 1000.0*ev.unc_vtx_mass > m_Ap + (Mbin/2): continue
        if 1000.0*ev.unc_vtx_mass < m_Ap - (Mbin/2): continue
        dNdm += mcScale['wab']
        pass
    dNdm = dNdm/Mbin
    dNdms.append(dNdm)
    wabFile.Close()
    print("Background Rate: %f"%dNdm)

    #Next get flat tuple from anaVtx and fill eff_vtx numerator
    lowMass = m_vdF - 2.8*massRes(m_vdF)/2.0
    highMass = m_vdF + 2.8*massRes(m_vdF)/2.0
    print("lowmass: ",lowMass)
    print("highmass: ",highMass)
    #zCut = 17.7702 + 138.166*m_vdFGeV - 5363.29*m_vdFGeV*m_vdFGeV + 44532.4*pow(m_vdFGeV,3) - 120578*pow(m_vdFGeV,4)
    zCut = zcuts[m_vdF] 
    zCutVals.append(zCut)
    vdSelZ_h = r.TH1F("vdSelZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 200, -50.3, 149.7)
    vdSelNoZ_h = r.TH1F("vdSelNoZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 200, -50.3, 149.7)
    vdFilename = "/sdf/group/hps/users/alspellm/projects/simps_2016/reach_estimate/mc/simps/recon/simp_masses/%i/hadd_mass_%i_simp_ana.root"%(m_vdI,m_vdI)
    vdFile = r.TFile(vdFilename)
    vdTree = vdFile.Get("vtxana_Tight_simpSIG/vtxana_Tight_simpSIG_tree")
    print("Counting Signal")
    for ev in vdTree:
        if 1000.0*ev.unc_vtx_mass > highMass: continue
        if 1000.0*ev.unc_vtx_mass < lowMass: continue
        if ev.true_vtx_z > 135.0: continue
        vdSelNoZ_h.Fill(ev.true_vtx_z)
        if ev.unc_vtx_z < zCut: continue
        vdSelZ_h.Fill(ev.true_vtx_z)
        pass
    vdFile.Close()

    #Make the efficiencies
    vdEffVtxZ_gae = r.TGraphAsymmErrors(vdSelZ_h, vdSimZ_h, "shortest")
    vdEffVtxZ_gae.SetName("vdEffVtxZ%i_gae"%m_vdI)
    vdEffVtxZ_e = r.TEfficiency(vdSelZ_h, vdSimZ_h)
    vdEffVtxZ_e.SetName("vdEffVtxZ%i_e"%m_vdI)
    tmpoutfile = r.TFile("checkme.root","RECREATE")
    tmpoutfile.cd()
    vdSelZ_h.Draw()
    vdSelZ_h.Write()
    vdSimZ_h.Draw()
    vdSimZ_h.Write()
    vdEffVtxZ_e.Draw()
    vdEffVtxZ_e.Write()
    tmpoutfile.Close()

    outFile.cd()
    vdSelNoZ_h.Write()
    vdSelZ_h.Write()
    vdEffVtxZ_gae.Write()
    vdEffVtxZ_e.Write()
    effCalc_h = vdEffVtxZ_e

    prevRate = 0.0
    excThr = 2.3
    print("Calculate expected signal rate")

    epsilons = []
    effVtxs = []
    gctaus = []
    apProduced = []
    NSigs = []
    checks = []
    for logEps2 in range(-1400, -100):
        tot_apProd = 0.
        Nsig = 0.
        logEps2 = logEps2/100.0
        eps2 = pow(10, logEps2)
        eps = float(np.sqrt(eps2))
        epsilons.append(logEps2)
        #Sensitive to two possible vector mesons produced. Add rates of each together for final rate
        #for vector_meson in ["rho","phi"]:
        for vector_meson in ["rho","phi"]:
            rho = False
            phi = False
            if "rho" in vector_meson:
                rho = True
            if "phi" in vector_meson:
                phi = True

            ctau = getCtau(m_Ap,m_pi,m_vdF,eps,alpha_D,f_pi,m_l,rho)
            E_V = 1.35 #hardcoded based on selected V_D MC energy distribution...need to improve in future!
            gcTau = ctau * gamma(m_vdFGeV, E_V)
            gctaus.append(gcTau)

            effVtx = 0.0
            check = 0.0
            for zbin in range(1,201):
                zz = vdSelZ_h.GetBinCenter(zbin)
                zzLow = vdSelZ_h.GetBinLowEdge(zbin)
                zzHigh = zzLow + vdSelZ_h.GetBinWidth(zbin)
                #Restric to beyond or at 2016 target position 
                #if zz < -4.3: continue
                if zz < 0.0: continue
                effVtx += (r.TMath.Exp((-4.3-zz)/gcTau)/gcTau)*effCalc_h.GetEfficiency(zbin)*vdSelZ_h.GetBinWidth(zbin)
                checks.append((r.TMath.Exp((-4.3-zz)/gcTau)/gcTau))
                #print("effCalc_h: ", effCalc_h.GetEfficiency(zbin))
                #print("vdSelZ_h width: ", vdSelZ_h.GetBinWidth(zbin))
                pass
            effVtxs.append(effVtx)

            #total production of A's before detector acceptance/eff
            #tot_apProd = (3.*137/2.)*3.14159*(mass*eps2*radFrac(mass)*dNdm)/totRadAcc(mass)  #I THINK THIS IS WRONG! DONT USE VD MASS
            tot_apProd = (3.*137/2.)*3.14159*(m_Ap*eps2*radFrac(m_Ap)*dNdm)/totRadAcc(m_Ap)
            apProduced.append(tot_apProd)
            #branching ratios
            br_Vpi_val = br_Vpi(m_Ap,m_pi,m_vdF,alpha_D,f_pi,rho,phi)
            br_V_to_ee = 1.0 
            #expected A' signal given V_D decays
            Nsig = Nsig + tot_apProd*effVtx*br_V_to_ee*br_Vpi_val

        if isinstance(Nsig, complex):
            Nsig = 0.0
        NSigs.append(Nsig)

        apProd_hh.Fill(m_Ap, logEps2, tot_apProd)
        Nsig_hh.Fill(m_Ap, logEps2, Nsig)
        effVtx_hh.Fill(m_Ap, logEps2, effVtx)
        gcTau_hh.Fill(m_Ap, logEps2, gcTau)
        if prevRate < excThr and Nsig > excThr:
            downExContourMass.append(m_Ap)
            downExContourEps2.append(logEps2)
            pass
        if prevRate > excThr and Nsig < excThr:
            upExContourMass.append(m_Ap)
            upExContourEps2.append(logEps2)
            pass
        prevRate = Nsig
        pass
    #debug
    effVtxEps_g = r.TGraph(len(epsilons),np.array(epsilons),np.array(effVtxs))
    effVtxEps_g.SetName("effVtx_eps_%i"%(m_vdI))
    effVtxEps_g.Draw()
    effVtxEps_g.Write()

    gamctauEps_g = r.TGraph(len(epsilons),np.array(epsilons),np.array(gctaus))
    gamctauEps_g.SetName("gctau_eps_%i"%(m_vdI))
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
    pass

upExContourMass.reverse()
upExContourEps2.reverse()
exContourMass = upExContourMass + downExContourMass
exContourEps2 = upExContourEps2 + downExContourEps2
exContourEps = [math.sqrt(pow(10,x)) for x in exContourEps2]
if(len(exContourEps) > 0):
    contOutFile = open("excContour.txt","w")
    for i in range(len(exContourMass)):
        contOutFile.write("%f\t%E\n"%(exContourMass[i], exContourEps[i]))
        pass
    contOutFile.close()
    print("len exContourMass", len(exContourMass))
    print("len exContourEps", len(exContourEps))
    excContour_g = r.TGraph(len(exContourMass), np.array(exContourMass), np.array(exContourEps))
    excContour_g.SetName("excContour_g")
    excContour_g.Write()

zCuts_g = r.TGraph(len(invMasses), np.array([float(x) for x in invMasses]), np.array(zCutVals))
zCuts_g.SetName("zCuts_g")
zCuts_g.Write()

dNdm_g = r.TGraph(len(ap_invMasses), np.array([float(x) for x in ap_invMasses]), np.array(dNdms))
dNdm_g.SetName("dNdm_g")
dNdm_g.Write()

apProd_hh.Write()
Nsig_hh.Write()
gcTau_hh.Write()
effVtx_hh.Write()
outFile.Close()
