#/bin/env python
import math
import glob
import numpy as np
import ROOT as r
import utilities as utils
import copy
import sys
from optparse import OptionParser

def radFrac(mass):
    #radF = (-3.08682e-1  + 9.65869e-3*mass + -9.86239e-5*pow(mass,2) + 4.81221e-7*pow(mass,3) + -1.11369e-9*pow(mass,4) + 9.718e-13*pow(mass,5) ) #alic 2019 simps 
    radF = -4.02992e-01 + 1.25999e-02*mass + -1.28959e-04*math.pow(mass,2) + 6.30470e-07*math.pow(mass,3) + -1.46441e-9*math.pow(mass,4) + 1.28540e-12*math.pow(mass,5) #2019 simps
    return radF

def totRadAcc(mass):
    #zeta = ( -7.61944e-1 + 5.54282e-2*mass + -1.62957e-3*pow(mass,2) + 2.51874e-05*pow(mass,3) + -2.26089e-7*pow(mass,4) + 1.24559e-9*pow(mass,5) + -4.27798e-12*pow(mass,6) + 8.92790e-15*pow(mass,7) + -1.03515e-17*pow(mass,8) + 5.10985e-21*pow(mass,9) ) #alic 2019 simps 
    radAcc = -7.77721e-1 + 5.67200e-2*mass + -1.67381e-3*math.pow(mass,2) + 2.60149e-05*math.pow(mass,3) + -2.35378e-7*math.pow(mass,4) + 1.31046e-9*math.pow(mass,5) + -4.56049e-12*math.pow(mass,6) + 9.67101e-15*math.pow(mass,7) + -1.14284e-17*math.pow(mass,8) + 5.76861e-21*math.pow(mass,9)
    return radAcc

def vtxRes(mass):
    mass = mass/1000.0 #cnv MeV to GeV 
    res = ( 2.74363 - 2.282014e1*mass + 1.27987e2*pow(mass,2) + -2.05207e2*pow(mass,3)) # 2019 simps
    return res

def massRes(mass):
    #res = (1.29712  + 1.93768e-02*mass + -3.67914e-07*pow(mass,2) + 9.77287e-8*pow(mass,3)) #alic 2019 simps 
    massRes = 1.46696 + 1.50421e-02*mass + 3.79468e-05*math.pow(mass,2) + -3.06407e-08*math.pow(mass,3) #MeV 2019 Simps
    return massRes

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

Lumi = 110. #1/pb

utils.SetStyle()

parser = OptionParser()

parser.add_option("-i", "--inputFile", type="string", dest="inputFile",
        help="Name of file to run on.", metavar="inputFile", default="toys/toys.root")
parser.add_option("-o", "--outputFile", type="string", dest="outputFile",
        help="Specify the output filename.", metavar="outputFile", default="expSigRate.root")
parser.add_option("-z", "--zcutFile", type="string", dest="zcutFile",
        help="Name of file containing zcut values.", metavar="zcutFile", default="/data/src/hpstr/run/reach/zcut/zcuts.dat")

(options, args) = parser.parse_args()
outFile = r.TFile(options.outputFile,"RECREATE")

#VD masses
invMasses = [30,35,40,45,50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200]
#AP masses
ap_invMasses = [round(x*(3./1.8),1) for x in invMasses]
nMasses = len(invMasses)
lowM = float(ap_invMasses[0] - 5.0)
highM = float(ap_invMasses[-1] + 5.0)

#Build list of totalRadiativeAcceptance for all A' masses
#Find largest A' mass where totRadAcc < 0, and set all lower masses to 0
totalRadiativeAcceptance={}
isNeg = False
for mass in ap_invMasses[::-1]:
    radacc = totRadAcc(mass)
    if radacc < 0:
        isNeg = True
    if isNeg:
        radacc = 0.0
    totalRadiativeAcceptance[mass] = radacc

print("totRadAcc: ", totalRadiativeAcceptance)

#Lumis = [0.0001,0.0005,0.001,0.01,0.1,0.2,0.5,1.0,2.0,10.0,50.,110.,200.]
Lumis = [110.]
for q, Lumi in enumerate(Lumis):
    
    print("LUMI: ", Lumi)

    #MC Scaling
    mcScale = {}
    mcScale['tritrig'] = 4.566e08*Lumi/(10000*1425) #pb2019
    mcScale['wab'] = 4.715e10*Lumi/(10000*9944) #pb2019

    #zcuts from dat file
    #zcuts = {}
    #zcutFile = open(options.zcutFile,"r")
    #for line in zcutFile:
    #    lineList = line.split()
    #    zcuts[float(lineList[0])] = float(lineList[1])
    #    pass
    #print(zcuts)
    
    #Plots
    apProd_hh = r.TH2D("apProd_Lumi_%s_hh"%str(Lumi), "apProd_Lumi_%s;m_{A'} [MeV];log_{10}(#epsilon^{2})"%str(Lumi), nMasses, lowM, highM, 620, -10.005, -3.905)
    Nsig_hh = r.TH2D("Nsig_Lumi_%s_hh"%str(Lumi), "Nsig_Lumi_%s;m_{A'} [MeV];log_{10}(#epsilon^{2})"%str(Lumi), nMasses, lowM, highM, 620, -10.005, -3.905)
    effVtx_hh = r.TH2D("effVtx_Lumi_%s_hh"%str(Lumi), "effVtx_Lumi_%s;m_{A'} [MeV];log_{10}(#epsilon^{2})"%str(Lumi), nMasses, lowM, highM, 620, -10.005, -3.905)
    gcTau_hh = r.TH2D("gcTau_Lumi_%s_hh"%str(Lumi), "gcTau_Lumi_%s;m_{A'} [MeV];log_{10}(#epsilon^{2})"%str(Lumi), nMasses, lowM, highM, 620, -10.005, -3.905)

    #Exclusion Contours
    upExContourMass = []
    upExContourEps2 = []
    downExContourMass = []
    downExContourEps2 = []

    #zcut values
    zCutVals = []
    #Background rate counted
    dNdms = []

    #Looping over A' masses, NOT VECTOR MASSES
    for m_Ap in ap_invMasses:

        #low stats causes failure when totalRadiativeAcceptance = 0. Skip those masses
        if totalRadiativeAcceptance[m_Ap] == 0.0:
            continue

        print("Running A' mass = %i MeV"%m_Ap)

        #SIMP Params in MeV units
        m_ApGeV = m_Ap/1000.0
        m_pi = m_Ap/3.0
        m_vdI = int(round(m_Ap*(1.8/3.0),0))
        m_vdF = float(m_vdI)
        m_vdFGeV = m_vdF/1000.0
        alpha_D = 0.01
        #f_pi = m_pi/3.
        f_pi = m_pi/(4*math.pi)
        m_l = 0.511
        #Scale VD lifetime
        lt_factor = 1.0

        #First grab the pretrigger vtx z distribution
        #vdSimFilename = "/sdf/group/hps/users/alspellm/projects/simps_2019/mc/simps/gen/slic/tuple_ana/hadd_mass_%i_simp_slic_ana.root"%(m_vdI)
        vdSimFilename = "/sdf/group/hps/users/alspellm/projects/simps_2019/mc/simps/gen/slic/tuple_ana/hadd_mass_%i_simp_400bins_ana.root"%(m_vdI)
        vdSimFile = r.TFile(vdSimFilename)
        vdSimZ_hcp = copy.deepcopy(vdSimFile.Get("mcAna/mcAna_mc625Z_h") )
        vdSimFile.Close()
        vdSimZ_hcp.SetName("vdSimZ%i_hcp"%m_vdI)
        #vdSimZ_h = r.TH1F("vdSimZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 200, -47.51, 152.49)
        vdSimZ_h = r.TH1F("vdSimZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 400, -47.5, 152.5)
        #vdSimZ_h = r.TH1F("vdSimZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 200, -50., 150.)
        for i in range(401):
            vdSimZ_h.SetBinContent(i, vdSimZ_hcp.GetBinContent(i))
            pass
        outFile.cd()
        if(q < 1):
            vdSimZ_h.Write()

        #Next count the differential background rate in 1 MeV bin
        dNdm = 0.0
        #tritrig
        ttFile = r.TFile("/sdf/group/hps/users/alspellm/projects/simps_2019/mc/tritrig_beam/hadd_tritrig_beam_ana_CR.root")
        ttTree = ttFile.Get("vtxana_kf_Tight_2019_simpCR/vtxana_kf_Tight_2019_simpCR_tree")
        ttTree.SetName("tritrig_Tight_tree")
        print("Counting background rate")
        Mbin = 30.0
        for ev in ttTree:
            if 1000.0*ev.unc_vtx_mass > m_Ap + (Mbin/2): continue
            if 1000.0*ev.unc_vtx_mass < m_Ap - (Mbin/2): continue
            dNdm += mcScale['tritrig']
            pass
        ttFile.Close()

        #WAB
        wabFile = r.TFile("/sdf/group/hps/users/alspellm/projects/simps_2019/mc/wab_beam/hadd_wab_beam_ana_CR.root")
        wabTree = wabFile.Get("vtxana_kf_Tight_2019_simpCR/vtxana_kf_Tight_2019_simpCR_tree")
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
        #zCut = zcuts[m_vdF] 
        #zCut = 0.727084 + -0.0221251*m_vdF + -0.000210165*pow(m_vdF,2) + 1.31807e-06*pow(m_vdF,3) + -1.89423e-09*pow(m_vdF,4) #handmade conservative zcut 2019 simps (OLD DONT USE)
        zCut = 9.71425 + -0.140865*m_vdF + 0.000441817*math.pow(m_vdF,2) + -4.73974e-07*math.pow(m_vdF,3) #More conservative zcut 2019 SIMPS 
        zCutVals.append(zCut)
        #vdSelZ_h = r.TH1F("vdSelZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 200, -47.51, 152.49)
        #vdSelNoZ_h = r.TH1F("vdSelNoZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 200, -47.51, 152.49)

        #vdSelZ_h = r.TH1F("vdSelZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 200, -50., 150.)
        #vdSelNoZ_h = r.TH1F("vdSelNoZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 200, -50., 150.)

        vdSelZ_h = r.TH1F("vdSelZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 400, -47.5, 152.5)
        vdSelNoZ_h = r.TH1F("vdSelNoZ%i_h"%m_vdI, ";true z_{vtx} [mm];MC Events", 400, -47.5, 152.5)

        vdFilename = "/sdf/group/hps/users/alspellm/projects/simps_2019/mc/simps/gen/recon/tuple_ana/hadd_mass_%i_simp_recon_ana.root"%(m_vdI)
        vdFile = r.TFile(vdFilename)
        vdTree = vdFile.Get("vtxana_kf_Tight_2019_simpSR/vtxana_kf_Tight_2019_simpSR_tree")
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

        outFile.cd()
        if(q<1):
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
                E_V = 2.3 #hardcoded based on selected V_D MC energy distribution...need to improve in future!
                gcTau = lt_factor * ctau * gamma(m_vdFGeV, E_V)
                gctaus.append(gcTau)

                effVtx = 0.0
                for zbin in range(1,401):
                    zz = vdSelZ_h.GetBinCenter(zbin)
                    zzLow = vdSelZ_h.GetBinLowEdge(zbin)
                    zzHigh = zzLow + vdSelZ_h.GetBinWidth(zbin)
                    if zz < -7.5: continue
                    effVtx += (r.TMath.Exp((-7.5-zz)/gcTau)/gcTau)*effCalc_h.GetEfficiency(zbin)*vdSelZ_h.GetBinWidth(zbin)
                    pass
                effVtxs.append(effVtx)

                #total production of A's before detector acceptance/eff
                tot_apProd = (3.*137/2.)*3.14159*(m_Ap*eps2*radFrac(m_Ap)*dNdm)/totalRadiativeAcceptance[m_Ap]

                #branching ratios
                br_Vpi_val = br_Vpi(m_Ap,m_pi,m_vdF,alpha_D,f_pi,rho,phi)
                #br_Vpi_val = br_Vpi(m_Ap,m_pi,m_vdF,alpha_D,m_pi/(4*math.pi),rho,phi) #maximize BR for ratio=3
                br_V_to_ee = 1.0 
                Nsig = Nsig + tot_apProd*effVtx*br_V_to_ee*br_Vpi_val
                if tot_apProd < 0.0:
                    print("Nsig: ", Nsig)

            apProd_hh.Fill(m_Ap, logEps2, tot_apProd)
            Nsig_hh.Fill(m_Ap, logEps2, Nsig)
            effVtx_hh.Fill(m_Ap, logEps2, effVtx)
            gcTau_hh.Fill(m_Ap, logEps2, gcTau)

            #Check if NSig > threshold
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

    outFile.cd()
    #Contour Plots
    upExContourMass.reverse()
    upExContourEps2.reverse()
    exContourMass = upExContourMass + downExContourMass
    exContourEps2 = upExContourEps2 + downExContourEps2
    exContourEps = [math.sqrt(pow(10,x)) for x in exContourEps2]
    if(len(exContourEps) > 0):
        #contOutFile = open("excContour.txt","w")
        #for i in range(len(exContourMass)):
        #    contOutFile.write("%f\t%E\n"%(exContourMass[i], exContourEps[i]))
        #    pass
        #contOutFile.close()
        print("len exContourMass", len(exContourMass))
        print("len exContourEps", len(exContourEps))
        exContourEps.append(exContourEps[0])
        exContourMass.append(exContourMass[0])
        excContour_g = r.TGraph(len(exContourMass), np.array(exContourMass), np.array(exContourEps))
        excContour_g.SetName("excContour_Lumi_%s_g"%str(Lumi))
        excContour_g.Write()

    zCuts_g = r.TGraph(len(invMasses), np.array([float(x) for x in invMasses]), np.array(zCutVals))
    zCuts_g.SetName("zCuts_g")

    dNdm_g = r.TGraph(len(ap_invMasses), np.array([float(x) for x in ap_invMasses]), np.array(dNdms))
    dNdm_g.SetName("dNdm_Lumi_%s_g"%str(Lumi))
    dNdm_g.Write()

    if(q<1):
        zCuts_g.Write()

    effVtx_hh.Write()
    Nsig_hh.Write()
    apProd_hh.Write()
    gcTau_hh.Write()

outFile.Close()
