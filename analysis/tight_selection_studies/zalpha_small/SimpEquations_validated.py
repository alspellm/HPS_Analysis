import math
import ROOT as r
class SimpEquations:

    def __init__(self, year = 2016, alpha_dark = 0.01, mass_ratio_Ap_to_Vd = 1.66, mass_ratio_Ap_to_Pid = 3.0, 
        ratio_mPi_to_fPi = 12.566, lepton_mass = 0.511):
        self.year = year
        self.alpha_dark = alpha_dark
        self.mass_ratio_Ap_to_Vd = mass_ratio_Ap_to_Vd
        self.mass_ratio_Ap_to_Pid = mass_ratio_Ap_to_Pid
        self.ratio_mPi_to_fPi = ratio_mPi_to_fPi
        self.lepton_mass = lepton_mass

    @staticmethod
    def rate_Ap_ee(m_Ap, eps):
        ml = 0.511
        r = ml/m_Ap
        coeff1 = ((1.0/137.0)*eps**2)/3.0
        coeff2 = (1.0 - 4.0*(r**2))**(0.5)
        coeff3 = (1.0 + 2.0*(r**2))*m_Ap
        return coeff1*coeff2*coeff3

    @staticmethod
    def rate_2pi(m_Ap, m_pi, m_V, alpha_dark):
        coeff = (2.0 * alpha_dark / 3.0) * m_Ap
        pow1 = math.pow((1 - (4 * m_pi * m_pi / (m_Ap * m_Ap))), 3 / 2.0)
        pow2 = math.pow(((m_V * m_V) / ((m_Ap * m_Ap) - (m_V * m_V))), 2)
        return coeff * pow1 * pow2

    @staticmethod
    def rate_Vrho_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi):
        x = m_pi / m_Ap
        y = m_V / m_Ap
        Tv = 3.0/4.0
        coeff = alpha_dark * Tv / (192.0 * math.pow(math.pi, 4))
        return coeff * math.pow((m_Ap / m_pi), 2) * math.pow(m_V / m_pi, 2) * math.pow((m_pi / f_pi), 4) * m_Ap * math.pow(SimpEquations.Beta(x, y), 3 / 2.0)

    def rate_Vphi_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi):
        x = m_pi / m_Ap
        y = m_V / m_Ap
        Tv = 3.0/2.0
        coeff = alpha_dark * Tv / (192.0 * math.pow(math.pi, 4))
        return coeff * math.pow((m_Ap / m_pi), 2) * math.pow(m_V / m_pi, 2) * math.pow((m_pi / f_pi), 4) * m_Ap * math.pow(SimpEquations.Beta(x, y), 3 / 2.0)

    def rate_Vcharged_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi):
        x = m_pi / m_Ap
        y = m_V / m_Ap
        Tv = 18.0 - ((3.0/2.0)+(3.0/4.0))
        coeff = alpha_dark * Tv / (192.0 * math.pow(math.pi, 4))
        return coeff * math.pow((m_Ap / m_pi), 2) * math.pow(m_V / m_pi, 2) * math.pow((m_pi / f_pi), 4) * m_Ap * math.pow(SimpEquations.Beta(x, y), 3 / 2.0)

    @staticmethod
    def br_2pi(m_Ap, m_pi, m_V, alpha_dark, f_pi):
        total_rate = SimpEquations.rate_Vrho_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) +SimpEquations.rate_Vphi_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) + SimpEquations.rate_Vcharged_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) + SimpEquations.rate_2pi(m_Ap, m_pi, m_V, alpha_dark)
        if m_Ap > 2.0*m_V:
            total_rate = total_rate + SimpEquations.rate_2V(m_Ap, m_V, alpha_dark)
        return SimpEquations.rate_2pi(m_Ap, m_pi, m_V, alpha_dark)/total_rate

    @staticmethod
    def br_Vrho_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi):
        total_rate = SimpEquations.rate_Vrho_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) +SimpEquations.rate_Vphi_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) + SimpEquations.rate_Vcharged_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) + SimpEquations.rate_2pi(m_Ap, m_pi, m_V, alpha_dark)
        if m_Ap > 2.0*m_V:
            total_rate = total_rate + SimpEquations.rate_2V(m_Ap, m_V, alpha_dark)
        return SimpEquations.rate_Vrho_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) / total_rate

    @staticmethod
    def br_Vphi_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi):
        total_rate = SimpEquations.rate_Vrho_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) +SimpEquations.rate_Vphi_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) + SimpEquations.rate_Vcharged_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) + SimpEquations.rate_2pi(m_Ap, m_pi, m_V, alpha_dark)
        if m_Ap > 2.0*m_V:
            total_rate = total_rate + SimpEquations.rate_2V(m_Ap, m_V, alpha_dark)
        return SimpEquations.rate_Vphi_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) / total_rate

    @staticmethod
    def br_Vcharged_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi):
        total_rate = SimpEquations.rate_Vrho_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) +SimpEquations.rate_Vphi_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) + SimpEquations.rate_Vcharged_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) + SimpEquations.rate_2pi(m_Ap, m_pi, m_V, alpha_dark)
        if m_Ap > 2.0*m_V:
            total_rate = total_rate + SimpEquations.rate_2V(m_Ap, m_V, alpha_dark)
        return SimpEquations.rate_Vcharged_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) / total_rate

    @staticmethod
    def br_2V(m_Ap, m_pi, m_V, alpha_dark, f_pi):
        total_rate = SimpEquations.rate_Vrho_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) +SimpEquations.rate_Vphi_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) + SimpEquations.rate_Vcharged_pi(m_Ap, m_pi, m_V, alpha_dark, f_pi) + SimpEquations.rate_2pi(m_Ap, m_pi, m_V, alpha_dark)
        if m_Ap > 2.0*m_V:
            total_rate = total_rate + SimpEquations.rate_2V(m_Ap, m_V, alpha_dark)
        if 2 * m_V >= m_Ap:
            return 0.0
        return SimpEquations.rate_2V(m_Ap, m_V, alpha_dark) / total_rate

    @staticmethod
    def Tv(rho, phi):
        if rho:
            return 3.0 / 4.0
        elif phi:
            return 3.0 / 2.0
        else:
            return 18.0

    @staticmethod
    def Beta(x, y):
        return (1 + math.pow(y, 2) - math.pow(x, 2) - 2 * y) * (1 + math.pow(y, 2) - math.pow(x, 2) + 2 * y)

    @staticmethod
    def rate_2V(m_Ap, m_V, alpha_dark):
        r = m_V / m_Ap
        return alpha_dark / 6.0 * m_Ap * SimpEquations.f(r)

    @staticmethod
    def f(r):
        # Define your function f(r) here
        # Example: return some_expression
        pass

    @staticmethod
    def rate_2l(m_Ap, m_pi, m_V, eps, alpha_dark, f_pi, m_l, rho):
        alpha = 1.0 / 137.0
        coeff = (16 * math.pi * alpha_dark * alpha * eps**2 * f_pi**2) / (3 * m_V**2)
        term1 = (m_V**2 / (m_Ap**2 - m_V**2))**2
        term2 = (1 - (4 * m_l**2 / m_V**2))**0.5
        term3 = 1 + (2 * m_l**2 / m_V**2)
        constant = 1 if not rho else 2
        return coeff * term1 * term2 * term3 * m_V * constant

    @staticmethod
    def getCtau(m_Ap, m_pi, m_V, eps, alpha_dark, f_pi, m_l, rho):
        #c = 3.00e10  # cm/s <! -- BUG
        c = 3.00e11 #mm/s
        hbar = 6.58e-22  # MeV*sec
        rate = SimpEquations.rate_2l(m_Ap, m_pi, m_V, eps, alpha_dark, f_pi, m_l, rho)  # MeV
        tau = hbar / rate
        ctau = c * tau
        return ctau

    @staticmethod
    def gamma(m_V, E_V):
        gamma = E_V / m_V
        return gamma

    @staticmethod
    def totalApProductionRate(m_Ap, eps, radFrac, radAcc, dNdm):
        # Total A' Production Rate
        apProduction = (3. * 137 / 2.) * 3.14159 * (m_Ap * eps * eps * radFrac * dNdm)/radAcc
        return apProduction

    def expectedSignalCalculation(self, m_V, eps, rho, phi, E_V, effCalc_h, target_pos, zcut):
        # Signal mass dependent SIMP parameters
        m_Ap = m_V * self.mass_ratio_Ap_to_Vd
        m_pi = m_Ap / self.mass_ratio_Ap_to_Pid
        f_pi = m_pi / self.ratio_mPi_to_fPi

        # Mass in MeV
        ctau = self.getCtau(m_Ap, m_pi, m_V, eps, self.alpha_dark, f_pi, self.m_l, rho)
        gcTau = ctau * self.gamma(m_V / 1000.0, E_V)  # E_V in GeV

        # Calculate the Efficiency Vertex (Displaced VD Acceptance)
        effVtx = 0.0
        for zbin in range(effCalc_h.GetTotalHistogram().GetNbinsX() + 1):
            zz = effCalc_h.GetTotalHistogram().GetBinLowEdge(zbin)
            if zz < zcut:
                continue
            effVtx += (math.exp((target_pos - zz) / gcTau) / gcTau) * \
                (effCalc_h.GetEfficiency(zbin) - effCalc_h.GetEfficiencyErrorLow(zbin)) * \
                effCalc_h.GetTotalHistogram().GetBinWidth(zbin)

        # Total A' Production Rate
        apProduction = (3. * 137 / 2.) * 3.14159 * (m_Ap * eps * eps * self.radiativeFraction(m_Ap) * self.controlRegionBackgroundRate(m_Ap)) \
            / self.radiativeAcceptance(m_Ap)

        # A' -> V+Pi Branching Ratio
        br_VPi = self.br_Vpi(m_Ap, m_pi, m_V, self.alpha_dark, f_pi, rho, phi)

        # Vector to e+e- BR = 1
        br_V_ee = 1.0

        # Expected Signal
        expSignal = apProduction * effVtx * br_VPi * br_V_ee

        return expSignal

    def expectedSignalCalculation(self, m_V, eps, rho, phi, E_V, radFrac, radAcc, dNdm, effCalc_h, target_pos, zcut):
        # Signal mass dependent SIMP parameters
        m_Ap = m_V * self.mass_ratio_Ap_to_Vd
        m_pi = m_Ap / self.mass_ratio_Ap_to_Pid
        f_pi = m_pi / self.ratio_mPi_to_fPi

        # Mass in MeV
        ctau = self.getCtau(m_Ap, m_pi, m_V, eps, self.alpha_dark, f_pi, self.lepton_mass, rho)
        gcTau = ctau * self.gamma(m_V / 1000.0, E_V)  # E_V in GeV

        # Calculate the Efficiency Vertex (Displaced VD Acceptance)
        effVtx = 0.0
        for zbin in range(effCalc_h.GetTotalHistogram().GetNbinsX() + 1):
            zz = effCalc_h.GetTotalHistogram().GetBinLowEdge(zbin)
            if zz < zcut:
                continue
            effVtx += (math.exp((target_pos - zz) / gcTau) / gcTau) * \
                (effCalc_h.GetEfficiency(zbin) - effCalc_h.GetEfficiencyErrorLow(zbin)) * \
                effCalc_h.GetTotalHistogram().GetBinWidth(zbin)

        # Total A' Production Rate
        apProduction = (3. * 137 / 2.) * 3.14159 * (m_Ap * eps * eps * radFrac * dNdm) \
            / radAcc

        # A' -> V+Pi Branching Ratio
        br_VPi = 0.0
        if(rho):
            br_VPi = self.br_Vrho_pi(m_Ap, m_pi, m_V, self.alpha_dark, f_pi)
        else:
            br_VPi = self.br_Vphi_pi(m_Ap, m_pi, m_V, self.alpha_dark, f_pi)

        # Vector to e+e- BR = 1
        br_V_ee = 1.0

        # Expected Signal
        expSignal = apProduction * effVtx * br_VPi * br_V_ee

        results = {"signal":expSignal, "tot_apProd":apProduction, "effVtx":effVtx, "gcTau":gcTau, "br_VPi":br_VPi}

        return results

################# temporary tools for calculating signal as of 20230927 ############
    @staticmethod
    #Calculated in 'makeRadFrac.py'
    def radiativeFraction(mass):
        radF = -1.04206e-01 + 9.92547e-03*mass + -1.99437e-04*pow(mass,2) + 1.83534e-06*pow(mass,3) + -7.93138e-9*pow(mass,4) + 1.30456e-11*pow(mass,5) #alic 2016 simps kf 11/15/22
        return radF

    #Calculated in 'makeTotRadAcc.py'
    @staticmethod
    def radiativeAcceptance(mass):
        acc = ( -7.35934e-01 + 9.75402e-02*mass + -5.22599e-03*pow(mass,2) + 1.47226e-04*pow(mass,3) + -2.41435e-06*pow(mass,4) + 2.45015e-08*pow(mass,5) + -1.56938e-10*pow(mass,6) + 6.19494e-13*pow(mass,7) + -1.37780e-15*pow(mass,8) + 1.32155e-18*pow(mass,9) ) #alic 2016 simps kf 11/15/22 
        return acc

    #Calculated in 'makeMassRes.py'
    @staticmethod
    def massRes(mass):
        res = 1.06314 + 3.45955e-02*mass + -6.62113e-05*pow(mass,2) # 2016 simps kf 11/15/22
        return res

    @staticmethod
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

    @staticmethod
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

    
