#/bin/env python
import glob
import numpy as np
import ROOT as r
import utilities as utils
import copy
from optparse import OptionParser

def radFrac(mass):
    #mass = mass/1000.0
    #radF = ( -1.65979e-1 + 1.54816e1*mass - 3.09729e2*pow(mass,2) + 2.79836e3*pow(mass,3) - 1.18254e4*pow(mass,4) + 1.90196e4*pow(mass,5) )
    #radF = ( -2.61669e-1 + 8.16834*mass - 8.20396e1*pow(mass,2) + 3.90014e2*pow(mass,3) - 8.84654e2*pow(mass,4) + 7.65127e2*pow(mass,5) )
    #radF = ( -1.14400e-01 + 1.15665e01*mass - 2.36420e02*pow(mass,2) + 2.17845e03*pow(mass,3) - 9.38035e3*pow(mass,4) + 1.53550e4*pow(mass,5) ) #alic
    radF = ( 1.52474 + -1.60368e-01*mass + 7.22991e-03*pow(mass,2) + -1.78679e-04*pow(mass,3) + 2.69923e-06*pow(mass,4) + -2.60729e-08*pow(mass,5) + 1.62032e-10*pow(mass,6) + -6.27371e-13*pow(mass,7) + 1.37734e-15*pow(mass,8) + -1.30909e-18*pow(mass,9) ) #alic 2016 simps Vd mass = 60
    return radF

def totRadAcc(mass):
    zeta = ( -5.18069e-01 + 7.10979e-02*mass + -3.88854e-03*pow(mass,2) + 1.09598e-04*pow(mass,3) + -1.75343e-06*pow(mass,4) + 1.69916e-08*pow(mass,5) + -1.02208e-10*pow(mass,6) + 3.74354e-13*pow(mass,7) + -7.66468e-16*pow(mass,8) + 6.73775e-19*pow(mass,9) ) #alic 2016 simps Vd mass = 60
    return zeta

def vtxRes(mass):
    #res = ( 3.74129 - 63.6758*mass + 536.299*pow(mass,2) - 2101.23*pow(mass,3) + 3073.72*pow(mass,4) )
    res = ( 6.27069 - -8.99768e01*mass + 5.59386e02*pow(mass,2) - -1.17788e03*pow(mass,3)) # 2016 simps
    return res

def massRes(mass):
    #res = 3.40881 + 4.29919e-2*mass - 5.28684e-5*mass*mass
    #res = 1.91442 + 1.87619e-2*mass - 5.60051e-6*mass*mass
    #res = 2.53854 + 3.30200e-3*mass + 1.20876e-4*mass*mass #alic
    #res = 6.86902 - 9.88222e-2*mass + 1.34223e-3*mass*mass - 3.56501e-6*mass*mass*mass
    single_point = 1.968e-03 * 1000.
    return single_point

def calcLifetime(mass, eps2):
    return 8.0*(4.55/10.0)*(1e-8/eps2)*pow(100.0/mass, 2)
    #gamma = 0.95
    #hbar_c = 1.973e-13
    #ct = hbar_c*3.0/(mass*(1/137.036)*10**logEps2)
    #gammact = hbar_c*3.0*2.3*gamma/(massFGeV*massFGeV*(1/137.036)*10**logEps2)
    #print "logEps2: %f    gcTau: %f    gammact: %f    percent diff: %f"%(logEps2, gcTau, gammact, 1 - gammact/gcTau)
    #gcTau = gammact

def rate_2pi(m_Ap,m_pi,m_V,alpha_D):
    coeff = 2*alpha_D/3 * m_Ap
    pow1 = (1-(4*m_pi**2/(m_Ap**2)))**(3/2.)
    pow2 = (m_V**2/(m_Ap**2-m_V**2))**2
    return coeff * pow1 * pow2

def rate_Vpi(m_Ap,m_pi,m_V,alpha_D,f_pi,rho,phi):
    x = m_pi/m_Ap
    y = m_V/m_Ap
    pi = 3.14159
    coeff = alpha_D*Tv(rho,phi)/(192*pi**4)
    return coeff * 1/(x**2) * (y**2/(x**2)) * (m_pi/f_pi)**4 * m_Ap*Beta(x,y)**(3/2.)

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
    term2 = (1-4*m_l**2/m_V**2)**0.5
    term3 = 1+2*m_l**2/m_V**2
    const = 1
    if rho:
        const = 2
    return coeff * term1 * term2 * term3 * m_V * const

def getCtau(m_Ap,m_pi,m_V,eps,alpha_D,f_pi,m_l,rho):
    #hbar_c = 1.973e-14 #GeV*cm
    c = 3.00e10 #cm/s
    hbar = 6.58e-22 #MeV*sec
    rate = rate_2l(m_Ap,m_pi,m_V,eps,alpha_D,f_pi,m_l,rho)#MeV
    tau = hbar/rate  
    ctau = c*tau
    #ctau = hbar_c/rate
    return ctau

def Vdistribution(z,targZ,gammact):
    return np.exp(targZ/gammact-1/gammact*z)/gammact

def gamma(m_V, E_V):
    gamma = E_V/m_V
    return gamma

outfile = r.TFile("testout.root","RECREATE")

#Replicate Figure 2
ratios = np.linspace(1,12,120)
brs = []
for ratio in ratios:
    mpi = 33.3
    mv = mpi*1.8
    f_pi = mpi/ratio
    mAp = 3.*mpi

    brs.append(br_Vpi(mAp,mpi,mv,1.0,f_pi,True,True))

vpi_gr = r.TGraph(len(ratios),np.array(ratios),np.array(brs))
vpi_gr.SetName("brs")
vpi_gr.Draw()
vpi_gr.Write()

alpha_D = 1.e-2
eps = 1.e-3
m_l = 0.511
rho = False
phi = True

vd_masses = np.linspace(10.,210.,num=2000)
ctaus = []
for mass in vd_masses:
    m_Ap = mass*(3./1.8)
    m_pi = mass/1.8
    f_pi = m_pi/3.
    print("mAP: %f | m_v: %f | mpi: %f | fpi: %f"%(m_Ap, mass,m_pi, f_pi))
    ctau = getCtau(m_Ap,m_pi,mass,eps,alpha_D,f_pi,m_l,rho)
    ctaus.append(ctau)

ctau_gr = r.TGraph(len(vd_masses),np.array(vd_masses),np.array(ctaus))

outfile.cd()
ctau_gr.Draw()
ctau_gr.Write()

    
