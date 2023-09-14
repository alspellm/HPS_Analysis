#/bin/env python3
import ROOT as r
import numpy as np

outfile = r.TFile("plotZcutOverBackground.root","RECREATE")

ttfile = r.TFile("/sdf/group/hps/users/alspellm/projects/simps_2016/reach_estimate/mc/tritrig_beam/SR_ana/hadd_tritrigv2-beamv6_2500kBunches_HPS-PhysicsRun2016-Pass2_v4_5_0_pairs1_SR_ana.root","READ")
ttrig_hh = ttfile.Get("vtxana_Tight_simpSIG/vtxana_Tight_simpSIG_vtx_InvM_vtx_z_hh")

zcutfile = r.TFile("/sdf/group/hps/users/alspellm/projects/simps_2016/reach_estimate/mc/scripts/zcuts.root","READ")
zcut_g = zcutfile.Get("zcuts_g")
y = []
x= []
oldx = []
oldy = []
for i in range(zcut_g.GetN()):
    y.append(zcut_g.GetY()[i])
    x.append(zcut_g.GetX()[i]*0.001)

    mass = zcut_g.GetX()[i]*0.001
    zCut = 17.7702 + 138.166*mass - 5363.29*mass*mass + 44532.4*pow(mass,3) - 120578*pow(mass,4)
    oldx.append(mass)
    oldy.append(zCut)
        
zcutgev_fit_h = r.TGraph(len(y),np.array(x),np.array(y))
zcutgev_fit_h.SetName("zcut_fit_from_mc")

zcutgev_old_h = r.TGraph(len(oldy),np.array(oldx),np.array(oldy))
zcutgev_old_h.SetName("2016_displaced_A'_zcut")

outfile.cd()
ttrig_hh.Write()
zcut_g.Write()
zcutgev_fit_h.Write()
zcutgev_old_h.Write()

c = r.TCanvas("canv","canv",1800,1000)
c.cd()
ttrig_hh.Draw("colz")
#zcut_g.Draw("same")
zcutgev_fit_h.Draw("same")

c.Write()






