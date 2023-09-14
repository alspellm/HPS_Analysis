#!/usr/bin/python3
import ROOT as r

outfile = r.TFile("outfile.root","RECREATE")

csel_3ns = {}
csel_3125ns = {}

with open('/sdf/group/hps/users/alspellm/projects/fit_shape_params/2019/corrected/pulse_shape_fits/comparisons/compare_to_25ns_delay/hpssvt_010353_database_params_3ns_cut_apv127.dat') as f1:
    lines = f1.readlines()
    for line in lines:
        vals = line.split(",")
        svtid = vals[0]
        tau1 = vals[3]
        tau2 = vals[4].strip()
        #print(svtid, "", tau1," ", tau2)
        taus = [float(tau1),float(tau2)]
        csel_3ns[float(svtid)]=taus

with open('/sdf/group/hps/users/alspellm/projects/fit_shape_params/2019/corrected/pulse_shape_fits/comparisons/compare_to_25ns_delay/hpssvt_010353_pulsefits_3125ns_cutapv127.dat') as f2:
    lines = f2.readlines()
    for line in lines:
        vals = line.split(",")
        svtid = vals[0]
        tau1 = vals[3]
        tau2 = vals[4].strip()
        #print(svtid, "", tau1," ", tau2)
        taus = [float(tau1),float(tau2)]
        csel_3125ns[float(svtid)]=taus

tau1_v_tau2_hh = r.TH2F("ratio_3ns_to_3125ns","ratio_3ns_to_3125ns;tau1_ratio;tau2_ratio",2000,0,2,2000,0,2)

for svtid in csel_3ns:
    tau1a = csel_3ns[svtid][0]
    tau2a = csel_3ns[svtid][1]

    tau1b = csel_3125ns[svtid][0]
    tau2b = csel_3125ns[svtid][1]

    tau1ratio = tau1a/tau1b
    tau2ratio = tau2a/tau2b

    tau1_v_tau2_hh.Fill(tau1ratio,tau2ratio)

outfile.cd()
tau1_v_tau2_hh.Write()

