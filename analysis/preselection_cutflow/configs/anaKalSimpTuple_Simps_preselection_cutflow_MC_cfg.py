import HpstrConf
import sys
import os
import baseConfig as base

base.parser.add_argument("-w", "--tracking", type=str, dest="tracking",
                         help="Which tracking to use to make plots", metavar="tracking", default="KF")
base.parser.add_argument("-f", "--makeFlatTuple", type=int, dest="makeFlatTuple",
                         help="Make True to make vertex ana flat tuple", metavar="makeFlatTuple", default=0)
base.parser.add_argument("-r", "--isRadPDG", type=int, dest="isRadPDG",
                         help="Set radiative trident PDG ID", metavar="isRadPDG", default=625)
base.parser.add_argument("-TS", "--trackstate", type=str, dest="trackstate",
                         help="Specify Track State | 'AtECal' or 'AtTarget'. Default is origin (AtIP)", metavar="trackstate", default="")
base.parser.add_argument("-bpc", "--beamPosCorr", type=str, dest="beamPosCorr",
                         help="Load beam position corrections from json", metavar="beamPosCorr", default="")
options = base.parser.parse_args()

# Use the input file to set the output file name
infile = options.inFilename
outfile = options.outFilename

outfile = outfile.split(".root")[0]+"_"+options.tracking+".root"

print('Input file: %s' % infile)
print('Output file: %s' % outfile)

p = HpstrConf.Process()

p.run_mode = 1
#p.max_events = 1000

# Library containing processors
p.add_library("libprocessors")

###############################
#          Processors         #
###############################

recoana_kf = HpstrConf.Processor('vtxana_kf', 'VertexAnaProcessor')
recoana_gbl = HpstrConf.Processor('vtxana_gbl', 'VertexAnaProcessor')
###############################
#   Processor Configuration   #
###############################
#RecoHitAna
recoana_kf.parameters["anaName"] = "vtxana_kf"
recoana_kf.parameters["trkColl"] = "KalmanFullTracks%s"%(options.trackstate)
recoana_kf.parameters["tsColl"] = "TSData"
recoana_kf.parameters["vtxColl"] = "UnconstrainedV0Vertices_KF"
recoana_kf.parameters["mcColl"] = "MCParticle"
recoana_kf.parameters["hitColl"] = "SiClusters"
recoana_kf.parameters["ecalColl"] = "RecoEcalClusters"
recoana_kf.parameters["vtxSelectionjson"] = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/configs/vertexSelection_2016_simp_nocuts.json'
recoana_kf.parameters["histoCfg"] = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/configs/vtxAnalysis_2016_simp_preselection.json'
recoana_kf.parameters["mcHistoCfg"] = os.environ['HPSTR_BASE']+'/analysis/plotconfigs/mc/basicMC.json'
#####
recoana_kf.parameters["beamE"] = base.beamE[str(options.year)]
recoana_kf.parameters["isData"] = options.isData
recoana_kf.parameters["analysis"] = options.analysis
recoana_kf.parameters["debug"] = 0
recoana_kf.parameters["isRadPDG"] = options.isRadPDG
recoana_kf.parameters["makeFlatTuple"] = options.makeFlatTuple
#recoana_kf.parameters["beamPosCfg"] = options.beamPosCorr
#recoana_kf.parameters["beamPosCfg"] = os.environ['HPSTR_BASE']+'/analysis/data/beamspot_positions_2016.json'
recoana_kf.parameters["eleTrackTimeBias"] = -2.2 #MC
recoana_kf.parameters["posTrackTimeBias"] = -2.1 #MC


CalTimeOffset = -999

if (options.isData == 1):
    CalTimeOffset = 56.
    print("Running on data file: Setting CalTimeOffset %d" % CalTimeOffset)

elif (options.isData == 0):
    CalTimeOffset = 43.
    print("Running on MC file: Setting CalTimeOffset %d" % CalTimeOffset)
else:
    print("Specify which type of ntuple you are running on: -t 1 [for Data] / -t 0 [for MC]")

recoana_kf.parameters["CalTimeOffset"] = CalTimeOffset
#Region definitions
RegionPath = '/sdf/group/hps/users/alspellm/projects/THESIS/analysis/preselection_cutflow/selections/cutflow/'

recoana_kf.parameters["regionDefinitions"] = [RegionPath+'posTrkTime_lt_selection_inclusive.json',
                                              RegionPath+'posTrkCluTimeDiff_lt_selection_inclusive.json',
                                              RegionPath+'posTrkChi2Ndf_lt_selection_inclusive.json',
                                              RegionPath+'posN2Dhits_gt_selection_inclusive.json',
                                              RegionPath+'posMom_gt_selection_inclusive.json',
                                              RegionPath+'Pair1_eq_selection_inclusive.json',
                                              RegionPath+'eleTrkTime_lt_selection_inclusive.json',
                                              RegionPath+'eleTrkCluTimeDiff_lt_selection_inclusive.json',
                                              RegionPath+'eleTrkChi2Ndf_lt_selection_inclusive.json',
                                              RegionPath+'eleposCluTimeDiff_lt_selection_inclusive.json',
                                              RegionPath+'eleN2Dhits_gt_selection_inclusive.json',
                                              RegionPath+'eleMom_lt_selection_inclusive.json',
                                              RegionPath+'eleMom_gt_selection_inclusive.json',
                                              RegionPath+'maxVtxMom_lt_selection_inclusive.json',
                                              RegionPath+'chi2unc_lt_selection_inclusive.json']


# Sequence which the processors will run.
#p.sequence = [recoana_kf,recoana_gbl]
if (options.tracking == "KF"):
    print("Run KalmanFullTracks analysis")
    p.sequence = [recoana_kf]  # ,mcana]

p.skip_events = options.skip_events
p.max_events = options.nevents

p.input_files = infile
p.output_files = [outfile]

p.printProcess()
