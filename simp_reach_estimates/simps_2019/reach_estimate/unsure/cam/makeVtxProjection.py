import HpstrConf
import sys
import baseConfig

parser = baseConfig.parser

baseConfig.parser.add_option("-w","--tracking", type="string", dest="tracking",
                             help="Which tracking to use to make plots", metavar="tracking", default="KF")
(options, args) = parser.parse_args()

# Use the input file to set the output file name
histo_file = options.inFilename
out_file   = options.outFilename

print('Histo file: %s' % histo_file)
print('Out file: %s' % out_file)

p = HpstrConf.Process()

p.run_mode = 2
#p.max_events = 1000

# Library containing processors
p.add_library("libprocessors")

###############################
#          Processors         #
###############################

vtxPostProc = HpstrConf.Processor('vtxPostProc', 'VtxHistoProcessor')

###############################
#   Processor Configuration   #
###############################
#MCParticles
vtxPostProc.parameters["debug"] = 1
vtxPostProc.parameters["rebin"] = 4
#KF
#vtxPostProc.parameters["selection"] = "vtxana_kf_vtxSelection"
#GBL
if (options.tracking == "KF"):
    vtxPostProc.parameters["selections"] = ["vtxana_Tight_e2pt3"]
elif (options.tracking == "GBL"):
    vtxPostProc.parameters["selections"] = ["vtxana_Tight_e2pt3"]
else:
    print("Error specify -w KF or -w GBL")
    exit(1)


vtxPostProc.parameters["projections"] = ["ele_d0_vs_p_hh","pos_d0_vs_p_hh",
                                         "ele_z0_vs_p_hh","pos_z0_vs_p_hh",
                                         "vtx_InvM_vtx_svt_z_hh","vtx_p_svt_z_hh"]


# Sequence which the processors will run.
p.sequence = [vtxPostProc]

p.input_files=[histo_file]
p.output_files = [out_file]

p.printProcess()
