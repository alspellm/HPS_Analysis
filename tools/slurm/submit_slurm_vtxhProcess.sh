#!/bin/bash

# Reset OPTIND to its initial value (1)
OPTIND=1
while getopts r: flag
do
    case "${flag}" in
        r) rundir=${OPTARG};;
    esac
done

echo "LOOK $rundir"
export Rundir=$(readlink -f "$rundir")
echo "$Rundir"


# List of input files
input_files=(
    "sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_100_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_105_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_110_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_115_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_120_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_125_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_130_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_135_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_140_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_145_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_150_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_155_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_160_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_165_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_170_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_175_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_180_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_185_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_190_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_195_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_200_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.roo"t
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_40_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_45_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_50_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_55_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_60_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_65_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_70_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_75_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_80_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_85_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_90_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"
"/sdf/home/a/alspellm/work/projects/THESIS/mc/2016/simps/signal_beam/20230713_slic/20230713_readout/20230713_recon/20230721_tuples/20230724_ana/ana_files/hadd_simp_95_MeV_egsv6_HPS-PhysicsRun2016-Pass2_recon_targetTracks_NObeamPosCorr_ana.root"

)
# Path to your Python script
hpstr_config="/sdf/group/hps/users/alspellm/projects/THESIS/ana/analysis_scripts/reach_estimate_tenpct/vtxhProcess.py"

# Loop through input files and submit jobs
for input_file in "${input_files[@]}"; do
    # Extract filename without path and extension
    filename=$(basename "$input_file")
    filename_without_extension="${filename%.*}"
    job_name=${filename_without_extension}
    echo "JOB NAME is $job_name"

    # SLURM job submission command
    sbatch <<EOT
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=00:30:00  # Adjust the time limit as needed
#SBATCH --mem=4000M
#SBATCH --job-name=${Rundir}/${job_name}
#SBATCH --output=${Rundir}/${job_name}.out
#SBATCH --error=${Rundir}/${job_name}.err
#SBATCH --partition=hps

echo "Running on input file ${input_file}"
source /sdf/home/a/alspellm/.bashrc

hpstr ${hpstr_config} -i ${input_file} 

EOT
done
