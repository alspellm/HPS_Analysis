#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TH2D.h"
#include "TH1D.h"


void fsp_macro(std::string inputFileName, std::string outputFileName) {
    
    std::string outfilename = outputFileName;
    TFile* outfile = new TFile(outfilename.c_str(), "RECREATE");

    std::string filename = "/sdf/group/hps/users/alspellm/projects/THESIS/ana/ecal_trig_res/20230628/hadd_tritrig-beam_fsp_2kfiles.root";
    TFile* file = new TFile(inputFileName.c_str(), "READ");
    TTree* tree = (TTree*)file->Get("HPS_Event");

    std::map<std::string,TH1D*> histos1d;
    //Plots 1D
    TH1D* clustE_h = new TH1D("fee_cluster_energy","fee_cluster_energy;Energy [GeV];Entries",250,0.0,5.0);
    TH1D* fidclustE_h = new TH1D("fiducial_fee_cluster_energy","fiducial_fee_cluster_energy;Energy [GeV];Entries",250,0.0,5.0);
    histos1d["fee_cluster_energy"] = clustE_h;
    histos1d["fiducial_fee_cluster_energy"] = fidclustE_h;

    std::map<std::string,TH2D*> histos2d;
    //Plots 2D
    TH2D* clust_pos = new TH2D("cluster_positions", "cluster_positions;x[mm];y[mm]",700,-300.0,400.0,170,-85.0, 85.0);
    TH2D* fidclust_pos = new TH2D("fiducial_cluster_positions", "fiducial_cluster_positions;x[mm];y[mm]",700,-300.0,400.0,170.0,-85.0, 85.0);
    TH2D* track_pos = new TH2D("track_positions", "track_positions;x[mm];y[mm]",700,-300.0,400.0,170,-85.0, 85.0);
    TH2D* fid_track_pos = new TH2D("fiducial_track_positions", "fiducial_track_positions;x[mm];y[mm]",700,-300.0,400.0,170.0,-85.0, 85.0);
    histos2d["cluster_positions"] = clust_pos;
    histos2d["fiducial_cluster_positions"] = fidclust_pos;
    histos2d["track_positions"] = track_pos;
    histos2d["fiducial_track_positions"] = fid_track_pos;


    std::vector<Particle*>* particles = new std::vector<Particle*>();
    tree->SetBranchAddress("FinalStateParticles",&particles);

    //ECal geometry (approx)
    //Define first order Ecal geometry --> assumed square here,
    //smaller than actual beamgap. Improve x geometry 
    double ecalx1 = -270.0; //mm
    double ecalx2 = 355.0;
    double ecaly1 = 85.0;
    double ecaly2 = -85.0;
    double bgapup = 29.0;
    double bgapdown = -29.0;
    double eholex1 = -87.5;
    double eholex2 = 12.5;
    double eholey1 = 42.0;
    double eholey2 = -42.0;
    double crystal_width = 13.0;


    int n = tree->GetEntries();
    for (int i = 0; i < n; i++){
        if(i%100000 == 0)
            std::cout << "Event " << i << " out of " << n << std::endl;
        tree->GetEntry(i);
        for( int pi = 0; pi < particles->size(); pi++){
            Particle* particle = particles->at(pi);    
            CalCluster cluster = particle->getCluster();
            double clusterE = cluster.getEnergy();
            double cluster_x = cluster.getPosition().at(0);
            double cluster_y = cluster.getPosition().at(1);
            double cluster_z = cluster.getPosition().at(2);

            Track track = particle->getTrack();
            int nhits = track.getTrackerHitCount();
            int charge = track.getCharge();
            double track_x = track.getPositionAtEcal().at(0);
            double track_y = track.getPositionAtEcal().at(1);
            if(charge > 0)
                continue;
            if(nhits < 12)
                continue;

            double p = track.getP();
            if(p < 2.1)
                continue;

            if(clusterE/p < 0.9)
                continue;

           histos1d["fee_cluster_energy"]->Fill(clusterE);
           histos2d["cluster_positions"]->Fill(cluster_x, cluster_y);
           histos2d["track_positions"]->Fill(track_x, track_y);

           //Fiducial region
            if(cluster_x < ecalx1 + crystal_width  || cluster_x > ecalx2 - crystal_width)
                continue;
            if(cluster_y > ecaly1 - crystal_width || cluster_y < ecaly2 + crystal_width)
                continue;
            if(cluster_y < bgapup + crystal_width && cluster_y > bgapdown - crystal_width)
                continue;
            if((cluster_x > eholex1 - crystal_width && cluster_x < eholex2 + crystal_width) && ( (cluster_y < eholey1 + crystal_width) && (cluster_y > eholey2 - crystal_width)))
                continue;

           //Fiducial region
            if(track_x < ecalx1 + crystal_width  || track_x > ecalx2 - crystal_width)
                continue;
            if(track_y > ecaly1 - crystal_width || track_y < ecaly2 + crystal_width)
                continue;
            if(track_y < bgapup + crystal_width && track_y > bgapdown - crystal_width)
                continue;
            if((track_x > eholex1 - crystal_width && track_x < eholex2 + crystal_width) && ( (track_y < eholey1 + crystal_width) && (track_y > eholey2 - crystal_width)))
                continue;

           histos1d["fiducial_fee_cluster_energy"]->Fill(clusterE);
           histos2d["fiducial_cluster_positions"]->Fill(cluster_x, cluster_y);
           histos2d["fiducial_track_positions"]->Fill(track_x, track_y);
        }
    }

    outfile->cd();

    //1d plots
    for(std::map<std::string, TH1D*>::iterator it = histos1d.begin(); it != histos1d.end(); it++){
        if(!it->second)
            continue;
        std::cout << "Writing " << it->first << std::endl;
        it->second->Draw();
        it->second->Write();
    }
    //2d plots
    for(std::map<std::string, TH2D*>::iterator it = histos2d.begin(); it != histos2d.end(); it++){
        if(!it->second)
            continue;
        std::cout << "Writing " << it->first << std::endl;
        it->second->Draw("colz");
        it->second->Write();
    }
}


int main(int argc, char* argv[]){
    std::string inputFileName = argv[1];
    std::string outputFileName = argv[2];
    fsp_macro(inputFileName, outputFileName);
    return 0;
}
