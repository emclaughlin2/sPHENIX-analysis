#include <iostream>
#include <TH2D.h>
#include <TH1D.h>
#include <TChain.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <sstream> //std::ostringstsream
#include <fstream> //std::ifstream
#include <iostream> //std::cout, std::endl
#include <cmath>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <string.h>
#include <set>
#include <TVector3.h>
#include <map>
#include <vector>
#include <TDatabasePDG.h>
#include <tuple>
#include <TProfile.h>
#include <TProfile2D.h>
#include "TH1D.h"
#include <TH2D.h>

using namespace std;

void hadron_shower_profile() {

	// create outfile 
	//string outfilepath = "hadron_shower_profile.root";
	string outfilepath = "cluster_test.root";
	TFile *out = new TFile(outfilepath.c_str(),"RECREATE");
	
	// define constants 
	float deltaeta = 0.0916667;
	float deltaphi = 0.0981748;

	TH1F* h_vz = new TH1F("h_vz","",400, -100, 100);
    TH2D* h_truth_particles_2D = new TH2D("h_truth_particles_2D","",50,-1.5,1.5,50,-M_PI,M_PI);
    TH2D* h_clusters_2D = new TH2D("h_clusters_2D","",50,-1.5,1.5,50,-M_PI,M_PI);
    TH2D* h_shower_profile = new TH2D("h_shower_profile","",100,-3.5,3.5,100,-3.5,3.5);
    TH2D* h_emcal_shower_profile = new TH2D("h_emcal_shower_profile","",100,-3.5,3.5,100,-3.5,3.5);
    TH2D* h_ihcal_shower_profile = new TH2D("h_ihcal_shower_profile","",100,-3.5,3.5,100,-3.5,3.5);
    TH2D* h_ohcal_shower_profile = new TH2D("h_ohcal_shower_profile","",100,-3.5,3.5,100,-3.5,3.5);
    TH1D* h_ep = new TH1D("h_ep","",2000,-2,2);

  	// input datasets from ttrees
 	TChain chain("T");
 	string wildcardPath = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run21_jet10_3sigma_output_0.root";
	chain.Add(wildcardPath.c_str());

	int clsmult = 0;
	float cluster_e[2000] = {0.0};
	float cluster_eta[2000] = {0.0};
	float cluster_phi[2000] = {0.0};
	int cluster_ntowers[2000];
    int cluster_tower_calo[200][500];
    int cluster_tower_ieta[200][500];
    int cluster_tower_iphi[200][500];
    float cluster_tower_e[200][500];

    chain.SetBranchStatus("*",0);

	chain.SetBranchStatus("clsmult", 1);
	chain.SetBranchStatus("cluster_e", 1);
	chain.SetBranchStatus("cluster_eta", 1);
	chain.SetBranchStatus("cluster_phi", 1);
	chain.SetBranchStatus("cluster_ntowers", 1);
	chain.SetBranchStatus("cluster_tower_e", 1);
	chain.SetBranchStatus("cluster_tower_calo", 1);
	chain.SetBranchStatus("cluster_tower_ieta", 1);
	chain.SetBranchStatus("cluster_tower_iphi", 1);

	chain.SetBranchAddress("clsmult",&clsmult);
	chain.SetBranchAddress("cluster_e",cluster_e);
	chain.SetBranchAddress("cluster_eta",cluster_eta);
	chain.SetBranchAddress("cluster_phi",cluster_phi);
	chain.SetBranchAddress("cluster_ntowers",cluster_ntowers);
	chain.SetBranchAddress("cluster_tower_e",cluster_tower_e);
	chain.SetBranchAddress("cluster_tower_calo",cluster_tower_calo);
	chain.SetBranchAddress("cluster_tower_ieta",cluster_tower_ieta);
	chain.SetBranchAddress("cluster_tower_iphi",cluster_tower_iphi);

	
    float zvtx; chain.SetBranchStatus("zvtx", 1); chain.SetBranchAddress("zvtx", &zvtx);
	int truthpar_n = 0; chain.SetBranchStatus("truthpar_n", 1); chain.SetBranchAddress("truthpar_n",&truthpar_n);
    float truthpar_e[100000] = {0.0}; chain.SetBranchStatus("truthpar_e", 1); chain.SetBranchAddress("truthpar_e",truthpar_e);
    float truthpar_eta[100000] = {0.0}; chain.SetBranchStatus("truthpar_eta", 1); chain.SetBranchAddress("truthpar_eta",truthpar_eta);
    float truthpar_phi[100000] = {0.0}; chain.SetBranchStatus("truthpar_phi", 1); chain.SetBranchAddress("truthpar_phi",truthpar_phi);
    int truthpar_pid[100000] = {0}; chain.SetBranchStatus("truthpar_pid", 1); chain.SetBranchAddress("truthpar_pid",truthpar_pid);

    std::vector<float>* unsubjet_pt = nullptr; chain.SetBranchStatus("pt", 1); chain.SetBranchAddress("pt", &unsubjet_pt);
    std::vector<float>* unsubjet_eta = nullptr; chain.SetBranchStatus("eta", 1); chain.SetBranchAddress("eta", &unsubjet_eta);
    std::vector<float>* unsubjet_phi = nullptr; chain.SetBranchStatus("phi", 1); chain.SetBranchAddress("phi", &unsubjet_phi);

	int eventnumber = 0; // number of events in tree (used for showing iteration through ttree)
	int events = 0; // number of events that pass event cuts (used for event level normalization)
    
    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;

    // main event loop 
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 1000000; ++entry) {
        chain.GetEntry(entry);
    	if (eventnumber % 10000 == 0) cout << "event " << eventnumber << endl;
    	eventnumber++;

    	if (isnan(zvtx)) { continue; }
  		if (zvtx < -30 || zvtx > 30) { continue; }
    	h_vz->Fill(zvtx);

    	TVector3 truth, truth_test, cluster;
    	bool isolated = true;

    	
    	for (int i = 0; i < truthpar_n; i++) {
    		if (truthpar_e[i] < 1.0 || fabs(truthpar_eta[i]) > 1.1 || fabs(truthpar_pid[i]) != 211) { continue; }
    		isolated = true;
    		truth.SetPtEtaPhi(truthpar_e[i]/cosh(truthpar_eta[i]),truthpar_eta[i],truthpar_phi[i]);
    		//std::cout << "test truth particle: pt " << truth.Pt() << " eta " << truth.Eta() << " phi " << truth.Phi() << std::endl;
    		for (int j = 0; j < truthpar_n; j++) {
    			if (j == i || truthpar_e[j] < 0.2 || fabs(truthpar_eta[j]) > 1.1) { continue; }
    			truth_test.SetPtEtaPhi(truthpar_e[j]/cosh(truthpar_eta[j]),truthpar_eta[j],truthpar_phi[j]);
    			//std::cout << "test truth particle2: pt " << truth_test.Pt() << " eta " << truth_test.Eta() << " phi " << truth_test.Phi() << " dR " << truth.DeltaR(truth_test) << std::endl;
    			if (truth.DeltaR(truth_test) < 0.4) {
    				//std::cout << "not isolated dR = " << truth.DeltaR(truth_test) << std::endl;
    				isolated = false;
    				break;
    			}
    		}
    		if (isolated) {
    			//std::cout << "isolated particle: pt " << truth.Pt() << " eta " << truth.Eta() << " phi " << truth.Phi() << std::endl;
    			h_truth_particles_2D->Fill(truth.Eta(),truth.Phi(),truth.Pt());
    			int ic = -1;
    			float match_dR = 9999;
    			for (int j = 0; j < clsmult; j++) {
    				cluster.SetPtEtaPhi(cluster_e[j]/cosh(cluster_eta[j]),cluster_eta[j],cluster_phi[j]);
    				//std::cout << "test cluster: pt " << cluster.Pt() << " eta " << cluster.Eta() << " phi " << cluster.Phi() << " dR " << truth.DeltaR(cluster) << std::endl;
    				if (truth.DeltaR(cluster) < 0.2 && truth.DeltaR(cluster) < match_dR) {
    					ic = j;
    					match_dR = truth.DeltaR(cluster);
    					//std::cout << "cluster match index " << ic << " match dR " << match_dR << " e " << cluster_e[j] << " eta " << cluster_eta[j] << " phi " << cluster_phi[j] << std::endl;
    				}
    			}
    			if (ic >= 0) {
    				h_clusters_2D->Fill(cluster_eta[ic], cluster_phi[ic], cluster_e[ic]/cosh(cluster_eta[ic]));
    				//std::cout << "cluster match index " << ic << " e " << cluster_e[ic] << " eta " << cluster_eta[ic] << " phi " << cluster_phi[ic] << std::endl;
				    
				    float eta_cm = 0;
				    float phi_cm = 0;
				    float e_tot = 0;
	    			for (int j = 0; j < cluster_ntowers[ic]; j++) {
	    				if (cluster_tower_calo[ic][j] == 1) {
	    					eta_cm += cluster_tower_e[ic][j] * cluster_tower_ieta[ic][j]/4.0;
		    				phi_cm += cluster_tower_e[ic][j] * cluster_tower_iphi[ic][j]/4.0;
		    				e_tot += cluster_tower_e[ic][j];
	    				} else {
		    				eta_cm += cluster_tower_e[ic][j] * cluster_tower_ieta[ic][j];
		    				phi_cm += cluster_tower_e[ic][j] * cluster_tower_iphi[ic][j];
		    				e_tot += cluster_tower_e[ic][j];
		    			}
	    			}
	    			eta_cm /= e_tot;
		    		phi_cm /= e_tot;
		    		//std::cout << "cluster ntowers " << cluster_ntowers[ic] << " eta " << cluster_eta[ic] << " phi " << cluster_phi[ic] << " eta_cm*deltaeta " << eta_cm*deltaeta - 1.1 << " phi_cm*deltaphi " << phi_cm*deltaphi << " eta_cm " << eta_cm << " phi_cm " << phi_cm << std::endl;
		    		h_ep->Fill(cluster_e[ic]/(truth.Pt()*cosh(cluster_eta[ic])));
		    		for (int j = 0; j < cluster_ntowers[ic]; j++) {
	    				if (cluster_tower_calo[ic][j] == 1) {
	    					//std::cout << "cluster emcal tower ieta " << cluster_tower_ieta[ic][j]/4.0 - eta_cm << " iphi " << cluster_tower_iphi[ic][j]/4.0 - phi_cm << " e " << cluster_tower_e[ic][j] << std::endl;
	    					h_shower_profile->Fill(cluster_tower_ieta[ic][j]/4.0 - eta_cm, cluster_tower_iphi[ic][j]/4.0 - phi_cm, cluster_tower_e[ic][j]/(truth.Pt()*cosh(cluster_eta[ic])));
	    					h_emcal_shower_profile->Fill(cluster_tower_ieta[ic][j]/4.0 - eta_cm, cluster_tower_iphi[ic][j]/4.0 - phi_cm, cluster_tower_e[ic][j]/(truth.Pt()*cosh(cluster_eta[ic])));
	    				} else {
		    				h_shower_profile->Fill(cluster_tower_ieta[ic][j] - eta_cm, cluster_tower_iphi[ic][j] - phi_cm, cluster_tower_e[ic][j]/(truth.Pt()*cosh(cluster_eta[ic])));
		    				if (cluster_tower_calo[ic][j] == 2) {
		    					//std::cout << "cluster ohcal tower ieta " << cluster_tower_ieta[ic][j] - eta_cm << " iphi " << cluster_tower_iphi[ic][j] - phi_cm << " e " << cluster_tower_e[ic][j] << std::endl;
		    					h_ohcal_shower_profile->Fill(cluster_tower_ieta[ic][j] - eta_cm, cluster_tower_iphi[ic][j] - phi_cm, cluster_tower_e[ic][j]/(truth.Pt()*cosh(cluster_eta[ic])));
		    				} else {
		    					//std::cout << "cluster ihcal tower ieta " << cluster_tower_ieta[ic][j] - eta_cm << " iphi " << cluster_tower_iphi[ic][j] - phi_cm << " e " << cluster_tower_e[ic][j] << std::endl;
		    					h_ihcal_shower_profile->Fill(cluster_tower_ieta[ic][j] - eta_cm, cluster_tower_iphi[ic][j] - phi_cm, cluster_tower_e[ic][j]/(truth.Pt()*cosh(cluster_eta[ic])));
		    				}
		    			}
	    			}
	    			events++;
    			}
    		}
    	}
  	}

  	out->Write();
  	out->Close();

}