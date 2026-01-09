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

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//																																			//
//	Note: 11.26.24 																															//
// 	This macro creates dijet QA plots, energy density QA plots with options for using towers, topoclusters or emcal clusters,				//
//	and plots of mean energy density vs jet variables (leading jet pT and dijet xj)															//
//																																			//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<int> run_numbers;
std::vector<double> collision_rates;
std::vector<double> pileup_rates;
std::vector<int> low_pu_runs;

void get_pileup_rate() {

    std::ifstream infile("pileup/0mrad_collision_rates.txt");
    if (!infile) {
        std::cerr << "Error: Unable to open input file!" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        int run_number;
        double collision_rate;  
        if (!(iss >> run_number >> collision_rate)) {
            std::cerr << "Error: Malformed line in input file!" << std::endl;
            continue;
        }
        std::cout << "run_number " << run_number << " collision_rate " << collision_rate << std::endl;
        run_numbers.push_back(run_number);
        collision_rates.push_back(collision_rate);

        // Calculate lambda
        double denominator = 78000.0 * 111.0;
        double lambda = collision_rate / denominator;
        double exp_neg_lambda = exp(-lambda);
        
        // Calculate pk
        double pk = 1 - lambda * exp_neg_lambda - exp_neg_lambda;
        
        // Calculate pileup rate
        double pileup_rate = pk / (pk + lambda * exp_neg_lambda);
        pileup_rates.push_back(pileup_rate);
    }
    infile.close();

    for (int i = 0; i < pileup_rates.size(); i++) {
    	if (pileup_rates[i] < 0.55) { low_pu_runs.push_back(run_numbers[i]); }
    }
    
}

void minbias_analysis(int runnumber = 51274, string outfilename = "dijet_calo_analysis_run21_jet10_topocluster.root", bool sim = false, bool clusters = true, const char* type = "1.5mrad") {

	// create outfile 
	std::cout << "Creating " << outfilename << std::endl;
	string outfilepath = "minbias_analysis_output/" + outfilename;
	TFile *out = new TFile(outfilepath.c_str(),"RECREATE");
	
	// define constants 
	float deltaeta = 0.0916667;
	float deltaphi = 0.0981748;
	float secteta = 2.2;
	float sectphi = 2.0*M_PI;
	float ptbins[] = {5,6,7,8,9,10,12,14,16,18,20,22,24,27,30,35,40,50,60,80};
    int nptbins = sizeof(ptbins) / sizeof(ptbins[0]) - 1;
    float etmin = -20;
    float etmax = 80;
    int netbins = int((10*(etmax - etmin)));
    float etbins[netbins+1];
    for (int i = 0; i <= netbins; i++) {
    	etbins[i] = etmin + i*(etmax - etmin) / netbins;
    }
 
  	// create UE histograms 
  	TH1F* h_vz = new TH1F("h_vz","",400, -100, 100);
  	TH2F* h_ue_2D_total = new TH2F("h_ue_2D_total","",24,-1.1,1.1,32,0,M_PI);
  	TH1F* h_ue_total = new TH1F("h_ue_total","",netbins, etbins);
  	TH1F* h_et_total_pos = new TH1F("h_et_total_pos","",800,0,80);
  	TH1F* h_et_total = new TH1F("h_et_total","",netbins, etbins);

  	// create topocluster histograms 
  	int topo_thresholds[] = {-9999,0,100,200,300,500,1000,2000};

  	TH1F* h_ntopo_total[8];
  	TH1F* h_topo_total[8];
  	TH2F* h_2D_topo_total[8];
  	TH1F* h_sume_topo_total[8];

  	for (int i = 0; i < 8; i++) {
  		h_ntopo_total[i] = new TH1F(Form("h_ntopo%d_total",topo_thresholds[i]),"",200,0,200);
  		h_topo_total[i] = new TH1F(Form("h_topo%d_total",topo_thresholds[i]),"",netbins, etbins);
  		h_2D_topo_total[i] = new TH2F(Form("h_2D_topo%d_total",topo_thresholds[i]),"",24,-1.1,1.1,32,0,M_PI);
  		h_sume_topo_total[i] = new TH1F(Form("h_sume_topo%d_total",topo_thresholds[i]),"",netbins, etbins);
  	}

  	// input datasets from ttrees
 	TChain chain("T");
 	string inputDirectory = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/";
 	string infilename;
 	if (sim) {
 		infilename = "topo_valid_3sigma_sim_run22_mb_output_*.root";
 	} else if (!strcmp(type,"1.5mrad")) {
 		infilename = "output_1.5mrad_ana468_calofitting_mb_events_" + to_string(runnumber) + "_*.root";
 	} else if (!strcmp(type,"0mrad")) {
 		infilename = "output_0mrad_ana468_" + to_string(runnumber) + "_*.root";
 	} else {
 		std::cout << "unknown type" << std::endl;
 		return;
 	}
	string wildcardPath = inputDirectory + infilename;
	chain.Add(wildcardPath.c_str());

	// define ttree branch variables 
	int m_event;
	int nJet;
	float zvtx;

	vector<int> *triggerVector = nullptr;
	vector<float> *eta = nullptr;
	vector<float> *phi = nullptr;
	vector<float> *e = nullptr;
	vector<float> *pt = nullptr;

	int emcaln = 0;
	float emcale[24576] = {0.0};
	float emcaleta[24576] = {0.0};
	float emcalphi[24576] = {0.0};
	int emcalieta[24576] = {0};
	int emcaliphi[24576] = {0};
	int emcalstatus[24576] = {0};

	int ihcaln = 0;
	float ihcale[1536] = {0.0};
	float ihcaleta[1536] = {0.0};
	float ihcalphi[1536] = {0.0};
	int ihcalieta[1536] = {0};
	int ihcaliphi[1536] = {0};
	int ihcalstatus[1536] = {0};

	int ohcaln = 0;
	float ohcale[1536] = {0.0};
	float ohcaleta[1536] = {0.0};
	float ohcalphi[1536] = {0.0};
	int ohcalieta[1536] = {0};
	int ohcaliphi[1536] = {0};
	int ohcalstatus[1536] = {0};

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

	chain.SetBranchStatus("m_event",1);
	chain.SetBranchStatus("nJet",1);
	chain.SetBranchStatus("zvtx",1);
	chain.SetBranchStatus("triggerVector",1);
	chain.SetBranchStatus("eta",1);
	chain.SetBranchStatus("phi",1);
	chain.SetBranchStatus("e",1);
	chain.SetBranchStatus("pt",1);

	chain.SetBranchAddress("m_event",&m_event);
	chain.SetBranchAddress("nJet",&nJet);
	chain.SetBranchAddress("zvtx",&zvtx);
	chain.SetBranchAddress("triggerVector",&triggerVector);
	chain.SetBranchAddress("eta",&eta);
	chain.SetBranchAddress("phi",&phi);
	chain.SetBranchAddress("e",&e);
	chain.SetBranchAddress("pt",&pt);

	if (!clusters) {
		chain.SetBranchStatus("emcaln", 1);
		chain.SetBranchStatus("emcale",1);
		chain.SetBranchStatus("emcaleta", 1);
		chain.SetBranchStatus("emcalphi", 1);
		chain.SetBranchStatus("emcalieta", 1);
		chain.SetBranchStatus("emcaliphi", 1);

		chain.SetBranchStatus("ihcaln", 1);
		chain.SetBranchStatus("ihcale", 1);
		chain.SetBranchStatus("ihcaleta", 1);
		chain.SetBranchStatus("ihcalphi", 1);
		chain.SetBranchStatus("ihcalieta", 1);
		chain.SetBranchStatus("ihcaliphi", 1);

		chain.SetBranchStatus("ohcaln", 1);
		chain.SetBranchStatus("ohcale", 1);
		chain.SetBranchStatus("ohcaleta", 1);
		chain.SetBranchStatus("ohcalphi", 1);
		chain.SetBranchStatus("ohcalieta", 1);
		chain.SetBranchStatus("ohcaliphi", 1);

		chain.SetBranchAddress("emcaln",&emcaln);
		chain.SetBranchAddress("emcale",emcale);
		chain.SetBranchAddress("emcaleta",emcaleta);
		chain.SetBranchAddress("emcalphi",emcalphi);
		chain.SetBranchAddress("emcalieta",emcalieta);
		chain.SetBranchAddress("emcaliphi",emcaliphi);
		chain.SetBranchAddress("emcalstatus",emcalstatus);

		chain.SetBranchAddress("ihcaln",&ihcaln);
		chain.SetBranchAddress("ihcale",ihcale);
		chain.SetBranchAddress("ihcaleta",ihcaleta);
		chain.SetBranchAddress("ihcalphi",ihcalphi);
		chain.SetBranchAddress("ihcalieta",ihcalieta);
		chain.SetBranchAddress("ihcaliphi",ihcaliphi);
		chain.SetBranchAddress("ihcalstatus",ihcalstatus);

		chain.SetBranchAddress("ohcaln",&ohcaln);
		chain.SetBranchAddress("ohcale",ohcale);
		chain.SetBranchAddress("ohcaleta",ohcaleta);
		chain.SetBranchAddress("ohcalphi",ohcalphi);
		chain.SetBranchAddress("ohcalieta",ohcalieta);
		chain.SetBranchAddress("ohcaliphi",ohcaliphi);
		chain.SetBranchAddress("ohcalstatus",ohcalstatus);
	} else {
		chain.SetBranchStatus("clsmult", 1);
		chain.SetBranchStatus("cluster_e", 1);
		chain.SetBranchStatus("cluster_eta", 1);
		chain.SetBranchStatus("cluster_phi", 1);
		if (!sim) {
			chain.SetBranchStatus("cluster_ntowers", 1);
			chain.SetBranchStatus("cluster_tower_e", 1);
			chain.SetBranchStatus("cluster_tower_calo", 1);
			chain.SetBranchStatus("cluster_tower_ieta", 1);
			chain.SetBranchStatus("cluster_tower_iphi", 1);
		}

		chain.SetBranchAddress("clsmult",&clsmult);
		chain.SetBranchAddress("cluster_e",cluster_e);
		chain.SetBranchAddress("cluster_eta",cluster_eta);
		chain.SetBranchAddress("cluster_phi",cluster_phi);
		if (!sim) {
			chain.SetBranchAddress("cluster_ntowers",cluster_ntowers);
			chain.SetBranchAddress("cluster_tower_e",cluster_tower_e);
			chain.SetBranchAddress("cluster_tower_calo",cluster_tower_calo);
			chain.SetBranchAddress("cluster_tower_ieta",cluster_tower_ieta);
			chain.SetBranchAddress("cluster_tower_iphi",cluster_tower_iphi);
		}
	}

	int eventnumber = 0; // number of events in tree (used for showing iteration through ttree)
	int events = 0; // number of events that pass event cuts (used for event level normalization)
    
    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;

    // main event loop 
    //for (Long64_t entry = 0; entry < nEntries; ++entry) {
    for (Long64_t entry = 0; entry < 1000000; ++entry) {
        chain.GetEntry(entry);
    	if (eventnumber % 10000 == 0) cout << "event " << eventnumber << endl;
    	eventnumber++;

    	// require jet trigger in data
  		bool jettrig = false;
  		for (auto t : *triggerVector) {
  			if (t == 12 || t == 10) {
  				jettrig = true;
  				break;
  			}
  		}

  		// require at least 2 jets in event and z vertex < 30 cm 
  		if (!jettrig && !sim) { continue; }
  		if (isnan(zvtx)) { continue; }
  		if (zvtx < -10 || zvtx > 10) { continue; }

    	h_vz->Fill(zvtx);
		float et_total = 0;
		int ntopo_total[] = {0,0,0,0,0,0,0,0};
		float sume_topo_total[] = {0,0,0,0,0,0,0,0};

		if (!clusters) { // using towers to find total energy in towards, transverse and away regions 
			for (int i = 0; i < emcaln; i++) {
				if (emcalstatus[i] & 1 || emcalstatus[i] & 4) { continue; }
				et_total += emcale[i]/cosh(emcaleta[i]);
				h_ue_2D_total->Fill(emcaleta[i],emcalphi[i],emcale[i]/cosh(emcaleta[i]));
			}
			for (int i = 0; i < ihcaln; i++) {
				if (ihcalstatus[i] & 1 || ihcalstatus[i] & 4) { continue; }
				et_total += ihcale[i]/cosh(ihcaleta[i]);
				h_ue_2D_total->Fill(ihcaleta[i],ihcalphi[i],ihcale[i]/cosh(ihcaleta[i]));
			}
			for (int i = 0; i < ohcaln; i++) {
				if (ohcalstatus[i] & 1 || ohcalstatus[i] & 4) { continue; }
				et_total += ohcale[i]/cosh(ohcaleta[i]);
				h_ue_2D_total->Fill(ohcaleta[i],ohcalphi[i],ohcale[i]/cosh(ohcaleta[i]));
			}
		} else { // using clusters to find total energy in towards, transverse and away regions 
			for (int i = 0; i < clsmult; i++) {
				float cluster_energy = cluster_e[i];
				if (!sim) {
					for (int n = 0; n < cluster_ntowers[i]; n++) {
						if (cluster_tower_e[i][n] < -1.4) {
							cluster_energy -= cluster_tower_e[i][n];
							std::cout << "cluster with very negative energy tower: calo " << cluster_tower_calo[i][n] << " ieta " << cluster_tower_ieta[i][n] << " iphi " << cluster_tower_iphi[i][n] << " energy " << cluster_tower_e[i][n] << std::endl;
							
						}
					}
				}
				h_ue_2D_total->Fill(cluster_eta[i],cluster_phi[i],cluster_energy/cosh(cluster_eta[i]));
				et_total += cluster_energy/cosh(cluster_eta[i]);
				for (int j = 0; j < 8; j++) {
					if (cluster_energy > float(topo_thresholds[j]/1000.0)) {
						ntopo_total[j] += 1;
						sume_topo_total[j] += cluster_energy/cosh(cluster_eta[i]);
						h_topo_total[j]->Fill(cluster_energy/cosh(cluster_eta[i]));
						h_2D_topo_total[j]->Fill(cluster_eta[i],cluster_phi[i],cluster_energy/cosh(cluster_eta[i]));
					}
				}
			}

  			h_et_total->Fill(et_total);
  			h_et_total_pos->Fill(et_total);
  			h_ue_total->Fill(et_total/(secteta*sectphi));

  			for (int i = 0; i < 8; i++) {
  				h_ntopo_total[i]->Fill(ntopo_total[i]);
  				h_sume_topo_total[i]->Fill(sume_topo_total[i]);
  			}

  			events++;
  		}

  	}

  	h_ue_2D_total->Scale(1.0/(events*deltaeta*deltaphi));

  	for (int i = 0; i < 6; i++) {
  		h_2D_topo_total[i]->Scale(1.0/(events*deltaeta*deltaphi));
  	}

  	out->Write();
  	out->Close();

}