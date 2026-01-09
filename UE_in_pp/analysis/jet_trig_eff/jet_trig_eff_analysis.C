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
#include <string>
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

//////////////////////////////////////////////////
//												//
//	Note: 5.26.25 								//
// 	Jet trigger numbers:						//
// 	Trigger | Name 								//
//  10      | MBD  								//
//  17      | Jet 8 + MBD 						//
//  18      | Jet 10 + MBD 						//
//  19      | Jet 12 + MBD 						//
//  12      | MBD + zvtx < 10 cm 				//
//  33      | Jet 8 + MBD + zvtx < 10 cm 		//
//  34      | Jet 10 + MBD + zvtx < 10 cm 		//
//  35      | Jet 12 + MBD + zvtx < 10 cm 		//
//  21      | Jet 8  							//
//  22      | Jet 10   							//
//  23      | Jet 12   							//
//												//
//////////////////////////////////////////////////

void jet_trig_eff_analysis(int runnumber = 51274) {

	// create outfile 
	TFile *out = new TFile(Form("analysis_output/output_%d.root", runnumber),"RECREATE");

    vector<string> vz_name = {"h_vz_mb","h_vz_jet8","h_vz_jet10","h_vz_jet12","h_vz_zcut_mb","h_vz_zcut_jet8","h_vz_zcut_jet10","h_vz_zcut_jet12"};
    vector<string> jet_name = {"h_leadjet_mb","h_leadjet_jet8","h_leadjet_jet10","h_leadjet_jet12","h_leadjet_zcut_mb","h_leadjet_zcut_jet8","h_leadjet_zcut_jet10","h_leadjet_zcut_jet12"};
    vector<string> full_jet_name = {"h_full_leadjet_mb","h_full_leadjet_jet8","h_full_leadjet_jet10","h_full_leadjet_jet12","h_full_leadjet_zcut_mb","h_full_leadjet_zcut_jet8","h_full_leadjet_zcut_jet10","h_full_leadjet_zcut_jet12"};
 
	// create event and jet histograms 
	TH1F h_vz[8];
	TH1F h_leadjet[8];
	TH1F h_leadjet_full[8];
	TEfficiency* leadjet_eff[6] = {nullptr};
	for (int i = 0; i < 8; i++) {
		h_vz[i] = TH1F(vz_name[i].c_str(),"",400,-100,100);
		h_leadjet[i] = TH1F(jet_name[i].c_str(),"",100,0,100);
		h_leadjet_full[i] = TH1F(full_jet_name[i].c_str(),"",100,0,100);
	}

	/////////////// Read Files ///////////////
	const char* baseDirJet = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput";
	TChain chain("T");
	chain.Add(Form("%s/output_1.5mrad_ana468_calofitting_jettrigeff_%d_*.root", baseDirJet, runnumber));
	chain.SetBranchStatus("*", 0);

	// define ttree branch variables 
	int nJet;
	float zvtx;
	vector<int> *triggerVector = nullptr;
	vector<float> *pt = nullptr;
	chain.SetBranchStatus("zvtx",1);
	chain.SetBranchStatus("triggerVector",1);
	chain.SetBranchStatus("pt",1);
	chain.SetBranchAddress("zvtx",&zvtx);
	chain.SetBranchAddress("triggerVector",&triggerVector);
	chain.SetBranchAddress("pt",&pt);

	int eventnumber = 0; // number of events in tree (used for showing iteration through ttree)
	
    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;

    // main event loop 
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 1000; ++entry) {
        chain.GetEntry(entry);
    	if (eventnumber % 10000 == 0) cout << "event " << eventnumber << endl;
    	eventnumber++;

    	// find the leading and subleading jets of the event
  		int ind_lead = -1;
  		float temp_lead = -1;
  		for (int i = 0; i < pt->size(); i++) {
  			if ((*pt)[i] > temp_lead) {
  				temp_lead = (*pt)[i];
  				ind_lead = i;
  			}
  		}

  		if (temp_lead < 0 || fabs(zvtx) > 30) { continue; }

  		bool isMB = false;
  		bool isMBzcut = false;

  		for (auto t : *triggerVector) {
  			if (t == 10) {
  				isMB = true;
  				h_leadjet_full[0].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[0].Fill(zvtx);
  					h_leadjet[0].Fill(temp_lead);
  				}
  			}
  			if (t == 17 && isMB) {
  				h_leadjet_full[1].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[1].Fill(zvtx);
  					h_leadjet[1].Fill(temp_lead);
  				}
  			} 
  			if (t == 18 && isMB) {
  				h_leadjet_full[2].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[2].Fill(zvtx);
  					h_leadjet[2].Fill(temp_lead);
  				}
  			} 
  			if (t == 19 && isMB) {
  				h_leadjet_full[3].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[3].Fill(zvtx);
  					h_leadjet[3].Fill(temp_lead);
  				}
  			} 
  			if (t == 12 && fabs(zvtx) <= 10) {
  				isMBzcut = true;
  				h_leadjet_full[4].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[4].Fill(zvtx);
  					h_leadjet[4].Fill(temp_lead);
  				}
  			}
			if (t == 33 && isMBzcut) {
  				h_leadjet_full[5].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[5].Fill(zvtx);
  					h_leadjet[5].Fill(temp_lead);
  				}
  			} 
  			if (t == 34 && isMBzcut) {
  				h_leadjet_full[6].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[6].Fill(zvtx);
  					h_leadjet[6].Fill(temp_lead);
  				}
  			} 
  			if (t == 35 && isMBzcut) {
  				h_leadjet_full[7].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[7].Fill(zvtx);
  					h_leadjet[7].Fill(temp_lead);
  				}
  			}
  		}
  	}
  	
  	for (int i = 1; i < 4; i++) {
  		if (TEfficiency::CheckConsistency(h_leadjet[i], h_leadjet[0])) {
  			leadjet_eff[i-1] = new TEfficiency(h_leadjet[i], h_leadjet[0]);
  			leadjet_eff[i-1]->SetName(Form("eff_%s", jet_name[i].c_str()));
  		}
  	}
  	
  	for (int i = 5; i < 8; i++) {
  		if (TEfficiency::CheckConsistency(h_leadjet[i], h_leadjet[4])) {
  			leadjet_eff[i-2] = new TEfficiency(h_leadjet[i], h_leadjet[4]);
  			leadjet_eff[i-2]->SetName(Form("eff_%s", jet_name[i].c_str()));
  		}
  	}

  	out->cd();
  	
  	for (int i = 0; i < 6; i++) {
  		if (leadjet_eff[i]) {
  			leadjet_eff[i]->Write();
  		} else {
  			std::cout << "leadjet_eff_" << i << " is null" << std::endl;
  		}
  	}

  	std::cout << "trigger 10 entries " << h_leadjet[0].GetEntries() << " trigger 12 entries " << h_leadjet[4].GetEntries() << std::endl;
  	std::cout << "trigger 17 entries " << h_leadjet[1].GetEntries() << " trigger 33 entries " << h_leadjet[5].GetEntries() << std::endl;
  	std::cout << "trigger 18 entries " << h_leadjet[2].GetEntries() << " trigger 34 entries " << h_leadjet[6].GetEntries() << std::endl;
  	std::cout << "trigger 19 entries " << h_leadjet[3].GetEntries() << " trigger 35 entries " << h_leadjet[7].GetEntries() << std::endl;

  	out->Write();
  	out->Close();

}