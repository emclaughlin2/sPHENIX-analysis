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

////////////////////////////////////////// Functions //////////////////////////////////////////
float get_dphi(float phi1, float phi2) {
  float dphi1 = phi1 - phi2;
  float dphi2 = phi1 - phi2 + 2*TMath::Pi();
  float dphi3 = phi1 - phi2 - 2*TMath::Pi();
  if (fabs(dphi1) > fabs(dphi2)) {
    dphi1 = dphi2;
  }
  if (fabs(dphi1) > fabs(dphi3)) {
    dphi1 = dphi3;
  }
  return dphi1;
}

void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et) {
  leadingjet_index = -1;
  subleadingjet_index = -1;
  float leadingjet_et = -9999;
  float subleadingjet_et = -9999;
  for (int ij = 0; ij < jet_et->size(); ++ij) {
    float jetet = jet_et->at(ij);
    if (jetet > leadingjet_et) {
      subleadingjet_et = leadingjet_et;
      subleadingjet_index = leadingjet_index;
      leadingjet_et = jetet;
      leadingjet_index = ij;
    } else if (jetet > subleadingjet_et) {
      subleadingjet_et = jetet;
      subleadingjet_index = ij;
    }
  }
}

bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi) {
  float dijet_min_phi = 3*TMath::Pi()/4.;
  //float dijet_min_phi = 7*TMath::Pi()/8.;
  float dphi = fabs(get_dphi(leadingjet_phi, subleadingjet_phi));
  return dphi > dijet_min_phi;
}

void jet_trig_eff_analysis(int runnumber = 51274) {

	// create outfile 
	TFile *out = new TFile(Form("ana509_analysis_output/dijet_and_mbdtimecut_output_%d.root", runnumber),"RECREATE");

    vector<string> vz_name = {"h_vz_mb","h_vz_mb_wzcut","h_vz_mb_wzcut_online","h_vz_jet8","h_vz_jet10",
    							"h_vz_jet12","h_vz_trig22","h_vz_trig22_wzcut","h_vz_trig22_wzcut_online","h_vz_zcut_mb",
    							"h_vz_zcut_jet8","h_vz_zcut_jet10","h_vz_zcut_jet12","h_vz_zcut_trig22"};
    vector<string> jet_name = {"h_leadjet_mb","h_leadjet_mb_wzcut","h_leadjet_mb_wzcut_online","h_leadjet_jet8","h_leadjet_jet10",
    							"h_leadjet_jet12","h_leadjet_trig22","h_leadjet_trig22_wzcut","h_leadjet_trig22_wzcut_online","h_leadjet_zcut_mb",
    							"h_leadjet_zcut_jet8","h_leadjet_zcut_jet10","h_leadjet_zcut_jet12","h_leadjet_zcut_trig22"};
    vector<string> full_jet_name = {"h_full_leadjet_mb","h_full_leadjet_mb_wzcut","h_full_leadjet_mb_wzcut_online","h_full_leadjet_jet8","h_full_leadjet_jet10",
    							"h_full_leadjet_jet12","h_full_leadjet_trig22","h_full_leadjet_trig22_wzcut","h_full_leadjet_trig22_wzcut_online","h_full_leadjet_zcut_mb",
    							"h_full_leadjet_zcut_jet8","h_full_leadjet_zcut_jet10","h_full_leadjet_zcut_jet12","h_full_leadjet_zcut_trig22"};
 	vector<string> jet_quad_name = {"h_leadjet_mb_1","h_leadjet_mb_2","h_leadjet_mb_3","h_leadjet_mb_4",
 								"h_leadjet_jet8_1","h_leadjet_jet8_2","h_leadjet_jet8_3","h_leadjet_jet8_4",
 								"h_leadjet_jet10_1","h_leadjet_jet10_2","h_leadjet_jet10_3","h_leadjet_jet10_4",
 								"h_leadjet_jet12_1","h_leadjet_jet12_2","h_leadjet_jet12_3","h_leadjet_jet12_4",
    							"h_leadjet_zcut_mb_1","h_leadjet_zcut_mb_2","h_leadjet_zcut_mb_3","h_leadjet_zcut_mb_4",
    							"h_leadjet_zcut_jet8_1","h_leadjet_zcut_jet8_2","h_leadjet_zcut_jet8_3","h_leadjet_zcut_jet8_4",
    							"h_leadjet_zcut_jet10_1","h_leadjet_zcut_jet10_2","h_leadjet_zcut_jet10_3","h_leadjet_zcut_jet10_4",
    							"h_leadjet_zcut_jet12_1","h_leadjet_zcut_jet12_2","h_leadjet_zcut_jet12_3","h_leadjet_zcut_jet12_4"};
	// create event and jet histograms 
	TH1F h_vz[14];
	TH1F h_leadjet[14];
	TH1F h_leadjet_full[14];
	TEfficiency* leadjet_eff[10] = {nullptr};
	for (int i = 0; i < 14; i++) {
		h_vz[i] = TH1F(vz_name[i].c_str(),"",400,-100,100);
		h_leadjet[i] = TH1F(jet_name[i].c_str(),"",100,0,100);
		h_leadjet_full[i] = TH1F(full_jet_name[i].c_str(),"",100,0,100);
	}
	// 
	TH1F h_leadjet_quad[32];
	TEfficiency* leadjet_quad_eff[24] = {nullptr};
	for (int i = 0; i < 32; i++) {
		h_leadjet_quad[i] = TH1F(jet_quad_name[i].c_str(),"",100,0,100);
	}

	/////////////// Read Files ///////////////
	const char* baseDirJet = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput";
	TChain chain("T");
	chain.Add(Form("%s/output_ana509_jet_trig_%d_*.root", baseDirJet, runnumber));
	chain.SetBranchStatus("*", 0);

	// define ttree branch variables 
	int nJet;
	float zvtx;
	vector<int> *triggerVector = nullptr;
	vector<float> *pt = nullptr;
	vector<float> *phi = nullptr;
	vector<float> *e = nullptr;
	std::vector<float>* jettime = nullptr;
	vector<int> *liveTriggerVector = nullptr;
	chain.SetBranchStatus("zvtx",1);
	chain.SetBranchStatus("triggerVector",1);
	chain.SetBranchStatus("liveTriggerVector",1);
	chain.SetBranchStatus("pt",1);
	chain.SetBranchStatus("e",1);
	chain.SetBranchStatus("phi",1);
	chain.SetBranchStatus("jettime", 1);
	chain.SetBranchAddress("zvtx",&zvtx);
	chain.SetBranchAddress("triggerVector",&triggerVector);
	chain.SetBranchAddress("liveTriggerVector",&liveTriggerVector);
	chain.SetBranchAddress("pt",&pt);
	chain.SetBranchAddress("phi",&phi);
	chain.SetBranchAddress("e",&e);
	chain.SetBranchAddress("jettime", &jettime);

	std::vector<float>* unsubjet_emcal_calo_e = nullptr; chain.SetBranchStatus("jetEmcalE", 1); chain.SetBranchAddress("jetEmcalE", &unsubjet_emcal_calo_e);
  	std::vector<float>* unsubjet_ihcal_calo_e = nullptr; chain.SetBranchStatus("jetIhcalE", 1); chain.SetBranchAddress("jetIhcalE", &unsubjet_ihcal_calo_e);
  	std::vector<float>* unsubjet_ohcal_calo_e = nullptr; chain.SetBranchStatus("jetOhcalE", 1); chain.SetBranchAddress("jetOhcalE", &unsubjet_ohcal_calo_e);
  	float mbd_t0; chain.SetBranchStatus("mbd_t0", 1); chain.SetBranchAddress("mbd_t0", &mbd_t0);

	// jet phi from [-pi,pi]

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
    	/*
  		int ind_lead = -1;
  		float temp_lead = -1;
  		for (int i = 0; i < pt->size(); i++) {
  			if ((*pt)[i] > temp_lead) {
  				temp_lead = (*pt)[i];
  				ind_lead = i;
  			}
  		}
  		if (temp_lead < 0 || fabs(zvtx) > 60) { continue; }
  		*/

    	// hanpu's jet background cut code
  		int ind_lead = -1;
	    int ind_sub = -1;
	    get_leading_subleading_jet(ind_lead, ind_sub, e);
	    if (ind_lead < 0) continue;
	    // hanpu's jet background cut code

  		bool pass_bkg = false;
	    bool pass_timing = false;
	    bool pass_njet = false;

        //float jete = e->at(ind_lead);
        //if (unsubjet_emcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_emcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ihcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ohcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_ohcal_calo_e->at(ind_lead)/jete < 0.9) {
        //    pass_bkg = true; 
        //} 
        //float lead_deltat = jettime->at(ind_lead)*17.6 - mbd_t0;
	    //if (fabs(jettime->at(ind_lead)*17.6+2) < 6 && lead_deltat >= -5.0 && lead_deltat <= 0.0) {
	    //	pass_timing = true;
	    //}

    	// hanpu's jet background cut code
       	if (ind_sub >= 0 && e->at(ind_sub) / (float) e->at(ind_lead) > 0.3) {
	      	pass_bkg = match_leading_subleading_jet(phi->at(ind_lead), phi->at(ind_sub));
	    }
	    if (ind_sub >= 0 && fabs(jettime->at(ind_lead)*17.6+2) < 6 && fabs(jettime->at(ind_lead)*17.6 - jettime->at(ind_sub)*17.6) < 3) {
	    	float lead_deltat = jettime->at(ind_lead)*17.6 - mbd_t0;
	    	if (lead_deltat >= -5.0 && lead_deltat <= 0.0) {
	    		pass_timing = true;
	    	}
		    //pass_timing = true;
		}
		// hanpu's jet background cut code

	    int n_jets = 0;
	    for (int ij = 0; ij < e->size(); ++ij) {
	    	if (e->at(ij) > 5.0) {
	    		n_jets++;
	    	}
	    }
	    if (n_jets <= 9) {
	    	pass_njet = true;
	    }

	    if (!pass_bkg || !pass_timing || !pass_njet) { 
	    	//std::cout << "pass bkg " << pass_bkg << " timing " << pass_timing << " njet " << pass_njet << std::endl;
	    	//if (!pass_bkg) {
	    	//	std::cout << "\t EMCal frac: " << unsubjet_emcal_calo_e->at(ind_lead)/jete << " IHCal frac: " << unsubjet_ihcal_calo_e->at(ind_lead)/jete << " OHCal frac: " << unsubjet_ohcal_calo_e->at(ind_lead)/jete << std::endl;
        	//}
	    	continue; 
	    }

    	// hanpu's jet background cut code
	    get_leading_subleading_jet(ind_lead, ind_sub, pt);
    	float temp_lead = -1;
    	if(ind_lead >= 0) {
      		temp_lead = pt->at(ind_lead);
    	}
    	if (temp_lead < 0) { continue; }
    	// hanpu's jet background cut code

  		int phi_index = -1;
  		if (phi->at(ind_lead) >= M_PI/2.0) {
  			phi_index = 1;
  		} else if (phi->at(ind_lead) >= 0) {
  			phi_index = 0;
  		} else if (phi->at(ind_lead) >= -M_PI/2.0) {
  			phi_index = 2;
  		} else if (phi->at(ind_lead) >= -M_PI) {
  			phi_index = 3;
  		} else {
  			std::cout << "Phi value out of bounds. Phi = " << phi->at(ind_lead) << std::endl;
  			continue;
  		}

  		bool isMB = false;
  		bool isMBzcut = false;

  		for (auto t : *triggerVector) {
  			if (t == 10) {
  				isMB = true;
  				h_leadjet_full[0].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[0].Fill(zvtx);
  					h_leadjet[0].Fill(temp_lead);
  					h_leadjet_quad[0+phi_index].Fill(temp_lead);
  				}
  				if (fabs(zvtx) <= 10) {
  					h_leadjet_full[1].Fill(temp_lead);
  					if (temp_lead >= 7) {
  						h_vz[1].Fill(zvtx);
  						h_leadjet[1].Fill(temp_lead);
  					}
  				}
  				if (find(liveTriggerVector->begin(), liveTriggerVector->end(), 12) != liveTriggerVector->end()) {
  					h_leadjet_full[2].Fill(temp_lead);
  					if (temp_lead >= 7) {
  						h_vz[2].Fill(zvtx);
  						h_leadjet[2].Fill(temp_lead);
  					}
  				}
  			}
  			if (t == 17 && isMB) {
  				h_leadjet_full[3].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[3].Fill(zvtx);
  					h_leadjet[3].Fill(temp_lead);
  					h_leadjet_quad[4+phi_index].Fill(temp_lead);
  				}
  			} 
  			if (t == 18 && isMB) {
  				h_leadjet_full[4].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[4].Fill(zvtx);
  					h_leadjet[4].Fill(temp_lead);
  					h_leadjet_quad[8+phi_index].Fill(temp_lead);
  				}
  			} 
  			if (t == 19 && isMB) {
  				h_leadjet_full[5].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[5].Fill(zvtx);
  					h_leadjet[5].Fill(temp_lead);
  					h_leadjet_quad[12+phi_index].Fill(temp_lead);
  				}
  			} 
  			if (t == 22 && isMB) {
  				h_leadjet_full[6].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[6].Fill(zvtx);
  					h_leadjet[6].Fill(temp_lead);
  				}
  				if (fabs(zvtx) <= 10) {
	  				h_leadjet_full[7].Fill(temp_lead);
	  				if (temp_lead >= 7) {
	  					h_vz[7].Fill(zvtx);
	  					h_leadjet[7].Fill(temp_lead);
	  				}
	  			}
  				if (find(liveTriggerVector->begin(), liveTriggerVector->end(), 12) != liveTriggerVector->end()) {
  					h_leadjet_full[8].Fill(temp_lead);
  					if (temp_lead >= 7) {
  						h_vz[8].Fill(zvtx);
  						h_leadjet[8].Fill(temp_lead);
  					}
  				}
  			}
  			if (t == 12 && fabs(zvtx) <= 10) {
  				isMBzcut = true;
  				h_leadjet_full[9].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[9].Fill(zvtx);
  					h_leadjet[9].Fill(temp_lead);
  					h_leadjet_quad[16+phi_index].Fill(temp_lead);
  				}
  			}
			if (t == 33 && isMBzcut) {
  				h_leadjet_full[10].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[10].Fill(zvtx);
  					h_leadjet[10].Fill(temp_lead);
  					h_leadjet_quad[20+phi_index].Fill(temp_lead);
  				}
  			} 
  			if (t == 34 && isMBzcut) {
  				h_leadjet_full[11].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[11].Fill(zvtx);
  					h_leadjet[11].Fill(temp_lead);
  					h_leadjet_quad[24+phi_index].Fill(temp_lead);
  				}
  			} 
  			if (t == 35 && isMBzcut) {
  				h_leadjet_full[12].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[12].Fill(zvtx);
  					h_leadjet[12].Fill(temp_lead);
  					h_leadjet_quad[28+phi_index].Fill(temp_lead);
  				}
  			}
  			if (t == 22 && isMBzcut) {
  				h_leadjet_full[13].Fill(temp_lead);
  				if (temp_lead >= 7) {
  					h_vz[13].Fill(zvtx);
  					h_leadjet[13].Fill(temp_lead);
  				}
  			}
  		}
  	}

  	// trig 10
  	for (int i = 3; i < 7; i++) {
  		if (TEfficiency::CheckConsistency(h_leadjet[i], h_leadjet[0])) {
  			leadjet_eff[i-3] = new TEfficiency(h_leadjet[i], h_leadjet[0]);
  			leadjet_eff[i-3]->SetName(Form("eff_%s", jet_name[i].c_str()));
  		}
  	}

  	// trig 10 w offline zvertex cut
  	if (TEfficiency::CheckConsistency(h_leadjet[7], h_leadjet[1])) {
  		leadjet_eff[4] = new TEfficiency(h_leadjet[7], h_leadjet[1]);
  		leadjet_eff[4]->SetName(Form("eff_%s", jet_name[7].c_str()));
  	}

  	// trig 10 w online zvertex cut
  	if (TEfficiency::CheckConsistency(h_leadjet[8], h_leadjet[2])) {
  		leadjet_eff[5] = new TEfficiency(h_leadjet[8], h_leadjet[2]);
  		leadjet_eff[5]->SetName(Form("eff_%s", jet_name[8].c_str()));
  	}


  	// trig 12
  	for (int i = 10; i < 14; i++) {
  		if (TEfficiency::CheckConsistency(h_leadjet[i], h_leadjet[9])) {
  			leadjet_eff[i-4] = new TEfficiency(h_leadjet[i], h_leadjet[9]);
  			leadjet_eff[i-4]->SetName(Form("eff_%s", jet_name[i].c_str()));
  		}
  	}

  	// quadrant study
  	for (int i = 0; i < 4; i++) {
  		if (TEfficiency::CheckConsistency(h_leadjet_quad[4+i], h_leadjet_quad[i])) {
  			leadjet_quad_eff[i] = new TEfficiency(h_leadjet_quad[4+i], h_leadjet_quad[i]);
  			leadjet_quad_eff[i]->SetName(Form("eff_%s", jet_quad_name[4+i].c_str()));
  		}
  		if (TEfficiency::CheckConsistency(h_leadjet_quad[8+i], h_leadjet_quad[i])) {
  			leadjet_quad_eff[i+4] = new TEfficiency(h_leadjet_quad[8+i], h_leadjet_quad[i]);
  			leadjet_quad_eff[i+4]->SetName(Form("eff_%s", jet_quad_name[8+i].c_str()));
  		}
  		if (TEfficiency::CheckConsistency(h_leadjet_quad[12+i], h_leadjet_quad[i])) {
  			leadjet_quad_eff[i+8] = new TEfficiency(h_leadjet_quad[12+i], h_leadjet_quad[i]);
  			leadjet_quad_eff[i+8]->SetName(Form("eff_%s", jet_quad_name[12+i].c_str()));
  		}
  		if (TEfficiency::CheckConsistency(h_leadjet_quad[20+i], h_leadjet_quad[16+i])) {
  			leadjet_quad_eff[i+12] = new TEfficiency(h_leadjet_quad[20+i], h_leadjet_quad[16+i]);
  			leadjet_quad_eff[i+12]->SetName(Form("eff_%s", jet_quad_name[20+i].c_str()));
  		}
  		if (TEfficiency::CheckConsistency(h_leadjet_quad[24+i], h_leadjet_quad[16+i])) {
  			leadjet_quad_eff[i+16] = new TEfficiency(h_leadjet_quad[24+i], h_leadjet_quad[16+i]);
  			leadjet_quad_eff[i+16]->SetName(Form("eff_%s", jet_quad_name[24+i].c_str()));
  		}
  		if (TEfficiency::CheckConsistency(h_leadjet_quad[28+i], h_leadjet_quad[16+i])) {
  			leadjet_quad_eff[i+20] = new TEfficiency(h_leadjet_quad[28+i], h_leadjet_quad[16+i]);
  			leadjet_quad_eff[i+20]->SetName(Form("eff_%s", jet_quad_name[28+i].c_str()));
  		}
  	}

  	out->cd();
  	
  	for (int i = 0; i < 10; i++) {
  		if (leadjet_eff[i]) {
  			leadjet_eff[i]->Write();
  		} else {
  			std::cout << "leadjet_eff_" << i << " is null" << std::endl;
  		}
  	}

  	for (int i = 0; i < 24; i++) {
  		if (leadjet_quad_eff[i]) {
  			leadjet_quad_eff[i]->Write();
  		} else {
  			std::cout << "leadjet_quad_eff_" << i << " is null" << std::endl;
  		}
  	}

  	std::cout << "trigger 10 entries " << h_leadjet[0].GetEntries() << " trigger 12 entries " << h_leadjet[7].GetEntries() << std::endl;
  	std::cout << "trigger 17 entries " << h_leadjet[1].GetEntries() << " trigger 33 entries " << h_leadjet[8].GetEntries() << std::endl;
  	std::cout << "trigger 18 entries " << h_leadjet[2].GetEntries() << " trigger 34 entries " << h_leadjet[9].GetEntries() << std::endl;
  	std::cout << "trigger 19 entries " << h_leadjet[3].GetEntries() << " trigger 35 entries " << h_leadjet[10].GetEntries() << std::endl;

  	out->Write();
  	out->Close();

}