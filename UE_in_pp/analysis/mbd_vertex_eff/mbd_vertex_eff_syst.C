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
#include "TRandom.h"
#include "TRandom3.h"
#include "../unfolding/unfold_Def.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include <TH2D.h>

using namespace std;

bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled);
void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et);
void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_et, std::vector<float>* jet_eta);
bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi);
void get_calibjet(float& calibjet_pt, float& calibjet_eta, float& calibjet_phi, float jet_pt, float jet_eta, float jet_phi, bool reco_bkg_cut, TF1* f_corr, float jes_para, float jer_para);
void get_truthjet(float& goodtruthjet_pt, float& goodtruthjet_eta, float& goodtruthjet_phi, float jet_pt, float jet_eta, float jet_phi);
void match_meas_truth(float meas_eta, float meas_phi, bool& matched, float truth_eta, float truth_phi, float jet_radius);
void find_match_meas_truth(float meas_eta, float meas_phi, bool& matched, float& deltaR, float truth_eta, float truth_phi, float jet_radius);

TRandom3 randGen(1234);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                       //
//  Note: 10.14.25                                                                                                                       //
//  This macro creates efficiency plots for jet background cuts for 2D unfolding of leading jet pT and transverse region energy density  //
//                                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mbd_vertex_eff_syst(int sim = 0, int runnumber = 51274, std::string bkg_cut = "dijet", float zvertex_cut = 20.0, std::string runtype = "mb", int start_seg = 0, int end_seg = 200)
{
  /////////////// General Set up ///////////////
  string outfilename;
  if (sim) {
    outfilename = "analysis_syst_output/output_sim_zvtx_lt_" + to_string(floor(zvertex_cut)) + "_" + bkg_cut + "_" + runtype + "_" + to_string(start_seg) + "_" + to_string(end_seg) + ".root";
  } else {
    outfilename = "analysis_syst_output/output_zvtx_lt_" + to_string(floor(zvertex_cut)) + "_" + bkg_cut + "_" + to_string(runnumber) + ".root";
  }
  TFile* f_out = new TFile(outfilename.c_str(), "RECREATE");

  /////////////// Read Files ///////////////
  const char* baseDirJet = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput";
  TChain chain("T");
  if (sim) {
    //for (int i = start_seg; i < end_seg; ++i) {
    //  chain.Add(Form("%s/mbd_eff_sim_run28_%s_output_%d.root", baseDirJet, runtype.c_str(), i));
    //}
    chain.Add(Form("mbd_eff_sim_run28_%s_output_2.root", runtype.c_str()));
  } else { 
    //chain.Add(Form("%s/mbd_eff_output_ana509_%d_*.root", baseDirJet, runnumber));
    chain.Add("mbd_eff_output_ana509.root");
  }
  chain.SetBranchStatus("*", 0);

  std::vector<int>* gl1_trigger_vector_scaled = nullptr;
  chain.SetBranchStatus("triggerVector", 1);
  chain.SetBranchAddress("triggerVector", &gl1_trigger_vector_scaled);
  float zvertex;
  chain.SetBranchStatus("zvtx", 1);
  chain.SetBranchAddress("zvtx", &zvertex);
  float truthvertex;
  if (sim) {
    chain.SetBranchStatus("ztruthvtx", 1);
    chain.SetBranchAddress("ztruthvtx", &truthvertex);
  }
  std::vector<float>* unsubjet_e = nullptr;
  chain.SetBranchStatus("e", 1);
  chain.SetBranchAddress("e", &unsubjet_e);
  std::vector<float>* unsubjet_pt = nullptr;
  chain.SetBranchStatus("pt", 1);
  chain.SetBranchAddress("pt", &unsubjet_pt);
  std::vector<float>* calibjet_pt = nullptr;
  chain.SetBranchStatus("calibpt", 1);
  chain.SetBranchAddress("calibpt", &calibjet_pt);
  std::vector<float>* unsubjet_eta = nullptr;
  chain.SetBranchStatus("eta", 1);
  chain.SetBranchAddress("eta", &unsubjet_eta);
  std::vector<float>* unsubjet_phi = nullptr;
  chain.SetBranchStatus("phi", 1);
  chain.SetBranchAddress("phi", &unsubjet_phi);
  std::vector<float>* unsubjet_emcal_calo_e = nullptr;
  chain.SetBranchStatus("jetEmcalE", 1);
  chain.SetBranchAddress("jetEmcalE", &unsubjet_emcal_calo_e);
  std::vector<float>* unsubjet_ihcal_calo_e = nullptr;
  chain.SetBranchStatus("jetIhcalE", 1);
  chain.SetBranchAddress("jetIhcalE", &unsubjet_ihcal_calo_e);
  std::vector<float>* unsubjet_ohcal_calo_e = nullptr;
  chain.SetBranchStatus("jetOhcalE", 1);
  chain.SetBranchAddress("jetOhcalE", &unsubjet_ohcal_calo_e);
  std::vector<float>* unsubjet_time = nullptr;
  chain.SetBranchStatus("jettime", 1);
  chain.SetBranchAddress("jettime", &unsubjet_time);
  float mbd_t0;
  chain.SetBranchStatus("mbd_t0", 1);
  chain.SetBranchAddress("mbd_t0", &mbd_t0);

  std::vector<float>* truthjet_pt = nullptr;
  std::vector<float>* truthjet_eta = nullptr;
  std::vector<float>* truthjet_phi = nullptr;
  if (sim) {
      chain.SetBranchStatus("truthPt", 1);
      chain.SetBranchAddress("truthPt", &truthjet_pt);
      chain.SetBranchStatus("truthEta", 1);
      chain.SetBranchAddress("truthEta", &truthjet_eta);
      chain.SetBranchStatus("truthPhi", 1);
      chain.SetBranchAddress("truthPhi", &truthjet_phi);
    }

  int clsmult = 0;
  float cluster_e[2000] = {0.0};
  float cluster_eta[2000] = {0.0};
  float cluster_phi[2000] = {0.0};

  chain.SetBranchStatus("clsmult", 1);
  chain.SetBranchAddress("clsmult", &clsmult);
  chain.SetBranchStatus("cluster_e", 1);
  chain.SetBranchAddress("cluster_e", cluster_e);
  chain.SetBranchStatus("cluster_eta", 1);
  chain.SetBranchAddress("cluster_eta", cluster_eta);
  chain.SetBranchStatus("cluster_phi", 1);
  chain.SetBranchAddress("cluster_phi", cluster_phi);

  ////////// JES func //////////
  TFile* corrFile = new TFile("../unfolding/JES_Calib_Default.root", "READ");
  if (!corrFile)
  {
    std::cout << "Error: cannot open JES_Calib_Default.root" << std::endl;
    return;
  }
  TF1* f_corr = (TF1*) corrFile->Get("JES_Calib_Default_Func");
  if (!f_corr)
  {
    std::cout << "Error: cannot open f_corr" << std::endl;
    return;
  }

    TFile *f_zvertex;
    TH1D* h_vertex_weight;
    string zvertex_filename = "output_vertex_reweight_run28_dijet_bkg_cut_iter_1.root";
    f_zvertex = new TFile(zvertex_filename.c_str(),"READ");
    TH1D* htemp = (TH1D*)f_zvertex->Get("vertex_weight");
    if (!htemp) {
        std::cerr << "ERROR: Could not find histogram 'vertex_weight'!" << std::endl;
        exit(1);
    }
    h_vertex_weight = (TH1D*)htemp->Clone("vertex_weight_clone");
    h_vertex_weight->SetDirectory(0);  // detach from file
    f_zvertex->Close();

  ////////////// Jet Background Timing Cut Efficiency ////////////
  double f_timingcut_eff = 0.95;
  double f_timingcut_eff_up = 0.99;
  double timingcut_scale = 1.0 / f_timingcut_eff;
  double timingcut_scale_up = 1.0 / f_timingcut_eff_up;
  double lead_time_cut[2] = {-8.0, 4.0};     // leading time cut
  double lead_time_var[2] = {-9.0, 5.0};     // leading time variation
  double deltat_mbd_cut[2] = {-5.0, 1.0};    // delta time cut
  double deltat_dijet_cut[2] = {-3.0, 3.0};  // delta time cut
  double deltat_mbd_var[2] = {-6.0, 2.0};    // delta time variation
  double deltat_dijet_var[2] = {-4.0, 4.0};  // delta timme variation

  double weight_scale = 1.0, truthjet_pt_min = 0, truthjet_pt_max = 3000;
  if (sim) {
    if (runtype == "jet12") {
        weight_scale = Jet12GeV_scale;
        truthjet_pt_min = 14;
        truthjet_pt_max = 21;
    } else if (runtype == "jet20") {
        weight_scale = Jet20GeV_scale;
        truthjet_pt_min = 21;
        truthjet_pt_max = 32;
    } else if (runtype == "jet30") {
        weight_scale = Jet30GeV_scale;
        truthjet_pt_min = 32;
        truthjet_pt_max = 42;
    } else if (runtype == "jet40") {
        weight_scale = Jet40GeV_scale;
        truthjet_pt_min = 42;
        truthjet_pt_max = 52;
    } else if (runtype == "jet50") {
        weight_scale = Jet50GeV_scale;
        truthjet_pt_min = 52;
        truthjet_pt_max = 62;
    } else if (runtype == "jet60") {
        weight_scale = Jet60GeV_scale;
        truthjet_pt_min = 62;
        truthjet_pt_max = 3000;
    }
  }

  ////////// Histograms //////////
  TH1D* h_zvertex = new TH1D("h_zvertex", ";Truth z-vertex [cm]; Weight", 400, -200, 200);
  TH1D* h_zvertex_noreweight = new TH1D("h_zvertex_noreweight", ";Truth z-vertex [cm]; Weight", 400, -200, 200);
  
  TH1D* h_lead_truth_spectra_record = new TH1D("h_lead_truth_spectra_record", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D* h_lead_spectra_record = new TH1D("h_lead_spectra_record", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D* h_et_transverse_record = new TH1D("h_et_transverse_record", ";#SigmaE_{T} [GeV]", calibnet, calibetbins);
  TH1D* h_nw_et_transverse_record = new TH1D("h_nw_et_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
  TH2D* h_ue_pt_transverse_record = new TH2D("h_ue_pt_transverse_record", "", calibnpt, calibptbins, calibnet, calibetbins);

  TH2D* h_total_measure_mbdvertex = new TH2D("h_total_measure_mbdvertex", "", calibnpt, calibptbins, calibnet, calibetbins);
  TH2D* h_pass_cut_measure_mbdvertex = new TH2D("h_pass_cut_measure_mbdvertex", "", calibnpt, calibptbins, calibnet, calibetbins);

  ////////// Event Loop //////////
  std::cout << "Data analysis started." << std::endl;
  Long64_t nEntries = chain.GetEntries();
  std::cout << "Total number of events: " << nEntries << std::endl;

  // Event variables setup.
  bool reco_bkg_cut = false;
  float calibjet_pt_dijet, calibjet_eta_dijet, calibjet_phi_dijet;
  bool calibjet_matched_dijet, qa_matched;

  for (Long64_t entry = 0; entry < nEntries; ++entry)
  {
    //for (Long64_t entry = 0; entry < 2000; ++entry) {
        if (entry % 10000 == 0) cout << "event " << entry << endl;
        chain.GetEntry(entry);

        if (!sim && !check_bad_trigger(gl1_trigger_vector_scaled)) continue;
        if (sim) {
            if (fabs(truthvertex) > zvertex_cut) continue;
        }
        //////////////////////////// SETUP JET VARIABLES FOR UNFOLDING ////////////////////////////
        // Fill z-vertex histogram.
        double truthvertex_weight = 1.0;
        int ind_truth_lead = -1;
        if (sim) {
            h_zvertex_noreweight->Fill(truthvertex, weight_scale);
            int truthvertex_bin = h_vertex_weight->FindBin(truthvertex);
            truthvertex_weight = h_vertex_weight->GetBinContent(truthvertex_bin);
            h_zvertex->Fill(truthvertex, truthvertex_weight*weight_scale);
            if (truthjet_pt->size() < 1) { continue; }
            get_leading_jet(ind_truth_lead, truthjet_pt, truthjet_eta);
            if (truthjet_pt->at(ind_truth_lead) < truthjet_pt_min || truthjet_pt->at(ind_truth_lead) > truthjet_pt_max) { continue; } 
        } else {
            if (fabs(zvertex) < 200) {
                h_zvertex->Fill(zvertex, weight_scale);
            }
        }

        // indices to find leading and subleading jets 
        int ind_lead = -1; int ind_sub = -1;
        float lead_e = 0;
        reco_bkg_cut = false;

        std::vector<float> recoe_new, recopt_new, recoeta_new, recophi_new, recoemcal_new, recoihcal_new, recoohcal_new;

        int nJetReq = 1;
        if (bkg_cut == "dijet" || bkg_cut == "all") { nJetReq = 2; }
        if (unsubjet_pt->size() < nJetReq) { continue; }

        int Njet = 0;
        for (size_t i = 0; i < unsubjet_pt->size(); i++) {
          if (unsubjet_pt->at(i) >= 5.0) {
            Njet++;
          }
        }

        for (size_t i = 0; i < unsubjet_eta->size(); ++i) {
            float vertex = fabs(zvertex) < 200.0 ? zvertex : 0.0;
            if (!check_bad_jet_eta(unsubjet_eta->at(i), vertex, jet_radius) && fabs(unsubjet_eta->at(i)) < 0.7) {
                recoe_new.push_back(unsubjet_e->at(i));
                recopt_new.push_back(unsubjet_pt->at(i));
                recoeta_new.push_back(unsubjet_eta->at(i));
                recophi_new.push_back(unsubjet_phi->at(i));
                recoemcal_new.push_back(unsubjet_emcal_calo_e->at(i));
                recoihcal_new.push_back(unsubjet_ihcal_calo_e->at(i));
                recoohcal_new.push_back(unsubjet_ohcal_calo_e->at(i));
            }
        }
        // Replace original contents
        *unsubjet_e = std::move(recoe_new);
        *unsubjet_pt = std::move(recopt_new);
        *unsubjet_eta = std::move(recoeta_new);
        *unsubjet_phi = std::move(recophi_new);
        *unsubjet_emcal_calo_e = std::move(recoemcal_new);
        *unsubjet_ihcal_calo_e = std::move(recoihcal_new);
        *unsubjet_ohcal_calo_e = std::move(recoohcal_new);

        if (unsubjet_pt->size() < nJetReq) { continue; }
        
        // find reco leading and subleading jets if nreco jets >= nJetReq
        TVector3 lead;
        if (unsubjet_pt->size() >= nJetReq) {
            if (bkg_cut == "dijet") {
                get_leading_subleading_jet(ind_lead, ind_sub, unsubjet_pt); 
                //std::cout << " sub/lead: " << unsubjet_e->at(ind_sub)/unsubjet_e->at(ind_lead) << " dPhi: " << get_dphi(unsubjet_phi->at(ind_lead), unsubjet_phi->at(ind_sub)) << std::endl; 
                if (unsubjet_e->at(ind_sub)/unsubjet_e->at(ind_lead) > 0.3 && match_leading_subleading_jet(unsubjet_phi->at(ind_lead), unsubjet_phi->at(ind_sub))) {
                    reco_bkg_cut = true;
                } else {
                    reco_bkg_cut = false;
                }
            } else if (bkg_cut == "efrac") {
                get_leading_jet(ind_lead, unsubjet_pt, unsubjet_eta);
                float jete = unsubjet_e->at(ind_lead);
                //std::cout << " EMCal frac: " << unsubjet_emcal_calo_e->at(ind_lead)/jete << " IHCal frac: " << unsubjet_ihcal_calo_e->at(ind_lead)/jete << " OHCal frac: " << unsubjet_ohcal_calo_e->at(ind_lead)/jete << std::endl;
                if (unsubjet_emcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_emcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ihcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ohcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_ohcal_calo_e->at(ind_lead)/jete < 0.9) {
                    reco_bkg_cut = true; 
                } else {
                    reco_bkg_cut = false;
                }
            } else if (bkg_cut == "all") {
                get_leading_subleading_jet(ind_lead, ind_sub, unsubjet_pt);
                float jete = unsubjet_e->at(ind_lead);
                if (unsubjet_e->at(ind_sub)/unsubjet_e->at(ind_lead) > 0.3 && match_leading_subleading_jet(unsubjet_phi->at(ind_lead), unsubjet_phi->at(ind_sub))) {
                    if (unsubjet_emcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_emcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ihcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ohcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_ohcal_calo_e->at(ind_lead)/jete < 0.9) {
                        
                        reco_bkg_cut = true;
                    }
                }
                if (Njet >= 9) {
                    reco_bkg_cut = false; 
                }
                if (!sim) {
                    float lead_time = unsubjet_time->at(ind_lead)*17.6;
                    if (lead_time < lead_time_cut[0] || lead_time > lead_time_cut[1]) {
                    reco_bkg_cut = false;
                    }
                    float lead_deltat;
                    lead_deltat = (unsubjet_time->at(ind_lead) - unsubjet_time->at(ind_sub))*17.6;
                    if (lead_deltat < deltat_dijet_cut[0] || lead_deltat > deltat_dijet_cut[1]) {
                        reco_bkg_cut = false;
                    } 
                    lead_deltat = unsubjet_time->at(ind_lead)*17.6 - mbd_t0;
                    if (lead_deltat < deltat_mbd_cut[0] || lead_deltat > deltat_mbd_cut[1]) {
                        reco_bkg_cut = false;
                    }
                }
            } else {
                reco_bkg_cut = true;
            }
            if (unsubjet_pt->size() >= 1) {
                get_leading_jet(ind_lead, unsubjet_pt, unsubjet_eta);
                lead.SetPtEtaPhi(calibjet_pt->at(ind_lead), unsubjet_eta->at(ind_lead), unsubjet_phi->at(ind_lead));
                lead_e = unsubjet_e->at(ind_lead);
            } else {
                lead.SetPtEtaPhi(0,0,0);
                lead_e = 0;
                reco_bkg_cut = false;
            }
        } else {
            lead.SetPtEtaPhi(0,0,0);
            lead_e = 0;
            reco_bkg_cut = false;
        }

        TVector3 caliblead;
        if (!sim) {
            caliblead.SetPtEtaPhi(lead.Pt(), lead.Eta(), lead.Phi());
        } else {
            caliblead.SetPtEtaPhi(lead.Pt() * (1 + randGen.Gaus(0.0, 0.1)) * 1, lead.Eta(), lead.Phi());
        }
        //std::cout << "Good reco lead: pt " << calibjet_pt_dijet << " eta " << calibjet_eta_dijet << " phi " << calibjet_phi_dijet << " bkg_cut " << reco_bkg_cut << std::endl;
    
        //////////////////////////// SETUP UE VARIABLES FOR UNFOLDING //////////////////////////////

        // find reco ET information 
        float et_transverse = 0;
        float total_et = 0;
        float radius_IH = 127.503;
        for (int i = 0; i < clsmult; i++) {
            float eta = cluster_eta[i];
            total_et += cluster_e[i]/cosh(cluster_eta[i]);
            float dphi = get_dphi(lead.Phi(),cluster_phi[i]);
            if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += cluster_e[i]/cosh(cluster_eta[i]); }
        }

        //////////////////////////// RECORD QA PLOTS FOR ALL EVENTS IN UNFOLDING PROCEDURE ////////////////////////////
        if (reco_bkg_cut) {
            if (sim) h_lead_truth_spectra_record->Fill(truthjet_pt->at(ind_truth_lead), weight_scale*truthvertex_weight);
            h_lead_spectra_record->Fill(caliblead.Pt(), weight_scale*truthvertex_weight);
            if (caliblead.Pt() > calibptbins[0]) {
                h_nw_et_transverse_record->Fill(et_transverse, weight_scale*truthvertex_weight);
                h_et_transverse_record->Fill(et_transverse, weight_scale*truthvertex_weight);
                h_ue_pt_transverse_record->Fill(caliblead.Pt(),et_transverse, weight_scale*truthvertex_weight);
            }
        }

        if (caliblead.Pt() >= calibptbins[0] && caliblead.Pt() < calibptbins[calibnpt] && reco_bkg_cut) {
            h_total_measure_mbdvertex->Fill(caliblead.Pt(), et_transverse, weight_scale*truthvertex_weight);
            if (fabs(zvertex) < zvertex_cut) {
                h_pass_cut_measure_mbdvertex->Fill(caliblead.Pt(), et_transverse, weight_scale*truthvertex_weight);
            }
        }
    }

    std::cout << "Writing histograms..." << std::endl;
    f_out->cd();
    h_zvertex->Write(); h_zvertex_noreweight->Write(); h_lead_truth_spectra_record->Write();
    h_lead_spectra_record->Write(); h_et_transverse_record->Write(); h_nw_et_transverse_record->Write(); h_ue_pt_transverse_record->Write();
    h_total_measure_mbdvertex->Write(); h_pass_cut_measure_mbdvertex->Write();
    std::cout << "All done!" << std::endl;
}

////////////////////////////////////////// Functions //////////////////////////////////////////
bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled) {
  if (std::find(gl1_trigger_vector_scaled->begin(), gl1_trigger_vector_scaled->end(), 22) != gl1_trigger_vector_scaled->end()) { return true; }
  return false;
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

void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_et, std::vector<float>* jet_eta) {
  leadingjet_index = -1;
  float leadingjet_et = -9999;
  for (int ij = 0; ij < jet_et->size(); ++ij) {
    float jetet = jet_et->at(ij);
    if (jetet > leadingjet_et) {
      leadingjet_et = jetet;
      leadingjet_index = ij;
    }
  }
}

bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi) {
  float dijet_min_phi = 7*TMath::Pi()/8.;
  float dphi = get_dphi(leadingjet_phi, subleadingjet_phi);
  return fabs(dphi) > dijet_min_phi;
}

void get_calibjet(float& calibjet_pt, float& calibjet_eta, float& calibjet_phi, float jet_pt, float jet_eta, float jet_phi, bool reco_bkg_cut, TF1* f_corr, float jes_para, float jer_para) {
    calibjet_pt = -9999;
    calibjet_eta = -9999;
    calibjet_phi = -9999;
    //if (!reco_bkg_cut) return;
    double calib_pt = f_corr->Eval(jet_pt) * (1 + randGen.Gaus(0.0, jer_para)) * jes_para;
    if (calib_pt < calibptbins[0] || calib_pt > calibptbins[calibnpt]) return;
    calibjet_pt = calib_pt;
    calibjet_eta = jet_eta;
    calibjet_phi = jet_phi;
}

void get_truthjet(float& goodtruthjet_pt, float& goodtruthjet_eta, float& goodtruthjet_phi, float jet_pt, float jet_eta, float jet_phi) {
    goodtruthjet_pt = -9999;
    goodtruthjet_eta = -9999;
    goodtruthjet_phi = -9999;
    if (jet_pt < truthptbins[0] || jet_pt > truthptbins[truthnpt]) return;
    goodtruthjet_pt = jet_pt;
    goodtruthjet_eta = jet_eta;
    goodtruthjet_phi = jet_phi;
}

void match_meas_truth(float meas_eta, float meas_phi, bool& matched, float truth_eta, float truth_phi, float jet_radius) {
    if (meas_eta == -9999 || meas_phi == -9999 || truth_eta == -9999 || truth_phi == -9999) {
        matched = false;
        return; 
    } else {
        float dR = get_dR(meas_eta, meas_phi, truth_eta, truth_phi);
        if (dR < jet_radius * 0.75) {
            //std::cout << "dR = " << dR << std::endl;
            matched = true;
            return;
        } else {
            //std::cout << "dR = " << dR << std::endl;
            matched = false;
            return;
        }
    }
}

void find_match_meas_truth(float meas_eta, float meas_phi, bool& matched, float& deltaR, float truth_eta, float truth_phi, float jet_radius) {
    if (meas_eta == -9999 || meas_phi == -9999 || truth_eta == -9999 || truth_phi == -9999) {
        matched = false;
        return; 
    } else {
        float dR = get_dR(meas_eta, meas_phi, truth_eta, truth_phi);
        deltaR = dR;
        if (dR < jet_radius * 0.75) {
            //std::cout << "dR = " << dR << std::endl;
            matched = true;
            return;
        } else {
            //std::cout << "dR = " << dR << std::endl;
            matched = false;
            return;
        }
    }
}