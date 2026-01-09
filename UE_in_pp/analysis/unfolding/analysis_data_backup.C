#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "unfold_Def.h"

bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled);
void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et);
bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi);
void get_jet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius);
void filter_jet(bool& jet_filter, float jet_e, float jet_pt, float jet_eta, float zvertex, float jet_radius);
void get_calibjet(float& calibjet_pt, float& calibjet_eta, float& calibjet_phi, bool& calibjet_background_dijet, bool& calibjet_background_frac, bool jet_filter, float jet_pt, float jet_eta, float jet_phi, bool background_dijet, bool background_frac, TF1* f_corr);

////////////////////////////////////////// Main Function //////////////////////////////////////////
void analysis_data(int runnumber = 51274, bool clusters = true, bool emcal_clusters = false)  {
  /////////////// General Set up ///////////////
  TFile *f_out = new TFile(Form("analysis_data_output/output_%d.root", runnumber),"RECREATE");

  /////////////// Read Files ///////////////
  const char* baseDirJet = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput";
  TChain chain("T");
  if (runnumber < 51274) {
    chain.Add(Form("%s/output_0mrad_ana468_%d_*.root", baseDirJet, runnumber));
  } else {
    chain.Add(Form("%s/output1.5mrad_ana468_%d_*.root", baseDirJet, runnumber));
  }
  chain.SetBranchStatus("*", 0);

  float zvertex; chain.SetBranchStatus("zvtx", 1); chain.SetBranchAddress("zvtx", &zvertex);
  std::vector<int>* gl1_trigger_vector_scaled = nullptr; chain.SetBranchStatus("triggerVector", 1); chain.SetBranchAddress("triggerVector", &gl1_trigger_vector_scaled);
  std::vector<float>* unsubjet_e = nullptr; chain.SetBranchStatus("e", 1); chain.SetBranchAddress("e", &unsubjet_e);
  std::vector<float>* unsubjet_pt = nullptr; chain.SetBranchStatus("pt", 1); chain.SetBranchAddress("pt", &unsubjet_pt);
  std::vector<float>* unsubjet_eta = nullptr; chain.SetBranchStatus("eta", 1); chain.SetBranchAddress("eta", &unsubjet_eta);
  std::vector<float>* unsubjet_phi = nullptr; chain.SetBranchStatus("phi", 1); chain.SetBranchAddress("phi", &unsubjet_phi);
  std::vector<float>* unsubjet_emcal_calo_e = nullptr; chain.SetBranchStatus("jetEmcalE", 1); chain.SetBranchAddress("jetEmcalE", &unsubjet_emcal_calo_e);
  std::vector<float>* unsubjet_ihcal_calo_e = nullptr; chain.SetBranchStatus("jetIhcalE", 1); chain.SetBranchAddress("jetIhcalE", &unsubjet_ihcal_calo_e);
  std::vector<float>* unsubjet_ohcal_calo_e = nullptr; chain.SetBranchStatus("jetOhcalE", 1); chain.SetBranchAddress("jetOhcalE", &unsubjet_ohcal_calo_e);

  int emcaln = 0; float emcale[24576] = {0.0}; float emcaleta[24576] = {0.0}; float emcalphi[24576] = {0.0};
  int ihcaln = 0; float ihcale[1536] = {0.0}; float ihcaleta[1536] = {0.0}; float ihcalphi[1536] = {0.0};
  int ohcaln = 0; float ohcale[1536] = {0.0}; float ohcaleta[1536] = {0.0}; float ohcalphi[1536] = {0.0};
  int clsmult = 0; float cluster_e[10000] = {0.0}; float cluster_eta[10000] = {0.0}; float cluster_phi[10000] = {0.0};
  if (!clusters) {
      chain.SetBranchStatus("emcaln", 1); chain.SetBranchAddress("emcaln",&emcaln);
      chain.SetBranchStatus("emcale", 1); chain.SetBranchAddress("emcale",emcale);
      chain.SetBranchStatus("emcaleta", 1); chain.SetBranchAddress("emcaleta",emcaleta);
      chain.SetBranchStatus("emcalphi", 1); chain.SetBranchAddress("emcalphi",emcalphi);
      chain.SetBranchStatus("ihcaln", 1); chain.SetBranchAddress("ihcaln",&ihcaln);
      chain.SetBranchStatus("ihcale", 1); chain.SetBranchAddress("ihcale",ihcale);
      chain.SetBranchStatus("ihcaleta", 1); chain.SetBranchAddress("ihcaleta",ihcaleta);
      chain.SetBranchStatus("ihcalphi", 1); chain.SetBranchAddress("ihcalphi",ihcalphi);
      chain.SetBranchStatus("ohcaln", 1); chain.SetBranchAddress("ohcaln",&ohcaln);
      chain.SetBranchStatus("ohcale", 1); chain.SetBranchAddress("ohcale",ohcale);
      chain.SetBranchStatus("ohcaleta", 1); chain.SetBranchAddress("ohcaleta",ohcaleta);
      chain.SetBranchStatus("ohcalphi", 1); chain.SetBranchAddress("ohcalphi",ohcalphi);
  } else if (clusters && !emcal_clusters) {
      chain.SetBranchStatus("clsmult", 1); chain.SetBranchAddress("clsmult",&clsmult);
      chain.SetBranchStatus("cluster_e", 1); chain.SetBranchAddress("cluster_e",cluster_e);
      chain.SetBranchStatus("cluster_eta", 1); chain.SetBranchAddress("cluster_eta",cluster_eta);
      chain.SetBranchStatus("cluster_phi", 1); chain.SetBranchAddress("cluster_phi",cluster_phi);
  } else {
      chain.SetBranchStatus("emcal_clsmult", 1); chain.SetBranchAddress("emcal_clsmult",&clsmult);
      chain.SetBranchStatus("emcal_cluster_e", 1); chain.SetBranchAddress("emcal_cluster_e",cluster_e);
      chain.SetBranchStatus("emcal_cluster_eta", 1); chain.SetBranchAddress("emcal_cluster_eta",cluster_eta);
      chain.SetBranchStatus("emcal_cluster_phi", 1); chain.SetBranchAddress("emcal_cluster_phi",cluster_phi);
  }

  /////////////// JES func ///////////////
  TFile *corrFile = new TFile("JES_Calib_Default.root", "READ");
  if (!corrFile) {
    std::cout << "Error: cannot open JES_Calib_Default.root" << std::endl;
    return;
  }
  TF1 *f_corr = (TF1*)corrFile->Get("JES_Calib_Default_Func");
  if (!f_corr) {
    std::cout << "Error: cannot open f_corr" << std::endl;
    return;
  }

  /////////////// Trigger Efficiency ///////////////
  TFile *f_trigger = new TFile("output_jettrigeff.root", "READ");
  if (!f_trigger) {
    std::cout << "Error: cannot open output_jettrigeff.root" << std::endl;
    return;
  }
  TF1* f_turnon = (TF1*)f_trigger->Get("f_fit");
  TF1* f_turnon_down = (TF1*)f_trigger->Get("f_fit_shiftdown"); // for trigger efficiency uncertainty
  TF1* f_turnon_up = (TF1*)f_trigger->Get("f_fit_shiftup"); // for trigger efficiency uncertainty
  if (!f_turnon || !f_turnon_down || !f_turnon_up) {
    std::cout << "Error: cannot open trigger efficiency functions" << std::endl;
    return;
  }

  /////////////// Histograms ///////////////
  TH1D *h_zvertex = new TH1D("h_zvertex", ";Z-vertex [cm]", 60, -30, 30);
  TH1D *h_recojet_pt_record = new TH1D("h_recojet_pt_record", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_dijet = new TH1D("h_recojet_pt_record_dijet", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_recojet_pt_record_frac = new TH1D("h_recojet_pt_record_frac", ";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_calibjet_pt_record = new TH1D("h_calibjet_pt_record", ";p_{T} [GeV]", 1000, 0, 100);
  TH2D *h_calibjet_pt = new TH2D("h_calibjet_pt", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, calibptbins, calibnet, calibetbins); // raw calib jet pT spectrum
  TH2D *h_calibjet_pt_dijet_eff = new TH2D("h_calibjet_pt_dijet_eff", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, calibptbins, calibnet, calibetbins); // with trigger efficiency + beam background efficiency correction applied
  TH2D *h_calibjet_pt_dijet_effdown = new TH2D("h_calibjet_pt_dijet_effdown", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, calibptbins, calibnet, calibetbins); // for trigger efficiency uncertainty
  TH2D *h_calibjet_pt_dijet_effup = new TH2D("h_calibjet_pt_dijet_effup", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, calibptbins, calibnet, calibetbins); // for trigger efficiency uncertainty
  TH2D *h_calibjet_pt_frac_eff = new TH2D("h_calibjet_pt_frac_eff", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, calibptbins, calibnet, calibetbins); // with trigger efficiency + beam background efficiency correction applied
  TH2D *h_calibjet_pt_frac_effdown = new TH2D("h_calibjet_pt_frac_effdown", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, calibptbins, calibnet, calibetbins); // for trigger efficiency uncertainty
  TH2D *h_calibjet_pt_frac_effup = new TH2D("h_calibjet_pt_frac_effup", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, calibptbins, calibnet, calibetbins); // for trigger efficiency uncertainty

  /////////////// Event Loop ///////////////
  std::cout << "Data analysis started." << std::endl;
  int n_events = chain.GetEntries();
  std::cout << "Total number of events: " << n_events << std::endl;
  // Event variables setup.
  bool background_dijet = false; bool background_frac = false;
  bool jet_filter = false;
  float calibjet_pt, calibjet_eta, calibjet_phi;
  bool calibjet_background_dijet, calibjet_background_frac;
  for (int ie = 0; ie < n_events; ++ie) { // event loop start
    // Load event.
    if (ie % 1000 == 0) {
      std::cout << "Processing event " << ie << "..." << std::endl;
    }
    chain.GetEntry(ie);

    // Trigger and Z-vertex cut.
    if (check_bad_trigger(gl1_trigger_vector_scaled)) continue;
    if (isnan(zvertex)) { continue; }
    if (fabs(zvertex) > 30) continue;

    //////////////////////////// SETUP JET VARIABLES ////////////////////////////

    // Get leading jet and subleading jet. Do basic jet cuts.
    int leadingunsubjet_index = -1;
    int subleadingunsubjet_index = -1;
    get_leading_subleading_jet(leadingunsubjet_index, subleadingunsubjet_index, unsubjet_pt);
    if (leadingunsubjet_index < 0) continue;

    // Fill z-vertex histogram.
    h_zvertex->Fill(zvertex);

    // Beam background cut. Event level recording.
    background_dijet = true;
    background_frac = true;
    
    // Dijet cut.
    bool match_dijet = false;
    if (subleadingunsubjet_index >= 0 && unsubjet_e->at(subleadingunsubjet_index) / (float) unsubjet_e->at(leadingunsubjet_index) > 0.3) {
      match_dijet = match_leading_subleading_jet(unsubjet_phi->at(leadingunsubjet_index), unsubjet_phi->at(subleadingunsubjet_index));
    }
    if (match_dijet) { background_dijet = false; }
    
    // Fraction cut.
    double emfrac = unsubjet_emcal_calo_e->at(leadingunsubjet_index) / (double)(unsubjet_e->at(leadingunsubjet_index));
    double ihfrac = unsubjet_ihcal_calo_e->at(leadingunsubjet_index) / (double)(unsubjet_e->at(leadingunsubjet_index));
    double ohfrac = unsubjet_ohcal_calo_e->at(leadingunsubjet_index) / (double)(unsubjet_e->at(leadingunsubjet_index));
    if (emfrac > 0.1 || emfrac < 0.9 || ohfrac > 0.1 || ohfrac < 0.9 || ihfrac < 0.9) { background_frac = false; }

    float lead_pt = unsubjet_pt->at(leadingunsubjet_index);
    float lead_e = unsubjet_e->at(leadingunsubjet_index);
    float lead_eta = unsubjet_eta->at(leadingunsubjet_index);
    float lead_phi = unsubjet_phi->at(leadingunsubjet_index);

    // Trigger efficiency.
    double jettrigeff = f_turnon->Eval(lead_pt);
    if (jettrigeff < 0.01) jettrigeff = 0.01;
    double jettrig_scale = 1.0 / jettrigeff;

    double jettrigeff_down = f_turnon_down->Eval(lead_pt);
    if (jettrigeff_down < 0.01) jettrigeff_down = 0.01;
    double jettrig_scale_down = 1.0 / jettrigeff_down;

    double jettrigeff_up = f_turnon_up->Eval(lead_pt);
    if (jettrigeff_up < 0.01) jettrigeff_up = 0.01;
    double jettrig_scale_up = 1.0 / jettrigeff_up;

    jet_filter = false;
    filter_jet(jet_filter, lead_e, lead_pt, lead_eta, zvertex, jet_radius);
    if (!jet_filter) h_recojet_pt_record->Fill(lead_pt);
    if (!jet_filter && !background_dijet) h_recojet_pt_record_dijet->Fill(lead_pt);
    if (!jet_filter && !background_frac) h_recojet_pt_record_frac->Fill(lead_pt);
    get_calibjet(calibjet_pt, calibjet_eta, calibjet_phi, calibjet_background_dijet, calibjet_background_frac, jet_filter, lead_pt, lead_eta, lead_phi, background_dijet, background_frac, f_corr);
    
    //////////////////////////// SETUP UE VARIABLES //////////////////////////////

    // find reco ET information 
    float et_transverse = 0;
    if (!clusters) {
        for (int i = 0; i < emcaln; i++) {
            float dphi = get_dphi(lead_phi,emcalphi[i]);
            if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += emcale[i]/cosh(emcaleta[i]); } 
        }
        for (int i = 0; i < ihcaln; i++) {
            float dphi = get_dphi(lead_phi,ihcalphi[i]);
            if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += ihcale[i]/cosh(ihcaleta[i]); } 
        }
        for (int i = 0; i < ohcaln; i++) {
            float dphi = get_dphi(lead_phi,ohcalphi[i]);
            if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += ohcale[i]/cosh(ohcaleta[i]); } 
        }
    } else {
        for (int i = 0; i < clsmult; i++) {
            float dphi = get_dphi(lead_phi,cluster_phi[i]);
            if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += cluster_e[i]/cosh(cluster_eta[i]); }
        }
    }

    //////////////////////////// FILL HISTOGRAMS FOR UNFOLDING ////////////////////////////
    if (calibjet_pt >= calibptbins[0] && et_transverse >= calibetbins[0] && et_transverse <= calibetbins[calibnet]) {
      h_calibjet_pt->Fill(calibjet_pt, et_transverse);
      h_calibjet_pt_record->Fill(calibjet_pt);
      if (!calibjet_background_dijet) {
        h_calibjet_pt_dijet_eff->Fill(calibjet_pt, et_transverse, jettrig_scale);
        h_calibjet_pt_dijet_effdown->Fill(calibjet_pt, et_transverse, jettrig_scale_down);
        h_calibjet_pt_dijet_effup->Fill(calibjet_pt, et_transverse, jettrig_scale_up);
      }
      if (!calibjet_background_frac) {
        h_calibjet_pt_frac_eff->Fill(calibjet_pt, et_transverse, jettrig_scale);
        h_calibjet_pt_frac_effdown->Fill(calibjet_pt, et_transverse, jettrig_scale_down);
        h_calibjet_pt_frac_effup->Fill(calibjet_pt, et_transverse, jettrig_scale_up);
      }
    }
  } // event loop end

  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();
  h_zvertex->Write();
  h_recojet_pt_record->Write();
  h_recojet_pt_record_dijet->Write();
  h_recojet_pt_record_frac->Write();
  h_calibjet_pt->Write();
  h_calibjet_pt_record->Write();
  h_calibjet_pt_dijet_eff->Write();
  h_calibjet_pt_dijet_effdown->Write();
  h_calibjet_pt_dijet_effup->Write();
  h_calibjet_pt_frac_eff->Write();
  h_calibjet_pt_frac_effdown->Write();
  h_calibjet_pt_frac_effup->Write();
  f_out->Close();
  std::cout << "All done!" << std::endl;
}

////////////////////////////////////////// Functions //////////////////////////////////////////
bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled) {
  if (std::find(gl1_trigger_vector_scaled->begin(), gl1_trigger_vector_scaled->end(), 18) != gl1_trigger_vector_scaled->end()) { return true; }
  if (std::find(gl1_trigger_vector_scaled->begin(), gl1_trigger_vector_scaled->end(), 34) != gl1_trigger_vector_scaled->end()) { return true; }
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

bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi) {
  float dijet_min_phi = 3*TMath::Pi()/4.;
  float dphi = get_dphi(leadingjet_phi, subleadingjet_phi);
  return dphi > dijet_min_phi;
}

void get_jet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius) {
  jet_filter.clear();
  int njet = jet_e->size();
  for (int ij = 0; ij < njet; ++ij) {
    jet_filter.push_back(jet_e->at(ij) < 0 || check_bad_jet_eta(jet_eta->at(ij), zvertex, jet_radius) || jet_pt->at(ij) < 1);
  }
}

void filter_jet(bool& jet_filter, float jet_e, float jet_pt, float jet_eta, float zvertex, float jet_radius) {
    jet_filter = (jet_e < 0 || check_bad_jet_eta(jet_eta, zvertex, jet_radius));
}

void get_calibjet(float& calibjet_pt, float& calibjet_eta, float& calibjet_phi, bool& calibjet_background_dijet, bool& calibjet_background_frac, bool jet_filter, float jet_pt, float jet_eta, float jet_phi, bool background_dijet, bool background_frac, TF1* f_corr) {
    calibjet_pt = -9999;
    calibjet_eta = -9999;
    calibjet_phi = -9999;
    calibjet_background_dijet = true;
    calibjet_background_frac = true;
    if (jet_filter) return;
    double calib_pt = f_corr->Eval(jet_pt);
    if (calib_pt < calibptbins[0] || calib_pt > calibptbins[calibnpt]) return;
    calibjet_pt = calib_pt;
    calibjet_eta = jet_eta;
    calibjet_phi = jet_phi;
    calibjet_background_dijet = background_dijet;
    calibjet_background_frac = background_frac;
}