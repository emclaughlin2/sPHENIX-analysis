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
#include "unfold_Def_hanpuversion.h"

void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et);
bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi);
void get_jet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius);
////////////////////////////////////////// Main Function //////////////////////////////////////////
void analysis_data_backupemmattree(int runnumber, string trig, int nseg, int iseg, double jet_radius)  {

  ////////// General Set up //////////
  double weight_scale = 1.0, truthjet_pt_min = 0, truthjet_pt_max = 3000;
  if (strstr(trig.c_str(), "10")) {
    truthjet_pt_min = 14;
    truthjet_pt_max = 17;
  } else if (strstr(trig.c_str(), "15")) {
    truthjet_pt_min = 17;
    truthjet_pt_max = 22;
  } else if (strstr(trig.c_str(), "20")) {
    truthjet_pt_min = 22;
    truthjet_pt_max = 35;
  } else if (strstr(trig.c_str(), "30")) {
    truthjet_pt_min = 35;
    truthjet_pt_max = 52;
  } else if (strstr(trig.c_str(), "50")) {
    truthjet_pt_min = 52;
    truthjet_pt_max = 71;
  } else if (strstr(trig.c_str(), "70")) {
    truthjet_pt_min = 71;
    truthjet_pt_max = 100;
  }

  /////////////// General Set up ///////////////
  int jet_radius_index = (int)(10 * jet_radius);
  TFile *f_out = new TFile(Form("output/emma_%s_output_r0%d_%d_%d_%d.root", trig.c_str(), jet_radius_index, runnumber, nseg, iseg),"RECREATE");
  //TFile *f_out = new TFile(Form("emma_%s_output_r0%d_%d_%d_%d.root", trig.c_str(), jet_radius_index, runnumber, nseg, iseg),"RECREATE");

  /////////////// Read Files ///////////////
  const char* baseDirJet = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput";
  TChain chain("T");
  //chain.Add(Form("%s/output_ana509_v2_%d_*.root", baseDirJet, runnumber));
  for (int i = nseg; i < iseg; i++) {
    chain.Add(Form("%s/sim_run%d_%s_3sigma_output_%d.root", baseDirJet, runnumber, trig.c_str(), i));
  }
  chain.SetBranchStatus("*", 0);
  float zvertex; chain.SetBranchStatus("zvtx", 1); chain.SetBranchAddress("zvtx", &zvertex);
  std::vector<int>* gl1_trigger_vector_scaled = nullptr; chain.SetBranchStatus("triggerVector", 1); chain.SetBranchAddress("triggerVector", &gl1_trigger_vector_scaled);
  std::vector<float>* unsubjet_e = nullptr; chain.SetBranchStatus("e", 1); chain.SetBranchAddress("e", &unsubjet_e);
  std::vector<float>* unsubjet_pt = nullptr; chain.SetBranchStatus("pt", 1); chain.SetBranchAddress("pt", &unsubjet_pt);
  std::vector<float>* unsubjet_eta = nullptr; chain.SetBranchStatus("eta", 1); chain.SetBranchAddress("eta", &unsubjet_eta);
  std::vector<float>* unsubjet_phi = nullptr; chain.SetBranchStatus("phi", 1); chain.SetBranchAddress("phi", &unsubjet_phi);

  std::vector<float>* truthjet_e = nullptr; chain.SetBranchStatus("truthE", 1); chain.SetBranchAddress("truthE", &truthjet_e);
  std::vector<float>* truthjet_pt = nullptr; chain.SetBranchStatus("truthPt", 1); chain.SetBranchAddress("truthPt", &truthjet_pt);
  std::vector<float>* truthjet_eta = nullptr; chain.SetBranchStatus("truthEta", 1); chain.SetBranchAddress("truthEta", &truthjet_eta);
  std::vector<float>* truthjet_phi = nullptr; chain.SetBranchStatus("truthPhi", 1); chain.SetBranchAddress("truthPhi", &truthjet_phi);
  
  /////////////// Histograms ///////////////
  TH1D* h_leadingunsubjet_pt = new TH1D("h_leadingunsubjet_pt", ";Leading Uncalibrated Jet p_{T} [GeV]", 1000, 0, 100);
  TH1D* h_leadingtruthjet_pt = new TH1D("h_leadingtruthjet_pt", ";Leading Truth Jet p_{T} [GeV]", 1000, 0, 100);
  TH1D* h_zvertex = new TH1D("h_zvertex",";z_vertex [cm]", 400, -100, 100);

  /////////////// Event Loop ///////////////
  int event = 0;
  std::cout << "Data analysis started." << std::endl;
  Long64_t n_events = chain.GetEntries();
  std::cout << "Total number of events: " << n_events << std::endl;
  std::vector<bool> jet_filter, truthjet_filter;
  for (Long64_t ie = 0; ie < n_events; ++ie) { // event loop start
    chain.GetEntry(ie);

    // Get Z-vertex and trigger check.
    if (runnumber > 30) {
      if (std::find(gl1_trigger_vector_scaled->begin(), gl1_trigger_vector_scaled->end(), 22) == gl1_trigger_vector_scaled->end()) { continue; } 
    }
    if (zvertex > 60 || zvertex < -60) continue; 

    // Get leading and subleading jets.
    int leadjet_index = -1;
    float leadjet_e = -9999;
    for (int ij = 0; ij < unsubjet_pt->size(); ++ij) {
      float jete = unsubjet_pt->at(ij);
      if (fabs(unsubjet_eta->at(ij)) > 2.5) { continue; }
      if (jete > leadjet_e) {
        leadjet_e = jete;
        leadjet_index = ij;
      }
    }

    if (leadjet_index < 0) continue;
    if (unsubjet_pt->at(leadjet_index) < 5.0) { continue; }

    // Get leading and subleading jets.
    int leadingtruthjet_index = -1;
    float leadingtruthjet_e = -9999;
    for (int ij = 0; ij < truthjet_e->size(); ++ij) {
      float jete = truthjet_e->at(ij);
      if (fabs(truthjet_eta->at(ij)) > 2.0) { continue; }
      if (truthjet_pt->at(ij) <= 5) { continue; }
      if (jete > leadingtruthjet_e) {
        leadingtruthjet_e = jete;
        leadingtruthjet_index = ij;
      }
    }
    if (leadingtruthjet_index >= 0) {
      if (truthjet_pt->at(leadingtruthjet_index) >= truthjet_pt_min && truthjet_pt->at(leadingtruthjet_index) <= truthjet_pt_max) {
        get_jet_filter(truthjet_filter, truthjet_e, truthjet_pt, truthjet_eta, zvertex, jet_radius);
        // Get leading and subleading jets.
        leadingtruthjet_index = -1;
        leadingtruthjet_e = -9999;
        for (int ij = 0; ij < truthjet_e->size(); ++ij) {
          float jete = truthjet_e->at(ij); 
          if (truthjet_filter.at(ij)) { continue; }
          if (truthjet_pt->at(ij) <= 5) { continue; }
          if (jete > leadingtruthjet_e) {
            leadingtruthjet_e = jete;
            leadingtruthjet_index = ij;
          }
        }
        if (leadingtruthjet_index >= 0 && truthjet_pt->at(leadingtruthjet_index) > 8) { 
          h_leadingtruthjet_pt->Fill(truthjet_pt->at(leadingtruthjet_index)); 
          std::cout << "zvertex: " << zvertex << " truth leading jet pt: " << truthjet_pt->at(leadingtruthjet_index) << std::endl;
          // fill in with debug if necessary
        }
      }
    }

    // 0.7 eta cut + acceptance cut.
    get_jet_filter(jet_filter, unsubjet_e, unsubjet_pt, unsubjet_eta, zvertex, jet_radius);

    // Get leading and subleading jets.
    int leadingunsubjet_index = -1;
    int subleadingunsubjet_index = -1;
    float leadingunsubjet_e = -9999;
    float subleadingunsubjet_e = -9999;
    for (int ij = 0; ij < unsubjet_e->size(); ++ij) {
      if (jet_filter.at(ij)) continue;
      if (unsubjet_e->at(ij) < 0) { continue; }
      if (unsubjet_pt->at(ij) <= 2) { continue; }
      float jete = unsubjet_e->at(ij);
      if (jete > leadingunsubjet_e) {
        subleadingunsubjet_e = leadingunsubjet_e;
        subleadingunsubjet_index = leadingunsubjet_index;
        leadingunsubjet_e = jete;
        leadingunsubjet_index = ij;
      } else if (jete > subleadingunsubjet_e) {
        subleadingunsubjet_e = jete;
        subleadingunsubjet_index = ij;
      }
    }
    
    if (leadingunsubjet_index < 0 || subleadingunsubjet_index < 0) continue;
    if (unsubjet_pt->at(leadingunsubjet_index) < 8) { continue; }

    // Dijet cut
    bool match_dijet = false;
    if (unsubjet_e->at(subleadingunsubjet_index) / (float) unsubjet_e->at(leadingunsubjet_index) > 0.3) {
      match_dijet = match_leading_subleading_jet(unsubjet_phi->at(leadingunsubjet_index), unsubjet_phi->at(subleadingunsubjet_index));
    }
    if (!match_dijet) continue;
    event++;
    std::cout << "event: " << event << " zvertex: " << zvertex << " leading jet pt: " << unsubjet_pt->at(leadingunsubjet_index) << std::endl;
    h_leadingunsubjet_pt->Fill(unsubjet_pt->at(leadingunsubjet_index));
    h_zvertex->Fill(zvertex);
  } // event loop end

  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();
  h_leadingunsubjet_pt->Write();
  h_leadingtruthjet_pt->Write();
  h_zvertex->Write();
  f_out->Close();
  std::cout << "All done!" << std::endl;
}

////////////////////////////////////////// Helper Functions //////////////////////////////////////////
float get_deta(float eta1, float eta2) {
  return eta1 - eta2;
}

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

float get_dR(float eta1, float phi1, float eta2, float phi2) {
  float deta = get_deta(eta1, eta2);
  float dphi = get_dphi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}

int get_hcal_up(int iphi) {
  if (iphi == 63) {
    return 0;
  } else {
    return iphi + 1;
  }
}

int get_hcal_bottom(int iphi) {
  if (iphi == 0) {
    return 63;
  } else {
    return iphi - 1;
  }
}

float get_emcal_mineta_zcorrected(float zvertex) {
  float minz_EM = -130.23;
  float radius_EM = 93.5;
  float z = minz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_emcal_maxeta_zcorrected(float zvertex) {
  float maxz_EM = 130.23;
  float radius_EM = 93.5;
  float z = maxz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_ihcal_mineta_zcorrected(float zvertex) {
  float minz_IH = -170.299;
  float radius_IH = 127.503;
  float z = minz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ihcal_maxeta_zcorrected(float zvertex) {
  float maxz_IH = 170.299;
  float radius_IH = 127.503;
  float z = maxz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ohcal_mineta_zcorrected(float zvertex) {
  float minz_OH = -301.683;
  float radius_OH = 225.87;
  float z = minz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

float get_ohcal_maxeta_zcorrected(float zvertex) {
  float maxz_OH = 301.683;
  float radius_OH = 225.87;
  float z = maxz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

bool check_bad_jet_eta(float jet_eta, float zvertex, float jet_radius) {
  float emcal_mineta = get_emcal_mineta_zcorrected(zvertex);
  float emcal_maxeta = get_emcal_maxeta_zcorrected(zvertex);
  float ihcal_mineta = get_ihcal_mineta_zcorrected(zvertex);
  float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zvertex);
  float ohcal_mineta = get_ohcal_mineta_zcorrected(zvertex);
  float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zvertex);
  float minlimit = emcal_mineta;
  if (ihcal_mineta > minlimit) minlimit = ihcal_mineta;
  if (ohcal_mineta > minlimit) minlimit = ohcal_mineta;
  float maxlimit = emcal_maxeta;
  if (ihcal_maxeta < maxlimit) maxlimit = ihcal_maxeta;
  if (ohcal_maxeta < maxlimit) maxlimit = ohcal_maxeta;
  minlimit += jet_radius;
  maxlimit -= jet_radius;
  return jet_eta < minlimit || jet_eta > maxlimit;
}

////////////////////////////////////////// Functions //////////////////////////////////////////
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
  float dphi = fabs(get_dphi(leadingjet_phi, subleadingjet_phi));
  return dphi > dijet_min_phi;
}

void get_jet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius) {
  jet_filter.clear();
  int njet = jet_e->size();
  for (int ij = 0; ij < njet; ++ij) {
    jet_filter.push_back(check_bad_jet_eta(jet_eta->at(ij), zvertex, jet_radius) || jet_eta->at(ij) > 1.1 - jet_radius || jet_eta->at(ij) < -1.1 + jet_radius);
  }
}