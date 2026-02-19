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

void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et);
bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi);
float doCalibration(const std::vector<std::vector<std::vector<TF1 *>>> &f_corr, float jet_radius, float jetPt, float zvrtx, float eta);

////////////////////////////////////////// Main Function //////////////////////////////////////////
void ppg09_jettrig_analysis(int runnumber)  {
  /////////////// General Set up ///////////////
  TFile *f_out = new TFile(Form("test_analysis_output/trig10_zvtx_lt_60_output_%d.root", runnumber),"RECREATE");

  /////////////// Read Files ///////////////
  const char* baseDirJet = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput";
  TChain chain("T");
  chain.Add(Form("%s/output_ana509_jet_trig_%d_*.root", baseDirJet, runnumber));
  chain.SetBranchStatus("*", 0);

  float zvertex; chain.SetBranchStatus("zvtx", 1); chain.SetBranchAddress("zvtx", &zvertex);
  std::vector<float>* unsubjet04_e = nullptr; chain.SetBranchStatus("e", 1); chain.SetBranchAddress("e", &unsubjet04_e);
  std::vector<float>* unsubjet04_pt = nullptr; chain.SetBranchStatus("pt", 1); chain.SetBranchAddress("pt", &unsubjet04_pt);
  std::vector<float>* unsubjet04_eta = nullptr; chain.SetBranchStatus("eta", 1); chain.SetBranchAddress("eta", &unsubjet04_eta);
  std::vector<float>* unsubjet04_phi = nullptr; chain.SetBranchStatus("phi", 1); chain.SetBranchAddress("phi", &unsubjet04_phi);
  std::vector<float>* unsubjet04_time = nullptr; chain.SetBranchStatus("jettime", 1); chain.SetBranchAddress("jettime", &unsubjet04_time);
  vector<int>* triggerVector = nullptr; chain.SetBranchStatus("triggerVector", 1); chain.SetBranchAddress("triggerVector",&triggerVector);
  vector<int>* liveTriggerVector = nullptr; chain.SetBranchStatus("liveTriggerVector", 1); chain.SetBranchAddress("liveTriggerVector",&liveTriggerVector);
  
  /////////////// Histograms ///////////////
  const int n_trig = 64;
  std::vector<int> check_trig_index = {10, 18, 22};
  TH1F *h_leadingjet04pt_all[n_trig], *h_leadingjet04e_all[n_trig], *h_leadingjet04calibpt_all[n_trig];
  TH1F *h_jet04pt_all[n_trig], *h_jet04e_all[n_trig], *h_jet04calibpt_all[n_trig];

  for (int it = 0; it < check_trig_index.size(); ++it) {
    int trig_index = check_trig_index[it];
    h_leadingjet04pt_all[trig_index] = new TH1F(Form("h_leadingjet04pt_all_%d", trig_index), "; p_{T} (GeV/c); Counts", 1000, 0, 100);
    h_leadingjet04e_all[trig_index] = new TH1F(Form("h_leadingjet04e_all_%d", trig_index), "; E (GeV); Counts", 1000, 0, 100);
    h_jet04pt_all[trig_index] = new TH1F(Form("h_jet04pt_all_%d", trig_index), "; p_{T} (GeV/c); Counts", 1000, 0, 100);
    h_jet04e_all[trig_index] = new TH1F(Form("h_jet04e_all_%d", trig_index), "; E (GeV); Counts", 1000, 0, 100);
  }

  /////////////// Event Loop ///////////////
  std::cout << "Data analysis started." << std::endl;
  int n_events = chain.GetEntries();
  std::cout << "Total number of events: " << n_events << std::endl;
  // Event variables setup.
  for (int ie = 0; ie < n_events; ++ie) { // event loop start
    // Load event.
    if (ie % 1000 == 0) {
      std::cout << "Processing event " << ie << "..." << std::endl;
    }
    chain.GetEntry(ie);
    // MBD >= 1 trigger check.
    if (find(triggerVector->begin(), triggerVector->end(), 10) == triggerVector->end()) { continue; }
    if (fabs(zvertex) > 60) { continue; }

    bool pass_dijet = false;
    bool pass_timing = false;
    bool pass_njet = false;

    int leadingunsubjet04_index = -1;
    int subleadingunsubjet04_index = -1;
    get_leading_subleading_jet(leadingunsubjet04_index, subleadingunsubjet04_index, unsubjet04_e);
    if (leadingunsubjet04_index < 0) continue;
    if (subleadingunsubjet04_index >= 0 && unsubjet04_e->at(subleadingunsubjet04_index) / (float) unsubjet04_e->at(leadingunsubjet04_index) > 0.3) {
      pass_dijet = match_leading_subleading_jet(unsubjet04_phi->at(leadingunsubjet04_index), unsubjet04_phi->at(subleadingunsubjet04_index));
    }
    if (subleadingunsubjet04_index >= 0 && fabs(unsubjet04_time->at(leadingunsubjet04_index)*17.6+2) < 6 && fabs(unsubjet04_time->at(leadingunsubjet04_index)*17.6 - unsubjet04_time->at(subleadingunsubjet04_index)*17.6) < 3) {
      pass_timing = true;
    }
    int n_jets = 0;
    for (int ij = 0; ij < unsubjet04_e->size(); ++ij) {
      if (unsubjet04_e->at(ij) > 5.0) {
        n_jets++;
      }
    }
    if (n_jets <= 9) {
      pass_njet = true;
    }

    float leading_jet04_e = 0; if(leadingunsubjet04_index >= 0) leading_jet04_e = unsubjet04_e->at(leadingunsubjet04_index);
    get_leading_subleading_jet(leadingunsubjet04_index, subleadingunsubjet04_index, unsubjet04_pt);

    float leading_jet04_pt = 0, leading_jet04_calibpt = 0;
    if(leadingunsubjet04_index >= 0) {
      leading_jet04_pt = unsubjet04_pt->at(leadingunsubjet04_index);
    }

    // Fill histograms for each trigger.
    for (int it = 0; it < check_trig_index.size(); ++it) {
      int trig_index = check_trig_index[it];
      if (find(liveTriggerVector->begin(), liveTriggerVector->end(), trig_index) == liveTriggerVector->end()) { continue; }

      if (pass_dijet && pass_timing && pass_njet) {
        h_leadingjet04pt_all[trig_index]->Fill(leading_jet04_pt);
        h_leadingjet04e_all[trig_index]->Fill(leading_jet04_e);
        for (int ij = 0; ij < unsubjet04_e->size(); ++ij) {
          h_jet04pt_all[trig_index]->Fill(unsubjet04_pt->at(ij));
          h_jet04e_all[trig_index]->Fill(unsubjet04_e->at(ij));
        }
      }
    }
  } // event loop end

  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();
  for (int it = 0; it < check_trig_index.size(); ++it) {
    int trig_index = check_trig_index[it];
    h_leadingjet04pt_all[trig_index]->Write();
    h_leadingjet04e_all[trig_index]->Write();
    h_jet04pt_all[trig_index]->Write();
    h_jet04e_all[trig_index]->Write();
  }
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
  //float dijet_min_phi = 7*TMath::Pi()/8.;
  float dphi = fabs(get_dphi(leadingjet_phi, subleadingjet_phi));
  return dphi > dijet_min_phi;
}

int getZvrtxBin(float zvrtx) {
  if (zvrtx >= -60.0 && zvrtx < -30.0) {
    return 1;  // -60 to -30
  } else if (zvrtx >= -30.0 && zvrtx < 30.0) {
    return 0;  // -30 to 30
  } else if (zvrtx >= 30.0 && zvrtx < 60) {
    return 2;  // 30 to 60
  } else if ((zvrtx < -60.0 && zvrtx >= -900) || (zvrtx >= 60.0 && zvrtx < 900)) {
    return 3;  // -inf to -60 or 60 to inf
  } else if (zvrtx < -900) {
    return 4;  // -999 no zvertex
  }
  return 0;  // Default case, should not happen
}

int getEtaBin(int zvrtxbin, float eta, float jet_radius) {
  float eta_low = 0;
  float eta_high = 0;
  if (zvrtxbin == 0) { // -30 < zvrtx < 30
    eta_low = -1.2 + jet_radius;
    eta_high = 1.2 - jet_radius;
  } else if (zvrtxbin == 1) { // -60 < zvrtx < -30
    eta_low = -0.95 + jet_radius;
    eta_high = 1.25 - jet_radius;
  } else if (zvrtxbin == 2) { // 30 < zvrtx < inf
    eta_low = -1.25 + jet_radius;
    eta_high = 0.95 - jet_radius;
  } else if (zvrtxbin == 3 || zvrtxbin == 4) {
    return 0;
  }
  float threshold1 = eta_low + ((eta_high - eta_low) / 4.0);
  float threshold2 = eta_low + ((eta_high - eta_low) / 2.0);
  float threshold3 = eta_low + (3.0 * (eta_high - eta_low) / 4.0);

  if (eta < threshold1) return 0;  // -inf to threshold1
  if (eta < threshold2) return 1;  // threshold1 to threshold2
  if (eta < threshold3) return 2;  // threshold2 to threshold3
  return 3;  // threshold3 to inf
}

float doCalibration(const std::vector<std::vector<std::vector<TF1 *>>> &f_corr, float jet_radius, float jetPt, float zvrtx, float eta) {
  float calib = jetPt;
  int zvrtxbin = 0;
  int etabin = 0;
  zvrtxbin = getZvrtxBin(zvrtx);
  if (zvrtxbin < 3) {
    etabin = getEtaBin(zvrtxbin, eta, jet_radius);
  } else {
    etabin = 0;
  }
  int radiusbin = (static_cast<int>(jet_radius * 10)) - 2;
  calib = f_corr[radiusbin][zvrtxbin][etabin]->Eval(jetPt);
  return calib;
}
