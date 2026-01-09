#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_pt);
void get_jet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius);
void get_recojet(std::vector<float>& recojet_pt, std::vector<float>& recojet_eta, std::vector<float>& recojet_phi, std::vector<bool>& jet_filter, std::vector<float>* jet_pt, std::vector<float>* jet_eta, std::vector<float>* jet_phi, float pt_lowthre, float pt_highthre);
void get_calibjet(std::vector<float>& calibjet_pt, std::vector<float>& calibjet_eta, std::vector<float>& calibjet_phi, std::vector<bool>& jet_filter, std::vector<float>* jet_pt, std::vector<float>* jet_eta, std::vector<float>* jet_phi, TF1* f_corr, float pt_lowthre, float pt_highthre);
void get_truthjet(std::vector<float>& goodtruthjet_pt, std::vector<float>& goodtruthjet_eta, std::vector<float>& goodtruthjet_phi, std::vector<bool>& truthjet_filter, std::vector<float>* jet_pt, std::vector<float>* jet_eta, std::vector<float>* jet_phi, float pt_lowthre, float pt_highthre);
void match_meas_truth(std::vector<float>& meas_eta, std::vector<float>& meas_phi, std::vector<float>& meas_matched, std::vector<float>& truth_eta, std::vector<float>& truth_phi, std::vector<float>& truth_matched, float jet_radius);

////////////////////////////////////////// Main Function //////////////////////////////////////////
void analysis_sim(std::string runtype, int nseg, int iseg)  {
  ////////// General Set up //////////
  const float PI = TMath::Pi();
  const float jet_radius = 0.4;

  //int n_MB_events = 10;
  //int n_Jet10GeV_events = 1;
  //int n_Jet30GeV_events = 1;
  double n_MB_events = 21.885464;
  double n_Jet10GeV_events = 2.664378;
  double n_Jet30GeV_events = 2.194879;

  double MB_cross_section = 4.197e-2;
  double Jet10GeV_cross_section = 3.646e-6;
  double Jet30GeV_cross_section = 2.505e-9;

  double MB_scale = MB_cross_section/(double)n_MB_events;
  double Jet10GeV_scale = Jet10GeV_cross_section/(double)n_Jet10GeV_events;
  double Jet30GeV_scale = Jet30GeV_cross_section/(double)n_Jet30GeV_events;

  double weight_scale = 1.0, truthjet_pt_min = 0, truthjet_pt_max = 3000;
  if (runtype == "MB") {
    weight_scale = MB_scale;
    truthjet_pt_min = 0;
    truthjet_pt_max = 14;
  } else if (runtype == "Jet10GeV") {
    weight_scale = Jet10GeV_scale;
    truthjet_pt_min = 14;
    truthjet_pt_max = 30;
  } else if (runtype == "Jet30GeV") {
    weight_scale = Jet30GeV_scale;
    truthjet_pt_min = 30;
    truthjet_pt_max = 3000;
  } else {
    std::cout << "Unknown runtype" << std::endl;
    return;
  }
 
  ////////// Files //////////
  TFile *f_out = new TFile(Form("output_sim/output_%s_%d_%d.root", runtype.c_str(), iseg, iseg+nseg), "RECREATE");
  TChain chain("ttree");
  for (int i = iseg; i < iseg + nseg; ++i) {
    chain.Add(Form("/sphenix/tg/tg01/jets/hanpuj/JES_MC_run21/%s/OutDir%d/output_sim.root", runtype.c_str(), i));
  }
  chain.SetBranchStatus("*", 0);

  float zvertex; chain.SetBranchStatus("z_vertex", 1); chain.SetBranchAddress("z_vertex", &zvertex);

  std::vector<float>* unsubjet_e = nullptr; chain.SetBranchStatus("unsubjet04_e", 1); chain.SetBranchAddress("unsubjet04_e", &unsubjet_e);
  std::vector<float>* unsubjet_pt = nullptr; chain.SetBranchStatus("unsubjet04_pt", 1); chain.SetBranchAddress("unsubjet04_pt", &unsubjet_pt);
  std::vector<float>* unsubjet_eta = nullptr; chain.SetBranchStatus("unsubjet04_eta", 1); chain.SetBranchAddress("unsubjet04_eta", &unsubjet_eta);
  std::vector<float>* unsubjet_phi = nullptr; chain.SetBranchStatus("unsubjet04_phi", 1); chain.SetBranchAddress("unsubjet04_phi", &unsubjet_phi);

  std::vector<float>* truthjet_e = nullptr; chain.SetBranchStatus("truthjet04_e", 1); chain.SetBranchAddress("truthjet04_e", &truthjet_e);
  std::vector<float>* truthjet_pt = nullptr; chain.SetBranchStatus("truthjet04_pt", 1); chain.SetBranchAddress("truthjet04_pt", &truthjet_pt);
  std::vector<float>* truthjet_eta = nullptr; chain.SetBranchStatus("truthjet04_eta", 1); chain.SetBranchAddress("truthjet04_eta", &truthjet_eta);
  std::vector<float>* truthjet_phi = nullptr; chain.SetBranchStatus("truthjet04_phi", 1); chain.SetBranchAddress("truthjet04_phi", &truthjet_phi);

  ////////// JES func //////////
  TFile *corrFile = new TFile("/sphenix/user/hanpuj/JES_MC_Calibration/offline/JES_Calib_Default.root", "READ");
  if (!corrFile) {
    std::cout << "Error: cannot open JES_Calib_Default.root" << std::endl;
    return;
  }
  TF1 *f_corr = (TF1*)corrFile->Get("JES_Calib_Default_Func");
  if (!f_corr) {
    std::cout << "Error: cannot open f_corr" << std::endl;
    return;
  }

  ////////// Histograms //////////
  double recoptbins[] = {12, 15, 18, 21, 24, 28, 32, 36, 40, 45, 50, 55, 60};
  double truthptbins[] = {6, 9, 12, 15, 18, 21, 24, 28, 32, 36, 40, 45, 50, 55, 60, 66, 72};

  int reconpt = sizeof(recoptbins) / sizeof(recoptbins[0]) - 1;
  int truthnpt = sizeof(truthptbins) / sizeof(truthptbins[0]) - 1;

  // Test histograms
  TH2D* h_fake_reco_eta = new TH2D("h_fake_reco_eta", "Reco Jet Fake; p_{T}; #eta; N_{jet}", 11, 6, 72, 120, -1.2, 1.2);
  TH2D* h_fake_calib_eta = new TH2D("h_fake_calib_eta", "Calib Jet Fake; p_{T}; #eta; N_{jet}", 11, 6, 72, 120, -1.2, 1.2);
  TH2D* h_miss_reco_eta = new TH2D("h_miss_reco_eta", "Truth Jet Miss; p_{T}; #eta; N_{jet}", 11, 6, 72, 120, -1.2, 1.2);
  TH2D* h_miss_calib_eta = new TH2D("h_miss_calib_eta", "Truth Jet Miss; p_{T}; #eta; N_{jet}", 11, 6, 72, 120, -1.2, 1.2);

  // General histograms
  TH1D* h_truth_reco_all_set = new TH1D("h_truth_reco_all_set", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH1D* h_measure_reco_all_set = new TH1D("h_measure_reco_all_set", ";p_{T}^{Reco jet} [GeV]", reconpt, recoptbins);
  TH2D* h_respmatrix_reco_all_set = new TH2D("h_respmatrix_reco_all_set", ";p_{T}^{Reco jet} [GeV];p_{T}^{Truth jet} [GeV]", reconpt, recoptbins, truthnpt, truthptbins);
  RooUnfoldResponse* resp_reco_all = new RooUnfoldResponse(h_measure_reco_all_set, h_truth_reco_all_set, h_respmatrix_reco_all_set, "resp_reco_all", "Response matrix for reco jet (all events)");
  TH1D* h_truth_reco_all_record = new TH1D("h_truth_reco_all_record", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH1D* h_measure_reco_all_record = new TH1D("h_measure_reco_all_record", ";p_{T}^{Reco jet} [GeV]", reconpt, recoptbins);
  TH2D* h_respmatrix_reco_all_record = new TH2D("h_respmatrix_reco_all_record", ";p_{T}^{Reco jet} [GeV];p_{T}^{Truth jet} [GeV]", reconpt, recoptbins, truthnpt, truthptbins);
  TH1D* h_fake_reco_all_record = new TH1D("h_fake_reco_all_record", "p_{T}^{Reco jet} [GeV]", reconpt, recoptbins);
  TH1D* h_miss_reco_all_record = new TH1D("h_miss_reco_all_record", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);

  TH1D* h_truth_calib_all_set = new TH1D("h_truth_calib_all_set", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH1D* h_measure_calib_all_set = new TH1D("h_measure_calib_all_set", ";p_{T}^{Calib jet} [GeV]", reconpt, recoptbins);
  TH2D* h_respmatrix_calib_all_set = new TH2D("h_respmatrix_calib_all_set", ";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", reconpt, recoptbins, truthnpt, truthptbins);
  RooUnfoldResponse* resp_calib_all = new RooUnfoldResponse(h_measure_calib_all_set, h_truth_calib_all_set, h_respmatrix_calib_all_set, "resp_calib_all", "Response matrix for calib jet (all events)");
  TH1D* h_truth_calib_all_record = new TH1D("h_truth_calib_all_record", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH1D* h_measure_calib_all_record = new TH1D("h_measure_calib_all_record", ";p_{T}^{Calib jet} [GeV]", reconpt, recoptbins);
  TH2D* h_respmatrix_calib_all_record = new TH2D("h_respmatrix_calib_all_record", ";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", reconpt, recoptbins, truthnpt, truthptbins);
  TH1D* h_fake_calib_all_record = new TH1D("h_fake_calib_all_record", "p_{T}^{Calib jet} [GeV]", reconpt, recoptbins);
  TH1D* h_miss_calib_all_record = new TH1D("h_miss_calib_all_record", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);

  TH1D* h_truth_reco_half_set = new TH1D("h_truth_reco_half_set", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH1D* h_measure_reco_half_set = new TH1D("h_measure_reco_half_set", ";p_{T}^{Reco jet} [GeV]", reconpt, recoptbins);
  TH2D* h_respmatrix_reco_half_set = new TH2D("h_respmatrix_reco_half_set", ";p_{T}^{Reco jet} [GeV];p_{T}^{Truth jet} [GeV]", reconpt, recoptbins, truthnpt, truthptbins);
  RooUnfoldResponse* resp_reco_half = new RooUnfoldResponse(h_measure_reco_half_set, h_truth_reco_half_set, h_respmatrix_reco_half_set, "resp_reco_half", "Response matrix for reco jet (half events)");
  TH1D* h_truth_reco_half_record = new TH1D("h_truth_reco_half_record", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins); // ie/10%2 == 1
  TH1D* h_measure_reco_half_record = new TH1D("h_measure_reco_half_record", ";p_{T}^{Reco jet} [GeV]", reconpt, recoptbins);
  TH2D* h_respmatrix_reco_half_record = new TH2D("h_respmatrix_reco_half_record", ";p_{T}^{Reco jet} [GeV];p_{T}^{Truth jet} [GeV]", reconpt, recoptbins, truthnpt, truthptbins);
  TH1D* h_fake_reco_half_record = new TH1D("h_fake_reco_half_record", "p_{T}^{Reco jet} [GeV]", reconpt, recoptbins);
  TH1D* h_miss_reco_half_record = new TH1D("h_miss_reco_half_record", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH1D* h_truth_reco_half_another = new TH1D("h_truth_reco_half_another", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins); // ie/10%2 == 0
  TH1D* h_measure_reco_half_another = new TH1D("h_measure_reco_half_another", ";p_{T}^{Reco jet} [GeV]", reconpt, recoptbins);

  TH1D* h_truth_calib_half_set = new TH1D("h_truth_calib_half_set", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH1D* h_measure_calib_half_set = new TH1D("h_measure_calib_half_set", ";p_{T}^{Calib jet} [GeV]", reconpt, recoptbins);
  TH2D* h_respmatrix_calib_half_set = new TH2D("h_respmatrix_calib_half_set", ";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", reconpt, recoptbins, truthnpt, truthptbins);
  RooUnfoldResponse* resp_calib_half = new RooUnfoldResponse(h_measure_calib_half_set, h_truth_calib_half_set, h_respmatrix_calib_half_set, "resp_calib_half", "Response matrix for calib jet (half events)");
  TH1D* h_truth_calib_half_record = new TH1D("h_truth_calib_half_record", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins); // ie/10%2 == 1
  TH1D* h_measure_calib_half_record = new TH1D("h_measure_calib_half_record", ";p_{T}^{Calib jet} [GeV]", reconpt, recoptbins);
  TH2D* h_respmatrix_calib_half_record = new TH2D("h_respmatrix_calib_half_record", ";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", reconpt, recoptbins, truthnpt, truthptbins);
  TH1D* h_fake_calib_half_record = new TH1D("h_fake_calib_half_record", "p_{T}^{Calib jet} [GeV]", reconpt, recoptbins);
  TH1D* h_miss_calib_half_record = new TH1D("h_miss_calib_half_record", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins);
  TH1D* h_truth_calib_half_another = new TH1D("h_truth_calib_half_another", "p_{T}^{Truth jet} [GeV]", truthnpt, truthptbins); // ie/10%2 == 0
  TH1D* h_measure_calib_half_another = new TH1D("h_measure_calib_half_another", ";p_{T}^{Calib jet} [GeV]", reconpt, recoptbins);

  std::cout << "xmin = " << resp_reco_all->Hresponse()->GetXaxis()->GetXmin() << std::endl;
  std::cout << "xmax = " << resp_reco_all->Hresponse()->GetXaxis()->GetXmax() << std::endl;
  std::cout << "ymin = " << resp_reco_all->Hresponse()->GetYaxis()->GetXmin() << std::endl;
  std::cout << "ymax = " << resp_reco_all->Hresponse()->GetYaxis()->GetXmax() << std::endl;

  ////////// Event Loop //////////
  std::cout << "Data analysis started." << std::endl;
  int n_events = chain.GetEntries();
  //n_events = 10;
  std::cout << "Total number of events: " << n_events << std::endl;
  // Event variables setup.
  std::vector<bool> jet_filter, truthjet_filter;
  std::vector<float> recojet_pt, recojet_eta, recojet_phi, recojet_matched;
  std::vector<float> calibjet_pt, calibjet_eta, calibjet_phi, calibjet_matched;
  std::vector<float> goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched;
  for (int ie = 0; ie < n_events; ++ie) { // event loop start
    // Load event.
    if (ie % 1000 == 0) {
      std::cout << "Processing event " << ie << "..." << std::endl;
    }
    chain.GetEntry(ie);

    // Z-vertex cut.
    if (fabs(zvertex) > 30) continue;

    int leadingtruthjet_index = -9999;
    get_leading_jet(leadingtruthjet_index, truthjet_pt);
    if (leadingtruthjet_index < 0) continue;
    float leadingtruthjet_pt = truthjet_pt->at(leadingtruthjet_index);
    if (leadingtruthjet_pt < truthjet_pt_min || leadingtruthjet_pt > truthjet_pt_max) continue;

    // Get good reco, calib, truth jet.
    get_jet_filter(jet_filter, unsubjet_e, unsubjet_pt, unsubjet_eta, zvertex, jet_radius);
    get_jet_filter(truthjet_filter, truthjet_e, truthjet_pt, truthjet_eta, zvertex, jet_radius);
    get_recojet(recojet_pt, recojet_eta, recojet_phi, jet_filter, unsubjet_pt, unsubjet_eta, unsubjet_phi, recoptbins[0], recoptbins[reconpt]);
    get_calibjet(calibjet_pt, calibjet_eta, calibjet_phi, jet_filter, unsubjet_pt, unsubjet_eta, unsubjet_phi, f_corr, recoptbins[0], recoptbins[reconpt]);
    get_truthjet(goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi, truthjet_filter, truthjet_pt, truthjet_eta, truthjet_phi, truthptbins[0], truthptbins[truthnpt]);

    // Fill record histograms.
    for (int ij = 0; ij < recojet_pt.size(); ++ij) {
      h_measure_reco_all_record->Fill(recojet_pt[ij], weight_scale);
      if (ie / 10 % 2 == 0) {
        h_measure_reco_half_another->Fill(recojet_pt[ij], weight_scale);
      } else {
        h_measure_reco_half_record->Fill(recojet_pt[ij], weight_scale);
      }
    }
    for (int ij = 0; ij < calibjet_pt.size(); ++ij) {
      h_measure_calib_all_record->Fill(calibjet_pt[ij], weight_scale);
      if (ie / 10 % 2 == 0) {
        h_measure_calib_half_another->Fill(calibjet_pt[ij], weight_scale);
      } else {
        h_measure_calib_half_record->Fill(calibjet_pt[ij], weight_scale);
      }
    }
    for (int ij = 0; ij < goodtruthjet_pt.size(); ++ij) {
      h_truth_reco_all_record->Fill(goodtruthjet_pt[ij], weight_scale);
      h_truth_calib_all_record->Fill(goodtruthjet_pt[ij], weight_scale);
      if (ie / 10 % 2 == 0) {
        h_truth_reco_half_another->Fill(goodtruthjet_pt[ij], weight_scale);
        h_truth_calib_half_another->Fill(goodtruthjet_pt[ij], weight_scale);
      } else {
        h_truth_reco_half_record->Fill(goodtruthjet_pt[ij], weight_scale);
        h_truth_calib_half_record->Fill(goodtruthjet_pt[ij], weight_scale);
      }
    }

    // Fill response matrix.
    match_meas_truth(recojet_eta, recojet_phi, recojet_matched, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_radius);
    for (int ij = 0; ij < recojet_pt.size(); ++ij) {
      if (recojet_matched[ij] < 0) {
        resp_reco_all->Fake(recojet_pt[ij], weight_scale);
        h_fake_reco_all_record->Fill(recojet_pt[ij], weight_scale);
        h_fake_reco_eta->Fill(recojet_pt[ij], recojet_eta[ij]);
        if (recojet_pt[ij] > 40 && recojet_pt[ij] < 45 && recojet_eta[ij] > -0.3 && recojet_eta[ij] < 0.3) {
          std::cout << "Event " << ie << std::endl;
          std::cout << "Z-vertex: " << zvertex << std::endl;
          std::cout << "Recojet_pt :" << std::endl;
          for (int ij = 0; ij < recojet_pt.size(); ++ij) {
            std::cout << ij << "  " << recojet_pt[ij] << "  " << recojet_eta[ij] << "  " << recojet_phi[ij] << "  " << recojet_matched[ij] << std::endl;
          }
          std::cout << "Truthjet_pt :" << std::endl;
          for (int ij = 0; ij < goodtruthjet_pt.size(); ++ij) {
            std::cout << ij << "  " << goodtruthjet_pt[ij] << "  " << goodtruthjet_eta[ij] << "  " << goodtruthjet_phi[ij] << "  " << goodtruthjet_matched[ij] << std::endl;
          }
          std::cout << std::endl;
        }
        if (ie / 10 % 2 == 1) {
          resp_reco_half->Fake(recojet_pt[ij], weight_scale);
          h_fake_reco_half_record->Fill(recojet_pt[ij], weight_scale);
        }
      } else {
        resp_reco_all->Fill(recojet_pt[ij], goodtruthjet_pt[recojet_matched[ij]], weight_scale);
        h_respmatrix_reco_all_record->Fill(recojet_pt[ij], goodtruthjet_pt[recojet_matched[ij]], weight_scale);
        if (ie / 10 % 2 == 1) {
          resp_reco_half->Fill(recojet_pt[ij], goodtruthjet_pt[recojet_matched[ij]], weight_scale);
          h_respmatrix_reco_half_record->Fill(recojet_pt[ij], goodtruthjet_pt[recojet_matched[ij]], weight_scale);
        }
      }
    }
    for (int ij = 0; ij < goodtruthjet_pt.size(); ++ij) {
      if (goodtruthjet_matched[ij] < 0) {
        resp_reco_all->Miss(goodtruthjet_pt[ij], weight_scale);
        h_miss_reco_all_record->Fill(goodtruthjet_pt[ij], weight_scale);
        h_miss_reco_eta->Fill(goodtruthjet_pt[ij], goodtruthjet_eta[ij]);
        if (ie / 10 % 2 == 1) {
          resp_reco_half->Miss(goodtruthjet_pt[ij], weight_scale);
          h_miss_reco_half_record->Fill(goodtruthjet_pt[ij], weight_scale);
        }
      }
    }

    match_meas_truth(calibjet_eta, calibjet_phi, calibjet_matched, goodtruthjet_eta, goodtruthjet_phi, goodtruthjet_matched, jet_radius);
    for (int ij = 0; ij < calibjet_pt.size(); ++ij) {
      if (calibjet_matched[ij] < 0) {
        resp_calib_all->Fake(calibjet_pt[ij], weight_scale);
        h_fake_calib_all_record->Fill(calibjet_pt[ij], weight_scale);
        h_fake_calib_eta->Fill(calibjet_pt[ij], calibjet_eta[ij]);
        if (ie / 10 % 2 == 1) {
          resp_calib_half->Fake(calibjet_pt[ij], weight_scale);
          h_fake_calib_half_record->Fill(calibjet_pt[ij], weight_scale);
        }
      } else {
        resp_calib_all->Fill(calibjet_pt[ij], goodtruthjet_pt[calibjet_matched[ij]], weight_scale);
        h_respmatrix_calib_all_record->Fill(calibjet_pt[ij], goodtruthjet_pt[calibjet_matched[ij]], weight_scale);
        if (ie / 10 % 2 == 1) {
          resp_calib_half->Fill(calibjet_pt[ij], goodtruthjet_pt[calibjet_matched[ij]], weight_scale);
          h_respmatrix_calib_half_record->Fill(calibjet_pt[ij], goodtruthjet_pt[calibjet_matched[ij]], weight_scale);
        }
      }
    }
    for (int ij = 0; ij < goodtruthjet_pt.size(); ++ij) {
      if (goodtruthjet_matched[ij] < 0) {
        resp_calib_all->Miss(goodtruthjet_pt[ij], weight_scale);
        h_miss_calib_all_record->Fill(goodtruthjet_pt[ij], weight_scale);
        h_miss_calib_eta->Fill(goodtruthjet_pt[ij], goodtruthjet_eta[ij]);
        if (ie / 10 % 2 == 1) {
          resp_calib_half->Miss(goodtruthjet_pt[ij], weight_scale);
          h_miss_calib_half_record->Fill(goodtruthjet_pt[ij], weight_scale);
        }
      }
    }
  } // event loop end

  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();

  h_fake_reco_eta->Write();
  h_fake_calib_eta->Write();
  h_miss_reco_eta->Write();
  h_miss_calib_eta->Write();

  h_truth_reco_all_set->Write();
  h_measure_reco_all_set->Write();
  h_respmatrix_reco_all_set->Write();
  resp_reco_all->Write();
  h_truth_reco_all_record->Write();
  h_measure_reco_all_record->Write();
  h_respmatrix_reco_all_record->Write();
  h_fake_reco_all_record->Write();
  h_miss_reco_all_record->Write();

  h_truth_calib_all_set->Write();
  h_measure_calib_all_set->Write();
  h_respmatrix_calib_all_set->Write();
  resp_calib_all->Write();
  h_truth_calib_all_record->Write();
  h_measure_calib_all_record->Write();
  h_respmatrix_calib_all_record->Write();
  h_fake_calib_all_record->Write();
  h_miss_calib_all_record->Write();
  
  h_truth_reco_half_set->Write();
  h_measure_reco_half_set->Write();
  h_respmatrix_reco_half_set->Write();
  resp_reco_half->Write();
  h_truth_reco_half_record->Write();
  h_measure_reco_half_record->Write();
  h_respmatrix_reco_half_record->Write();
  h_fake_reco_half_record->Write();
  h_miss_reco_half_record->Write();
  h_truth_reco_half_another->Write();
  h_measure_reco_half_another->Write();
  
  h_truth_calib_half_set->Write();
  h_measure_calib_half_set->Write();
  h_respmatrix_calib_half_set->Write();
  resp_calib_half->Write();
  h_truth_calib_half_record->Write();
  h_measure_calib_half_record->Write();
  h_respmatrix_calib_half_record->Write();
  h_fake_calib_half_record->Write();
  h_miss_calib_half_record->Write();
  h_truth_calib_half_another->Write();
  h_measure_calib_half_another->Write();

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

bool check_bad_jet_eta(float jet_eta, float zertex, float jet_radius) {
  float emcal_mineta = get_emcal_mineta_zcorrected(zertex);
  float emcal_maxeta = get_emcal_maxeta_zcorrected(zertex);
  float ihcal_mineta = get_ihcal_mineta_zcorrected(zertex);
  float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zertex);
  float ohcal_mineta = get_ohcal_mineta_zcorrected(zertex);
  float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zertex);
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
void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_pt) {
  leadingjet_index = -1;
  float leadingjet_pt = -9999;
  for (int ij = 0; ij < jet_pt->size(); ++ij) {
    float jetpt = jet_pt->at(ij);
    if (jetpt > leadingjet_pt) {
      leadingjet_pt = jetpt;
      leadingjet_index = ij;
    }
  }
}

void get_jet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius) {
  jet_filter.clear();
  int njet = jet_e->size();
  for (int ij = 0; ij < njet; ++ij) {
    jet_filter.push_back(jet_e->at(ij) < 0 || check_bad_jet_eta(jet_eta->at(ij), zvertex, jet_radius));
  }
}

void get_recojet(std::vector<float>& recojet_pt, std::vector<float>& recojet_eta, std::vector<float>& recojet_phi, std::vector<bool>& jet_filter, std::vector<float>* jet_pt, std::vector<float>* jet_eta, std::vector<float>* jet_phi, float pt_lowthre, float pt_highthre) {
  recojet_pt.clear();
  recojet_eta.clear();
  recojet_phi.clear();
  for (int ij = 0; ij < jet_filter.size(); ++ij) {
    if (jet_filter.at(ij)) continue;
    if (jet_pt->at(ij) < pt_lowthre || jet_pt->at(ij) > pt_highthre) continue;
    recojet_pt.push_back(jet_pt->at(ij));
    recojet_eta.push_back(jet_eta->at(ij));
    recojet_phi.push_back(jet_phi->at(ij));
  }
}

void get_calibjet(std::vector<float>& calibjet_pt, std::vector<float>& calibjet_eta, std::vector<float>& calibjet_phi, std::vector<bool>& jet_filter, std::vector<float>* jet_pt, std::vector<float>* jet_eta, std::vector<float>* jet_phi, TF1* f_corr, float pt_lowthre, float pt_highthre) {
  calibjet_pt.clear();
  calibjet_eta.clear();
  calibjet_phi.clear();
  for (int ij = 0; ij < jet_filter.size(); ++ij) {
    if (jet_filter.at(ij)) continue;
    float calib_pt = f_corr->Eval(jet_pt->at(ij));
    if (calib_pt < pt_lowthre || calib_pt > pt_highthre) continue;
    calibjet_pt.push_back(calib_pt);
    calibjet_eta.push_back(jet_eta->at(ij));
    calibjet_phi.push_back(jet_phi->at(ij));
  }
}

void get_truthjet(std::vector<float>& goodtruthjet_pt, std::vector<float>& goodtruthjet_eta, std::vector<float>& goodtruthjet_phi, std::vector<bool>& truthjet_filter, std::vector<float>* jet_pt, std::vector<float>* jet_eta, std::vector<float>* jet_phi, float pt_lowthre, float pt_highthre) {
  goodtruthjet_pt.clear();
  goodtruthjet_eta.clear();
  goodtruthjet_phi.clear();
  for (int ij = 0; ij < truthjet_filter.size(); ++ij) {
    if (truthjet_filter.at(ij)) continue;
    if (jet_pt->at(ij) < pt_lowthre || jet_pt->at(ij) > pt_highthre) continue;
    goodtruthjet_pt.push_back(jet_pt->at(ij));
    goodtruthjet_eta.push_back(jet_eta->at(ij));
    goodtruthjet_phi.push_back(jet_phi->at(ij));
  }
}

void match_meas_truth(std::vector<float>& meas_eta, std::vector<float>& meas_phi, std::vector<float>& meas_matched, std::vector<float>& truth_eta, std::vector<float>& truth_phi, std::vector<float>& truth_matched, float jet_radius) {
  meas_matched.assign(meas_eta.size(), -1);
  truth_matched.assign(truth_eta.size(), -1);
  float max_match_dR = jet_radius * 0.75;
  for (int im = 0; im < meas_eta.size(); ++im) {
    float min_dR = 100;
    int match_index = -9999;
    for (int it = 0; it < truth_eta.size(); ++it) {
      float dR = get_dR(meas_eta[im], meas_phi[im], truth_eta[it], truth_phi[it]);
      if (dR < min_dR) {
        match_index = it;
        min_dR = dR;
      }
    }
    if (min_dR < max_match_dR) {
      if (truth_matched[match_index] == -1) {
        meas_matched[im] = match_index;
        truth_matched[match_index] = im;
      } else {
        float dR1 = get_dR(meas_eta[truth_matched[match_index]], meas_phi[truth_matched[match_index]], truth_eta[match_index], truth_phi[match_index]);
        if (min_dR < dR1) {
          meas_matched[truth_matched[match_index]] = -1;
          meas_matched[im] = match_index;
          truth_matched[match_index] = im;
        }
      }
    }
  }
}
