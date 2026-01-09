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
#include "/sphenix/user/hanpuj/CaloDataAna24_skimmed/src/draw_template.C"               

void do_unfolding_closure() {
  SetAtlasStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  //********** General Set up **********//
  const float PI = TMath::Pi();
  const float jet_radius = 0.4;

  std::vector<TH1F*> h_input;
  std::vector<int> color;
  std::vector<int> markerstyle;
  std::vector<std::string> text;
  std::vector<std::string> legend;
 
  //********** Files **********//
  TFile *f_sim = new TFile("output_sim.root", "READ");
  TFile *f_purityefficiency = new TFile("output_purityefficiency.root", "READ");
  TFile *f_out = new TFile("output_closure.root", "RECREATE");

  TH2D* h_respmatrix_calib_all_record = (TH2D*)f_sim->Get("h_respmatrix_calib_frac");
  TH1D* h_truth_calib_all_record = (TH1D*)f_sim->Get("h_truth_calib_frac");
  TH1D* h_measure_calib_all_record = (TH1D*)f_sim->Get("h_measure_calib_frac");

  TH2D* h_respmatrix_calib_frac_half1_record = (TH2D*)f_sim->Get("h_respmatrix_calib_frac_half1");
  TH1D* h_truth_calib_frac_half1_record = (TH1D*)f_sim->Get("h_truth_calib_frac_half1");
  TH1D* h_measure_calib_frac_half1_record = (TH1D*)f_sim->Get("h_measure_calib_frac_half1");
  TH2D* h_respmatrix_calib_frac_half2_record = (TH2D*)f_sim->Get("h_respmatrix_calib_frac_half2");
  TH1D* h_truth_calib_frac_half2_record = (TH1D*)f_sim->Get("h_truth_calib_frac_half2");
  TH1D* h_measure_calib_frac_half2_record = (TH1D*)f_sim->Get("h_measure_calib_frac_half2");

  TH1D* h_purity_calib_all = (TH1D*)f_purityefficiency->Get("h_purity_calib_frac");
  TH1D* h_efficiency_calib_all = (TH1D*)f_purityefficiency->Get("h_efficiency_calib_frac");
  TH1D* h_purity_calib_frac_half1 = (TH1D*)f_purityefficiency->Get("h_purity_calib_frac_half1");
  TH1D* h_efficiency_calib_frac_half1 = (TH1D*)f_purityefficiency->Get("h_efficiency_calib_frac_half1");

  // Full closure test on calib //
  TH1D* h_meas_calib_all = (TH1D*)h_respmatrix_calib_all_record->ProjectionX("h_meas_calib_all");
  TH1D* h_truth_calib_all = (TH1D*)h_respmatrix_calib_all_record->ProjectionY("h_truth_calib_all");
  RooUnfoldResponse* resp_calib_all = new RooUnfoldResponse(h_meas_calib_all, h_truth_calib_all, h_respmatrix_calib_all_record, "resp_calib_all", "Response matrix for calib jet (all events)");
  RooUnfoldBayes unfold_calib_all_fulltest_1(resp_calib_all, h_meas_calib_all, 1);
  TH1D* h_unfold_calib_all_fulltest_1 = (TH1D*)unfold_calib_all_fulltest_1.Hunfold(RooUnfold::kErrors); h_unfold_calib_all_fulltest_1->SetName("h_unfold_calib_all_fulltest_1");
  RooUnfoldBayes unfold_calib_all_fulltest_2(resp_calib_all, h_meas_calib_all, 2);
  TH1D* h_unfold_calib_all_fulltest_2 = (TH1D*)unfold_calib_all_fulltest_2.Hunfold(RooUnfold::kErrors); h_unfold_calib_all_fulltest_2->SetName("h_unfold_calib_all_fulltest_2");
  for (int ib = 1; ib <= h_purity_calib_all->GetNbinsX(); ++ib) {
    h_measure_calib_all_record->SetBinContent(ib, h_measure_calib_all_record->GetBinContent(ib) * h_purity_calib_all->GetBinContent(ib));
    h_measure_calib_all_record->SetBinError(ib, h_measure_calib_all_record->GetBinError(ib) * h_purity_calib_all->GetBinContent(ib));
  }
  for (int ib = 1; ib <= h_efficiency_calib_all->GetNbinsX(); ++ib) {
    h_truth_calib_all_record->SetBinContent(ib, h_truth_calib_all_record->GetBinContent(ib) * h_efficiency_calib_all->GetBinContent(ib));
    h_truth_calib_all_record->SetBinError(ib, h_truth_calib_all_record->GetBinError(ib) * h_efficiency_calib_all->GetBinContent(ib));
  }
  RooUnfoldBayes unfold_calib_all_fulltest_extra_1(resp_calib_all, h_measure_calib_all_record, 1);
  TH1D* h_unfold_calib_all_fulltest_extra_1 = (TH1D*)unfold_calib_all_fulltest_extra_1.Hunfold(RooUnfold::kErrors); h_unfold_calib_all_fulltest_extra_1->SetName("h_unfold_calib_all_fulltest_extra_1");
  RooUnfoldBayes unfold_calib_all_fulltest_extra_2(resp_calib_all, h_measure_calib_all_record, 2);
  TH1D* h_unfold_calib_all_fulltest_extra_2 = (TH1D*)unfold_calib_all_fulltest_extra_2.Hunfold(RooUnfold::kErrors); h_unfold_calib_all_fulltest_extra_2->SetName("h_unfold_calib_all_fulltest_extra_2");

  // Half closure test on calib //
  TH1D* h_meas_calib_frac_half1 = (TH1D*)h_respmatrix_calib_frac_half1_record->ProjectionX("h_meas_calib_frac_half1");
  TH1D* h_truth_calib_frac_half1 = (TH1D*)h_respmatrix_calib_frac_half1_record->ProjectionY("h_truth_calib_frac_half1");
  RooUnfoldResponse* resp_calib_frac_half1 = new RooUnfoldResponse(h_meas_calib_frac_half1, h_truth_calib_frac_half1, h_respmatrix_calib_frac_half1_record, "resp_calib_frac_half1", "Response matrix for calib jet (half1)");
  TH1D* h_meas_calib_frac_half2 = (TH1D*)h_respmatrix_calib_frac_half2_record->ProjectionX("h_meas_calib_frac_half2");
  TH1D* h_truth_calib_frac_half2 = (TH1D*)h_respmatrix_calib_frac_half2_record->ProjectionY("h_truth_calib_frac_half2");
  RooUnfoldBayes unfold_calib_frac_half_1(resp_calib_frac_half1, h_meas_calib_frac_half2, 1);
  TH1D* h_unfold_calib_frac_half_1 = (TH1D*)unfold_calib_frac_half_1.Hunfold(RooUnfold::kErrors); h_unfold_calib_frac_half_1->SetName("h_unfold_calib_frac_half_1");

  // Draw histograms.
  h_input.push_back((TH1F*)h_truth_calib_all);
  h_input.push_back((TH1F*)h_meas_calib_all);
  h_input.push_back((TH1F*)h_unfold_calib_all_fulltest_1);
  color.push_back(kRed);
  color.push_back(kBlack);
  color.push_back(kGreen);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("|z_{vertex}| < 30 cm");
  text.push_back("#eta^{jet} within detector range");
  text.push_back("E^{jet} > 0 GeV");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Reco jet spectrum");
  legend.push_back("Unfolded jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              true, 1e-9, 2e-1, true,
                              true, 0., 2.,
                              true, "p_{T}^{jet} [GeV]", "Arbitrary Unit", "Ratio", 0,
                              false, "Run 21 MB + Jet10GeV + Jet30GeV simulation data set",
                              true, text, 0.55, 0.81, 0.04,
                              true, legend, 0.58, 0.66, 0.04,
                              "figure/unfold_calib_fulltest_1.png");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_truth_calib_all_record);
  h_input.push_back((TH1F*)h_measure_calib_all_record);
  h_input.push_back((TH1F*)h_unfold_calib_all_fulltest_extra_1);
  color.push_back(kRed);
  color.push_back(kBlack);
  color.push_back(kGreen);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("|z_{vertex}| < 30 cm");
  text.push_back("#eta^{jet} within detector range");
  text.push_back("E^{jet} > 0 GeV");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Reco jet spectrum");
  legend.push_back("Unfolded jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              true, 1e-9, 2e-1, true,
                              true, 0., 2.,
                              true, "p_{T}^{jet} [GeV]", "Arbitrary Unit", "Ratio", 0,
                              false, "Run 21 MB + Jet10GeV + Jet30GeV simulation data set",
                              true, text, 0.55, 0.81, 0.04,
                              true, legend, 0.58, 0.66, 0.04,
                              "figure/unfold_calib_fulltest_extra_1.png");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_truth_calib_frac_half1);
  h_input.push_back((TH1F*)h_meas_calib_frac_half2);
  h_input.push_back((TH1F*)h_unfold_calib_frac_half_1);
  color.push_back(kRed);
  color.push_back(kBlack);
  color.push_back(kGreen);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("|z_{vertex}| < 30 cm");
  text.push_back("#eta^{jet} within detector range");
  text.push_back("E^{jet} > 0 GeV");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Reco jet spectrum");
  legend.push_back("Unfolded jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              true, 1e-9, 2e-1, true,
                              true, 0., 2.,
                              true, "p_{T}^{jet} [GeV]", "Arbitrary Unit", "Ratio", 0,
                              false, "Run 21 MB + Jet10GeV + Jet30GeV simulation data set",
                              true, text, 0.55, 0.81, 0.04,
                              true, legend, 0.58, 0.66, 0.04,
                              "figure/unfold_calib_halftest_1.png");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();
  h_unfold_calib_all_fulltest_1->Write();
  h_unfold_calib_all_fulltest_extra_1->Write();
  h_unfold_calib_frac_half_1->Write();
  f_out->Close();
  std::cout << "All done!" << std::endl;
}