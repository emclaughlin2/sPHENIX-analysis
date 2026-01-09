#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <vector>
#include "TLatex.h"
#include <string>
#include "unfold_Def.h"
#include "/sphenix/user/hanpuj/CaloDataAna24_skimmed/src/draw_template.C"               

ifstream f_luminosity("luminosity.txt");
double luminosity;
f_luminosity >> luminosity;
double delta_eta = 1.4;
int n_underflow_bin = 2;
int n_overflow_bin = 1;
std::cout << "Luminosity: " << luminosity << " mb^-1" << std::endl;
std::cout << "Delta eta: " << delta_eta << std::endl;

void draw_1D_multiple_plot_uncertainty(std::vector<TH1F*> h_input, float x_min, float x_max, float ratio_min, float ratio_max, std::string xtitle, std::string ytitle, std::vector<std::string> legend, float xstart_legend, float ystart_legend, std::vector<int> color, std::vector<int> markerstyle, std::string output_name);
void set_tpad(TPad* pad, TH2D* h, TBox* box1, TBox* box2);

void get_finalspectrum() {

  std::vector<TH1F*> h_input;
  std::vector<TGraphAsymmErrors*> g_input;
  std::vector<int> color;
  std::vector<int> markerstyle;
  std::vector<std::string> text;
  std::vector<std::string> legend;

  TFile *f_data = new TFile("output_data.root", "READ");
  TH1D* h_recojet_pt_record = (TH1D*)f_data->Get("h_recojet_pt_record");
  h_recojet_pt_record->Rebin(10);
  TH1D* h_recojet_pt_record_dijet = (TH1D*)f_data->Get("h_recojet_pt_record_dijet");
  h_recojet_pt_record_dijet->Rebin(10);
  TH1D* h_recojet_pt_record_frac = (TH1D*)f_data->Get("h_recojet_pt_record_frac");
  h_recojet_pt_record_frac->Rebin(10);
  
  h_input.push_back((TH1F*)h_recojet_pt_record);
  h_input.push_back((TH1F*)h_recojet_pt_record_dijet);
  h_input.push_back((TH1F*)h_recojet_pt_record_frac);
  color.push_back(kBlack);
  color.push_back(kAzure+1);
  color.push_back(kRed+1);
  markerstyle.push_back(20);
  markerstyle.push_back(21);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Preliminary");
  text.push_back("p+p#sqrt{s} = 200 GeV,#kern[+0.7]{#scale[0.5]{#int}}#kern[+0.2]{#it{L}}dt = 16.6 pb^{-1}");
  text.push_back("anti-k_{t} #kern[-0.3]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  legend.push_back("Reco jet spectrum");
  legend.push_back("After dijet cut");
  legend.push_back("After fraction cut");
  draw_1D_multiple_plot(h_input, color, markerstyle,
                        false, 10, true,
                        true, 10, 100, false,
                        true, 0.3, 4e+06, true,
                        true, "p_{T}^{jet} [GeV]", "dN_{jet} / dp_{T} [1/GeV]",
                        true, text, 0.32, 0.9, 0.05,
                        true, legend, 0.48, 0.65, 0.05,
                        "figure/raw_reco_jet_spectrum.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  f_data->Close();
  std::cout << "check0" << std::endl;

  TFile* f_spectrum = new TFile("output_unfolded.root", "READ");

  TH1D* h_unfold_calib_dijet_reweighted_1 = new TH1D("h_spectra_calib_dijet_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_reweighted_1");
  TH1D* h_unfold_calib_dijet_effdown_reweighted_1 = new TH1D("h_spectra_calib_dijet_effdown_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_effdown_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_effdown_reweighted_1");
  TH1D* h_unfold_calib_dijet_effup_reweighted_1 = new TH1D("h_spectra_calib_dijet_effup_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_effup_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_effup_reweighted_1");
  TH1D* h_unfold_calib_dijet_jesdown_reweighted_1 = new TH1D("h_spectra_calib_dijet_jesdown_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_jesdown_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_jesdown_reweighted_1");
  TH1D* h_unfold_calib_dijet_jesup_reweighted_1 = new TH1D("h_spectra_calib_dijet_jesup_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_jesup_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_jesup_reweighted_1");
  TH1D* h_unfold_calib_dijet_jerdown_reweighted_1 = new TH1D("h_spectra_calib_dijet_jerdown_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_jerdown_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_jerdown_reweighted_1");
  TH1D* h_unfold_calib_dijet_jerup_reweighted_1 = new TH1D("h_spectra_calib_dijet_jerup_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_jerup_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_jerup_reweighted_1");

  TH1D* h_unfold_calib_frac_reweighted_1 = new TH1D("h_spectra_calib_frac_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_reweighted_1");
  TH1D* h_unfold_calib_frac_effdown_reweighted_1 = new TH1D("h_spectra_calib_frac_effdown_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_effdown_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_effdown_reweighted_1");
  TH1D* h_unfold_calib_frac_effup_reweighted_1 = new TH1D("h_spectra_calib_frac_effup_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_effup_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_effup_reweighted_1");
  TH1D* h_unfold_calib_frac_jesdown_reweighted_1 = new TH1D("h_spectra_calib_frac_jesdown_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_jesdown_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_jesdown_reweighted_1");
  TH1D* h_unfold_calib_frac_jesup_reweighted_1 = new TH1D("h_spectra_calib_frac_jesup_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_jesup_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_jesup_reweighted_1");
  TH1D* h_unfold_calib_frac_jerdown_reweighted_1 = new TH1D("h_spectra_calib_frac_jerdown_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_jerdown_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_jerdown_reweighted_1");
  TH1D* h_unfold_calib_frac_jerup_reweighted_1 = new TH1D("h_spectra_calib_frac_jerup_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_jerup_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_jerup_reweighted_1");

  TH1D* h_unfold_calib_dijet_reweighted_2 = new TH1D("h_spectra_calib_dijet_reweighted_2", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_reweighted_2_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_reweighted_2");
  TH1D* h_unfold_calib_dijet_reweighted_3 = new TH1D("h_spectra_calib_dijet_reweighted_3", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_reweighted_3_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_reweighted_3");
  TH1D* h_unfold_calib_dijet_mbddown_reweighted_1 = new TH1D("h_spectra_calib_dijet_mbddown_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_mbddown_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_mbddown_reweighted_1");
  TH1D* h_unfold_calib_dijet_mbdup_reweighted_1 = new TH1D("h_spectra_calib_dijet_mbdup_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_mbdup_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_mbdup_reweighted_1");

  TH1D* h_unfold_calib_frac_reweighted_2 = new TH1D("h_spectra_calib_frac_reweighted_2", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_reweighted_2_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_reweighted_2");
  TH1D* h_unfold_calib_frac_reweighted_3 = new TH1D("h_spectra_calib_frac_reweighted_3", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_reweighted_3_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_reweighted_3");
  TH1D* h_unfold_calib_frac_mbddown_reweighted_1 = new TH1D("h_spectra_calib_frac_mbddown_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_mbddown_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_mbddown_reweighted_1");
  TH1D* h_unfold_calib_frac_mbdup_reweighted_1 = new TH1D("h_spectra_calib_frac_mbdup_reweighted_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_mbdup_reweighted_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_mbdup_reweighted_1");

  TH1D* h_unfold_calib_dijet_1 = new TH1D("h_spectra_calib_dijet_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_dijet_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_dijet_1");
  TH1D* h_unfold_calib_frac_1 = new TH1D("h_spectra_calib_frac_1", "", calibnpt, calibptbins); TH1D* h_unfold_calib_frac_1_temp = (TH1D*)f_spectrum->Get("h_unfold_calib_frac_1");

  for (int ib = 1; ib <= h_unfold_calib_dijet_reweighted_1->GetNbinsX(); ++ib) {
    h_unfold_calib_dijet_reweighted_1->SetBinContent(ib, h_unfold_calib_dijet_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_reweighted_1->SetBinError(ib, h_unfold_calib_dijet_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_dijet_effdown_reweighted_1->SetBinContent(ib, h_unfold_calib_dijet_effdown_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_effdown_reweighted_1->SetBinError(ib, h_unfold_calib_dijet_effdown_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_dijet_effup_reweighted_1->SetBinContent(ib, h_unfold_calib_dijet_effup_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_effup_reweighted_1->SetBinError(ib, h_unfold_calib_dijet_effup_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_dijet_jesdown_reweighted_1->SetBinContent(ib, h_unfold_calib_dijet_jesdown_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_jesdown_reweighted_1->SetBinError(ib, h_unfold_calib_dijet_jesdown_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_dijet_jesup_reweighted_1->SetBinContent(ib, h_unfold_calib_dijet_jesup_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_jesup_reweighted_1->SetBinError(ib, h_unfold_calib_dijet_jesup_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_dijet_jerdown_reweighted_1->SetBinContent(ib, h_unfold_calib_dijet_jerdown_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_jerdown_reweighted_1->SetBinError(ib, h_unfold_calib_dijet_jerdown_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_dijet_jerup_reweighted_1->SetBinContent(ib, h_unfold_calib_dijet_jerup_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_jerup_reweighted_1->SetBinError(ib, h_unfold_calib_dijet_jerup_reweighted_1_temp->GetBinError(ib + n_underflow_bin));

    h_unfold_calib_frac_reweighted_1->SetBinContent(ib, h_unfold_calib_frac_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_reweighted_1->SetBinError(ib, h_unfold_calib_frac_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_frac_effdown_reweighted_1->SetBinContent(ib, h_unfold_calib_frac_effdown_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_effdown_reweighted_1->SetBinError(ib, h_unfold_calib_frac_effdown_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_frac_effup_reweighted_1->SetBinContent(ib, h_unfold_calib_frac_effup_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_effup_reweighted_1->SetBinError(ib, h_unfold_calib_frac_effup_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_frac_jesdown_reweighted_1->SetBinContent(ib, h_unfold_calib_frac_jesdown_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_jesdown_reweighted_1->SetBinError(ib, h_unfold_calib_frac_jesdown_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_frac_jesup_reweighted_1->SetBinContent(ib, h_unfold_calib_frac_jesup_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_jesup_reweighted_1->SetBinError(ib, h_unfold_calib_frac_jesup_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_frac_jerdown_reweighted_1->SetBinContent(ib, h_unfold_calib_frac_jerdown_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_jerdown_reweighted_1->SetBinError(ib, h_unfold_calib_frac_jerdown_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_frac_jerup_reweighted_1->SetBinContent(ib, h_unfold_calib_frac_jerup_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_jerup_reweighted_1->SetBinError(ib, h_unfold_calib_frac_jerup_reweighted_1_temp->GetBinError(ib + n_underflow_bin));

    h_unfold_calib_dijet_reweighted_2->SetBinContent(ib, h_unfold_calib_dijet_reweighted_2_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_reweighted_2->SetBinError(ib, h_unfold_calib_dijet_reweighted_2_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_dijet_reweighted_3->SetBinContent(ib, h_unfold_calib_dijet_reweighted_3_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_reweighted_3->SetBinError(ib, h_unfold_calib_dijet_reweighted_3_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_dijet_mbddown_reweighted_1->SetBinContent(ib, h_unfold_calib_dijet_mbddown_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_mbddown_reweighted_1->SetBinError(ib, h_unfold_calib_dijet_mbddown_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_dijet_mbdup_reweighted_1->SetBinContent(ib, h_unfold_calib_dijet_mbdup_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_mbdup_reweighted_1->SetBinError(ib, h_unfold_calib_dijet_mbdup_reweighted_1_temp->GetBinError(ib + n_underflow_bin));

    h_unfold_calib_frac_reweighted_2->SetBinContent(ib, h_unfold_calib_frac_reweighted_2_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_reweighted_2->SetBinError(ib, h_unfold_calib_frac_reweighted_2_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_frac_reweighted_3->SetBinContent(ib, h_unfold_calib_frac_reweighted_3_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_reweighted_3->SetBinError(ib, h_unfold_calib_frac_reweighted_3_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_frac_mbddown_reweighted_1->SetBinContent(ib, h_unfold_calib_frac_mbddown_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_mbddown_reweighted_1->SetBinError(ib, h_unfold_calib_frac_mbddown_reweighted_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_frac_mbdup_reweighted_1->SetBinContent(ib, h_unfold_calib_frac_mbdup_reweighted_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_mbdup_reweighted_1->SetBinError(ib, h_unfold_calib_frac_mbdup_reweighted_1_temp->GetBinError(ib + n_underflow_bin));

    h_unfold_calib_dijet_1->SetBinContent(ib, h_unfold_calib_dijet_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_dijet_1->SetBinError(ib, h_unfold_calib_dijet_1_temp->GetBinError(ib + n_underflow_bin));
    h_unfold_calib_frac_1->SetBinContent(ib, h_unfold_calib_frac_1_temp->GetBinContent(ib + n_underflow_bin)); h_unfold_calib_frac_1->SetBinError(ib, h_unfold_calib_frac_1_temp->GetBinError(ib + n_underflow_bin));
  }

  // Get Uncertainty
  double uncertainty_unfold_upper[9] = {0};
  double uncertainty_unfold_lower[9] = {0};
  double uncertainty_cut_upper[9] = {0};
  double uncertainty_cut_lower[9] = {0};
  double uncertainty_mbd_upper[9] = {0}; // in ratio
  double uncertainty_mbd_lower[9] = {0}; // in ratio
  double uncertainty_other_upper[9] = {0};
  double uncertainty_other_lower[9] = {0};
  double uncertainty_final_upper[9] = {0};
  double uncertainty_final_lower[9] = {0};

  double uncertainty_unfold_upper_record[9] = {0}, uncertainty_unfold_lower_record[9] = {0};
  double uncertainty_cut_upper_record[9] = {0}, uncertainty_cut_lower_record[9] = {0};
  double uncertainty_jes_upper_record[9] = {0}, uncertainty_jes_lower_record[9] = {0};
  double uncertainty_jer_upper_record[9] = {0}, uncertainty_jer_lower_record[9] = {0};
  double uncertainty_jeteff_upper_record[9] = {0}, uncertainty_jeteff_lower_record[9] = {0};
  double uncertainty_mbdeff_upper_record[9] = {0}, uncertainty_mbdeff_lower_record[9] = {0};
  double uncertainty_mbdcs_upper_record[9] = {0}, uncertainty_mbdcs_lower_record[9] = {0};
  double uncertainty_stat_upper_record[9] = {0}, uncertainty_stat_lower_record[9] = {0};
  double uncertainty_total_upper_record[9] = {0}, uncertainty_total_lower_record[9] = {0};

  TH1D* h_unfold_calib_nominal = new TH1D("h_unfold_calib_nominal", "", calibnpt, calibptbins);
  for (int ib = 1; ib <= h_unfold_calib_nominal->GetNbinsX(); ++ib) {
    h_unfold_calib_nominal->SetBinContent(ib, (h_unfold_calib_dijet_reweighted_1->GetBinContent(ib) + h_unfold_calib_frac_reweighted_1->GetBinContent(ib)) / 2.);
    h_unfold_calib_nominal->SetBinError(ib, TMath::Sqrt(TMath::Power(h_unfold_calib_dijet_reweighted_1->GetBinError(ib), 2) + TMath::Power(h_unfold_calib_frac_reweighted_1->GetBinError(ib), 2)) / 2.);
  }

  TH1D* h_unfold_calib_nominal_unfold = new TH1D("h_unfold_calib_nominal_unfold", "", calibnpt, calibptbins);
  for (int ib = 1; ib <= h_unfold_calib_nominal_unfold->GetNbinsX(); ++ib) {
    h_unfold_calib_nominal_unfold->SetBinContent(ib, (h_unfold_calib_dijet_1->GetBinContent(ib) + h_unfold_calib_frac_1->GetBinContent(ib)) / 2.);
    h_unfold_calib_nominal_unfold->SetBinError(ib, TMath::Sqrt(TMath::Power(h_unfold_calib_dijet_1->GetBinError(ib), 2) + TMath::Power(h_unfold_calib_frac_1->GetBinError(ib), 2)) / 2.);
  }

  // MBD Cross Section Uncertainty
  for (int ib = 0; ib < 9; ++ib) {
    uncertainty_mbd_upper[ib] = 2.3 / 25.2;
    uncertainty_mbd_lower[ib] = 1.7 / 25.2;
    uncertainty_mbdcs_upper_record[ib] = uncertainty_mbd_upper[ib];
    uncertainty_mbdcs_lower_record[ib] = uncertainty_mbd_lower[ib];
  }

  // Unfolding Uncertainty
  h_input.push_back((TH1F*)h_unfold_calib_nominal);
  h_input.push_back((TH1F*)h_unfold_calib_nominal_unfold);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  legend.push_back("Nominal distribution");
  legend.push_back("Nominal distribution without reweighting");
  draw_1D_multiple_plot_uncertainty(h_input, calibptbins[0], calibptbins[calibnpt], 0.7, 1.3, "p_{T}^{jet} [GeV]", "d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]", legend, 0.22, 0.2, color, markerstyle, "figure/specturm_uncertainty_unfold.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  legend.clear();
  std::cout << std::endl << "Uncertainty from unfolding: " << std::endl << "up\t\tdown" << std::endl;
  for (int ib = 1; ib <= h_unfold_calib_nominal->GetNbinsX(); ++ib) {
    double unc = h_unfold_calib_nominal_unfold->GetBinContent(ib) - h_unfold_calib_nominal->GetBinContent(ib);
    double unc_up = 0;
    double unc_down = 0;
    unc_up = TMath::Sqrt(unc_up * unc_up + unc * unc);
    unc_down = TMath::Sqrt(unc_down * unc_down + unc * unc);

    uncertainty_unfold_upper[ib-1] = TMath::Sqrt(uncertainty_unfold_upper[ib-1] * uncertainty_unfold_upper[ib-1] + unc * unc);
    uncertainty_unfold_lower[ib-1] = TMath::Sqrt(uncertainty_unfold_lower[ib-1] * uncertainty_unfold_lower[ib-1] + unc * unc);

    uncertainty_unfold_upper_record[ib-1] = unc_up / h_unfold_calib_nominal->GetBinContent(ib);
    uncertainty_unfold_lower_record[ib-1] = unc_down / h_unfold_calib_nominal->GetBinContent(ib);
    std::cout << unc_up / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%      \t" << unc_down / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%" << std::endl;
  }

  // Beam Background Cut Uncertainty
  h_input.push_back((TH1F*)h_unfold_calib_nominal);
  h_input.push_back((TH1F*)h_unfold_calib_dijet_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_frac_reweighted_1);
  color.push_back(kBlack);
  color.push_back(kAzure+1);
  color.push_back(kRed+1);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  legend.push_back("Nominal distribution");
  legend.push_back("With dijet cut");
  legend.push_back("With fraction cut");
  draw_1D_multiple_plot_uncertainty(h_input, calibptbins[0], calibptbins[calibnpt], 0.8, 1.2, "p_{T}^{jet} [GeV]", "d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]", legend, 0.22, 0.26, color, markerstyle, "figure/specturm_uncertainty_cut.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  legend.clear();
  std::cout << std::endl << "Uncertainty from cut: " << std::endl << "up\t\tdown" << std::endl;
  for (int ib = 1; ib <= h_unfold_calib_nominal->GetNbinsX(); ++ib) {
    double unc1 = h_unfold_calib_dijet_reweighted_1->GetBinContent(ib) - h_unfold_calib_nominal->GetBinContent(ib);
    double unc2 = h_unfold_calib_frac_reweighted_1->GetBinContent(ib) - h_unfold_calib_nominal->GetBinContent(ib);
    double unc_up = 0;
    double unc_down = 0;
    if (unc1 > 0) {
      uncertainty_cut_upper[ib-1] = TMath::Sqrt(uncertainty_cut_upper[ib-1] * uncertainty_cut_upper[ib-1] + unc1 * unc1);
      unc_up = TMath::Sqrt(unc_up * unc_up + unc1 * unc1);
    } else {
      uncertainty_cut_lower[ib-1] = TMath::Sqrt(uncertainty_cut_lower[ib-1] * uncertainty_cut_lower[ib-1] + unc1 * unc1);
      unc_down = TMath::Sqrt(unc_down * unc_down + unc1 * unc1);
    }
    if (unc2 > 0) {
      uncertainty_cut_upper[ib-1] = TMath::Sqrt(uncertainty_cut_upper[ib-1] * uncertainty_cut_upper[ib-1] + unc2 * unc2);
      unc_up = TMath::Sqrt(unc_up * unc_up + unc2 * unc2);
    } else {
      uncertainty_cut_lower[ib-1] = TMath::Sqrt(uncertainty_cut_lower[ib-1] * uncertainty_cut_lower[ib-1] + unc2 * unc2);
      unc_down = TMath::Sqrt(unc_down * unc_down + unc2 * unc2);
    }
    uncertainty_cut_upper_record[ib-1] = unc_up / h_unfold_calib_nominal->GetBinContent(ib);
    uncertainty_cut_lower_record[ib-1] = unc_down / h_unfold_calib_nominal->GetBinContent(ib);
    std::cout << unc_up / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%      \t" << unc_down / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%" << std::endl;
  }

  // Jet Trigger Efficiency Uncertainty
  h_input.push_back((TH1F*)h_unfold_calib_dijet_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_dijet_effdown_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_dijet_effup_reweighted_1);
  color.push_back(kBlack);
  color.push_back(kRed);
  color.push_back(kBlue);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  legend.push_back("Unfolded spectrum (Dijet cut)");
  legend.push_back("With jet trigger efficiency shift down");
  legend.push_back("With jet trigger efficiency shift up");
  draw_1D_multiple_plot_uncertainty(h_input, calibptbins[0], calibptbins[calibnpt], 0.8, 1.2, "p_{T}^{jet} [GeV]", "d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]", legend, 0.22, 0.26, color, markerstyle, "figure/specturm_uncertainty_eff_dijet.pdf");
  h_input.clear();
  legend.clear();
  h_input.push_back((TH1F*)h_unfold_calib_frac_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_frac_effdown_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_frac_effup_reweighted_1);
  legend.push_back("Unfolded spectrum (Fraction cut)");
  legend.push_back("With jet trigger efficiency shift down");
  legend.push_back("With jet trigger efficiency shift up");
  draw_1D_multiple_plot_uncertainty(h_input, calibptbins[0], calibptbins[calibnpt], 0.8, 1.2, "p_{T}^{jet} [GeV]", "d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]", legend, 0.22, 0.26, color, markerstyle, "figure/specturm_uncertainty_eff_frac.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  legend.clear();
  std::cout << std::endl << "Uncertainty from jet trigger efficiency: " << std::endl << "up\t\t\tdown\t\t\tunc1\t\t\tunc2\t\t\tdijet1\t\t\tdijet2\t\t\tfrac1\t\t\tfrac2" << std::endl;
  for (int ib = 1; ib <= h_unfold_calib_nominal->GetNbinsX(); ++ib) {
    double unc1_dijet = h_unfold_calib_dijet_effdown_reweighted_1->GetBinContent(ib) - h_unfold_calib_dijet_reweighted_1->GetBinContent(ib);
    double unc2_dijet = h_unfold_calib_dijet_effup_reweighted_1->GetBinContent(ib) - h_unfold_calib_dijet_reweighted_1->GetBinContent(ib);
    double unc1_frac = h_unfold_calib_frac_effdown_reweighted_1->GetBinContent(ib) - h_unfold_calib_frac_reweighted_1->GetBinContent(ib);
    double unc2_frac = h_unfold_calib_frac_effup_reweighted_1->GetBinContent(ib) - h_unfold_calib_frac_reweighted_1->GetBinContent(ib);
    double unc1 = (unc1_dijet + unc1_frac) / 2.;
    double unc2 = (unc2_dijet + unc2_frac) / 2.;
    if (unc1 > 0) {
      uncertainty_other_upper[ib-1] = TMath::Sqrt(uncertainty_other_upper[ib-1] * uncertainty_other_upper[ib-1] + unc1 * unc1);
      uncertainty_jeteff_upper_record[ib-1] = TMath::Sqrt(uncertainty_jeteff_upper_record[ib-1] * uncertainty_jeteff_upper_record[ib-1] + unc1 * unc1);
    } else {
      uncertainty_other_lower[ib-1] = TMath::Sqrt(uncertainty_other_lower[ib-1] * uncertainty_other_lower[ib-1] + unc1 * unc1);
      uncertainty_jeteff_lower_record[ib-1] = TMath::Sqrt(uncertainty_jeteff_lower_record[ib-1] * uncertainty_jeteff_lower_record[ib-1] + unc1 * unc1);
    }
    if (unc2 > 0) {
      uncertainty_other_upper[ib-1] = TMath::Sqrt(uncertainty_other_upper[ib-1] * uncertainty_other_upper[ib-1] + unc2 * unc2);
      uncertainty_jeteff_upper_record[ib-1] = TMath::Sqrt(uncertainty_jeteff_upper_record[ib-1] * uncertainty_jeteff_upper_record[ib-1] + unc2 * unc2);
    } else {
      uncertainty_other_lower[ib-1] = TMath::Sqrt(uncertainty_other_lower[ib-1] * uncertainty_other_lower[ib-1] + unc2 * unc2);
      uncertainty_jeteff_lower_record[ib-1] = TMath::Sqrt(uncertainty_jeteff_lower_record[ib-1] * uncertainty_jeteff_lower_record[ib-1] + unc2 * unc2);
    }
    uncertainty_jeteff_upper_record[ib-1] = uncertainty_jeteff_upper_record[ib-1] / h_unfold_calib_nominal->GetBinContent(ib);
    uncertainty_jeteff_lower_record[ib-1] = uncertainty_jeteff_lower_record[ib-1] / h_unfold_calib_nominal->GetBinContent(ib);
    std::cout << uncertainty_other_upper[ib-1] << " (" << uncertainty_other_upper[ib-1] / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t" << uncertainty_other_lower[ib-1] << " (" << uncertainty_other_lower[ib-1] / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1 << " (" << unc1 / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t" << unc2 << " (" << unc2 / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1_dijet << " (" << unc1_dijet / h_unfold_calib_dijet_reweighted_1->GetBinContent(ib) * 100 << "\%) \t" << unc2_dijet << " (" << unc2_dijet / h_unfold_calib_dijet_reweighted_1->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1_frac << " (" << unc1_frac / h_unfold_calib_frac_reweighted_1->GetBinContent(ib) * 100 << "\%) \t" << unc2_frac << " (" << unc2_frac / h_unfold_calib_frac_reweighted_1->GetBinContent(ib) * 100 << "\%)" << std::endl;
  }

  // JES Uncertainty
  h_input.push_back((TH1F*)h_unfold_calib_dijet_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_dijet_jesdown_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_dijet_jesup_reweighted_1);
  color.push_back(kBlack);
  color.push_back(kRed);
  color.push_back(kBlue);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  legend.push_back("Unfolded spectrum (Dijet cut)");
  legend.push_back("With JES shift down");
  legend.push_back("With JES shift up");
  draw_1D_multiple_plot_uncertainty(h_input, calibptbins[0], calibptbins[calibnpt], 0, 2.43, "p_{T}^{jet} [GeV]", "d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]", legend, 0.22, 0.26, color, markerstyle, "figure/specturm_uncertainty_jes_dijet.pdf");
  h_input.clear();
  legend.clear();
  h_input.push_back((TH1F*)h_unfold_calib_frac_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_frac_jesdown_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_frac_jesup_reweighted_1);
  legend.push_back("Unfolded spectrum (Fraction cut)");
  legend.push_back("With JES shift down");
  legend.push_back("With JES shift up");
  draw_1D_multiple_plot_uncertainty(h_input, calibptbins[0], calibptbins[calibnpt], 0, 2.43, "p_{T}^{jet} [GeV]", "d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]", legend, 0.22, 0.26, color, markerstyle, "figure/specturm_uncertainty_jes_frac.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  legend.clear();
  std::cout << std::endl << "Uncertainty from JES: " << std::endl << "up\t\t\tdown\t\t\tunc1\t\t\tunc2\t\t\tdijet1\t\t\tdijet2\t\t\tfrac1\t\t\tfrac2" << std::endl;
  for (int ib = 1; ib <= h_unfold_calib_nominal->GetNbinsX(); ++ib) {
    double unc1_dijet = h_unfold_calib_dijet_jesdown_reweighted_1->GetBinContent(ib) - h_unfold_calib_dijet_reweighted_1->GetBinContent(ib);
    double unc2_dijet = h_unfold_calib_dijet_jesup_reweighted_1->GetBinContent(ib) - h_unfold_calib_dijet_reweighted_1->GetBinContent(ib);
    double unc1_frac = h_unfold_calib_frac_jesdown_reweighted_1->GetBinContent(ib) - h_unfold_calib_frac_reweighted_1->GetBinContent(ib);
    double unc2_frac = h_unfold_calib_frac_jesup_reweighted_1->GetBinContent(ib) - h_unfold_calib_frac_reweighted_1->GetBinContent(ib);
    double unc1 = (unc1_dijet + unc1_frac) / 2.;
    double unc2 = (unc2_dijet + unc2_frac) / 2.;
    if (unc1 > 0) {
      uncertainty_other_upper[ib-1] = TMath::Sqrt(uncertainty_other_upper[ib-1] * uncertainty_other_upper[ib-1] + unc1 * unc1);
      uncertainty_jes_upper_record[ib-1] = TMath::Sqrt(uncertainty_jes_upper_record[ib-1] * uncertainty_jes_upper_record[ib-1] + unc1 * unc1);
    } else {
      uncertainty_other_lower[ib-1] = TMath::Sqrt(uncertainty_other_lower[ib-1] * uncertainty_other_lower[ib-1] + unc1 * unc1);
      uncertainty_jes_lower_record[ib-1] = TMath::Sqrt(uncertainty_jes_lower_record[ib-1] * uncertainty_jes_lower_record[ib-1] + unc1 * unc1);
    }
    if (unc2 > 0) {
      uncertainty_other_upper[ib-1] = TMath::Sqrt(uncertainty_other_upper[ib-1] * uncertainty_other_upper[ib-1] + unc2 * unc2);
      uncertainty_jes_upper_record[ib-1] = TMath::Sqrt(uncertainty_jes_upper_record[ib-1] * uncertainty_jes_upper_record[ib-1] + unc2 * unc2);
    } else {
      uncertainty_other_lower[ib-1] = TMath::Sqrt(uncertainty_other_lower[ib-1] * uncertainty_other_lower[ib-1] + unc2 * unc2);
      uncertainty_jes_lower_record[ib-1] = TMath::Sqrt(uncertainty_jes_lower_record[ib-1] * uncertainty_jes_lower_record[ib-1] + unc2 * unc2);
    }
    uncertainty_jes_upper_record[ib-1] = uncertainty_jes_upper_record[ib-1] / h_unfold_calib_nominal->GetBinContent(ib);
    uncertainty_jes_lower_record[ib-1] = uncertainty_jes_lower_record[ib-1] / h_unfold_calib_nominal->GetBinContent(ib);
    std::cout << uncertainty_other_upper[ib-1] << " (" << uncertainty_other_upper[ib-1] / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t" << uncertainty_other_lower[ib-1] << " (" << uncertainty_other_lower[ib-1] / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1 << " (" << unc1 / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t" << unc2 << " (" << unc2 / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1_dijet << " (" << unc1_dijet / h_unfold_calib_dijet_reweighted_1->GetBinContent(ib) * 100 << "\%) \t" << unc2_dijet << " (" << unc2_dijet / h_unfold_calib_dijet_reweighted_1->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1_frac << " (" << unc1_frac / h_unfold_calib_frac_reweighted_1->GetBinContent(ib) * 100 << "\%) \t" << unc2_frac << " (" << unc2_frac / h_unfold_calib_frac_reweighted_1->GetBinContent(ib) * 100 << "\%)" << std::endl;
  }

  // JER Uncertainty
  h_input.push_back((TH1F*)h_unfold_calib_dijet_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_dijet_jerdown_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_dijet_jerup_reweighted_1);
  color.push_back(kBlack);
  color.push_back(kRed);
  color.push_back(kBlue);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  legend.push_back("Unfolded spectrum (Dijet cut)");
  legend.push_back("With JER shift down");
  legend.push_back("With JER shift up");
  draw_1D_multiple_plot_uncertainty(h_input, calibptbins[0], calibptbins[calibnpt], 0, 2, "p_{T}^{jet} [GeV]", "d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]", legend, 0.22, 0.26, color, markerstyle, "figure/specturm_uncertainty_jer_dijet.pdf");
  h_input.clear();
  legend.clear();
  h_input.push_back((TH1F*)h_unfold_calib_frac_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_frac_jerdown_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_frac_jerup_reweighted_1);
  legend.push_back("Unfolded spectrum (Fraction cut)");
  legend.push_back("With JER shift down");
  legend.push_back("With JER shift up");
  draw_1D_multiple_plot_uncertainty(h_input, calibptbins[0], calibptbins[calibnpt], 0, 2, "p_{T}^{jet} [GeV]", "d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]", legend, 0.22, 0.26, color, markerstyle, "figure/specturm_uncertainty_jer_frac.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  legend.clear();
  std::cout << std::endl << "Uncertainty from JER: " << std::endl << "up\t\t\tdown\t\t\tunc1\t\t\tunc2\t\t\tdijet1\t\t\tdijet2\t\t\tfrac1\t\t\tfrac2" << std::endl;
  for (int ib = 1; ib <= h_unfold_calib_nominal->GetNbinsX(); ++ib) {
    double unc1_dijet = h_unfold_calib_dijet_jerdown_reweighted_1->GetBinContent(ib) - h_unfold_calib_dijet_reweighted_1->GetBinContent(ib);
    double unc2_dijet = h_unfold_calib_dijet_jerup_reweighted_1->GetBinContent(ib) - h_unfold_calib_dijet_reweighted_1->GetBinContent(ib);
    double unc1_frac = h_unfold_calib_frac_jerdown_reweighted_1->GetBinContent(ib) - h_unfold_calib_frac_reweighted_1->GetBinContent(ib);
    double unc2_frac = h_unfold_calib_frac_jerup_reweighted_1->GetBinContent(ib) - h_unfold_calib_frac_reweighted_1->GetBinContent(ib);
    double unc1 = (unc1_dijet + unc1_frac) / 2.;
    double unc2 = (unc2_dijet + unc2_frac) / 2.;
    if (unc1 > 0) {
      uncertainty_other_upper[ib-1] = TMath::Sqrt(uncertainty_other_upper[ib-1] * uncertainty_other_upper[ib-1] + unc1 * unc1);
      uncertainty_jer_upper_record[ib-1] = TMath::Sqrt(uncertainty_jer_upper_record[ib-1] * uncertainty_jer_upper_record[ib-1] + unc1 * unc1);
    } else {
      uncertainty_other_lower[ib-1] = TMath::Sqrt(uncertainty_other_lower[ib-1] * uncertainty_other_lower[ib-1] + unc1 * unc1);
      uncertainty_jer_lower_record[ib-1] = TMath::Sqrt(uncertainty_jer_lower_record[ib-1] * uncertainty_jer_lower_record[ib-1] + unc1 * unc1);
    }
    if (unc2 > 0) {
      uncertainty_other_upper[ib-1] = TMath::Sqrt(uncertainty_other_upper[ib-1] * uncertainty_other_upper[ib-1] + unc2 * unc2);
      uncertainty_jer_upper_record[ib-1] = TMath::Sqrt(uncertainty_jer_upper_record[ib-1] * uncertainty_jer_upper_record[ib-1] + unc2 * unc2);
    } else {
      uncertainty_other_lower[ib-1] = TMath::Sqrt(uncertainty_other_lower[ib-1] * uncertainty_other_lower[ib-1] + unc2 * unc2);
      uncertainty_jer_lower_record[ib-1] = TMath::Sqrt(uncertainty_jer_lower_record[ib-1] * uncertainty_jer_lower_record[ib-1] + unc2 * unc2);
    }
    uncertainty_jer_upper_record[ib-1] = uncertainty_jer_upper_record[ib-1] / h_unfold_calib_nominal->GetBinContent(ib);
    uncertainty_jer_lower_record[ib-1] = uncertainty_jer_lower_record[ib-1] / h_unfold_calib_nominal->GetBinContent(ib);
    std::cout << uncertainty_other_upper[ib-1] << " (" << uncertainty_other_upper[ib-1] / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t" << uncertainty_other_lower[ib-1] << " (" << uncertainty_other_lower[ib-1] / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1 << " (" << unc1 / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t" << unc2 << " (" << unc2 / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1_dijet << " (" << unc1_dijet / h_unfold_calib_dijet_reweighted_1->GetBinContent(ib) * 100 << "\%) \t" << unc2_dijet << " (" << unc2_dijet / h_unfold_calib_dijet_reweighted_1->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1_frac << " (" << unc1_frac / h_unfold_calib_frac_reweighted_1->GetBinContent(ib) * 100 << "\%) \t" << unc2_frac << " (" << unc2_frac / h_unfold_calib_frac_reweighted_1->GetBinContent(ib) * 100 << "\%)" << std::endl;
  }

  // MBD Trigger Efficiency Uncertainty
  h_input.push_back((TH1F*)h_unfold_calib_dijet_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_dijet_mbddown_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_dijet_mbdup_reweighted_1);
  color.push_back(kBlack);
  color.push_back(kRed);
  color.push_back(kBlue);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  legend.push_back("Unfolded spectrum (Dijet cut)");
  legend.push_back("With MBD trigger efficiency shift down");
  legend.push_back("With MBD trigger efficiency shift up");
  draw_1D_multiple_plot_uncertainty(h_input, calibptbins[0], calibptbins[calibnpt], 0.8, 1.2, "p_{T}^{jet} [GeV]", "d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]", legend, 0.22, 0.26, color, markerstyle, "figure/specturm_uncertainty_mbd_dijet.pdf");
  h_input.clear();
  legend.clear();
  h_input.push_back((TH1F*)h_unfold_calib_frac_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_frac_mbddown_reweighted_1);
  h_input.push_back((TH1F*)h_unfold_calib_frac_mbdup_reweighted_1);
  legend.push_back("Unfolded spectrum (Fraction cut)");
  legend.push_back("With MBD trigger efficiency shift down");
  legend.push_back("With MBD trigger efficiency shift up");
  draw_1D_multiple_plot_uncertainty(h_input, calibptbins[0], calibptbins[calibnpt], 0.8, 1.2, "p_{T}^{jet} [GeV]", "d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]", legend, 0.22, 0.26, color, markerstyle, "figure/specturm_uncertainty_mbd_frac.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  legend.clear();
  std::cout << std::endl << "Uncertainty from MBD trigger efficiency: " << std::endl << "up\t\t\tdown\t\t\tunc1\t\t\tunc2\t\t\tdijet1\t\t\tdijet2\t\t\tfrac1\t\t\tfrac2" << std::endl;
  for (int ib = 1; ib <= h_unfold_calib_nominal->GetNbinsX(); ++ib) {
    double unc1_dijet = h_unfold_calib_dijet_mbddown_reweighted_1->GetBinContent(ib) - h_unfold_calib_dijet_reweighted_1->GetBinContent(ib);
    double unc2_dijet = h_unfold_calib_dijet_mbdup_reweighted_1->GetBinContent(ib) - h_unfold_calib_dijet_reweighted_1->GetBinContent(ib);
    double unc1_frac = h_unfold_calib_frac_mbddown_reweighted_1->GetBinContent(ib) - h_unfold_calib_frac_reweighted_1->GetBinContent(ib);
    double unc2_frac = h_unfold_calib_frac_mbdup_reweighted_1->GetBinContent(ib) - h_unfold_calib_frac_reweighted_1->GetBinContent(ib);
    double unc1 = (unc1_dijet + unc1_frac) / 2.;
    double unc2 = (unc2_dijet + unc2_frac) / 2.;
    if (unc1 > 0) {
      uncertainty_other_upper[ib-1] = TMath::Sqrt(uncertainty_other_upper[ib-1] * uncertainty_other_upper[ib-1] + unc1 * unc1);
      uncertainty_mbdeff_upper_record[ib-1] = TMath::Sqrt(uncertainty_mbdeff_upper_record[ib-1] * uncertainty_mbdeff_upper_record[ib-1] + unc1 * unc1);
    } else {
      uncertainty_other_lower[ib-1] = TMath::Sqrt(uncertainty_other_lower[ib-1] * uncertainty_other_lower[ib-1] + unc1 * unc1);
      uncertainty_mbdeff_lower_record[ib-1] = TMath::Sqrt(uncertainty_mbdeff_lower_record[ib-1] * uncertainty_mbdeff_lower_record[ib-1] + unc1 * unc1);
    }
    if (unc2 > 0) {
      uncertainty_other_upper[ib-1] = TMath::Sqrt(uncertainty_other_upper[ib-1] * uncertainty_other_upper[ib-1] + unc2 * unc2);
      uncertainty_mbdeff_upper_record[ib-1] = TMath::Sqrt(uncertainty_mbdeff_upper_record[ib-1] * uncertainty_mbdeff_upper_record[ib-1] + unc2 * unc2);
    } else {
      uncertainty_other_lower[ib-1] = TMath::Sqrt(uncertainty_other_lower[ib-1] * uncertainty_other_lower[ib-1] + unc2 * unc2);
      uncertainty_mbdeff_lower_record[ib-1] = TMath::Sqrt(uncertainty_mbdeff_lower_record[ib-1] * uncertainty_mbdeff_lower_record[ib-1] + unc2 * unc2);
    }
    uncertainty_mbdeff_upper_record[ib-1] = uncertainty_mbdeff_upper_record[ib-1] / h_unfold_calib_nominal->GetBinContent(ib);
    uncertainty_mbdeff_lower_record[ib-1] = uncertainty_mbdeff_lower_record[ib-1] / h_unfold_calib_nominal->GetBinContent(ib);
    std::cout << uncertainty_other_upper[ib-1] << " (" << uncertainty_other_upper[ib-1] / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t" << uncertainty_other_lower[ib-1] << " (" << uncertainty_other_lower[ib-1] / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1 << " (" << unc1 / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t" << unc2 << " (" << unc2 / h_unfold_calib_nominal->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1_dijet << " (" << unc1_dijet / h_unfold_calib_dijet_reweighted_1->GetBinContent(ib) * 100 << "\%) \t" << unc2_dijet << " (" << unc2_dijet / h_unfold_calib_dijet_reweighted_1->GetBinContent(ib) * 100 << "\%) \t";
    std::cout << unc1_frac << " (" << unc1_frac / h_unfold_calib_frac_reweighted_1->GetBinContent(ib) * 100 << "\%) \t" << unc2_frac << " (" << unc2_frac / h_unfold_calib_frac_reweighted_1->GetBinContent(ib) * 100 << "\%)" << std::endl;
  }

  // Get total uncertainty
  std::cout << std::endl << "Total uncertainty with background cut and unfolding uncertainty: " << std::endl << "nomiY\t\t\tnoup\t\t\tnodown" << std::endl;
  for (int ib = 0; ib < calibnpt; ++ib) {
    uncertainty_final_upper[ib] = TMath::Sqrt(uncertainty_other_upper[ib] * uncertainty_other_upper[ib] + uncertainty_cut_upper[ib] * uncertainty_cut_upper[ib] + uncertainty_unfold_upper[ib] * uncertainty_unfold_upper[ib]);
    uncertainty_final_lower[ib] = TMath::Sqrt(uncertainty_other_lower[ib] * uncertainty_other_lower[ib] + uncertainty_cut_lower[ib] * uncertainty_cut_lower[ib] + uncertainty_unfold_lower[ib] * uncertainty_unfold_lower[ib]);
    std::cout << h_unfold_calib_nominal->GetBinContent(ib+1) << "  \t\t" << uncertainty_final_upper[ib] << " (" << uncertainty_final_upper[ib] / h_unfold_calib_nominal->GetBinContent(ib+1) * 100 << "\%)\t" << uncertainty_final_lower[ib] << " (" << uncertainty_final_lower[ib] / h_unfold_calib_nominal->GetBinContent(ib+1) * 100 << "\%)" << std::endl;
  }

  std::cout << std::endl << std::endl;
  std::cout << "pT range\td^{2}#sigma/(d#eta dp_{T})\t#delta#sigma_{stat}\t\t#delta#sigma_{sys}^{low}\t#delta#sigma_{sys}^{high}" << std::endl;
  std::cout << "[GeV]\t\t[pb/GeV]\t\t\t[pb/GeV]\t\t\t[pb/GeV]\t\t\t[pb/GeV]" << std::endl;
  double sampled_N = luminosity * 25.2;
  for (int ib = 0; ib < calibnpt; ++ib) {
    h_unfold_calib_nominal->SetBinContent(ib+1, h_unfold_calib_nominal->GetBinContent(ib+1) / (float)(sampled_N * delta_eta * h_unfold_calib_nominal->GetBinWidth(ib+1)));
    h_unfold_calib_nominal->SetBinError(ib+1, h_unfold_calib_nominal->GetBinError(ib+1) / (float)(sampled_N * delta_eta * h_unfold_calib_nominal->GetBinWidth(ib+1)));
    uncertainty_final_upper[ib] = uncertainty_final_upper[ib] / (float)(sampled_N * delta_eta * h_unfold_calib_nominal->GetBinWidth(ib+1)) / (float)h_unfold_calib_nominal->GetBinContent(ib+1);
    uncertainty_final_lower[ib] = uncertainty_final_lower[ib] / (float)(sampled_N * delta_eta * h_unfold_calib_nominal->GetBinWidth(ib+1)) / (float)h_unfold_calib_nominal->GetBinContent(ib+1);
    h_unfold_calib_nominal->SetBinContent(ib+1, h_unfold_calib_nominal->GetBinContent(ib+1) * 25.2);
    h_unfold_calib_nominal->SetBinError(ib+1, h_unfold_calib_nominal->GetBinError(ib+1) * 25.2);
    uncertainty_final_upper[ib] = h_unfold_calib_nominal->GetBinContent(ib+1) * TMath::Sqrt(uncertainty_final_upper[ib] * uncertainty_final_upper[ib] + uncertainty_mbd_upper[ib] * uncertainty_mbd_upper[ib]);
    uncertainty_final_lower[ib] = h_unfold_calib_nominal->GetBinContent(ib+1) * TMath::Sqrt(uncertainty_final_lower[ib] * uncertainty_final_lower[ib] + uncertainty_mbd_lower[ib] * uncertainty_mbd_lower[ib]);
    std::cout << calibptbins[ib] << " - " << calibptbins[ib+1] << "  \t" << h_unfold_calib_nominal->GetBinContent(ib+1) << "  \t\t\t" << h_unfold_calib_nominal->GetBinError(ib+1) << " (" << h_unfold_calib_nominal->GetBinError(ib+1) / h_unfold_calib_nominal->GetBinContent(ib+1) * 100 << "\%)\t\t" << uncertainty_final_lower[ib] << " (" << uncertainty_final_lower[ib] / h_unfold_calib_nominal->GetBinContent(ib+1) * 100 << "\%)\t\t" << uncertainty_final_upper[ib] << " (" << uncertainty_final_upper[ib] / h_unfold_calib_nominal->GetBinContent(ib+1) * 100 << "\%)" << std::endl;
  }

  // Draw the final spectrum
  // sPHENIX result
  double g_sphenix_x[calibnpt], g_sphenix_xerrdown[calibnpt], g_sphenix_xerrup[calibnpt], g_sphenix_y[calibnpt], g_sphenix_yerrdown[calibnpt], g_sphenix_yerrup[calibnpt], g_sphenix_yerrstat[calibnpt];
  for (int ib = 0; ib < calibnpt; ++ ib) {
    g_sphenix_x[ib] = (calibptbins[ib] + calibptbins[ib+1]) / 2.;
    g_sphenix_xerrdown[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_sphenix_xerrup[ib] = (calibptbins[ib+1] - calibptbins[ib]) / 2.;
    g_sphenix_y[ib] = h_unfold_calib_nominal->GetBinContent(ib+1);
    g_sphenix_yerrstat[ib] = h_unfold_calib_nominal->GetBinError(ib+1);
    g_sphenix_yerrdown[ib] = TMath::Sqrt(g_sphenix_yerrstat[ib]*g_sphenix_yerrstat[ib] + uncertainty_final_lower[ib]*uncertainty_final_lower[ib]);
    g_sphenix_yerrup[ib] = TMath::Sqrt(g_sphenix_yerrstat[ib]*g_sphenix_yerrstat[ib] + uncertainty_final_upper[ib]*uncertainty_final_upper[ib]);

    uncertainty_stat_upper_record[ib] = g_sphenix_yerrstat[ib] / g_sphenix_y[ib];
    uncertainty_stat_lower_record[ib] = g_sphenix_yerrstat[ib] / g_sphenix_y[ib];
    uncertainty_total_upper_record[ib] = g_sphenix_yerrup[ib] / g_sphenix_y[ib];
    uncertainty_total_lower_record[ib] = g_sphenix_yerrdown[ib] / g_sphenix_y[ib];
  }
  TGraphAsymmErrors* g_sphenix_result = new TGraphAsymmErrors(calibnpt, g_sphenix_x, g_sphenix_y, g_sphenix_xerrdown, g_sphenix_xerrup, g_sphenix_yerrdown, g_sphenix_yerrup);
  g_sphenix_result->SetMarkerStyle(20);
  g_sphenix_result->SetMarkerSize(1.5);
  g_sphenix_result->SetMarkerColor(kAzure+2);
  g_sphenix_result->SetLineColor(kAzure+2);
  g_sphenix_result->SetFillColorAlpha(kAzure + 1, 0.5);
  ofstream fout_sphenix("sphenix_result_point.txt");
  for (int ib = 0; ib < calibnpt; ++ib) {
    fout_sphenix << g_sphenix_x[ib] << " " << g_sphenix_xerrdown[ib] << " " << g_sphenix_xerrup[ib] << " " << g_sphenix_y[ib] << " " << g_sphenix_yerrdown[ib] << " " << g_sphenix_yerrup[ib] << " " << g_sphenix_yerrstat[ib] << std::endl;
  }
  TGraphErrors* g_sphenix_result_point = new TGraphErrors(calibnpt, g_sphenix_x, g_sphenix_y, g_sphenix_xerrdown, g_sphenix_yerrstat);
  g_sphenix_result_point->SetMarkerStyle(20);
  g_sphenix_result_point->SetMarkerSize(1.5);
  g_sphenix_result_point->SetMarkerColor(kAzure+2);
  g_sphenix_result_point->SetLineColor(kAzure+2);
  g_sphenix_result_point->SetLineWidth(2);

  // Pythia result
  TFile *f_sim_specturm = new TFile("output_sim.root", "READ");
  TH1D* h_pythia_spectrum_dijet = (TH1D*)f_sim_specturm->Get("h_truth_calib_dijet");
  TH1D* h_pythia_spectrum_frac = (TH1D*)f_sim_specturm->Get("h_truth_calib_frac");
  TH1D* h_pythia_spectrum = (TH1D*)h_pythia_spectrum_dijet->Clone("h_pythia_spectrum");
  h_pythia_spectrum->Add(h_pythia_spectrum_frac);
  h_pythia_spectrum->Scale(1. / 2.);
  h_pythia_spectrum->Scale(1e+06); // N_event is scaled down by 1e+06 in analysis sim.
  for (int ib = 1; ib <= h_pythia_spectrum->GetNbinsX(); ++ib) {
    h_pythia_spectrum->SetBinContent(ib, h_pythia_spectrum->GetBinContent(ib)/ (float)(delta_eta * h_pythia_spectrum->GetBinWidth(ib)));
    h_pythia_spectrum->SetBinError(ib, h_pythia_spectrum->GetBinError(ib)/ (float)(delta_eta * h_pythia_spectrum->GetBinWidth(ib)));
  }
  double g_pythia_x[truthnpt], g_pythia_xerrdown[truthnpt], g_pythia_xerrup[truthnpt], g_pythia_y[truthnpt], g_pythia_yerrdown[truthnpt], g_pythia_yerrup[truthnpt];
  for (int ib = 0; ib < truthnpt; ++ib) {
    g_pythia_x[ib] = (truthptbins[ib] + truthptbins[ib+1]) / 2;
    g_pythia_xerrdown[ib] = (truthptbins[ib+1] - truthptbins[ib]) / 2;
    g_pythia_xerrup[ib] = (truthptbins[ib+1] - truthptbins[ib]) / 2;
    g_pythia_y[ib] = h_pythia_spectrum->GetBinContent(ib+1);
    g_pythia_yerrdown[ib] = h_pythia_spectrum->GetBinError(ib+1);
    g_pythia_yerrup[ib] = h_pythia_spectrum->GetBinError(ib+1);
  }
  TGraphAsymmErrors* g_pythia_result = new TGraphAsymmErrors(truthnpt, g_pythia_x, g_pythia_y, g_pythia_xerrdown, g_pythia_xerrup, g_pythia_yerrdown, g_pythia_yerrup);
  g_pythia_result->SetMarkerStyle(21);
  g_pythia_result->SetMarkerSize(1.5);
  g_pythia_result->SetMarkerColor(kMagenta+2);
  g_pythia_result->SetLineColor(kMagenta+2);
  ofstream fout_pythia("pythia_result_point.txt");
  for (int ib = 0; ib < truthnpt; ++ib) {
    fout_pythia << g_pythia_x[ib] << " " << g_pythia_xerrdown[ib] << " " << g_pythia_xerrup[ib] << " " << g_pythia_y[ib] << " " << g_pythia_yerrdown[ib] << " " << g_pythia_yerrup[ib] << std::endl;
  }

  // Pythia ratio
  TH1D* h_pythia_ratio = (TH1D*)h_unfold_calib_nominal->Clone("h_pythia_ratio");
  for (int ib = 1; ib <= h_pythia_ratio->GetNbinsX(); ++ib) {
    h_pythia_ratio->SetBinContent(ib, h_pythia_spectrum->GetBinContent(ib+n_underflow_bin) / h_pythia_ratio->GetBinContent(ib));
    h_pythia_ratio->SetBinError(ib, h_pythia_ratio->GetBinContent(ib) * TMath::Sqrt(TMath::Power(h_pythia_spectrum->GetBinError(ib+n_underflow_bin) / h_pythia_spectrum->GetBinContent(ib+n_underflow_bin), 2) + TMath::Power(h_unfold_calib_nominal->GetBinError(ib) / h_unfold_calib_nominal->GetBinContent(ib), 2)));
  }

  // STAR old result
  double starpointxold[] = {8.3, 10.3, 12.6, 15.5, 19.0, 23.4, 28.7, 35.3, 43.3};
  const static int starnpointold = sizeof(starpointxold) / sizeof(starpointxold[0]);
  double starpointxerrdownold[starnpointold] = {0.7, 1.0, 1.2, 1.4, 1.7, 2.1, 2.5, 3.1, 3.7};
  double starpointxerruppld[starnpointold] = {1.0, 1.1, 1.5, 1.8, 2.3, 2.8, 3.5, 4.3, 5.4};
  double starpointyold[] = {6.40e+05*2*TMath::Pi(), 2.40e+05*2*TMath::Pi(), 5.90e+04*2*TMath::Pi(), 1.19e+04*2*TMath::Pi(), 3.50e+03*2*TMath::Pi(), 7.80e+02*2*TMath::Pi(), 1.41e+02*2*TMath::Pi(), 2.20e+01*2*TMath::Pi(), 2.60e+00*2*TMath::Pi()};
  double starpointyerrdownold[] = {0.90e+05, 0.30e+05, 0.60e+04, 0.10e+04, 0.20e+03, 0.40e+02, 0.09e+02, 0.30e+01, 0.60e+00};
  double starpointyerrupold[] = {3.10e+05, 1.10e+05, 2.80e+04, 0.57e+04, 1.70e+03, 0.37e+02, 0.68e+02, 1.10e+01, 1.30e+00};
  TGraphAsymmErrors* g_star_result_old = new TGraphAsymmErrors(starnpointold, starpointxold, starpointyold, starpointxerrdownold, starpointxerruppld, starpointyerrdownold, starpointyerrupold);
  g_star_result_old->SetMarkerStyle(20);
  g_star_result_old->SetMarkerColor(kGreen+3);
  g_star_result_old->SetLineColor(kGreen+3);

  // STAR result
  double starpointx[] = {7.8, 10.4, 13.85, 18.45, 24.6, 32.8, 43.75};
  const static int starnpoint = sizeof(starpointx) / sizeof(starpointx[0]);
  double starpointxerrdown[starnpoint] = {1.1, 1.5, 1.95, 2.65, 3.5, 4.7, 6.25};
  double starpointxerrup[starnpoint] = {1.1, 1.5, 1.95, 2.65, 3.5, 4.7, 6.25};
  double starpointy[starnpoint] = {5.655596e+06, 9.684281e+05, 1.603961e+05, 2.641261e+04, 3.391784e+03, 2.871360e+02, 2.178426e+01};
  double starpointyerrdown[starnpoint] = {7.270251e+05, 1.409256e+05, 1.897396e+04, 3.829357e+03, 5.537599e+02, 5.478192e+01, 5.218297e+00};
  double starpointyerrup[starnpoint] = {7.486261e+05, 1.452015e+05, 1.897113e+04, 4.021414e+03, 5.756504e+02, 6.336199e+01, 5.765260e+00};
  double starpointunkownstat[starnpoint] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  TGraphAsymmErrors* g_star_result = new TGraphAsymmErrors(starnpoint, starpointx, starpointy, starpointxerrdown, starpointxerrup, starpointyerrdown, starpointyerrup);
  g_star_result->SetMarkerStyle(20);
  g_star_result->SetMarkerSize(1.5);
  g_star_result->SetMarkerColor(kGreen+3);
  g_star_result->SetLineColor(kGreen+3);
  g_star_result->SetFillColorAlpha(kGreen+2, 0.5);
  TGraphErrors* g_star_result_point = new TGraphErrors(starnpoint, starpointx, starpointy, starpointxerrdown, starpointunkownstat);
  g_star_result_point->SetMarkerStyle(20);
  g_star_result_point->SetMarkerSize(1.5);
  g_star_result_point->SetMarkerColor(kGreen+3);
  g_star_result_point->SetLineColor(kGreen+3);
  g_star_result_point->SetLineWidth(2);

  // PHENIX result
  //double phenixpointx[] = {8.45355, 9.45319, 10.8252, 13.0123, 15.7055, 18.739, 22.0908, 26.2576, 31.1771, 37.4739};
  //const static int phenixnpoint = sizeof(phenixpointx) / sizeof(phenixpointx[0]);
  //double phenixpointxerrdown[phenixnpoint] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
  //double phenixpointxerrup[phenixnpoint] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
  //double phenixpointy[phenixnpoint] = {5.86972e+05, 2.99829e+05, 1.3335e+05, 4.07571e+04, 1.24647e+04, 4.28294e+03, 1.39962e+03, 4.11469e+02, 1.01631e+02, 2.49519e+01};
  //double phenixpointysysdown[phenixnpoint] = {2.32513e+05, 1.08448e+05, 4.80741e+04, 1.44942e+04, 4.19428e+03, 1.21361e+03, 3.87209e+02, 9.41715e+01, 1.20866e+01, 1.06009};
  //double phenixpointysysup[phenixnpoint] = {1.90074e+05, 9.64104e+04, 4.43763e+04, 1.42428e+04, 4.28854e+03, 1.30758e+03, 4.68572e+02, 1.15733e+02, 1.5672e+01, 2.62781};
  //double phenixpointystat[phenixnpoint] = {3.5247e+03, 2.16461e+03, 1.15278e+03, 4.6706e+02, 1.86848e+02, 1.04529e+02, 4.38652e+01, 2.17658e+01, 9.3723, 4.46448};
  double phenixpointx[] = {13.01, 15.7, 18.74, 22.08, 26.26, 31.19, 37.47};
  const static int phenixnpoint = sizeof(phenixpointx) / sizeof(phenixpointx[0]);
  double phenixpointxerrdown[phenixnpoint] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
  double phenixpointxerrup[phenixnpoint] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
  double phenixpointy[] = {7.78e+04, 2.46e+04, 8.8e+03, 2.95e+03, 8e+02, 1.9e+02, 4.1e+01};
  double phenixpointysysdown[] = {2.27e+04, 6.9e+03, 3e+03, 1.21e+03, 3.8e+02, 8e+01, 2.1e+01};
  double phenixpointysysup[] = {2.3e+04, 7e+03, 3e+03, 1.21e+03, 3.8e+02, 8e+01, 2.1e+01};
  double phenixpointystat[] = {3.2e+03, 1.3e+03, 7e+02, 3e+02, 1.2e+02, 5e+01, 2.1e+01};
  double phenixpointyerrdown[phenixnpoint];
  double phenixpointyerrup[phenixnpoint];
  for (int i = 0; i < phenixnpoint; ++i) {
    phenixpointyerrdown[i] = TMath::Sqrt(phenixpointystat[i]*phenixpointystat[i] + phenixpointysysdown[i]*phenixpointysysdown[i]);
    phenixpointyerrup[i] = TMath::Sqrt(phenixpointystat[i]*phenixpointystat[i] + phenixpointysysup[i]*phenixpointysysup[i]);
  }
  TGraphAsymmErrors* g_phenix_result = new TGraphAsymmErrors(phenixnpoint, phenixpointx, phenixpointy, phenixpointxerrdown, phenixpointxerrup, phenixpointyerrdown, phenixpointyerrup);
  g_phenix_result->SetMarkerStyle(20);
  g_phenix_result->SetMarkerSize(1.5);
  g_phenix_result->SetMarkerColor(kGray+2);
  g_phenix_result->SetLineColor(kGray+2);
  g_phenix_result->SetFillColorAlpha(kGray+1, 0.5);
  ofstream fout_phenix("phenix_result_point.txt");
  for (int ib = 0; ib < phenixnpoint; ++ib) {
    fout_phenix << phenixpointx[ib] << " " << phenixpointxerrdown[ib] << " " << phenixpointxerrup[ib] << " " << phenixpointy[ib] << " " << phenixpointyerrdown[ib] << " " << phenixpointyerrup[ib] << " " << phenixpointystat[ib] << std::endl;
  }
  TGraphErrors* g_phenix_result_point = new TGraphErrors(phenixnpoint, phenixpointx, phenixpointy, phenixpointxerrdown, phenixpointystat);
  g_phenix_result_point->SetMarkerStyle(20);
  g_phenix_result_point->SetMarkerSize(1.5);
  g_phenix_result_point->SetMarkerColor(kGray+2);
  g_phenix_result_point->SetLineColor(kGray+2);
  g_phenix_result_point->SetLineWidth(2);

  // NLO result
  ifstream myfile1, myfile05, myfile2;
  TGraph *tjet = new TGraph();
  TGraph *tjet05 = new TGraph();
  TGraph *tjet2 = new TGraph();
  double scalefactor = 2.0 * 3.14159;
  myfile1.open("sphenix_nlo/sPhenix-sc1.dat");
  myfile05.open("sphenix_nlo/sPhenix-sc05.dat");
  myfile2.open("sphenix_nlo/sPhenix-sc2.dat");
  double drop1, drop2, drop3, pt1, yield1, pt05, yield05, pt2, yield2;
  std::vector<double> xpoint, ypoint05, ypoint2;
  for (int index=0; index<66; index++) {
    myfile1 >> drop1 >> drop2 >> drop3 >> pt1 >> yield1;
    yield1 = yield1 * scalefactor * pt1;
    tjet->SetPoint(index, pt1, yield1);
    xpoint.push_back(pt1);

    myfile05 >> drop1 >> drop2 >> drop3 >> pt05 >> yield05;
    yield05 = yield05 * scalefactor * pt05;
    tjet05->SetPoint(index, pt05, yield05);
    ypoint05.push_back(yield05);

    myfile2 >> drop1 >> drop2 >> drop3 >> pt2 >> yield2;
    yield2 = yield2 * scalefactor * pt2;
    tjet2->SetPoint(index, pt2, yield2);
    ypoint2.push_back(yield2);
  }
  std::vector <double> xpoint_shadow, ypoint_shadow;
  for (int i = 0; i < xpoint.size(); i++) {
    xpoint_shadow.push_back(xpoint[i]);
    ypoint_shadow.push_back(ypoint05[i]);
  }
  for (int i = xpoint.size(); i > 0; i--) {
    xpoint_shadow.push_back(xpoint[i-1]);
    ypoint_shadow.push_back(ypoint2[i-1]);
  }
  TGraph *tjet_shadow = new TGraph(2*xpoint.size(), &xpoint_shadow[0], &ypoint_shadow[0]);
  tjet_shadow->SetMarkerStyle(20);
  tjet_shadow->SetMarkerSize(1.5);
  tjet_shadow->SetMarkerColor(kRed-4);
  tjet_shadow->SetLineColor(kRed);
  tjet_shadow->SetFillColorAlpha(kRed-4, 0.5);

  myfile1.close();
  myfile05.close();
  myfile2.close();
  tjet->SetMarkerStyle(21);
  tjet->SetMarkerColor(2);
  tjet->SetLineColor(2);
  tjet05->SetMarkerStyle(21);
  tjet05->SetMarkerColor(2);
  tjet05->SetLineColor(2);
  tjet2->SetMarkerStyle(21);
  tjet2->SetMarkerColor(2);
  tjet2->SetLineColor(2);

  TImage *img = TImage::Open("L_int.png");
  if (!img) {
     printf("Could not load image.\n");
     return;
  }
  int imgWidth = img->GetWidth();
  int imgHeight = img->GetHeight();
  double aspectRatio = static_cast<double>(imgWidth) / imgHeight;
  TFile *f_L_int = new TFile("L_int.root", "READ");
  TH2D* h_L_int = (TH2D*)f_L_int->Get("h_L_int");
  double desiredHeight;
  TBox *coverTop = new TBox(0, imgHeight-1, imgWidth, imgHeight);
  coverTop->SetFillColor(kWhite);
  coverTop->SetLineColor(kWhite);
  coverTop->SetLineWidth(0);
  TBox *coverRight = new TBox(imgWidth-1, 0, imgWidth, imgHeight);
  coverRight->SetFillColor(kWhite);
  coverRight->SetLineColor(kWhite);
  coverRight->SetLineWidth(0);

  // Draw the final spectrum
  TCanvas *can_nlo = new TCanvas("can_nlo", "", 850, 800);
  gStyle->SetPalette(57);
  can_nlo->SetTopMargin(0.03);
  can_nlo->SetLeftMargin(0.15);
  can_nlo->SetBottomMargin(0.12);
  can_nlo->SetRightMargin(0.05);
  can_nlo->SetLogy();
  TH2F *frame = new TH2F("frame", "", 10, 15, 72, 10, 1e-03, 2e+05);
  frame->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame->GetXaxis()->SetTitleOffset(0.98);
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetXaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->Draw();
  g_sphenix_result->SetMarkerStyle(20);
  g_sphenix_result->SetMarkerColor(kAzure+2);
  g_sphenix_result->SetLineColor(kAzure+2);
  g_sphenix_result->SetFillColorAlpha(kAzure + 1, 0.5);
  g_sphenix_result->Draw("2 same");
  g_sphenix_result_point->SetMarkerStyle(20);
  g_sphenix_result_point->SetMarkerColor(kAzure+2);
  g_sphenix_result_point->SetLineColor(kAzure+2);
  g_sphenix_result_point->SetLineWidth(2);
  g_sphenix_result_point->Draw("P same");
  myText(0.4, 0.9, 1, "#bf{#it{sPHENIX}} Preliminary", 0.04);
  myText(0.4, 0.85, 1, "p+p#sqrt{s} = 200 GeV,#kern[+0.7]{#scale[0.5]{#int}}#kern[+0.2]{#it{L}}dt = 16.6 pb^{-1}", 0.04);
  myText(0.4, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  myText(0.4, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  TLegend *leg = new TLegend(0.2, 0.18, 0.9, 0.29);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->AddEntry(g_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg->AddEntry(tjet, "NLO pQCD (No Hadronization)", "pl");
  leg->AddEntry(tjet, "Werner Vogelsang", "");
  leg->Draw();
  tjet->Draw("pl same");
  tjet05->Draw("l same");
  tjet2->Draw("l same");
  can_nlo->SaveAs("figure/spectrum_calib_nlo.pdf");

  TCanvas *can_nlo_shadow = new TCanvas("can_nlo_shadow", "", 850, 800);
  gStyle->SetPalette(57);
  can_nlo_shadow->SetTopMargin(0.03);
  can_nlo_shadow->SetLeftMargin(0.15);
  can_nlo_shadow->SetBottomMargin(0.12);
  can_nlo_shadow->SetRightMargin(0.05);
  can_nlo_shadow->SetLogy();
  TH2F *frame_shadow = new TH2F("frame_shadow", "", 10, 15, 72, 10, 1e-03, 2e+05);
  frame_shadow->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_shadow->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_shadow->GetXaxis()->SetTitleOffset(0.98);
  frame_shadow->GetYaxis()->SetTitleOffset(1.15);
  frame_shadow->GetXaxis()->SetLabelSize(0.045);
  frame_shadow->GetYaxis()->SetLabelSize(0.045);
  frame_shadow->GetXaxis()->CenterTitle();
  frame_shadow->GetYaxis()->CenterTitle();
  frame_shadow->Draw();
  g_sphenix_result->SetMarkerStyle(20);
  g_sphenix_result->SetMarkerColor(kAzure+2);
  g_sphenix_result->SetLineColor(kAzure+2);
  g_sphenix_result->SetFillColorAlpha(kAzure + 1, 0.5);
  g_sphenix_result->Draw("2 same");
  g_sphenix_result_point->SetMarkerStyle(20);
  g_sphenix_result_point->SetMarkerColor(kAzure+2);
  g_sphenix_result_point->SetLineColor(kAzure+2);
  g_sphenix_result_point->SetLineWidth(2);
  g_sphenix_result_point->Draw("P same");
  myText(0.4, 0.9, 1, "#bf{#it{sPHENIX}} Preliminary", 0.04);
  myText(0.4, 0.85, 1, "p+p#sqrt{s} = 200 GeV,#kern[+0.7]{#scale[0.5]{#int}}#kern[+0.2]{#it{L}}dt = 16.6 pb^{-1}", 0.04);
  myText(0.4, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  myText(0.4, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  TLegend *leg_shadow = new TLegend(0.2, 0.18, 0.9, 0.29);
  leg_shadow->SetBorderSize(0);
  leg_shadow->SetFillStyle(0);
  leg_shadow->SetTextSize(0.03);
  leg_shadow->SetTextFont(42);
  leg_shadow->AddEntry(g_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_shadow->AddEntry(tjet_shadow, "NLO pQCD (No Hadronization)", "f");
  leg_shadow->AddEntry(tjet_shadow, "Werner Vogelsang", "");
  leg_shadow->Draw();
  tjet_shadow->Draw("F same");
  can_nlo_shadow->SaveAs("figure/spectrum_calib_nlo_shadow.pdf");

  TCanvas *can_sphenix = new TCanvas("can_sphenix", "", 850, 800);
  gStyle->SetPalette(57);
  can_sphenix->SetTopMargin(0.03);
  can_sphenix->SetLeftMargin(0.15);
  can_sphenix->SetBottomMargin(0.12);
  can_sphenix->SetRightMargin(0.05);
  can_sphenix->SetLogy();
  TH2F *frame_sphenix = new TH2F("frame_sphenix", "", 10, 15, 72, 10, 1e-03, 2e+05);
  frame_sphenix->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_sphenix->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_sphenix->GetXaxis()->SetTitleOffset(0.98);
  frame_sphenix->GetYaxis()->SetTitleOffset(1.15);
  frame_sphenix->GetXaxis()->SetLabelSize(0.045);
  frame_sphenix->GetYaxis()->SetLabelSize(0.045);
  frame_sphenix->GetXaxis()->CenterTitle();
  frame_sphenix->GetYaxis()->CenterTitle();
  frame_sphenix->Draw();
  g_sphenix_result->Draw("2 same");
  g_sphenix_result_point->Draw("P same");
  myText(0.4, 0.9, 1, "#bf{#it{sPHENIX}} Preliminary", 0.04);
  myText(0.4, 0.85, 1, "p+p#sqrt{s} = 200 GeV,       = 16.6 pb^{-1}", 0.04);
  //desiredHeight = 0.04; TPad *pad_sphenix = new TPad("pad_sphenix", "", 0.68, 0.85 - desiredHeight/2. - 0.005, 0.68+desiredHeight*aspectRatio, 0.85+desiredHeight/2. - 0.005); pad_sphenix->Draw(); pad_sphenix->cd(); img->Draw("xxx"); pad_sphenix->Modified(); pad_sphenix->Update(); can_sphenix->cd(); can_sphenix->Modified(); can_sphenix->Update();
  myText(0.4, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  myText(0.4, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  TLegend *leg_sphenix = new TLegend(0.2, 0.18, 0.9, 0.21);
  leg_sphenix->SetBorderSize(0);
  leg_sphenix->SetFillStyle(0);
  leg_sphenix->SetTextSize(0.03);
  leg_sphenix->SetTextFont(42);
  leg_sphenix->AddEntry(g_sphenix_result, "sPHENIX Run 2024 data", "pf");
  //leg_sphenix->Draw();
  can_sphenix->Print("figure/spectrum_calib_sphenix.pdf");

  TCanvas *can_data = new TCanvas("can_data", "", 850, 800);
  gStyle->SetPalette(57);
  can_data->SetTopMargin(0.03);
  can_data->SetLeftMargin(0.15);
  can_data->SetBottomMargin(0.12);
  can_data->SetRightMargin(0.05);
  can_data->SetLogy();
  TH2F *frame_data = new TH2F("frame_data", "", 10, 15, 72, 10, 1e-03, 2e+05);
  frame_data->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_data->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_data->GetXaxis()->SetTitleOffset(0.98);
  frame_data->GetYaxis()->SetTitleOffset(1.15);
  frame_data->GetXaxis()->SetLabelSize(0.045);
  frame_data->GetYaxis()->SetLabelSize(0.045);
  frame_data->GetXaxis()->CenterTitle();
  frame_data->GetYaxis()->CenterTitle();
  frame_data->Draw();
  myText(0.4, 0.9, 1, "#bf{#it{sPHENIX}} Preliminary", 0.04);
  myText(0.4, 0.85, 1, "p+p#sqrt{s} = 200 GeV,         = 16.6 pb^{-1}", 0.04);
  //desiredHeight = 0.04; TPad *pad_data = new TPad("pad_data", "", 0.68, 0.85 - desiredHeight/2. - 0.005, 0.68+desiredHeight*aspectRatio, 0.85+desiredHeight/2. - 0.005); pad_data->Draw(); pad_data->cd(); img->Draw("xxx"); pad_data->Modified(); pad_data->Update(); can_data->cd(); can_data->Modified(); can_data->Update();
  myText(0.4, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  myText(0.4, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  TLegend *leg_data = new TLegend(0.15, 0.15, 0.85, 0.4);
  leg_data->SetBorderSize(0);
  leg_data->SetFillStyle(0);
  leg_data->SetTextSize(0.03);
  leg_data->SetTextFont(42);
  leg_data->AddEntry(g_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_data->AddEntry(g_star_result, "STAR Run 2012 Preliminary", "pf");
  leg_data->AddEntry(g_star_result, "|#eta| < 0.8, anti-k_{#kern[-0.9]{t}} #kern[-0.3]{#it{R}} = 0.6", "");
  leg_data->AddEntry(g_phenix_result, "PHENIX Run 2008", "pf");
  leg_data->AddEntry(g_phenix_result, "|#eta| < 0.15, anti-k_{#kern[-1]{t}} #kern[-0.3]{#it{R}} = 0.3", "");
  leg_data->Draw();
  g_star_result->Draw("2 same");
  g_phenix_result->Draw("2 same");
  g_sphenix_result->Draw("2 same");
  g_star_result_point->Draw("p same");
  g_phenix_result_point->Draw("p same");
  g_sphenix_result_point->Draw("P same");
  can_data->SaveAs("figure/spectrum_calib_datacompare.png");

  TCanvas *can_dataold = new TCanvas("can_dataold", "", 850, 800);
  gStyle->SetPalette(57);
  can_dataold->SetTopMargin(0.03);
  can_dataold->SetLeftMargin(0.15);
  can_dataold->SetBottomMargin(0.12);
  can_dataold->SetRightMargin(0.05);
  can_dataold->SetLogy();
  TH2F *frame_dataold = new TH2F("frame_dataold", "", 10, 15, 72, 10, 1e-03, 2e+05);
  frame_dataold->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_dataold->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_dataold->GetXaxis()->SetTitleOffset(0.98);
  frame_dataold->GetYaxis()->SetTitleOffset(1.15);
  frame_dataold->GetXaxis()->SetLabelSize(0.045);
  frame_dataold->GetYaxis()->SetLabelSize(0.045);
  frame_dataold->GetXaxis()->CenterTitle();
  frame_dataold->GetYaxis()->CenterTitle();
  frame_dataold->Draw();
  g_sphenix_result->Draw("2 same");
  g_sphenix_result_point->Draw("P same");
  myText(0.4, 0.9, 1, "#bf{#it{sPHENIX}} Preliminary", 0.04);
  myText(0.4, 0.85, 1, "p+p#sqrt{s} = 200 GeV,#kern[+0.7]{#scale[0.5]{#int}}#kern[+0.2]{#it{L}}dt = 16.6 pb^{-1}", 0.04);
  myText(0.4, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  myText(0.4, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  TLegend *leg_dataold = new TLegend(0.2, 0.18, 0.9, 0.29);
  leg_dataold->SetBorderSize(0);
  leg_dataold->SetFillStyle(0);
  leg_dataold->SetTextSize(0.03);
  leg_dataold->SetTextFont(42);
  leg_dataold->AddEntry(g_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_dataold->AddEntry(g_star_result_old, "STAR Midpoint Cone#it{R}=0.4 (Combined HT)", "pl");
  leg_dataold->AddEntry(g_phenix_result, "PHENIX Run12+15, anti-k_{t} R=0.3", "pl");
  leg_dataold->Draw();
  g_star_result_old->Draw("p same");
  g_phenix_result->Draw("p same");
  can_dataold->SaveAs("figure/spectrum_calib_datacompareold.pdf");

  TCanvas *can_sim = new TCanvas("can_sim", "", 850, 800);
  gStyle->SetPalette(57);
  can_sim->SetTopMargin(0.03);
  can_sim->SetLeftMargin(0.15);
  can_sim->SetBottomMargin(0.12);
  can_sim->SetRightMargin(0.05);
  can_sim->SetLogy();
  TH2F *frame_sim = new TH2F("frame_sim", "", 10, 15, 72, 10, 1e-03, 2e+05);
  frame_sim->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_sim->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_sim->GetXaxis()->SetTitleOffset(0.98);
  frame_sim->GetYaxis()->SetTitleOffset(1.15);
  frame_sim->GetXaxis()->SetLabelSize(0.045);
  frame_sim->GetYaxis()->SetLabelSize(0.045);
  frame_sim->GetXaxis()->CenterTitle();
  frame_sim->GetYaxis()->CenterTitle();
  frame_sim->Draw();
  g_sphenix_result->Draw("2 same");
  myText(0.4, 0.9, 1, "#bf{#it{sPHENIX}} Preliminary", 0.04);
  myText(0.4, 0.85, 1, "p+p#sqrt{s} = 200 GeV,#kern[+0.7]{#scale[0.5]{#int}}#kern[+0.2]{#it{L}}dt = 16.6 pb^{-1}", 0.04);
  myText(0.4, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  myText(0.4, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  TLegend *leg_sim = new TLegend(0.2, 0.18, 0.9, 0.25);
  leg_sim->SetBorderSize(0);
  leg_sim->SetFillStyle(0);
  leg_sim->SetTextSize(0.03);
  leg_sim->SetTextFont(42);
  leg_sim->AddEntry(g_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_sim->AddEntry(g_pythia_result, "PYTHIA8 truth jet spectrum", "pl");
  leg_sim->Draw();
  g_pythia_result->Draw("p same");
  g_sphenix_result_point->Draw("P same");
  can_sim->SaveAs("figure/spectrum_calib_simcompare.pdf");

  TCanvas *can_simnlo = new TCanvas("can_simnlo", "", 850, 800);
  gStyle->SetPalette(57);
  can_simnlo->SetTopMargin(0.03);
  can_simnlo->SetLeftMargin(0.15);
  can_simnlo->SetBottomMargin(0.12);
  can_simnlo->SetRightMargin(0.05);
  can_simnlo->SetLogy();
  TH2F *frame_simnlo = new TH2F("frame_simnlo", "", 10, 15, 72, 10, 1e-03, 2e+05);
  frame_simnlo->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_simnlo->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_simnlo->GetXaxis()->SetTitleOffset(0.98);
  frame_simnlo->GetYaxis()->SetTitleOffset(1.15);
  frame_simnlo->GetXaxis()->SetLabelSize(0.045);
  frame_simnlo->GetYaxis()->SetLabelSize(0.045);
  frame_simnlo->GetXaxis()->CenterTitle();
  frame_simnlo->GetYaxis()->CenterTitle();
  frame_simnlo->Draw();
  g_sphenix_result->Draw("2 same");
  myText(0.4, 0.9, 1, "#bf{#it{sPHENIX}} Preliminary", 0.04);
  myText(0.4, 0.85, 1, "p+p#sqrt{s} = 200 GeV,#kern[+0.7]{#scale[0.5]{#int}}#kern[+0.2]{#it{L}}dt = 16.6 pb^{-1}", 0.04);
  myText(0.4, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  myText(0.4, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  TLegend *leg_simnlo = new TLegend(0.2, 0.18, 0.9, 0.33);
  leg_simnlo->SetBorderSize(0);
  leg_simnlo->SetFillStyle(0);
  leg_simnlo->SetTextSize(0.03);
  leg_simnlo->SetTextFont(42);
  leg_simnlo->AddEntry(g_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_simnlo->AddEntry(g_pythia_result, "PYTHIA8 truth jet spectrum", "pl");
  leg_simnlo->AddEntry(tjet, "NLO pQCD (No Hadronization)", "pl");
  leg_simnlo->AddEntry(tjet, "Werner Vogelsang", "");
  leg_simnlo->Draw();
  g_pythia_result->Draw("p same");
  g_sphenix_result_point->Draw("P same");
  tjet->Draw("pl same");
  tjet05->Draw("l same");
  tjet2->Draw("l same");
  can_simnlo->SaveAs("figure/spectrum_calib_simnlocompare.pdf");

  TCanvas *can_simnlo_shadow = new TCanvas("can_simnlo_shadow", "", 2550, 2400);
  gStyle->SetPalette(57);
  can_simnlo_shadow->SetTopMargin(0.03);
  can_simnlo_shadow->SetLeftMargin(0.15);
  can_simnlo_shadow->SetBottomMargin(0.12);
  can_simnlo_shadow->SetRightMargin(0.05);
  can_simnlo_shadow->SetLogy();
  TH2F *frame_simnlo_shadow = new TH2F("frame_simnlo_shadow", "", 10, 15, 72, 10, 1e-03, 2e+05);
  frame_simnlo_shadow->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_simnlo_shadow->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_simnlo_shadow->GetXaxis()->SetTitleOffset(0.98);
  frame_simnlo_shadow->GetYaxis()->SetTitleOffset(1.15);
  frame_simnlo_shadow->GetXaxis()->SetLabelSize(0.045);
  frame_simnlo_shadow->GetYaxis()->SetLabelSize(0.045);
  frame_simnlo_shadow->GetXaxis()->CenterTitle();
  frame_simnlo_shadow->GetYaxis()->CenterTitle();
  frame_simnlo_shadow->Draw();
  g_sphenix_result->Draw("2 same");
  myText(0.4, 0.9, 1, "#bf{#it{sPHENIX}} Preliminary", 0.04);
  myText(0.4, 0.85, 1, "p+p#sqrt{s} = 200 GeV,       = 16.6 pb^{-1}", 0.04);
  //desiredHeight = 0.04; TPad* pad_simnlo_shadow = new TPad("pad_simnlo_shadow", "", 0.68, 0.85 - desiredHeight/2. - 0.002, 0.68+(desiredHeight-0.006)*aspectRatio, 0.85+desiredHeight/2. - 0.008); pad_simnlo_shadow->Draw(); set_tpad(pad_simnlo_shadow, h_L_int, coverTop, coverRight); can_simnlo_shadow->cd(); can_simnlo_shadow->Update(); 
  //TPad *pad_simnlo_shadow = new TPad("pad_simnlo_shadow", "", 0.68, 0.85 - desiredHeight/2. - 0.005, 0.68+desiredHeight*aspectRatio, 0.85+desiredHeight/2. - 0.005); pad_simnlo_shadow->Draw(); pad_simnlo_shadow->cd(); img->Draw("X"); can_simnlo_shadow->cd(); can_simnlo_shadow->Update();
  //myText(0.4, 0.85, 1, "p+p#sqrt{s} = 200 GeV,#kern[+0.7]{#scale[0.5]{#int}}#kern[+0.2]{#it{L}}dt = 16.6 pb^{-1}", 0.04);
  myText(0.4, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  myText(0.4, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  TLegend *leg_simnlo_shadow = new TLegend(0.2, 0.18, 0.9, 0.33);
  leg_simnlo_shadow->SetBorderSize(0);
  leg_simnlo_shadow->SetFillStyle(0);
  leg_simnlo_shadow->SetTextSize(0.03);
  leg_simnlo_shadow->SetTextFont(42);
  g_pythia_result->SetMarkerSize(4.5);
  g_sphenix_result->SetMarkerSize(4.5);
  g_sphenix_result_point->SetMarkerSize(4.5);
  leg_simnlo_shadow->AddEntry(g_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_simnlo_shadow->AddEntry(g_pythia_result, "PYTHIA8 truth jet spectrum", "pl");
  leg_simnlo_shadow->AddEntry(tjet_shadow, "NLO pQCD (No Hadronization)", "f");
  leg_simnlo_shadow->AddEntry(tjet_shadow, "Werner Vogelsang", "");
  leg_simnlo_shadow->Draw();
  g_pythia_result->Draw("p same");
  g_sphenix_result_point->Draw("P same");
  tjet_shadow->Draw("F same");
  can_simnlo_shadow->SaveAs("figure/spectrum_calib_simnlocompare_shadow.pdf");
  gStyle->SetPalette(55);
  g_pythia_result->SetMarkerSize(1.5);
  g_sphenix_result->SetMarkerSize(1.5);
  g_sphenix_result_point->SetMarkerSize(1.5);

  TCanvas *can_simratio = new TCanvas("can_simratio", "", 850, 1091);
  can_simratio->Divide(1, 2);
  gStyle->SetPalette(57);
  TPad *pad_1 = (TPad*)can_simratio->cd(1);
  pad_1->SetPad(0, 0.33333, 1, 1);
  pad_1->SetTopMargin(0.03);
  pad_1->SetLeftMargin(0.15);
  pad_1->SetBottomMargin(0.035);
  pad_1->SetRightMargin(0.05);
  pad_1->SetLogy();
  TH2F *frame_simratio = new TH2F("frame_simratio", "", 10, 15, 72, 10, 1e-03, 2e+05);
  frame_simratio->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_simratio->GetYaxis()->SetTitle("d^{2}#sigma/(d#eta dp_{T}) [pb/GeV]");
  frame_simratio->GetXaxis()->SetTitleOffset(0.98);
  frame_simratio->GetYaxis()->SetTitleOffset(1.15);
  frame_simratio->GetXaxis()->SetLabelSize(0.045);
  frame_simratio->GetYaxis()->SetLabelSize(0.045);
  frame_simratio->GetXaxis()->SetLabelOffset(2);
  frame_simratio->GetXaxis()->CenterTitle();
  frame_simratio->GetYaxis()->CenterTitle();
  frame_simratio->Draw();
  g_sphenix_result->Draw("2 same");
  g_sphenix_result_point->Draw("P same");
  myText(0.4, 0.9, 1, "#bf{#it{sPHENIX}} Preliminary", 0.04);
  myText(0.4, 0.85, 1, "p+p#sqrt{s} = 200 GeV,#kern[+0.7]{#scale[0.5]{#int}}#kern[+0.2]{#it{L}}dt = 16.6 pb^{-1}", 0.04);
  myText(0.4, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.04);
  myText(0.4, 0.75, 1, "|#eta^{jet}| < 0.7", 0.04);
  TLegend *leg_simratio = new TLegend(0.2, 0.18, 0.9, 0.25);
  leg_simratio->SetBorderSize(0);
  leg_simratio->SetFillStyle(0);
  leg_simratio->SetTextSize(0.03);
  leg_simratio->SetTextFont(42);
  leg_simratio->AddEntry(g_sphenix_result, "sPHENIX Run 2024 data", "pf");
  leg_simratio->AddEntry(g_pythia_result, "PYTHIA8 truth jet spectrum", "pl");
  leg_simratio->Draw();
  g_pythia_result->Draw("p same");
  TPad *pad_2 = (TPad*)can_simratio->cd(2);
  pad_2->SetPad(0, 0, 1, 0.333333);
  pad_2->SetTopMargin(0.02);
  pad_2->SetLeftMargin(0.15);
  pad_2->SetBottomMargin(0.25);
  pad_2->SetRightMargin(0.05);
  TH2F *frame_ratio = new TH2F("frame_ratio", "", 10, 15, 72, 2, 0.6, 1.9);
  frame_ratio->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_ratio->GetYaxis()->SetTitle("PYTHIA8/Data");
  frame_ratio->GetXaxis()->SetTitleOffset(0.98);
  frame_ratio->GetYaxis()->SetTitleOffset(0.5);
  frame_ratio->GetXaxis()->SetLabelSize(0.045);
  frame_ratio->GetYaxis()->SetLabelSize(0.045);
  frame_ratio->GetXaxis()->CenterTitle();
  frame_ratio->GetYaxis()->CenterTitle();
  frame_ratio->GetYaxis()->SetTitleOffset(frame_simratio->GetYaxis()->GetTitleOffset()*1/2.);
  frame_ratio->GetYaxis()->SetLabelOffset(frame_simratio->GetYaxis()->GetLabelOffset()*1/2.);
  frame_ratio->GetXaxis()->SetLabelSize(frame_simratio->GetXaxis()->GetLabelSize()*2);
  frame_ratio->GetYaxis()->SetLabelSize(frame_simratio->GetYaxis()->GetLabelSize()*2);
  frame_ratio->GetXaxis()->SetTitleSize(frame_simratio->GetXaxis()->GetTitleSize()*2);
  frame_ratio->GetYaxis()->SetTitleSize(frame_simratio->GetYaxis()->GetTitleSize()*2);
  frame_ratio->Draw();
  h_pythia_ratio->SetMarkerStyle(20);
  h_pythia_ratio->SetMarkerColor(kBlack);
  h_pythia_ratio->Draw("P same");
  TLine *line_ratio = new TLine(15, 1, 72, 1);
  line_ratio->SetLineStyle(3);
  line_ratio->Draw("same");
  can_simratio->SaveAs("figure/spectrum_calib_simratio.pdf");

  // Draw the uncertainty
  TH1D* h_uncertainty_unfold_upper = new TH1D("h_uncertainty_unfold_upper", "", calibnpt, calibptbins); TH1D* h_uncertainty_unfold_lower = new TH1D("h_uncertainty_unfold_lower", "", calibnpt, calibptbins);
  TH1D* h_uncertainty_cut_upper = new TH1D("h_uncertainty_cut_upper", "", calibnpt, calibptbins); TH1D* h_uncertainty_cut_lower = new TH1D("h_uncertainty_cut_lower", "", calibnpt, calibptbins);
  TH1D* h_uncertainty_jes_upper = new TH1D("h_uncertainty_jes_upper", "", calibnpt, calibptbins); TH1D* h_uncertainty_jes_lower = new TH1D("h_uncertainty_jes_lower", "", calibnpt, calibptbins);
  TH1D* h_uncertainty_jer_upper = new TH1D("h_uncertainty_jer_upper", "", calibnpt, calibptbins); TH1D* h_uncertainty_jer_lower = new TH1D("h_uncertainty_jer_lower", "", calibnpt, calibptbins);
  TH1D* h_uncertainty_jeteff_upper = new TH1D("h_uncertainty_jeteff_upper", "", calibnpt, calibptbins); TH1D* h_uncertainty_jeteff_lower = new TH1D("h_uncertainty_jeteff_lower", "", calibnpt, calibptbins);
  TH1D* h_uncertainty_mbdeff_upper = new TH1D("h_uncertainty_mbdeff_upper", "", calibnpt, calibptbins); TH1D* h_uncertainty_mbdeff_lower = new TH1D("h_uncertainty_mbdeff_lower", "", calibnpt, calibptbins);
  TH1D* h_uncertainty_mbdcs_upper = new TH1D("h_uncertainty_mbdcs_upper", "", calibnpt, calibptbins); TH1D* h_uncertainty_mbdcs_lower = new TH1D("h_uncertainty_mbdcs_lower", "", calibnpt, calibptbins);
  TH1D* h_uncertainty_stat_upper = new TH1D("h_uncertainty_stat_upper", "", calibnpt, calibptbins); TH1D* h_uncertainty_stat_lower = new TH1D("h_uncertainty_stat_lower", "", calibnpt, calibptbins);
  TH1D* h_uncertainty_total_upper = new TH1D("h_uncertainty_total_upper", "", calibnpt, calibptbins); TH1D* h_uncertainty_total_lower = new TH1D("h_uncertainty_total_lower", "", calibnpt, calibptbins);
  for (int ib = 1; ib <= h_uncertainty_cut_upper->GetNbinsX(); ++ib) {
    h_uncertainty_unfold_upper->SetBinContent(ib, uncertainty_unfold_upper_record[ib-1]);
    h_uncertainty_unfold_lower->SetBinContent(ib, -uncertainty_unfold_lower_record[ib-1]);
    h_uncertainty_cut_upper->SetBinContent(ib, uncertainty_cut_upper_record[ib-1]);
    h_uncertainty_cut_lower->SetBinContent(ib, -uncertainty_cut_lower_record[ib-1]);
    h_uncertainty_jes_upper->SetBinContent(ib, uncertainty_jes_upper_record[ib-1]);
    h_uncertainty_jes_lower->SetBinContent(ib, -uncertainty_jes_lower_record[ib-1]);
    h_uncertainty_jer_upper->SetBinContent(ib, uncertainty_jer_upper_record[ib-1]);
    h_uncertainty_jer_lower->SetBinContent(ib, -uncertainty_jer_lower_record[ib-1]);
    h_uncertainty_jeteff_upper->SetBinContent(ib, uncertainty_jeteff_upper_record[ib-1]);
    h_uncertainty_jeteff_lower->SetBinContent(ib, -uncertainty_jeteff_lower_record[ib-1]);
    h_uncertainty_mbdeff_upper->SetBinContent(ib, uncertainty_mbdeff_upper_record[ib-1]);
    h_uncertainty_mbdeff_lower->SetBinContent(ib, -uncertainty_mbdeff_lower_record[ib-1]);
    h_uncertainty_mbdcs_upper->SetBinContent(ib, uncertainty_mbdcs_upper_record[ib-1]);
    h_uncertainty_mbdcs_lower->SetBinContent(ib, -uncertainty_mbdcs_lower_record[ib-1]);
    h_uncertainty_stat_upper->SetBinContent(ib, uncertainty_stat_upper_record[ib-1]);
    h_uncertainty_stat_lower->SetBinContent(ib, -uncertainty_stat_lower_record[ib-1]);
    h_uncertainty_total_upper->SetBinContent(ib, uncertainty_total_upper_record[ib-1]);
    h_uncertainty_total_lower->SetBinContent(ib, -uncertainty_total_lower_record[ib-1]);
  }
  TCanvas *can_unc = new TCanvas("can_unc", "", 850, 800);
  gStyle->SetPalette(57);
  can_unc->SetTopMargin(0.03);
  can_unc->SetLeftMargin(0.15);
  can_unc->SetBottomMargin(0.12);
  can_unc->SetRightMargin(0.05);
  TH2F *frame_unc = new TH2F("frame_unc", "", 10, 15, 72, 2, -1.8, 2);
  frame_unc->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  frame_unc->GetYaxis()->SetTitle("Relative Uncertainty");
  frame_unc->GetXaxis()->SetTitleOffset(0.98);
  frame_unc->GetYaxis()->SetTitleOffset(1.17);
  frame_unc->GetXaxis()->SetLabelSize(0.045);
  frame_unc->GetYaxis()->SetLabelSize(0.045);
  frame_unc->GetXaxis()->CenterTitle();
  frame_unc->GetYaxis()->CenterTitle();
  frame_unc->Draw();
  h_uncertainty_unfold_upper->SetLineColor(kYellow+1); h_uncertainty_unfold_lower->SetLineColor(kYellow+1); h_uncertainty_unfold_upper->Draw("hist same"); h_uncertainty_unfold_lower->Draw("hist same");
  h_uncertainty_cut_upper->SetLineColor(kAzure+2); h_uncertainty_cut_lower->SetLineColor(kAzure+2); h_uncertainty_cut_upper->Draw("hist same"); h_uncertainty_cut_lower->Draw("hist same");
  h_uncertainty_jes_upper->SetLineColor(kGreen+2); h_uncertainty_jes_lower->SetLineColor(kGreen+2); h_uncertainty_jes_upper->Draw("hist same"); h_uncertainty_jes_lower->Draw("hist same");
  h_uncertainty_jer_upper->SetLineColor(kMagenta); h_uncertainty_jer_lower->SetLineColor(kMagenta); h_uncertainty_jer_upper->Draw("hist same"); h_uncertainty_jer_lower->Draw("hist same");
  h_uncertainty_jeteff_upper->SetLineColor(kOrange+1); h_uncertainty_jeteff_lower->SetLineColor(kOrange+1); h_uncertainty_jeteff_upper->Draw("hist same"); h_uncertainty_jeteff_lower->Draw("hist same");
  h_uncertainty_mbdeff_upper->SetLineColor(kRed); h_uncertainty_mbdeff_lower->SetLineColor(kRed); h_uncertainty_mbdeff_upper->Draw("hist same"); h_uncertainty_mbdeff_lower->Draw("hist same");
  h_uncertainty_mbdcs_upper->SetLineColor(kCyan); h_uncertainty_mbdcs_lower->SetLineColor(kCyan); h_uncertainty_mbdcs_upper->Draw("hist same"); h_uncertainty_mbdcs_lower->Draw("hist same");
  h_uncertainty_stat_upper->SetLineColor(kViolet-1); h_uncertainty_stat_lower->SetLineColor(kViolet-1); h_uncertainty_stat_upper->Draw("hist same"); h_uncertainty_stat_lower->Draw("hist same");
  h_uncertainty_total_upper->SetLineColor(kBlack); h_uncertainty_total_lower->SetLineColor(kBlack);
  h_uncertainty_total_upper->SetMarkerStyle(20); h_uncertainty_total_lower->SetMarkerStyle(20); 
  h_uncertainty_total_upper->Draw("p same"); h_uncertainty_total_lower->Draw("p same");
  myText(0.17, 0.92, 1, "#bf{#it{sPHENIX}} Preliminary", 0.05);
  myText(0.17, 0.86, 1, "p+p#sqrt{s} = 200 GeV,        = 16.6 pb^{-1}", 0.05);
  //desiredHeight = 0.05;
  //TPad *pad_uncertainty = new TPad("pad_sphenix", "", 0.52, 0.86 - desiredHeight/2. - 0.005, 0.52+desiredHeight*aspectRatio, 0.86+desiredHeight/2. - 0.005); pad_uncertainty->Draw(); pad_uncertainty->cd(); img->Draw("X"); can_unc->cd(); can_unc->Update();
  myText(0.17, 0.8, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.05);
  myText(0.17, 0.74, 1, "|#eta^{jet}| < 0.7", 0.05);
  //TLegend *leg_unc1 = new TLegend(0.155, 0.15, 0.5, 0.34);
  //leg_unc1->SetBorderSize(0);
  //leg_unc1->SetFillStyle(0);
  //leg_unc1->AddEntry(h_uncertainty_unfold_upper, "Unfolding", "l");
  //leg_unc1->AddEntry(h_uncertainty_cut_upper, "Background Cut", "l");
  //leg_unc1->AddEntry(h_uncertainty_jes_upper, "Jet Energy Scale", "l");
  //leg_unc1->AddEntry(h_uncertainty_jer_upper, "Jet Energy Resolution", "l");
  //leg_unc1->AddEntry(h_uncertainty_jeteff_upper, "Jet Efficiency", "l");
  //leg_unc1->Draw();
  //TLegend *leg_unc2 = new TLegend(0.505, 0.15, 0.85, 0.30);
  //leg_unc2->SetBorderSize(0);
  //leg_unc2->SetFillStyle(0);
  //leg_unc2->AddEntry(h_uncertainty_mbdeff_upper, "MBD Efficiency", "l");
  //leg_unc2->AddEntry(h_uncertainty_mbdcs_upper, "MBD Cross Section    ", "l");
  //leg_unc2->AddEntry(h_uncertainty_stat_upper, "Statistical", "l");
  //leg_unc2->AddEntry(h_uncertainty_total_upper, "Total Systematic", "p");
  //leg_unc2->Draw();
  TLegend *leg_unc1 = new TLegend(0.15, 0.13, 0.95, 0.37);
  leg_unc1->SetNColumns(2);
  leg_unc1->SetBorderSize(0);
  leg_unc1->SetFillStyle(0);
  leg_unc1->SetTextFont(42);
  leg_unc1->AddEntry(h_uncertainty_cut_upper, "Background Rejection", "l");
  leg_unc1->AddEntry(h_uncertainty_unfold_upper, "Unfolding", "l");
  leg_unc1->AddEntry(h_uncertainty_jes_upper, "Jet Energy Scale", "l");
  leg_unc1->AddEntry(h_uncertainty_jeteff_upper, "Jet Trigger Efficiency", "l");
  leg_unc1->AddEntry(h_uncertainty_jer_upper, "Jet Energy Resolution", "l");
  leg_unc1->AddEntry(h_uncertainty_mbdeff_upper, "MBD Trigger Efficiency", "l");
  leg_unc1->AddEntry(h_uncertainty_mbdcs_upper, "MBD Cross Section", "l");
  leg_unc1->AddEntry(h_uncertainty_stat_upper, "Statistical", "l");
  leg_unc1->AddEntry(h_uncertainty_total_upper, "Total Uncertainty", "p");
  leg_unc1->Draw();
  TLine *line = new TLine(15, 0, 72, 0);
  line->SetLineStyle(10);
  line->Draw("same");
  can_unc->SaveAs("figure/spectrum_uncertainty.pdf");
}

void draw_1D_multiple_plot_uncertainty(std::vector<TH1F*> h_input, float x_min, float x_max, float ratio_min, float ratio_max, std::string xtitle, std::string ytitle, std::vector<std::string> legend, float xstart_legend, float ystart_legend, std::vector<int> color, std::vector<int> markerstyle, std::string output_name) {
  if (h_input.size() <= 1) {
    std::cout << "Error: no enough input histograms for draw_1D_multiple_plot" << std::endl;
    return;
  }
  std::vector<TH1F*> h;
  for (int i = 0; i < h_input.size(); ++i) {
    h.push_back((TH1F*)h_input[i]->Clone(Form("h_%d", i)));
  }
  TCanvas *can = new TCanvas("can", "", 800, 963);
  can->Divide(1, 2);
  gStyle->SetPalette(57);
  TPad *pad_1 = (TPad*)can->cd(1);
  pad_1->SetPad(0, 0.4, 1, 1);
  pad_1->SetTopMargin(0.03);
  pad_1->SetLeftMargin(0.15);
  pad_1->SetBottomMargin(0.035);
  pad_1->SetRightMargin(0.08);
  pad_1->SetLogy();
  for (int ih = 0; ih < h.size(); ++ih) {
    h.at(ih)->SetMarkerStyle(markerstyle.at(ih));
    h.at(ih)->SetMarkerColor(color.at(ih));
    h.at(ih)->SetLineColor(color.at(ih));
    for (int ib = 1; ib <= h.at(1)->GetNbinsX(); ++ib) {
      h.at(ih)->SetBinContent(ib, h.at(ih)->GetBinContent(ib)/ (float)(/*2*TMath::Pi() * */luminosity * delta_eta * h.at(ih)->GetBinWidth(ib)));
      h.at(ih)->SetBinError(ib, h.at(ih)->GetBinError(ib)/ (float)(/*2*TMath::Pi() * */luminosity * delta_eta * h.at(ih)->GetBinWidth(ib)));
    }
  }
  h.at(0)->GetXaxis()->SetTitle(xtitle.c_str());
  h.at(0)->GetYaxis()->SetTitle(ytitle.c_str());
  h.at(0)->GetXaxis()->SetRangeUser(x_min, x_max);
  h.at(0)->GetXaxis()->SetTitleSize(0.065);
  h.at(0)->GetYaxis()->SetTitleSize(0.065);
  h.at(0)->GetXaxis()->SetTitleOffset(1.04);
  h.at(0)->GetYaxis()->SetTitleOffset(1.06);
  h.at(0)->GetXaxis()->SetLabelSize(0.055);
  h.at(0)->GetYaxis()->SetLabelSize(0.055);
  h.at(0)->GetXaxis()->SetLabelOffset(2);
  h.at(0)->GetXaxis()->CenterTitle();
  h.at(0)->GetYaxis()->CenterTitle();
  h.at(0)->Draw();
  for (int i = 0; i < h.size(); ++i) {
    if (i == 0) continue;
    h.at(i)->Draw("same");
  }

  myText(0.4, 0.9, 1, "#bf{#it{sPHENIX}} Preliminary", 0.05);
  myText(0.4, 0.84, 1, "p+p#sqrt{s} = 200 GeV,#kern[+0.7]{#scale[0.5]{#int}}#kern[+0.2]{#it{L}}dt = 16.6 pb^{-1}", 0.05);
  myText(0.4, 0.78, 1, "anti-k_{t} #kern[-0.3]{#it{R}} = 0.4", 0.05);
  myText(0.4, 0.72, 1, "|#eta^{jet}| < 0.7", 0.05);
  for (int i = 0; i < legend.size(); i++) {
    myMarkerLineText(xstart_legend, ystart_legend-i*0.06, 1, color.at(i), markerstyle.at(i), color.at(i), 1, legend.at(i).c_str(), 0.05, true);
  }

  TPad *pad_2 = (TPad*)can->cd(2);
  pad_2->SetPad(0, 0, 1, 0.4);
  pad_2->SetTopMargin(0.03);
  pad_2->SetLeftMargin(0.15);
  pad_2->SetBottomMargin(0.25);
  pad_2->SetRightMargin(0.08);
  std::vector<TH1F*> h_ratio;
  for (int i = 0; i < h.size()-1; ++i) {
    TH1F *h_temp = (TH1F*)h.at(i+1)->Clone(Form("h_ratio_%d", i));
    h_temp->Divide(h.at(0));
    h_ratio.push_back(h_temp);
    h_ratio.at(i)->SetMarkerStyle(markerstyle.at(i+1));
    h_ratio.at(i)->SetMarkerColor(color.at(i+1));
    h_ratio.at(i)->SetLineColor(color.at(i+1));
  }
  h_ratio.at(0)->GetXaxis()->SetTitle(xtitle.c_str());
  h_ratio.at(0)->GetYaxis()->SetTitle("Ratio");
  h_ratio.at(0)->GetXaxis()->SetRangeUser(x_min, x_max);
  h_ratio.at(0)->GetYaxis()->SetRangeUser(ratio_min, ratio_max);
  h_ratio.at(0)->GetXaxis()->CenterTitle();
  h_ratio.at(0)->GetYaxis()->CenterTitle();
  h_ratio.at(0)->GetXaxis()->SetTitleOffset(0.95);
  h_ratio.at(0)->GetYaxis()->SetTitleOffset(h.at(0)->GetYaxis()->GetTitleOffset()*4/6.);
  //h_ratio.at(0)->GetYaxis()->SetLabelOffset(h.at(0)->GetYaxis()->GetLabelOffset()*4/6.);
  h_ratio.at(0)->GetXaxis()->SetLabelSize(h.at(0)->GetXaxis()->GetLabelSize()*6/4.);
  h_ratio.at(0)->GetYaxis()->SetLabelSize(h.at(0)->GetYaxis()->GetLabelSize()*6/4.);
  h_ratio.at(0)->GetXaxis()->SetTitleSize(h.at(0)->GetXaxis()->GetTitleSize()*6/4.);
  h_ratio.at(0)->GetYaxis()->SetTitleSize(h.at(0)->GetYaxis()->GetTitleSize()*6/4.);
  h_ratio.at(0)->Draw();
  for (int i = 0; i < h_ratio.size(); ++i) {
    if (i == 0) continue;
    h_ratio.at(i)->Draw("same");
  }

  TLine *line;
  line = new TLine(x_min, 1, x_max, 1);
  line->SetLineColor(kBlack);
  line->SetLineStyle(3);
  line->Draw("same");

  can->SaveAs(output_name.c_str());
  delete can;
}

void set_tpad(TPad* pad, TH2D* h, TBox* box1, TBox* box2) {
  pad->SetFillColor(0);
  pad->SetBorderMode(0);
  pad->SetBorderSize(0);
  pad->SetTickx(0);
  pad->SetTicky(0);
  pad->SetTopMargin(0);
  pad->SetBottomMargin(0);
  pad->SetLeftMargin(0);
  pad->SetRightMargin(0);

  h->SetStats(false);
  h->SetMinimum(0);
  h->SetMaximum(1);
  h->GetXaxis()->SetLabelSize(0);
  h->GetXaxis()->SetTickLength(0);
  h->GetXaxis()->SetTitle("");
  h->GetYaxis()->SetLabelSize(0);
  h->GetYaxis()->SetTickLength(0);
  h->GetYaxis()->SetTitle("");
  h->GetXaxis()->SetAxisColor(0);
  h->GetYaxis()->SetAxisColor(0);
  h->SetLineColor(0);
  h->SetFillColor(0);

  Int_t black_white[2] = {kWhite, kBlack};
  gStyle->SetPalette(2, black_white);
  h->SetContour(2);

  pad->cd();
  h->Draw("col");
  box1->Draw("same");
  box2->Draw("same");
  pad->Update();

  //delete coverRight;
  //delete coverTop;
}
