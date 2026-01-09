#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include "/sphenix/user/hanpuj/CaloDataAna24_skimmed/src/draw_template.C"

void get_purityefficiency() {
  TFile *f_in = new TFile("output_sim.root", "READ");
  if (!f_in) {
    std::cout << "Error: cannot open output_sim.root" << std::endl;
    return;
  }

  TH1D* h_measure_reco_all_record = (TH1D*)f_in->Get("h_measure_reco_all_record");
  TH1D* h_fake_reco_all_record = (TH1D*)f_in->Get("h_fake_reco_all_record");
  TH1D* h_truth_reco_all_record = (TH1D*)f_in->Get("h_truth_reco_all_record");
  TH1D* h_miss_reco_all_record = (TH1D*)f_in->Get("h_miss_reco_all_record");
  TH2D* h_respmatrix_reco_all_record = (TH2D*)f_in->Get("h_respmatrix_reco_all_record");

  TH1D* h_measure_calib_all_record = (TH1D*)f_in->Get("h_measure_calib_all_record");
  TH1D* h_fake_calib_all_record = (TH1D*)f_in->Get("h_fake_calib_all_record");
  TH1D* h_truth_calib_all_record = (TH1D*)f_in->Get("h_truth_calib_all_record");
  TH1D* h_miss_calib_all_record = (TH1D*)f_in->Get("h_miss_calib_all_record");
  TH2D* h_respmatrix_calib_all_record = (TH2D*)f_in->Get("h_respmatrix_calib_all_record");

  double recoptbins[] = {12, 15, 18, 21, 24, 28, 32, 36, 40, 45, 50, 55, 60};
  double truthptbins[] = {6, 9, 12, 15, 18, 21, 24, 28, 32, 36, 40, 45, 50, 55, 60, 66, 72};

  int reconpt = sizeof(recoptbins) / sizeof(recoptbins[0]) - 1;
  int truthnpt = sizeof(truthptbins) / sizeof(truthptbins[0]) - 1;
  TH1D* h_purity_reco_all = new TH1D("h_purity_reco_all", "Reco Jet Purity; p_{T} [GeV]; Purity", reconpt, recoptbins);
  TH1D* h_purity_calib_all = new TH1D("h_purity_calib_all", "Calib Jet Purity; p_{T} [GeV]; Purity", reconpt, recoptbins);
  TH1D* h_efficiency_reco_all = new TH1D("h_efficiency_reco_all", "Reco Jet Efficiency; p_{T} [GeV]; Efficiency", truthnpt, truthptbins);
  TH1D* h_efficiency_calib_all = new TH1D("h_efficiency_calib_all", "Calib Jet Efficiency; p_{T} [GeV]; Efficiency", truthnpt, truthptbins);

  std::cout << std::endl;
  std::cout << "Purity:" << std::endl;
  std::cout << "Reco     \tCalib" << std::endl;
  for (int ib = 1; ib <= reconpt; ++ib) {
    double measure_reco_all = h_measure_reco_all_record->GetBinContent(ib);
    double fake_reco_all = h_fake_reco_all_record->GetBinContent(ib);
    double purity_reco_all = 1 - fake_reco_all / (double)measure_reco_all;
    std::cout << purity_reco_all << " \t";
    h_purity_reco_all->SetBinContent(ib, purity_reco_all);

    double measure_calib_all = h_measure_calib_all_record->GetBinContent(ib);
    double fake_calib_all = h_fake_calib_all_record->GetBinContent(ib);
    double purity_calib_all = 1 - fake_calib_all / (double)measure_calib_all;
    std::cout << purity_calib_all << std::endl;
    h_purity_calib_all->SetBinContent(ib, purity_calib_all);
  }

  std::cout << std::endl;
  std::cout << "Efficiency:" << std::endl;
  std::cout << "Reco     \tCalib" << std::endl;
  for (int ib = 1; ib <= truthnpt; ++ib) {
    double truth_reco_all = h_truth_reco_all_record->GetBinContent(ib);
    double miss_reco_all = h_miss_reco_all_record->GetBinContent(ib);
    double efficiency_reco_all = 1 - miss_reco_all / (double)truth_reco_all;
    std::cout << efficiency_reco_all << " \t";
    h_efficiency_reco_all->SetBinContent(ib, efficiency_reco_all);

    double truth_calib_all = h_truth_calib_all_record->GetBinContent(ib);
    double miss_calib_all = h_miss_calib_all_record->GetBinContent(ib);
    double efficiency_calib_all = 1 - miss_calib_all / (double)truth_calib_all;
    std::cout << efficiency_calib_all << std::endl;
    h_efficiency_calib_all->SetBinContent(ib, efficiency_calib_all);
  }
  std::cout << std::endl;

  TH1D* h_measure_reco_all_proj = (TH1D*)h_respmatrix_reco_all_record->ProjectionX("h_measure_reco_all_proj");
  TH1D* h_measure_calib_all_proj = (TH1D*)h_respmatrix_calib_all_record->ProjectionX("h_measure_calib_all_proj");
  TH1D* h_truth_reco_all_proj = (TH1D*)h_respmatrix_reco_all_record->ProjectionY("h_truth_reco_all_proj");
  TH1D* h_truth_calib_all_proj = (TH1D*)h_respmatrix_calib_all_record->ProjectionY("h_truth_calib_all_proj");

  std::cout << std::endl;
  std::cout << "Check measure:" << std::endl;
  std::cout << "Meas     \tFake     \tCalc     \tProj" << std::endl;
  for (int ib = 1; ib <= reconpt; ++ib) {
    std::cout << h_measure_reco_all_record->GetBinContent(ib) << " \t";
    std::cout << h_fake_reco_all_record->GetBinContent(ib) << " \t";
    std::cout << h_measure_reco_all_record->GetBinContent(ib) - h_fake_reco_all_record->GetBinContent(ib) << " \t";
    std::cout << h_measure_reco_all_proj->GetBinContent(ib) << std::endl;
  }

  std::cout << std::endl;
  std::cout << "Check truth:" << std::endl;
  std::cout << "Truth    \tMiss     \tCalc     \tProj" << std::endl;
  for (int ib = 1; ib <= truthnpt; ++ib) {
    std::cout << h_truth_reco_all_record->GetBinContent(ib) << " \t";
    std::cout << h_miss_reco_all_record->GetBinContent(ib) << " \t";
    std::cout << h_truth_reco_all_record->GetBinContent(ib) - h_miss_reco_all_record->GetBinContent(ib) << " \t";
    std::cout << h_truth_reco_all_proj->GetBinContent(ib) << std::endl;
  }

  std::vector<TH1F*> h_input;
  std::vector<int> color;
  std::vector<int> markerstyle;
  std::vector<std::string> text;
  std::vector<std::string> legend;

  h_input.push_back((TH1F*)h_measure_reco_all_record);
  h_input.push_back((TH1F*)h_measure_reco_all_proj);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("|z_{vertex}| < 30 cm");
  text.push_back("#eta^{jet} within detector range");
  text.push_back("E^{jet} > 0 GeV");
  text.push_back("1 to 1 match with dR < 0.3");
  legend.push_back("Reco jet spectrum");
  legend.push_back("Matched reco jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 12, 60, false,
                              false, 0, 0.5, true,
                              true, 0.75, 1.1,
                              true, "p_{T}^{reco jet} [GeV]", "Arbitrary Unit", "Purity", 1,
                              false, "Run 21 MB + Jet10GeV + Jet30GeV simulation data set",
                              true, text, 0.55, 0.81, 0.04,
                              true, legend, 0.58, 0.61, 0.04,
                              "figure/purity_reco.png");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_measure_calib_all_record);
  h_input.push_back((TH1F*)h_measure_calib_all_proj);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("|z_{vertex}| < 30 cm");
  text.push_back("#eta^{jet} within detector range");
  text.push_back("E^{jet} > 0 GeV");
  text.push_back("1 to 1 match with dR < 0.3");
  text.push_back("With JES calibration");
  legend.push_back("Reco jet spectrum");
  legend.push_back("Matched reco jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 12, 60, false,
                              false, 0, 0.5, true,
                              true, 0.75, 1.1,
                              true, "p_{T}^{reco jet} [GeV]", "Arbitrary Unit", "Purity", 1,
                              false, "Run 21 MB + Jet10GeV + Jet30GeV simulation data set",
                              true, text, 0.55, 0.81, 0.04,
                              true, legend, 0.58, 0.56, 0.04,
                              "figure/purity_calib.png");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_truth_reco_all_record);
  h_input.push_back((TH1F*)h_truth_reco_all_proj);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("|z_{vertex}| < 30 cm");
  text.push_back("#eta^{jet} within detector range");
  text.push_back("E^{jet} > 0 GeV");
  text.push_back("1 to 1 match with dR < 0.3");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Matched truth jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 12, 60, false,
                              false, 0, 0.5, true,
                              true, 0, 1.2,
                              true, "p_{T}^{truth jet} [GeV]", "Arbitrary Unit", "Efficiency", 1,
                              false, "Run 21 MB + Jet10GeV + Jet30GeV simulation data set",
                              true, text, 0.55, 0.81, 0.04,
                              true, legend, 0.58, 0.61, 0.04,
                              "figure/efficiency_reco.png");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_truth_calib_all_record);
  h_input.push_back((TH1F*)h_truth_calib_all_proj);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("|z_{vertex}| < 30 cm");
  text.push_back("#eta^{jet} within detector range");
  text.push_back("E^{jet} > 0 GeV");
  text.push_back("1 to 1 match with dR < 0.3");
  text.push_back("With JES calibration");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Matched truth jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 12, 60, false,
                              false, 0, 0.5, true,
                              true, 0, 1.2,
                              true, "p_{T}^{truth jet} [GeV]", "Arbitrary Unit", "Efficiency", 1,
                              false, "Run 21 MB + Jet10GeV + Jet30GeV simulation data set",
                              true, text, 0.55, 0.81, 0.04,
                              true, legend, 0.58, 0.56, 0.04,
                              "figure/efficiency_calib.png");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  TFile *f_out = new TFile("output_purityefficiency.root", "RECREATE");
  h_purity_reco_all->Write();
  h_purity_calib_all->Write();
  h_efficiency_reco_all->Write();
  h_efficiency_calib_all->Write();
  f_out->Close();
}