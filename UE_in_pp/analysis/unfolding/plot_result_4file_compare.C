#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TColor.h>
#include <TProfile.h>
#include <TPad.h>
#include <TLine.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <unfold_Def.h>
#include "RooUnfoldResponse.h"

struct CompareSample {
  std::string tag;
  std::string file_path;
  std::string reco_hist_name;
  std::string truth_hist_name;
  TFile* file = nullptr;
  TH2D* reco_raw = nullptr;
  TH2D* truth_raw = nullptr;
  TH2D* reco_var = nullptr;
  TH2D* truth_var = nullptr;
};

TH2D* make_variable_binned_copy(const TH2D* input, bool is_reco, const std::string& output_name)
{
  if (!input) return nullptr;

  const int nx = is_reco ? calibnpt : truthnpt;
  const int ny = is_reco ? calibnet : truthnet;
  const double* xbins = is_reco ? calibptbins : truthptbins;
  const double* ybins = is_reco ? calibetbins : truthetbins;

  TH2D* out = new TH2D(output_name.c_str(), "", nx, xbins, ny, ybins);
  out->SetDirectory(nullptr);

  const int copy_nx = std::min(nx, input->GetNbinsX());
  const int copy_ny = std::min(ny, input->GetNbinsY());
  if (copy_nx != nx || copy_ny != ny) {
    std::cout << "Warning: bin-count mismatch for " << output_name
              << " (source: " << input->GetNbinsX() << "x" << input->GetNbinsY()
              << ", target: " << nx << "x" << ny << "). Copying overlapping bins only."
              << std::endl;
  }

  for (int ix = 1; ix <= copy_nx; ++ix) {
    for (int iy = 1; iy <= copy_ny; ++iy) {
      out->SetBinContent(ix, iy, input->GetBinContent(ix, iy));
      out->SetBinError(ix, iy, input->GetBinError(ix, iy));
    }
  }
  return out;
}

void scale_to_integral(TH1* h, double target_integral)
{
  if (!h) return;
  const double current = h->Integral();
  if (current > 0.0 && target_integral > 0.0) {
    h->Scale(target_integral / current);
  }
}

void style_hist(TH1* h, int color, int marker)
{
  if (!h) return;
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetLineWidth(2);
  h->SetStats(0);
}

void draw_overlay(const std::vector<TH1*>& hists,
                  const std::vector<std::string>& labels,
                  const std::string& xtitle,
                  const std::string& ytitle,
                  const std::string& canvas_name,
                  const std::string& output_name,
                  bool is_profile = false)
{
  if (hists.empty()) return;

  TCanvas* c = new TCanvas(canvas_name.c_str(), "", 900, 700);
  c->cd();

  TPad* pad_top = new TPad(Form("%s_top", canvas_name.c_str()), "", 0.0, 0.40, 1.0, 1.0);
  TPad* pad_ratio = new TPad(Form("%s_ratio", canvas_name.c_str()), "", 0.0, 0.0, 1.0, 0.40);
  pad_top->SetBottomMargin(0.02);
  pad_ratio->SetTopMargin(0.03);
  pad_ratio->SetBottomMargin(0.40);
  pad_ratio->SetGridy();
  pad_top->Draw();
  pad_ratio->Draw();

  TLegend* leg = new TLegend(0.52, 0.67, 0.9, 0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.032);

  int ref_idx = -1;
  for (std::size_t i = 0; i < hists.size(); ++i) {
    if (hists[i]) {
      ref_idx = static_cast<int>(i);
      break;
    }
  }
  if (ref_idx < 0) return;

  TH1* ref_hist = hists[ref_idx];

  pad_top->cd();
  bool first_drawn = false;
  for (std::size_t i = 0; i < hists.size(); ++i) {
    if (!hists[i]) continue;
    hists[i]->GetXaxis()->SetTitle(xtitle.c_str());
    hists[i]->GetYaxis()->SetTitle(ytitle.c_str());
    hists[i]->GetXaxis()->SetLabelSize(0.0);
    hists[i]->GetXaxis()->SetTitleSize(0.0);
    hists[i]->SetTitle("");
    const char* draw_opt = (!first_drawn) ? (is_profile ? "E1" : "hist E1")
                                           : (is_profile ? "E1 same" : "hist E1 same");
    hists[i]->Draw(draw_opt);
    leg->AddEntry(hists[i], labels[i].c_str(), "lp");
    first_drawn = true;
  }

  leg->Draw();

  pad_ratio->cd();
  std::vector<TH1*> ratio_hists;
  ratio_hists.reserve(hists.size());
  double ratio_min = std::numeric_limits<double>::max();
  double ratio_max = std::numeric_limits<double>::lowest();

  for (std::size_t i = 0; i < hists.size(); ++i) {
    if (!hists[i]) continue;
    TH1* ratio = dynamic_cast<TH1*>(hists[i]->Clone(Form("%s_ratio_hist_%zu", canvas_name.c_str(), i)));
    if (!ratio) continue;
    ratio->SetDirectory(nullptr);
    ratio->Divide(ref_hist);
    ratio->GetYaxis()->SetTitle("Ratio");
    ratio->GetXaxis()->SetTitle(xtitle.c_str());
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetXaxis()->SetTitleSize(0.11);
    ratio->GetXaxis()->SetLabelSize(0.10);
    ratio->GetYaxis()->SetTitleSize(0.085);
    ratio->GetYaxis()->SetTitleOffset(0.55);
    ratio->GetYaxis()->SetLabelSize(0.09);
    ratio->SetLineWidth(2);
    ratio_hists.push_back(ratio);

    for (int ibin = 1; ibin <= ratio->GetNbinsX(); ++ibin) {
      const double den = ref_hist->GetBinContent(ibin);
      if (den == 0.0) continue;
      const double val = ratio->GetBinContent(ibin);
      const double err = ratio->GetBinError(ibin);
      ratio_min = std::min(ratio_min, val - err);
      ratio_max = std::max(ratio_max, val + err);
    }
  }

  if (!ratio_hists.empty()) {
    if (!(ratio_min < ratio_max)) {
      ratio_min = 0.5;
      ratio_max = 1.5;
    }
    const double pad = 0.15 * (ratio_max - ratio_min);
    ratio_min = std::max(0.0, ratio_min - pad);
    ratio_max += pad;

    ratio_hists[0]->SetMinimum(0.8);
    ratio_hists[0]->SetMaximum(1.2);
    ratio_hists[0]->Draw("E1");
    for (std::size_t i = 1; i < ratio_hists.size(); ++i) {
      ratio_hists[i]->Draw("E1 same");
    }

    TLine* one = new TLine(ratio_hists[0]->GetXaxis()->GetXmin(), 1.0,
                           ratio_hists[0]->GetXaxis()->GetXmax(), 1.0);
    one->SetLineStyle(2);
    one->SetLineWidth(2);
    one->Draw("same");
  }

  c->SaveAs(output_name.c_str());
}

void plot_result_4file_compare()
{
  gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
  gROOT->ProcessLine("SetsPhenixStyle()");
  gStyle->SetOptStat(0);

  std::vector<CompareSample> samples = {
      {"File1: sim run21 test iter 1",
       "analysis_sim_run21_test_output/output_dijet_bkg_cut_sim_iter_1.root",
       "calib_dijet",
       "h_truth_calib_dijet"},
      {"File2: sim run28 test iter 1",
       "analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_1.root",
       "calib_dijet",
       "h_truth_calib_dijet"},
       {"File3: sim run21 test iter 2",
       "analysis_sim_run21_test_output/output_dijet_bkg_cut_sim_iter_2.root",
       "calib_dijet_trim_10",
       "h_truth_calib_dijet_trim_10"},
       {"File4: sim run28 test iter 2",
       "analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_2.root",
       "calib_dijet_trim_10",
       "h_truth_calib_dijet_trim_10"},
       {"File5: sim run21 test iter 3",
       "analysis_sim_run21_test_output/output_dijet_bkg_cut_sim_iter_3.root",
       "calib_dijet_reweight_trim_10",
       "h_truth_calib_dijet_reweight_trim_10"},
       {"File6: sim run28 test iter 3",
       "analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_3.root",
       "calib_dijet_reweight_trim_10",
       "h_truth_calib_dijet_reweight_trim_10"},
      {"File7: data run21 test",
       "run21_test_output_files/output_mbd_correct_w_syst_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root",
       "h_calibjet_pt_dijet_eff",
       "h_unfold_calib_dijet_reweight_trim_10_2_etEffCorrected"},
      {"File8: data run28",
       "run28_output_files/output_mbd_correct_w_syst_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root",
       "h_calibjet_pt_dijet_eff",
       "h_unfold_calib_dijet_reweight_trim_10_2_etEffCorrected"}};

  const int colors[8] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kYellow + 1, kMagenta + 1, kCyan + 1, kOrange + 1};
  const int markers[8] = {24, 24, 24, 24, 24, 24, 20, 20};

  for (std::size_t i = 0; i < samples.size(); ++i) {
    samples[i].file = TFile::Open(samples[i].file_path.c_str(), "READ");
    if (!samples[i].file || samples[i].file->IsZombie()) {
      std::cout << "Error: failed to open " << samples[i].file_path << std::endl;
      return;
    }

    if (i < 6) {
      RooUnfoldResponse* sample_respmatrix = dynamic_cast<RooUnfoldResponse*>(samples[i].file->Get(Form("h_respmatrix_%s", samples[i].reco_hist_name.c_str())));
      samples[i].reco_raw = dynamic_cast<TH2D*>(sample_respmatrix->Hmeasured());
      samples[i].truth_raw = dynamic_cast<TH2D*>(sample_respmatrix->Htruth());
    } else {
      samples[i].reco_raw = dynamic_cast<TH2D*>(samples[i].file->Get(samples[i].reco_hist_name.c_str()));
      samples[i].truth_raw = dynamic_cast<TH2D*>(samples[i].file->Get(samples[i].truth_hist_name.c_str()));
    }

    if (!samples[i].reco_raw) {
      std::cout << "Error: missing reco TH2D " << samples[i].reco_hist_name
                << " in " << samples[i].file_path << std::endl;
      return;
    }
    if (!samples[i].truth_raw) {
      std::cout << "Error: missing truth TH2D " << samples[i].truth_hist_name
                << " in " << samples[i].file_path << std::endl;
      return;
    }

    samples[i].reco_var = make_variable_binned_copy(samples[i].reco_raw, true, Form("h_reco_var_%zu", i));
    samples[i].truth_var = make_variable_binned_copy(samples[i].truth_raw, false, Form("h_truth_var_%zu", i));
  }

  std::vector<TH1*> reco_projx;
  std::vector<TH1*> reco_projy;
  std::vector<TH1*> reco_profilex;
  std::vector<TH1*> truth_projx;
  std::vector<TH1*> truth_projy;
  std::vector<TH1*> truth_profilex;
  std::vector<std::string> labels;

  for (std::size_t i = 0; i < samples.size(); ++i) {
    labels.push_back(samples[i].tag);

    TH1D* rx = samples[i].reco_var->ProjectionX(Form("h_reco_projx_%zu", i));
    TH1D* ry = samples[i].reco_var->ProjectionY(Form("h_reco_projy_%zu", i));
    TProfile* rp = samples[i].reco_var->ProfileX(Form("h_reco_prof_%zu", i));
    TH1D* tx = samples[i].truth_var->ProjectionX(Form("h_truth_projx_%zu", i));
    TH1D* ty = samples[i].truth_var->ProjectionY(Form("h_truth_projy_%zu", i));
    TProfile* tp = samples[i].truth_var->ProfileX(Form("h_truth_prof_%zu", i));

    style_hist(rx, colors[i], markers[i]);
    style_hist(ry, colors[i], markers[i]);
    style_hist(rp, colors[i], markers[i]);
    style_hist(tx, colors[i], markers[i]);
    style_hist(ty, colors[i], markers[i]);
    style_hist(tp, colors[i], markers[i]);

    reco_projx.push_back(rx);
    reco_projy.push_back(ry);
    reco_profilex.push_back(rp);
    truth_projx.push_back(tx);
    truth_projy.push_back(ty);
    truth_profilex.push_back(tp);
  }

  // Additional fixed simulation scaling.
  const double sim1_scale = 3.997;
  const double sim2_scale = 1.4903 * 0.81608598;
  reco_projx[0]->Scale(sim1_scale);
  reco_projy[0]->Scale(sim1_scale);
  truth_projx[0]->Scale(sim1_scale);
  truth_projy[0]->Scale(sim1_scale);
  reco_projx[1]->Scale(sim2_scale);
  reco_projy[1]->Scale(sim2_scale);
  truth_projx[1]->Scale(sim2_scale);
  truth_projy[1]->Scale(sim2_scale);
  reco_projx[2]->Scale(sim1_scale);
  reco_projy[2]->Scale(sim1_scale);
  truth_projx[2]->Scale(sim1_scale);
  truth_projy[2]->Scale(sim1_scale);
  reco_projx[3]->Scale(sim2_scale);
  reco_projy[3]->Scale(sim2_scale);
  truth_projx[3]->Scale(sim2_scale);
  truth_projy[3]->Scale(sim2_scale);
  reco_projx[4]->Scale(sim1_scale);
  reco_projy[4]->Scale(sim1_scale);
  truth_projx[4]->Scale(sim1_scale);
  truth_projy[4]->Scale(sim1_scale);
  reco_projx[5]->Scale(sim2_scale);
  reco_projy[5]->Scale(sim2_scale);
  truth_projx[5]->Scale(sim2_scale);
  truth_projy[5]->Scale(sim2_scale);

  // Normalize file 3 and file 4 projections to the integral of scaled file 1.
  // (No normalization applied to ProfileX.)
  const double reco_x_ref = reco_projx[1]->Integral();
  const double reco_y_ref = reco_projy[1]->Integral();
  const double truth_x_ref = truth_projx[1]->Integral();
  const double truth_y_ref = truth_projy[1]->Integral();

  scale_to_integral(reco_projx[6], reco_x_ref);
  scale_to_integral(reco_projx[7], reco_x_ref);
  scale_to_integral(reco_projy[6], reco_y_ref);
  scale_to_integral(reco_projy[7], reco_y_ref);
  scale_to_integral(truth_projx[6], truth_x_ref);
  scale_to_integral(truth_projx[7], truth_x_ref);
  scale_to_integral(truth_projy[6], truth_y_ref);
  scale_to_integral(truth_projy[7], truth_y_ref);

  draw_overlay(reco_projx, labels, "Calibrated p_{T,lead} [GeV]", "Counts", "c_reco_projx",
               "run21_test_output_files/compare4_reco_projectionX.png");
  draw_overlay(reco_projy, labels, "Reco #SigmaE_{T} [GeV]", "Counts", "c_reco_projy",
               "run21_test_output_files/compare4_reco_projectionY.png");
  draw_overlay(reco_profilex, labels, "Calibrated p_{T,lead} [GeV]", "<Reco #SigmaE_{T}> [GeV]", "c_reco_profilex",
               "run21_test_output_files/compare4_reco_profileX.png", true);

  draw_overlay(truth_projx, labels, "Truthp_{T,lead} [GeV]", "Counts", "c_truth_projx",
               "run21_test_output_files/compare4_truth_projectionX.png");
  draw_overlay(truth_projy, labels, "Truth #SigmaE_{T} [GeV]", "Counts", "c_truth_projy",
               "run21_test_output_files/compare4_truth_projectionY.png");
  draw_overlay(truth_profilex, labels, "Truth p_{T,lead} [GeV]", "<Truth #SigmaE_{T}> [GeV]", "c_truth_profilex",
               "run21_test_output_files/compare4_truth_profileX.png", true);

  std::cout << "Saved output plots in run21_test_output_files/ (compare4_*). " << std::endl;
}
