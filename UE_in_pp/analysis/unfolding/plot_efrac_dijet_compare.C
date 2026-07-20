#include <TCanvas.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TProfile.h>
#include <TROOT.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <unfold_Def.h>

namespace
{

  std::vector<std::string> GetSystList()
  {
    return {
        "calib_dijet_reweight_trim_10_2_etEffCorrected",
        "calib_dijet_jesdown_reweight_trim_10_2_etEffCorrected",
        "calib_dijet_jesup_reweight_trim_10_2_etEffCorrected",
        "calib_dijet_jerdown_reweight_trim_10_2_etEffCorrected",
        "calib_dijet_jerup_reweight_trim_10_2_etEffCorrected",
        "calib_dijet_trim_10_2_etEffCorrected",
        "h_calibjet_pt_dijet_effup_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_dijet_effdown_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_dijet_pu_correct_et_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_timingeffup_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_emcal_scale_up_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_emcal_scale_down_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_ihcal_scale_up_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_ihcal_scale_down_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_ohcal_scale_up_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_ohcal_scale_down_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_had_resp_up_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_had_resp_down_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_phi_res_calib_dijet_reweight_trim_10_2_etEffCorrected",
        "calib_dijet_clus_smear_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_2sigma_noise_calib_dijet_2sigma_noise_reweight_trim_10_2_etEffCorrected",
        "h_calibjet_pt_4sigma_noise_calib_dijet_4sigma_noise_reweight_trim_10_2_etEffCorrected",
        "calib_dijet_ohcal_mc_data_var_reweight_trim_10_2_etEffCorrected",
        "calib_dijet_reweight_trim_10_2_var_up_etEffCorrected",
        "calib_dijet_reweight_trim_10_2_var_down_etEffCorrected",
        "calib_dijet_reweight_trim_10_2_etEffCorrected"};
  }

  TGraphAsymmErrors *BuildTotalSystBand(TH1D *h_nom, const std::vector<TH1D *> &h_vars, int color)
  {
    const int n = h_nom->GetNbinsX();
    std::vector<double> x(n), ex(n), y(n), eyl(n,0.0), eyh(n,0.0);

    for (int b = 1; b <= n; ++b)
    {
      x[b - 1] = h_nom->GetBinCenter(b);
      ex[b - 1] = 0.5 * h_nom->GetBinWidth(b);
      y[b - 1] = h_nom->GetBinContent(b);

      double up2 = 0.0, dn2 = 0.0;
      for (size_t i = 1; i < h_vars.size(); ++i)
      {
        const double delta = h_vars[i]->GetBinContent(b) - y[b-1];
        if (delta >= 0)
          up2 += delta * delta;
        else
          dn2 += delta * delta;
      }
      eyh[b - 1] = std::sqrt(up2);
      eyl[b - 1] = std::sqrt(dn2);
    }

    auto* g = new TGraphAsymmErrors(n, x.data(), y.data(), ex.data(), ex.data(), eyl.data(), eyh.data());
    g->SetFillColorAlpha(color, 0.50);
    g->SetFillStyle(1001);
    g->SetLineWidth(0);
    return g;
  };

  TH1D *BuildProfileHistogram(TFile *f, const std::string &syst_name, const std::string &out_name)
  {
    TH2D *h2_uniform = dynamic_cast<TH2D *>(f->Get(("h_unfold_" + syst_name).c_str()));
    if (!h2_uniform)
    {
      std::cerr << "Missing histogram h_unfold_" << syst_name << std::endl;
      return nullptr;
    }

    TH2D *h2_var = new TH2D((out_name + "_var").c_str(), "", truthnpt, truthptbins, truthnet, truthetbins);
    for (int ix = 1; ix <= h2_uniform->GetNbinsX(); ++ix)
    {
      for (int iy = 1; iy <= h2_uniform->GetNbinsY(); ++iy)
      {
        h2_var->SetBinContent(ix, iy, h2_uniform->GetBinContent(ix, iy));
        h2_var->SetBinError(ix, iy, h2_uniform->GetBinError(ix, iy));
      }
    }

    TProfile *prof = h2_var->ProfileX((out_name + "_prof").c_str(), 1, h2_var->GetNbinsY() - 1);
    const int nbins = prof->GetNbinsX();
    std::vector<double> edges(nbins + 1, 0.0);
    for (int b = 0; b < nbins; ++b)
      edges[b] = prof->GetXaxis()->GetBinLowEdge(b + 1);
    edges[nbins] = prof->GetXaxis()->GetBinUpEdge(nbins);

    TH1D *out = new TH1D(out_name.c_str(), "", nbins, edges.data());
    for (int b = 1; b <= nbins; ++b)
    {
      out->SetBinContent(b, prof->GetBinContent(b));
      out->SetBinError(b, prof->GetBinError(b));
    }
    return out;
  }

  bool Valid(double a, double b) { return (a > 0.0 && b > 0.0); }

}  // namespace

void plot_efrac_dijet_compare(
    const char *dijet_input = "sphenix_primary_run28_output_files/output_mbd_correct_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root",
    const char *efrac_input = "sphenix_primary_run28_output_files/output_mbd_correct_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root",
    const char *dijet_herwig_input = "sphenix_primary_run28_output_files/output_mbd_correct_herwig_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root",
    const char *efrac_herwig_input = "sphenix_primary_run28_output_files/output_mbd_correct_herwig_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root",
    const char *output_name = "sphenix_primary_run28_output_files/h_compare_efrac_over_dijet_shape_only_unc_figure.pdf")
{
  gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
  gROOT->ProcessLine("SetsPhenixStyle()");

  TFile *f_dijet = TFile::Open(dijet_input);
  TFile *f_efrac = TFile::Open(efrac_input);
  TFile *f_dijet_herwig = TFile::Open(dijet_herwig_input);
  TFile *f_efrac_herwig = TFile::Open(efrac_herwig_input);
  if (!f_dijet || f_dijet->IsZombie() || !f_efrac || f_efrac->IsZombie() || !f_dijet_herwig || f_dijet_herwig->IsZombie() || !f_efrac_herwig || f_efrac_herwig->IsZombie())
  {
    std::cerr << "Could not open input files" << std::endl;
    return;
  }

  const std::vector<std::string> syst = GetSystList();
  std::vector<TH1D *> h_dijet, h_efrac;
  h_dijet.reserve(syst.size());
  h_efrac.reserve(syst.size());

  for (size_t i = 0; i < syst.size()-1; ++i)
  {
    TH1D *hd = BuildProfileHistogram(f_dijet, syst[i], Form("h_dijet_%zu", i));
    TH1D *he = BuildProfileHistogram(f_efrac, syst[i], Form("h_efrac_%zu", i));
    if (!hd || !he)
    {
      std::cerr << "Failed at syst index " << i << std::endl;
      return;
    }
    h_dijet.push_back(hd);
    h_efrac.push_back(he);
  }
  TH1D *hd_herwig = BuildProfileHistogram(f_dijet_herwig, syst[syst.size()-1], Form("h_dijet_herwig_%zu", syst.size()-1));
  TH1D *he_herwig = BuildProfileHistogram(f_efrac_herwig, syst[syst.size()-1], Form("h_efrac_herwig_%zu", syst.size()-1));
  if (!hd_herwig || !he_herwig)
  {
    std::cerr << "Failed at syst index " << syst.size()-1 << std::endl;
    return;
  }
  h_dijet.push_back(hd_herwig);
  h_efrac.push_back(he_herwig);

  for (size_t i = 0; i < syst.size(); ++i)
  {
    h_dijet[i]->Scale(3.0 / (2.2*2*M_PI));
    h_efrac[i]->Scale(3.0 / (2.2*2*M_PI));
  }

  TH1D *h_dijet_nom = h_dijet[0];
  TH1D *h_efrac_nom = h_efrac[0];
  const int nbins = h_dijet_nom->GetNbinsX();

  TH1D *h_ratio = dynamic_cast<TH1D *>(h_efrac_nom->Clone("h_ratio_efrac_over_dijet"));
  for (int b = 1; b <= nbins; ++b)
  {
    const double e = h_efrac_nom->GetBinContent(b), d = h_dijet_nom->GetBinContent(b);
    const double e_stat = h_efrac_nom->GetBinError(b), d_stat = h_dijet_nom->GetBinError(b);
    if (Valid(e, d))
    {
      const double r = e / d;
      const double rel2 = (e_stat * e_stat) / (e * e) + (d_stat * d_stat) / (d * d);
      h_ratio->SetBinContent(b, r);
      h_ratio->SetBinError(b, r * std::sqrt(rel2));
    }
    else
    {
      h_ratio->SetBinContent(b, 0.0);
      h_ratio->SetBinError(b, 0.0);
    }
  }

  // shape-only band: remove per-source global mean shift in ratio.
  std::vector<double> x(nbins), ex(nbins), y(nbins), eyl2(nbins, 0.0), eyh2(nbins, 0.0);
  for (int b = 0; b < nbins; ++b)
  {
    x[b] = h_ratio->GetBinCenter(b + 1);
    ex[b] = 0.5 * h_ratio->GetBinWidth(b + 1);
    y[b] = h_ratio->GetBinContent(b + 1);
  }

  for (size_t i = 1; i < syst.size(); ++i)
  {
    std::vector<double> fshift(nbins, 0.0);
    double avg = 0.0;
    int n = 0;

    for (int b = 1; b <= nbins; ++b)
    {
      const double e_nom = h_efrac_nom->GetBinContent(b), d_nom = h_dijet_nom->GetBinContent(b);
      const double e_var = h_efrac[i]->GetBinContent(b), d_var = h_dijet[i]->GetBinContent(b);
      if (!Valid(e_nom, d_nom) || !Valid(e_var, d_var))
        continue;
      const double r_nom = e_nom / d_nom, r_var = e_var / d_var;
      if (r_nom <= 0.0)
        continue;
      fshift[b - 1] = (r_var / r_nom) - 1.0;
      avg += fshift[b - 1];
      ++n;
    }
    if (n == 0)
      continue;
    avg /= static_cast<double>(n);

    for (int b = 1; b <= nbins; ++b)
    {
      const double r_nom = h_ratio->GetBinContent(b);
      if (r_nom <= 0.0)
        continue;
      const double dshape = (fshift[b - 1] - avg) * r_nom;
      if (dshape >= 0.0)
        eyh2[b - 1] += dshape * dshape;
      else
        eyl2[b - 1] += dshape * dshape;
    }
  }

  std::vector<double> eyl(nbins), eyh(nbins);
  for (int b = 0; b < nbins; ++b)
  {
    eyl[b] = std::sqrt(eyl2[b]);
    eyh[b] = std::sqrt(eyh2[b]);
  }

  TGraphAsymmErrors *g_shape =
      new TGraphAsymmErrors(nbins, x.data(), y.data(), ex.data(), ex.data(), eyl.data(), eyh.data());
  g_shape->SetFillColorAlpha(kAzure - 9, 0.60);

  g_shape->SetFillStyle(1001);
  g_shape->SetLineWidth(0);

  TCanvas *c = new TCanvas("c_efrac_dijet_compare", "", 1000, 500);

  TPad *p_top = new TPad("p_top", "", 0, 0, 0.5, 1.0);
  //p_top->SetBottomMargin(0.02);
  p_top->Draw();
  p_top->cd();

  h_dijet_nom->SetStats(0);
  h_dijet_nom->SetLineColor(kBlack);
  h_dijet_nom->SetLineWidth(0);
  h_dijet_nom->SetMarkerColor(kBlack);
  h_dijet_nom->SetMarkerStyle(20);
  h_dijet_nom->GetXaxis()->SetTitle("p_{T,lead} [GeV]");
  h_dijet_nom->GetXaxis()->SetTitleSize(25);
  h_dijet_nom->GetXaxis()->SetTitleFont(43);
  h_dijet_nom->GetXaxis()->SetLabelFont(43);
  h_dijet_nom->GetXaxis()->SetLabelSize(25);
  h_dijet_nom->GetXaxis()->SetRangeUser(21, 63);
  h_dijet_nom->GetYaxis()->SetRangeUser(0.2, 0.7);
  h_dijet_nom->GetYaxis()->SetTitleSize(25);
  h_dijet_nom->GetYaxis()->SetTitleFont(43);
  h_dijet_nom->GetYaxis()->SetTitleOffset(1.5);
  h_dijet_nom->GetYaxis()->SetLabelFont(43);
  h_dijet_nom->GetYaxis()->SetLabelSize(25);
  h_dijet_nom->GetYaxis()->SetTitle("<#SigmaE_{T}/#delta#eta#delta#phi> [GeV]");

  h_efrac_nom->SetLineColor(kAzure + 2);
  h_efrac_nom->SetLineWidth(0);
  h_efrac_nom->SetMarkerColor(kAzure + 2);
  h_efrac_nom->SetMarkerStyle(20);

  auto *g_dijet_syst = BuildTotalSystBand(h_dijet_nom, h_dijet, kGray + 2);
  auto *g_efrac_syst = BuildTotalSystBand(h_efrac_nom, h_efrac, kAzure - 9);
  g_dijet_syst->SetMarkerColor(kBlack);
  g_efrac_syst->SetMarkerColor(kAzure + 2);
  g_dijet_syst->SetMarkerStyle(20);
  g_efrac_syst->SetMarkerStyle(20);
  g_efrac_syst->SetLineColor(kAzure + 2);
  g_dijet_syst->SetLineColor(kBlack);
  g_dijet_syst->SetLineWidth(2);
  g_efrac_syst->SetLineWidth(2);

  h_dijet_nom->Draw("E1");
  g_dijet_syst->Draw("EP same");
  h_dijet_nom->Draw("E1 same");
  g_efrac_syst->Draw("EP same");
  h_efrac_nom->Draw("E1 same");

  TLegend *leg_top = new TLegend(0.17, 0.72, 0.85, 0.92);
  leg_top->AddEntry("", "#bf{#it{sPHENIX}} Internal", "");
  leg_top->AddEntry("", "200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7", "");
  //leg_top->AddEntry(h_dijet_nom, "Dijet (nominal, stat. unc.)", "lp");
  //leg_top->AddEntry(h_efrac_nom, "Inclusive (nominal, stat. unc.)", "lp");
  leg_top->AddEntry(g_dijet_syst, "Exclusive Dijet", "ple");
  leg_top->AddEntry(g_efrac_syst, "Inclusive Jet", "ple");
  leg_top->SetTextSize(0.04);
  leg_top->Draw();

  c->cd();
  TPad *p_bot = new TPad("p_bot", "", 0.5, 0.0, 1, 1.0);
  //p_bot->SetTopMargin(0.02);
  //p_bot->SetBottomMargin(0.2);
  p_bot->Draw();
  p_bot->cd();

  h_ratio->SetStats(0);
  h_ratio->SetLineColor(kBlack);
  h_ratio->SetMarkerColor(kBlack);
  h_ratio->SetMarkerStyle(20);
  h_ratio->GetYaxis()->SetTitle("Inclusive / Dijet");
  h_ratio->GetYaxis()->SetNdivisions(510);
  h_ratio->GetYaxis()->SetRangeUser(0.8, 1.2);
  h_ratio->GetYaxis()->SetTitleSize(25);
  h_ratio->GetYaxis()->SetTitleFont(43);
  h_ratio->GetYaxis()->SetTitleOffset(1.5);
  h_ratio->GetYaxis()->SetLabelFont(43);
  h_ratio->GetYaxis()->SetLabelSize(25);
  h_ratio->GetXaxis()->SetTitle("p_{T,lead} [GeV]");
  h_ratio->GetXaxis()->SetTitleSize(25);
  h_ratio->GetXaxis()->SetTitleFont(43);
  h_ratio->GetXaxis()->SetLabelFont(43);
  h_ratio->GetXaxis()->SetLabelSize(25);
  h_ratio->GetXaxis()->SetRangeUser(21, 63);

  h_ratio->Draw("E1");
  g_shape->Draw("2 same");
  h_ratio->Draw("E1 same");

  TLine *l0 = new TLine(21, 1.0, 63, 1.0);
  l0->SetLineStyle(2);
  l0->Draw("same");

  TLegend *leg_bot = new TLegend(0.17, 0.78, 0.90, 0.92);
  leg_bot->AddEntry(h_ratio, "Inclusive/Dijet", "lp");
  leg_bot->AddEntry(g_shape, "Syst. Unc.", "f");
  leg_bot->SetTextSize(0.04);
  leg_bot->Draw();

  c->Update();
  c->SaveAs(output_name);
  std::cout << "Saved: " << output_name << std::endl;
}