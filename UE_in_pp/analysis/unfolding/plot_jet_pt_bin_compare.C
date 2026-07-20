#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
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
        "calib_dijet_reweight_trim_10_2_etEffCorrected",                                          // [0] nominal
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
        "calib_dijet_reweight_trim_10_2_etEffCorrected"};     // [last] Herwig, read from herwig file
  }

  // --------------------------------------------------------------------------
  // BuildProfileHistogram: unchanged from original
  // --------------------------------------------------------------------------
  TH1D *BuildProfileHistogram(TFile *f, const std::string &syst_name, const std::string &out_name)
  {
    TH2D *h2_uniform = dynamic_cast<TH2D *>(f->Get(("h_unfold_" + syst_name).c_str()));
    if (!h2_uniform)
    {
      std::cerr << "Missing histogram h_unfold_" << syst_name << std::endl;
      return nullptr;
    }

    TH2D *h2_var = new TH2D((out_name + "_var").c_str(), "",
                             truthnpt, truthptbins, truthnet, truthetbins);
    for (int ix = 1; ix <= h2_uniform->GetNbinsX(); ++ix)
      for (int iy = 1; iy <= h2_uniform->GetNbinsY(); ++iy)
      {
        h2_var->SetBinContent(ix, iy, h2_uniform->GetBinContent(ix, iy));
        h2_var->SetBinError  (ix, iy, h2_uniform->GetBinError  (ix, iy));
      }

    TProfile *prof = h2_var->ProfileX((out_name + "_prof").c_str(),
                                       1, h2_var->GetNbinsY() - 1);
    const int nbins = prof->GetNbinsX();
    std::vector<double> edges(nbins + 1, 0.0);
    for (int b = 0; b < nbins; ++b)
      edges[b] = prof->GetXaxis()->GetBinLowEdge(b + 1);
    edges[nbins] = prof->GetXaxis()->GetBinUpEdge(nbins);

    TH1D *out = new TH1D(out_name.c_str(), "", nbins, edges.data());
    for (int b = 1; b <= nbins; ++b)
    {
      out->SetBinContent(b, prof->GetBinContent(b));
      out->SetBinError  (b, prof->GetBinError  (b));
    }
    return out;
  }

  // --------------------------------------------------------------------------
  // FitLinear
  //
  // Fit h with f(pT) = p0 + p1*pT over [ptlo, pthi].
  // Returns the fitted TF1 (owned by caller).  p0 and p1 are set via reference.
  // --------------------------------------------------------------------------
  TF1 *FitLinear(TH1D *h, double ptlo, double pthi,
                 double &p0, double &p1,
                 const std::string &fname)
  {
    TF1 *f = new TF1(fname.c_str(), "[0] + [1]*x", ptlo, pthi);
    h->Fit(f, "RQN0");   // Range, Quiet, No-draw, no-store
    p0 = f->GetParameter(0);
    p1 = f->GetParameter(1);
    return f;
  }

  // --------------------------------------------------------------------------
  // BuildSlopeBand
  //
  // For a nominal histogram h_nom and its variation histograms h_vars
  // (index 0 = nominal, 1..N-1 = systematics), build a TGraphAsymmErrors
  // that represents the uncertainty on the LINEAR SLOPE only.
  //
  // Method per source i:
  //   1. Fit nominal  -> slope b_nom
  //   2. Fit variation i -> slope b_i
  //   3. The slope-only variation at bin centre pT_j is:
  //        delta_i(pT_j) = (b_i - b_nom) * (pT_j - ptlo)
  //      where ptlo is the fit range low edge, so the band pivots at ptlo
  //      (the intercept is held fixed = nominal intercept, i.e. scale removed).
  //   4. Accumulate delta_i^2 into up/dn quadrature sums asymmetrically.
  //
  // The returned graph is centred on the NOMINAL FIT LINE (not the data),
  // so it shows purely the slope uncertainty envelope around the fit.
  // --------------------------------------------------------------------------
  TGraphAsymmErrors *BuildSlopeBand(TH1D             *h_nom,
                                    const std::vector<TH1D *> &h_vars,
                                    double            ptlo,
                                    double            pthi,
                                    int               color,
                                    float             alpha,
                                    const std::string &tag)
  {
    // --- fit nominal ---
    double a_nom = 0., b_nom = 0.;
    TF1 *f_nom = FitLinear(h_nom, ptlo, pthi, a_nom, b_nom, tag + "_nom_fit");

    const int n = h_nom->GetNbinsX();
    std::vector<double> x(n), ex(n), y_fit(n);
    std::vector<double> up2(n, 0.), dn2(n, 0.);

    // centre the graph on the nominal fit line
    for (int b = 1; b <= n; ++b)
    {
      x[b-1]     = h_nom->GetBinCenter(b);
      ex[b-1]    = 0.5 * h_nom->GetBinWidth(b);
      y_fit[b-1] = f_nom->Eval(x[b-1]);
    }

    // --- loop over systematic variations ---
    for (size_t i = 1; i < h_vars.size(); ++i)
    {
      double a_var = 0., b_var = 0.;
      TF1 *f_var = FitLinear(h_vars[i], ptlo, pthi, a_var, b_var,
                              Form("%s_var%zu_fit", tag.c_str(), i));

      const double db = b_var - b_nom;   // slope difference

      for (int b = 0; b < n; ++b)
      {
        // slope-only shift: pivot at ptlo so the intercept (scale) is fixed
        const double delta = db * (x[b] - 33);
        if (delta >= 0.) up2[b] += delta * delta;
        else             dn2[b] += delta * delta;
      }
      delete f_var;
    }

    std::vector<double> eyl(n), eyh(n);
    for (int b = 0; b < n; ++b)
    {
      eyl[b] = std::sqrt(dn2[b]);
      eyh[b] = std::sqrt(up2[b]);
    }

    auto *g = new TGraphAsymmErrors(n,
                                    x.data(), y_fit.data(),
                                    ex.data(), ex.data(),
                                    eyl.data(), eyh.data());
    g->SetFillColorAlpha(color, alpha);
    g->SetFillStyle(1001);
    g->SetLineWidth(0);

    delete f_nom;
    return g;
  }

} // namespace

// ============================================================================
void plot_jet_pt_bin_compare(
    const char *dijet_input        = "sphenix_primary_run28_output_files/output_mbd_correct_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root",
    const char *efrac_input        = "sphenix_primary_run28_output_files/output_mbd_correct_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root",
    const char *dijet_herwig_input = "sphenix_primary_run28_output_files/output_mbd_correct_herwig_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root",
    const char *efrac_herwig_input = "sphenix_primary_run28_output_files/output_mbd_correct_herwig_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root",
    const char *output_name        = "sphenix_primary_run28_output_files/h_jet_pt_bin_unc_figure.pdf",
    double      ptlo               = 21.,
    double      pthi               = 63.)
{
  gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
  gROOT->ProcessLine("SetsPhenixStyle()");

  TFile *f_dijet        = TFile::Open(dijet_input);
  TFile *f_efrac        = TFile::Open(efrac_input);
  TFile *f_dijet_herwig = TFile::Open(dijet_herwig_input);
  TFile *f_efrac_herwig = TFile::Open(efrac_herwig_input);
  if (!f_dijet        || f_dijet->IsZombie()        ||
      !f_efrac        || f_efrac->IsZombie()        ||
      !f_dijet_herwig || f_dijet_herwig->IsZombie() ||
      !f_efrac_herwig || f_efrac_herwig->IsZombie())
  {
    std::cerr << "Could not open input files" << std::endl;
    return;
  }

  // --------------------------------------------------------------------------
  // Load profile histograms for all systematic variations
  // --------------------------------------------------------------------------
  const std::vector<std::string> syst = GetSystList();
  std::vector<TH1D *> h_dijet, h_efrac;
  h_dijet.reserve(syst.size());
  h_efrac.reserve(syst.size());

  // indices 0 .. N-2: read from dijet/efrac files
  for (size_t i = 0; i < syst.size() - 1; ++i)
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
  // last index: Herwig MC variation, read from herwig files
  TH1D *hd_hw = BuildProfileHistogram(f_dijet_herwig, syst.back(),
                                       Form("h_dijet_herwig_%zu", syst.size()-1));
  TH1D *he_hw = BuildProfileHistogram(f_efrac_herwig, syst.back(),
                                       Form("h_efrac_herwig_%zu", syst.size()-1));
  if (!hd_hw || !he_hw)
  {
    std::cerr << "Failed at Herwig syst index" << std::endl;
    return;
  }
  h_dijet.push_back(hd_hw);
  h_efrac.push_back(he_hw);

  // --------------------------------------------------------------------------
  // Common normalisation scale
  // --------------------------------------------------------------------------
  const double norm = 3.0 / (2.2 * 2. * M_PI);
  for (size_t i = 0; i < h_dijet.size(); ++i)
  {
    h_dijet[i]->Scale(norm);
    h_efrac[i]->Scale(norm);
  }

  TH1D *h_dijet_nom = h_dijet[0];
  TH1D *h_efrac_nom = h_efrac[0];

  // --------------------------------------------------------------------------
  // Fit nominals and get the nominal fit functions for drawing
  // --------------------------------------------------------------------------
  double a_dnom = 0., b_dnom = 0., a_enom = 0., b_enom = 0.;
  TF1 *f_dijet_nom = FitLinear(h_dijet_nom, ptlo, pthi, a_dnom, b_dnom, "f_dijet_nom");
  TF1 *f_efrac_nom = FitLinear(h_efrac_nom, ptlo, pthi, a_enom, b_enom, "f_efrac_nom");

  f_dijet_nom->SetLineColor(kBlack);
  f_dijet_nom->SetLineWidth(2);
  f_dijet_nom->SetLineStyle(1);

  f_efrac_nom->SetLineColor(kAzure + 2);
  f_efrac_nom->SetLineWidth(2);
  f_efrac_nom->SetLineStyle(1);

  // --------------------------------------------------------------------------
  // Build slope-only uncertainty bands (centred on nominal fit line)
  // --------------------------------------------------------------------------
  TGraphAsymmErrors *g_dijet_slope = BuildSlopeBand(h_dijet_nom, h_dijet,
                                                     ptlo, pthi,
                                                     kGray + 2, 0.50,
                                                     "dijet");

  TGraphAsymmErrors *g_efrac_slope = BuildSlopeBand(h_efrac_nom, h_efrac,
                                                     ptlo, pthi,
                                                     kAzure - 9, 0.5,
                                                     "efrac");

  // --------------------------------------------------------------------------
  // Canvas and single pad
  // --------------------------------------------------------------------------
  TCanvas *c = new TCanvas("c_slope_unc", "", 700, 700);
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.12);

  // --- axis histogram (use dijet nominal as template) ---
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
  h_dijet_nom->GetXaxis()->SetRangeUser(ptlo, pthi);
  h_dijet_nom->GetYaxis()->SetRangeUser(0.2, 0.7);
  h_dijet_nom->GetYaxis()->SetTitle("<#SigmaE_{T}/#delta#eta#delta#phi> [GeV]");
  h_dijet_nom->GetYaxis()->SetTitleSize(25);
  h_dijet_nom->GetYaxis()->SetTitleFont(43);
  h_dijet_nom->GetYaxis()->SetTitleOffset(1.5);
  h_dijet_nom->GetYaxis()->SetLabelFont(43);
  h_dijet_nom->GetYaxis()->SetLabelSize(25);

  h_efrac_nom->SetStats(0);
  h_efrac_nom->SetLineColor(kAzure + 2);
  h_efrac_nom->SetLineWidth(0);
  h_efrac_nom->SetMarkerColor(kAzure + 2);
  h_efrac_nom->SetMarkerStyle(20);

  g_dijet_slope->SetMarkerColor(kBlack);
  g_efrac_slope->SetMarkerColor(kAzure + 2);

  // --------------------------------------------------------------------------
  // Draw: slope bands first (behind), then fit lines, then data points on top
  // --------------------------------------------------------------------------
  h_dijet_nom->Draw("E1");          // sets axes, draws dijet stat errors

  g_dijet_slope->Draw("E2 same");    // dijet slope band
  g_efrac_slope->Draw("E2 same");    // efrac slope band

  f_dijet_nom->Draw("same");        // nominal linear fit — dijet
  f_efrac_nom->Draw("same");        // nominal linear fit — efrac

  h_dijet_nom->Draw("E1 same");     // dijet data points over band
  h_efrac_nom->Draw("E1 same");     // efrac data points over band

   // Also retrieve fit uncertainties from the TF1 objects
  const double a_dnom_err = f_dijet_nom->GetParError(0);
  const double b_dnom_err = f_dijet_nom->GetParError(1);
  const double a_enom_err = f_efrac_nom->GetParError(0);
  const double b_enom_err = f_efrac_nom->GetParError(1);

  // --------------------------------------------------------------------------
  // Legend
  // --------------------------------------------------------------------------
  TLegend *leg = new TLegend(0.15, 0.6, 0.83, 0.92);
  leg->SetTextFont(43);
  leg->SetTextSize(20);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry("",            "#bf{#it{sPHENIX}} Internal",                     "");
  leg->AddEntry("",            "200 GeV p+p  anti-k_{t} R=0.4  |#eta_{jet}|<0.7","");
  leg->AddEntry("", "Bin-Correlated Scale Uncertainty Removed","");
  leg->AddEntry(g_dijet_slope, "Excluded Dijet",                       "pf");
  leg->AddEntry(g_efrac_slope, "Inclusive Jet",                         "pf");
  leg->AddEntry(f_dijet_nom,Form("f_{dijet}(p_{T}) = %.3f #pm %.3f + (%.4f #pm %.4f) p_{T}",a_dnom, a_dnom_err, b_dnom, b_dnom_err),"l");
  leg->AddEntry(f_efrac_nom,Form("f_{jet}(p_{T}) = %.3f #pm %.3f + (%.4f #pm %.4f) p_{T}",a_enom, a_enom_err, b_enom, b_enom_err),"l");
  leg->Draw();

  c->Update();
  c->SaveAs(output_name);
  std::cout << "Saved: " << output_name << std::endl;
}
