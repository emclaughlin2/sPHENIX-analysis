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
#include <iomanip>
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
  // BuildProfileHistogramFromKey
  //
  // Reads a 2D histogram named exactly `hist_key` from file `f`, copies it onto
  // the truth binning, profiles over Y (bins 1 .. NbinsY-1), and returns the
  // resulting TH1D.  This is the generic core; BuildProfileHistogram wraps it
  // for the unfold histograms, and the truth histograms call it directly.
  // --------------------------------------------------------------------------
  TH1D *BuildProfileHistogramFromKey(TFile *f, const std::string &hist_key, const std::string &out_name)
  {
    TH2D *h2_uniform = dynamic_cast<TH2D *>(f->Get(hist_key.c_str()));
    if (!h2_uniform)
    {
      std::cerr << "Missing histogram " << hist_key << std::endl;
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
  // BuildProfileHistogram: unchanged behaviour — prefixes "h_unfold_".
  // --------------------------------------------------------------------------
  TH1D *BuildProfileHistogram(TFile *f, const std::string &syst_name, const std::string &out_name)
  {
    return BuildProfileHistogramFromKey(f, "h_unfold_" + syst_name, out_name);
  }

  // --------------------------------------------------------------------------
  // FitLinear
  //
  // Fit h with f(pT) = p0 + p1*pT over [ptlo, pthi].
  // Returns the fitted TF1 (owned by caller).  p0 and p1 are set via reference.
  // Still used by the driver for the nominal and truth fit lines that are
  // drawn on the figure.
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
  // SlopeSystError
  //
  // Systematic uncertainty on the linear-fit SLOPE, built from the variations
  // (not from the nominal fit's parameter error). For each variation i (1..N-1)
  // the line is refit, its slope b_i extracted, and the deviation from the
  // nominal slope (b_i - b_nom) accumulated in quadrature:
  //
  //     sigma_b^syst = sqrt( sum_i (b_i - b_nom)^2 ).
  //
  // The loop spans the same variation set as BuildScaleRemovedShapeBand, so the
  // quoted slope error and the band describe the same systematics. A fully-
  // correlated scale shifts the intercept, not the slope, so it contributes
  // negligibly here — consistent with the "scale removed" figure.
  //
  // NOTE: each variation is treated independently (up/down pairs both counted),
  // matching the band's convention. Group pairs first if your note requires it.
  // --------------------------------------------------------------------------
  double SlopeSystError(TH1D *h_nom, const std::vector<TH1D *> &h_vars,
                        double ptlo, double pthi, double b_nom,
                        const std::string &tag)
  {
    (void)h_nom;
    double sum2 = 0.;
    for (size_t i = 1; i < h_vars.size(); ++i)
    {
      double a_i = 0., b_i = 0.;
      TF1 *f = FitLinear(h_vars[i], ptlo, pthi, a_i, b_i,
                         Form("%s_sysfit_%zu", tag.c_str(), i));
      const double db = b_i - b_nom;
      sum2 += db * db;
      delete f;
    }
    return std::sqrt(sum2);
  }

  // --------------------------------------------------------------------------
  // BuildScaleRemovedShapeBand
  //
  // For a nominal histogram h_nom and its variation histograms h_vars
  // (index 0 = nominal, 1..N-1 = systematics), build a TGraphAsymmErrors
  // that represents the SHAPE uncertainty with the fully-correlated
  // (multiplicative scale) component removed.
  //
  // Method per source i:
  //   1. Form the bin-by-bin ratio  r_b = h_var_b / h_nom_b.
  //   2. Compute the inverse-variance-weighted mean ratio k = <r_b>_w, using
  //      weights w_b = 1/sigma_nom_b^2. This k is the single best-fit common
  //      scale factor, i.e. the fully-correlated ("scale") part of the
  //      variation.
  //   3. Divide k out and convert back to absolute units:
  //        delta_i(b) = (r_b - k) * h_nom_b
  //      This keeps the genuine bin-by-bin shape difference (including any
  //      curvature) while removing only the correlated scale mode. No pivot
  //      and no linearity assumption are used.
  //   4. Accumulate delta_i^2 into up/dn quadrature sums asymmetrically.
  //
  // The returned graph is centred on the NOMINAL BIN CONTENTS (not a fit line),
  // which is the consistent choice now that the band carries bin-by-bin shape.
  //
  // NOTE: this removes the bin-correlated component from EVERY source fed in,
  // which is the intended generalisation of the old slope-only band (the
  // figure shows "Bin-Correlated Scale Uncertainty Removed"). Feed the same
  // variation list as before.
  // --------------------------------------------------------------------------
  TGraphAsymmErrors *BuildScaleRemovedShapeBand(TH1D                      *h_nom,
                                                const std::vector<TH1D *> &h_vars,
                                                int                        color,
                                                float                      alpha,
                                                const std::string         &tag,
                                                const std::vector<std::string> &names = {})
  {
    const int n = h_nom->GetNbinsX();

    std::vector<double> x(n), ex(n), y_nom(n);
    std::vector<double> up2(n, 0.), dn2(n, 0.);
    std::vector<double> w(n);

    // --- per-bin weights from nominal stat errors, band centred on nominal ---
    for (int b = 0; b < n; ++b)
    {
      const double e = h_nom->GetBinError(b + 1);
      w[b]     = (e > 0.) ? 1.0 / (e * e) : 0.;
      x[b]     = h_nom->GetBinCenter(b + 1);
      ex[b]    = 0.5 * h_nom->GetBinWidth(b + 1);
      y_nom[b] = h_nom->GetBinContent(b + 1);
    }

    // --- header for the per-systematic, per-bin diagnostic ------------------
    std::cout << "\n=== [" << tag << "] scale-removed systematic shifts per pT "
              << "bin (GeV; + = up, - = down) ===\n";
    std::cout << std::setw(34) << std::left << "pT centre [GeV]:" << std::right;
    for (int b = 0; b < n; ++b) std::cout << std::setw(11) << x[b];
    std::cout << "\n";

    // --- loop over systematic variations ---
    for (size_t i = 1; i < h_vars.size(); ++i)
    {
      // weighted-mean ratio = the fully-correlated (scale) component
      double sw = 0., swr = 0.;
      for (int b = 0; b < n; ++b)
      {
        const double nom = y_nom[b];
        if (nom == 0.) continue;
        const double r = h_vars[i]->GetBinContent(b + 1) / nom;
        sw  += w[b];
        swr += w[b] * r;
      }
      const double k = (sw > 0.) ? swr / sw : 1.0;   // common scale factor, divided out

      // label: short name if provided, else index
      std::string label = "src " + std::to_string(i);
      if (i < names.size()) label += " " + names[i].substr(0, 24);

      std::cout << std::setw(34) << std::left << label << std::right;
      for (int b = 0; b < n; ++b)
      {
        const double nom = y_nom[b];
        const double r   = (nom != 0.) ? h_vars[i]->GetBinContent(b + 1) / nom : 1.0;
        const double delta = (r - k) * nom;          // scale-removed, kept bin-by-bin
        if (delta >= 0.) up2[b] += delta * delta;
        else             dn2[b] += delta * delta;
        std::cout << std::setw(11) << std::fixed << std::setprecision(5) << delta;
      }
      std::cout << "\n";
    }

    // --- asymmetric errors from quadrature sums ---
    std::vector<double> eyl(n), eyh(n);
    for (int b = 0; b < n; ++b)
    {
      eyl[b] = std::sqrt(dn2[b]);
      eyh[b] = std::sqrt(up2[b]);
    }

    // --- total band per bin ---
    std::cout << std::setw(34) << std::left << "TOTAL syst up (+):" << std::right;
    for (int b = 0; b < n; ++b)
      std::cout << std::setw(11) << std::fixed << std::setprecision(5) << eyh[b];
    std::cout << "\n";
    std::cout << std::setw(34) << std::left << "TOTAL syst down (-):" << std::right;
    for (int b = 0; b < n; ++b)
      std::cout << std::setw(11) << std::fixed << std::setprecision(5) << (-eyl[b]);
    std::cout << "\n";
    std::cout << std::setw(34) << std::left << "TOTAL syst up [%]:" << std::right;
    for (int b = 0; b < n; ++b)
      std::cout << std::setw(11) << std::fixed << std::setprecision(2)
                << (y_nom[b] != 0. ? 100. * eyh[b] / y_nom[b] : 0.);
    std::cout << "\n";
    std::cout << std::setw(34) << std::left << "TOTAL syst down [%]:" << std::right;
    for (int b = 0; b < n; ++b)
      std::cout << std::setw(11) << std::fixed << std::setprecision(2)
                << (y_nom[b] != 0. ? -100. * eyl[b] / y_nom[b] : 0.);
    std::cout << "\n" << std::defaultfloat;

    auto *g = new TGraphAsymmErrors(n,
                                    x.data(), y_nom.data(),
                                    ex.data(), ex.data(),
                                    eyl.data(), eyh.data());
    g->SetFillColorAlpha(color, alpha);
    g->SetFillStyle(1001);
    g->SetLineWidth(0);
    (void)tag;   // retained in signature for naming consistency with other bands
    return g;
  }

} // namespace

TGraphErrors *HistToGraph(TH1D *h)
{
  const int n = h->GetNbinsX();
  std::vector<double> x(n), y(n), ex(n, 0.), ey(n);
  for (int b = 0; b < n; ++b)
  {
    x[b]  = h->GetBinCenter(b + 1);
    y[b]  = h->GetBinContent(b + 1);
    ey[b] = h->GetBinError(b + 1);
  }
  auto *g = new TGraphErrors(n, x.data(), y.data(), ex.data(), ey.data());
  g->SetMarkerStyle(h->GetMarkerStyle());
  g->SetMarkerColor(h->GetMarkerColor());
  g->SetMarkerSize (h->GetMarkerSize());
  g->SetLineColor  (h->GetLineColor());
  g->SetLineWidth  (h->GetLineWidth());
  return g;
}


// ============================================================================
void figure6_plots_update(
    const char *dijet_input        = "output_mbd_correct_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root",
    const char *efrac_input        = "output_mbd_correct_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root",
    const char *dijet_herwig_input = "output_mbd_correct_herwig_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root",
    const char *efrac_herwig_input = "output_mbd_correct_herwig_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root",
    const char *output_name        = "h_jet_pt_bin_unc_figure_update.pdf",
    double      ptlo               = 21.,
    double      pthi               = 63.5)
{
  gROOT->LoadMacro("sPhenixStyle.C");
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
  // Load generator-level TRUTH profiles (no variations, central slope only).
  // Same profile->fit path as the unfold, but the input histogram is
  // "h_truth_calib_dijet" instead of "h_unfold_<syst>".
  //   PYTHIA truth -> dijet (Pythia) file
  //   HERWIG truth -> dijet Herwig file
  // --------------------------------------------------------------------------
  TH1D *h_truth_pythia = BuildProfileHistogramFromKey(f_dijet,        "h_truth_calib_dijet", "h_truth_pythia");
  TH1D *h_truth_herwig = BuildProfileHistogramFromKey(f_dijet_herwig, "h_truth_calib_dijet", "h_truth_herwig");
  if (!h_truth_pythia || !h_truth_herwig)
  {
    std::cerr << "Failed to load truth histograms" << std::endl;
    return;
  }

  // --------------------------------------------------------------------------
  // Common normalisation scale (apply to truth too, so it's comparable)
  // --------------------------------------------------------------------------
  const double norm = 3.0 / (2.2 * 2. * M_PI);
  for (size_t i = 0; i < h_dijet.size(); ++i)
  {
    h_dijet[i]->Scale(norm);
    h_efrac[i]->Scale(norm);
  }
  h_truth_pythia->Scale(norm);
  h_truth_herwig->Scale(norm);

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
  f_efrac_nom->SetMarkerStyle(33);
  f_efrac_nom->SetMarkerSize(1.4);

  // --------------------------------------------------------------------------
  // Fit TRUTH (central slope only — no systematic band)
  // --------------------------------------------------------------------------
  double a_tpy = 0., b_tpy = 0., a_thw = 0., b_thw = 0.;
  TF1 *f_truth_pythia = FitLinear(h_truth_pythia, ptlo, pthi, a_tpy, b_tpy, "f_truth_pythia");
  TF1 *f_truth_herwig = FitLinear(h_truth_herwig, ptlo, pthi, a_thw, b_thw, "f_truth_herwig");

  f_truth_pythia->SetLineColor(kRed + 1);
  f_truth_pythia->SetLineWidth(2);
  f_truth_pythia->SetLineStyle(7);   // dashed

  f_truth_herwig->SetLineColor(kSpring - 6);
  f_truth_herwig->SetLineWidth(2);
  f_truth_herwig->SetLineStyle(9);   // dashed

  h_truth_pythia->SetStats(0);
  //h_truth_pythia->SetLineColor(kRed + 1);
  h_truth_pythia->SetLineWidth(0);
  h_truth_pythia->SetMarkerColor(kRed + 1);
  h_truth_pythia->SetMarkerStyle(21);   // open circle
  //h_truth_pythia->SetMarkerSize(0.8);

  h_truth_herwig->SetStats(0);
  //h_truth_herwig->SetLineColor(kSpring - 6);
  h_truth_herwig->SetLineWidth(0);
  h_truth_herwig->SetMarkerColor(kSpring - 6);
  h_truth_herwig->SetMarkerStyle(22);   // open square
  //h_truth_herwig->SetMarkerSize(0.8);

  // --------------------------------------------------------------------------
  // Build scale-removed shape uncertainty bands (centred on nominal points).
  // The bin-correlated (scale) component is divided out of each variation;
  // the remaining bin-by-bin shape difference is summed in quadrature.
  // --------------------------------------------------------------------------
  TGraphAsymmErrors *g_dijet_band = BuildScaleRemovedShapeBand(h_dijet_nom, h_dijet,
                                                               kGray + 2, 0.60,
                                                               "dijet", syst);

  TGraphAsymmErrors *g_efrac_band = BuildScaleRemovedShapeBand(h_efrac_nom, h_efrac,
                                                               kAzure - 9, 0.6,
                                                               "efrac", syst);

  // --------------------------------------------------------------------------
  // Canvas and single pad
  // --------------------------------------------------------------------------
  TCanvas *c = new TCanvas("c_slope_unc", "", 700, 700);
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.12);

  // --- axis histogram (use dijet nominal as template) ---
  h_dijet_nom->SetStats(0);
  h_dijet_nom->SetLineColor(kBlack);
  h_dijet_nom->SetLineWidth(2);
  h_dijet_nom->SetMarkerColor(kBlack);
  h_dijet_nom->SetMarkerStyle(20);
  h_dijet_nom->GetXaxis()->SetTitle("p_{T,lead} [GeV]");
  h_dijet_nom->GetXaxis()->SetTitleSize(25);
  h_dijet_nom->GetXaxis()->SetTitleFont(43);
  h_dijet_nom->GetXaxis()->SetLabelFont(43);
  h_dijet_nom->GetXaxis()->SetLabelSize(25);
  h_dijet_nom->GetXaxis()->SetRangeUser(ptlo, pthi);
  h_dijet_nom->GetYaxis()->SetRangeUser(0.3, 0.7);
  h_dijet_nom->GetYaxis()->SetTitle("<#SigmaE_{T}/#delta#eta#delta#phi> [GeV]");
  h_dijet_nom->GetYaxis()->SetTitleSize(25);
  h_dijet_nom->GetYaxis()->SetTitleFont(43);
  h_dijet_nom->GetYaxis()->SetTitleOffset(1.5);
  h_dijet_nom->GetYaxis()->SetLabelFont(43);
  h_dijet_nom->GetYaxis()->SetLabelSize(25);

  h_efrac_nom->SetStats(0);
  h_efrac_nom->SetLineColor(kAzure + 2);
  h_efrac_nom->SetLineWidth(2);
  h_efrac_nom->SetMarkerColor(kAzure + 2);
  h_efrac_nom->SetMarkerStyle(33);
  h_efrac_nom->SetMarkerSize(1.4);

  g_dijet_band->SetMarkerColor(kBlack);
  g_dijet_band->SetMarkerStyle(20);
  g_efrac_band->SetMarkerColor(kAzure + 2);
  g_efrac_band->SetMarkerStyle(33);
  g_efrac_band->SetMarkerSize(1.4);

  TGraphErrors *g_dijet_pts = HistToGraph(h_dijet_nom);
  TGraphErrors *g_efrac_pts = HistToGraph(h_efrac_nom);

  h_dijet_nom->SetLineWidth(0);

  // --------------------------------------------------------------------------
  // Draw: shape bands first (behind), then fit lines, then data points on top,
  //       then truth fit lines + truth points
  // --------------------------------------------------------------------------
  h_dijet_nom->Draw("E1");          // sets axes, draws dijet stat errors

  g_dijet_band->Draw("E2 same");    // dijet scale-removed shape band
  g_efrac_band->Draw("E2 same");    // efrac scale-removed shape band

  f_dijet_nom->Draw("same");        // nominal linear fit — dijet
  f_efrac_nom->Draw("same");        // nominal linear fit — efrac

  f_truth_pythia->Draw("same");     // truth linear fit — Pythia
  f_truth_herwig->Draw("same");     // truth linear fit — Herwig
  h_truth_pythia->Draw("P same");  // truth points — Pythia
  h_truth_herwig->Draw("P same");  // truth points — Herwig

  g_dijet_pts->Draw("P same");   // data points carry STATISTICAL errors only
  g_efrac_pts->Draw("P same");

   // Also retrieve fit uncertainties from the TF1 objects
  const double a_dnom_err = f_dijet_nom->GetParError(0);
  const double b_dnom_err = f_dijet_nom->GetParError(1);
  const double a_enom_err = f_efrac_nom->GetParError(0);
  const double b_enom_err = f_efrac_nom->GetParError(1);
  const double b_tpy_err  = f_truth_pythia->GetParError(1);
  const double b_thw_err  = f_truth_herwig->GetParError(1);

  // Systematic uncertainty on the slope from the variations (replaces the
  // statistical parameter error in the legend below).
  const double b_dnom_syst = SlopeSystError(h_dijet_nom, h_dijet, ptlo, pthi, b_dnom, "dijet");
  const double b_enom_syst = SlopeSystError(h_efrac_nom, h_efrac, ptlo, pthi, b_enom, "efrac");
  std::cout << "dijet slope = " << b_dnom << "  stat = " << b_dnom_err
            << "  syst = " << b_dnom_syst << std::endl;
  std::cout << "efrac slope = " << b_enom << "  stat = " << b_enom_err
            << "  syst = " << b_enom_syst << std::endl;

  // --------------------------------------------------------------------------
  // Legend
  // --------------------------------------------------------------------------

  f_truth_pythia->SetMarkerStyle(21);
  f_truth_pythia->SetMarkerColor(kRed + 1);
  f_truth_herwig->SetMarkerStyle(22);
  f_truth_herwig->SetMarkerColor(kSpring - 6);
  TLegend *leg = new TLegend(0.15, 0.6, 0.63, 0.92);
  leg->SetTextFont(43);
  leg->SetTextSize(20);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry("",            "#bf{#it{sPHENIX}} Internal",                     "");
  leg->AddEntry("",            "200 GeV p+p  anti-k_{t} R=0.4  |#eta_{jet}|<0.7","");
  leg->AddEntry("", "Bin-Correlated Scale Uncertainty Removed","");
  leg->AddEntry(g_dijet_band, "Excluded Dijet",                       "pf");
  leg->AddEntry(g_efrac_band, "Inclusive Jet",                         "pf");
  leg->AddEntry(f_dijet_nom,Form("Ex. Dijet: slope = %.4f #kern[-0.15]{#pm %.4f}",b_dnom, b_dnom_syst),"l");
  leg->AddEntry(f_efrac_nom,Form("Inc. Jet: slope = %.4f #kern[-0.15]{#pm %.4f}",b_enom, b_enom_syst),"l");
  leg->AddEntry(f_truth_pythia,Form("PYTHIA truth: slope = %.4f", b_tpy),"lp");
  leg->AddEntry(f_truth_herwig,Form("HERWIG truth: slope = %.4f", b_thw),"lp");
  leg->Draw();

  c->Update();
  c->SaveAs(output_name);
  std::cout << "Saved: " << output_name << std::endl;
}