#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TROOT.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <unfold_Def.h>

namespace
{
  // ==========================================================================
  //  Systematic-source grouping.
  //
  //  Each entry is the list of indices (into GetSystList / h_vars) that form ONE
  //  systematic SOURCE. The Monte Carlo throws ONE Gaussian per source per toy,
  //  applied coherently across all four pT bins — this is what makes correlated
  //  systematics (e.g. a global scale) propagate into the INTERCEPT while leaving
  //  the SLOPE largely unaffected, with no explicit "scale removal" needed.
  //
  //    pair  {up, down} -> per-bin sigma_j = 0.5*|var_up_j - var_down_j|
  //    single {idx}     -> per-bin sigma_j = |var_idx_j - nom_j|   (symmetrised)
  //
  //  The Herwig entry (last index) is the GENERATOR comparison shown as its own
  //  curve; it is NOT included as a throw source here. Flip include_herwig if you
  //  want it treated as an unfolding/model systematic instead.
  //
  //  >>> VERIFY this grouping against your systematics note before trusting the
  //  >>> quoted parameter uncertainties. <<<
  // ==========================================================================
  std::vector<std::vector<int>> GetSystSources()
  {
    return {
        {1, 2},     // JES up/down
        {3, 4},     // JER up/down
        {5},        // reweight closure (unweighted trim)
        {6, 7},     // dijet efficiency up/down
        {8},        // pileup ET correction
        {9},        // timing efficiency
        {10, 11},   // EMCal scale up/down
        {12, 13},   // IHCal scale up/down
        {14, 15},   // OHCal scale up/down
        {16, 17},   // hadronic response up/down
        {18},       // phi resolution
        {19},       // cluster smearing
        {20, 21},   // noise threshold (2sigma/4sigma about nominal 3sigma)
        {22},       // OHCal MC/data variation
        {23, 24},   // ET reweight variation up/down
        {25}
    };
  }

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

  TH1D *BuildProfileHistogram(TFile *f, const std::string &syst_name, const std::string &out_name)
  {
    return BuildProfileHistogramFromKey(f, "h_unfold_" + syst_name, out_name);
  }

  // --------------------------------------------------------------------------
  // ProjectionXbar
  //
  // The X axis of the 2D histogram `key` in file `f` is the jet pT axis, so
  // ProjectionX gives the yield (dN/dpT) spectrum. For each analysis pT bin
  // [edges[j], edges[j+1]) we form the yield-weighted mean pT,
  //
  //     <pT>_j = sum_i ( N_i * pTcenter_i ) / sum_i N_i ,
  //
  // over the fine projection bins i whose centre falls in the analysis bin.
  // This is exactly the abscissa the linear-model folding requires.
  //
  // REQUIRES the projected X axis to be (a) physical pT and (b) finer than the
  // analysis bins. If only one fine bin lands in an analysis bin, <pT> collapses
  // to the bin centre and a warning is printed — that means the stored 2D is
  // already on the coarse binning and carries no sub-bin information.
  // --------------------------------------------------------------------------
  std::vector<double> ProjectionXbar(TFile *f, const std::string &key,
                                     const std::vector<double> &edges,
                                     const std::string &tag)
  {
    const int nb = (int)edges.size() - 1;
    std::vector<double> xbar(nb, 0.);

    TH2D *h2 = dynamic_cast<TH2D *>(f->Get(key.c_str()));
    if (!h2)
    {
      std::cerr << "ProjectionXbar: missing 2D histogram " << key
                << " (using bin centres)" << std::endl;
      for (int j = 0; j < nb; ++j) xbar[j] = 0.5 * (edges[j] + edges[j + 1]);
      return xbar;
    }

    TH1D *spec = h2->ProjectionX((tag + "_projx").c_str());
    for (int j = 0; j < nb; ++j)
    {
      const double lo = edges[j], hi = edges[j + 1];
      double sw = 0., swx = 0.; int nfine = 0;
      for (int i = 1; i <= spec->GetNbinsX(); ++i)
      {
        const double cx = spec->GetBinCenter(i);
        if (cx < lo || cx >= hi) continue;
        const double w = spec->GetBinContent(i);
        if (w > 0.) { sw += w; swx += w * cx; ++nfine; }
      }
      xbar[j] = (sw > 0.) ? swx / sw : 0.5 * (lo + hi);
      if (nfine <= 1)
        std::cerr << "[WARNING] " << tag << " pT bin [" << lo << "," << hi
                  << "): only " << nfine << " fine projection bin(s); <pT> ~ bin "
                  << "centre. Stored 2D is not finer than the analysis bins."
                  << std::endl;
    }
    return xbar;
  }

  // --------------------------------------------------------------------------
  // Closed-form weighted linear least squares for y = a + b*x.
  // Returns a, b and the 2x2 statistical covariance [[Va,Cab],[Cab,Vb]].
  // Only bins with sig > 0 and x in [xlo,xhi] are used.
  // --------------------------------------------------------------------------
  struct WLS { double a, b, Va, Vb, Cab; double chi2; int ndf; };

  WLS SolveWLS(const std::vector<double> &x,
               const std::vector<double> &y,
               const std::vector<double> &sig,
               double xlo, double xhi)
  {
    double S = 0, Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
    int nused = 0;
    for (size_t j = 0; j < x.size(); ++j)
    {
      if (sig[j] <= 0. || x[j] < xlo || x[j] > xhi) continue;
      const double w = 1.0 / (sig[j] * sig[j]);
      S += w; Sx += w * x[j]; Sy += w * y[j];
      Sxx += w * x[j] * x[j]; Sxy += w * x[j] * y[j];
      ++nused;
    }
    const double D = S * Sxx - Sx * Sx;
    WLS r;
    r.a = (Sxx * Sy - Sx * Sxy) / D;
    r.b = (S * Sxy - Sx * Sy) / D;
    r.Va = Sxx / D; r.Vb = S / D; r.Cab = -Sx / D;
    r.chi2 = 0.; 
    for (size_t j = 0; j < x.size(); ++j)
    {
      if (sig[j] <= 0. || x[j] < xlo || x[j] > xhi) continue;
      const double d = y[j] - (r.a + r.b * x[j]);
      r.chi2 += d * d / (sig[j] * sig[j]);
    }
    r.ndf = nused - 2;
    return r;
  }

  // --------------------------------------------------------------------------
  // FitResult: central params, total covariance, correlation, band-ready info.
  // --------------------------------------------------------------------------
  struct FitResult
  {
    double a, b;                  // central (best) fit to nominal points
    double Va, Vb, Cab;           // TOTAL covariance (stat + syst)
    double sa, sb, rho;           // total sigmas and correlation
    double sb_stat, sb_syst;      // slope error split (for the legend)
    double sa_stat, sa_syst;      // intercept error split
    double chi2; int ndf;
  };

  // h_nom + h_vars are the loaded profiles (index 0 = nominal). xbar holds the
  // yield-weighted <pT> per bin. Toys throw systematics ONLY; the statistical
  // covariance is taken from the nominal WLS fit and added afterwards.
  FitResult FitWithMCSystematics(TH1D *h_nom,
                                 const std::vector<TH1D *> &h_vars,
                                 const std::vector<std::vector<int>> &sources,
                                 const std::vector<double> &xbar,
                                 double xlo, double xhi,
                                 int ntoys, TRandom3 &rng,
                                 const std::string &tag = "")
  {
    const int n = h_nom->GetNbinsX();
    std::vector<double> y(n), sig(n);
    for (int j = 0; j < n; ++j)
    {
      y[j]   = h_nom->GetBinContent(j + 1);
      sig[j] = h_nom->GetBinError(j + 1);   // statistical error per bin
    }

    // central (nominal) fit -> central params + statistical covariance
    const WLS c = SolveWLS(xbar, y, sig, xlo, xhi);

    // per-source, per-bin 1-sigma shift vectors (symmetrised)
    std::vector<std::vector<double>> dsrc;   // [source][bin]
    dsrc.reserve(sources.size());
    for (const auto &grp : sources)
    {
      std::vector<double> d(n, 0.);
      if (grp.size() == 2)
      {
        TH1D *up = h_vars[grp[0]];
        TH1D *dn = h_vars[grp[1]];
        for (int j = 0; j < n; ++j)
          d[j] = 0.5 * (up->GetBinContent(j + 1) - dn->GetBinContent(j + 1));
      }
      else // single
      {
        TH1D *v = h_vars[grp[0]];
        for (int j = 0; j < n; ++j)
          d[j] = v->GetBinContent(j + 1) - y[j];
      }
      dsrc.push_back(std::move(d));
    }

    // ----- per-source diagnostic (exact analytic linear response) -----------
    // The fit (a,b) is linear in y, so a coherent shift v of the points moves
    // (a,b) by the WLS solution applied to v. Each source's 1-sigma slope/inter-
    // cept contribution is therefore exact; the quadrature sum equals the toy
    // systematic and reveals any single source that dominates the band.
    {
      double S = 0, Sx = 0, Sxx = 0;
      for (int j = 0; j < n; ++j)
      {
        if (sig[j] <= 0. || xbar[j] < xlo || xbar[j] > xhi) continue;
        const double w = 1.0 / (sig[j] * sig[j]);
        S += w; Sx += w * xbar[j]; Sxx += w * xbar[j] * xbar[j];
      }
      const double D = S * Sxx - Sx * Sx;
      double q_a2 = 0, q_b2 = 0;
      std::cout << "  [" << tag << "] per-source 1-sigma contributions "
                << "(da, db):" << std::endl;
      for (size_t s = 0; s < dsrc.size(); ++s)
      {
        double SyV = 0, SxyV = 0;
        for (int j = 0; j < n; ++j)
        {
          if (sig[j] <= 0. || xbar[j] < xlo || xbar[j] > xhi) continue;
          const double w = 1.0 / (sig[j] * sig[j]);
          SyV += w * dsrc[s][j]; SxyV += w * xbar[j] * dsrc[s][j];
        }
        const double da = (Sxx * SyV - Sx * SxyV) / D;
        const double db = (S * SxyV - Sx * SyV) / D;
        q_a2 += da * da; q_b2 += db * db;
        std::cout << "    src " << s << " idx{";
        for (size_t q = 0; q < sources[s].size(); ++q)
          std::cout << sources[s][q] << (q + 1 < sources[s].size() ? "," : "");
        std::cout << "}  da=" << da << "  db=" << db << std::endl;
      }
      std::cout << "    -> quad-sum syst:  da=" << std::sqrt(q_a2)
                << "  db=" << std::sqrt(q_b2) << std::endl;
    }

    // Monte Carlo: throw one Gaussian per source, refit, accumulate (a,b) moments
    double sa = 0, sb = 0, saa = 0, sbb = 0, sab = 0;
    std::vector<double> ytoy(n);
    for (int t = 0; t < ntoys; ++t)
    {
      ytoy = y;
      for (size_t s = 0; s < dsrc.size(); ++s)
      {
        const double g = rng.Gaus(0., 1.);   // ONE throw for the whole source
        for (int j = 0; j < n; ++j) ytoy[j] += g * dsrc[s][j];
      }
      const WLS w = SolveWLS(xbar, ytoy, sig, xlo, xhi);
      sa += w.a; sb += w.b;
      saa += w.a * w.a; sbb += w.b * w.b; sab += w.a * w.b;
    }
    const double inv = 1.0 / ntoys;
    const double ma = sa * inv, mb = sb * inv;
    const double Va_syst = saa * inv - ma * ma;
    const double Vb_syst = sbb * inv - mb * mb;
    const double Cab_syst = sab * inv - ma * mb;

    FitResult r;
    r.a = c.a; r.b = c.b;
    r.chi2 = c.chi2; r.ndf = c.ndf;
    r.Va = c.Va + Va_syst;
    r.Vb = c.Vb + Vb_syst;
    r.Cab = c.Cab + Cab_syst;
    r.sa = std::sqrt(r.Va);
    r.sb = std::sqrt(r.Vb);
    r.rho = r.Cab / (r.sa * r.sb);
    r.sa_stat = std::sqrt(c.Va);  r.sa_syst = std::sqrt(std::max(0., Va_syst));
    r.sb_stat = std::sqrt(c.Vb);  r.sb_syst = std::sqrt(std::max(0., Vb_syst));
    return r;
  }

  // --------------------------------------------------------------------------
  // Build the Delta-chi2 contour band as a filled TGraph.
  //
  //   band(x) = yhat(x) +/- sqrt( k * Var(yhat(x)) ),
  //   Var(yhat(x)) = Va + x^2 Vb + 2 x Cab   (TOTAL covariance)
  //
  // k = 2.30 is the 68% CL Delta-chi2 for 2 parameters (use 6.18 for 95%).
  // This is exactly the envelope of all lines inside the chi2_best + k contour.
  // --------------------------------------------------------------------------
  TGraph *BuildContourBand(const FitResult &r, double xlo, double xhi,
                           int color, float alpha, double k = 2.30, int npts = 100)
  {
    std::vector<double> xu, yu, xl, yl;
    for (int i = 0; i <= npts; ++i)
    {
      const double x = xlo + (xhi - xlo) * i / npts;
      const double var = r.Va + x * x * r.Vb + 2. * x * r.Cab;
      const double half = std::sqrt(std::max(0., k * var));
      const double yhat = r.a + r.b * x;
      xu.push_back(x); yu.push_back(yhat + half);
      xl.push_back(x); yl.push_back(yhat - half);
    }
    // closed polygon: upper edge left->right, lower edge right->left
    std::vector<double> px, py;
    px.insert(px.end(), xu.begin(), xu.end());
    py.insert(py.end(), yu.begin(), yu.end());
    for (int i = (int)xl.size() - 1; i >= 0; --i) { px.push_back(xl[i]); py.push_back(yl[i]); }
    px.push_back(xu.front()); py.push_back(yu.front());

    TGraph *g = new TGraph((int)px.size(), px.data(), py.data());
    g->SetFillColorAlpha(color, alpha);
    g->SetLineWidth(0);
    return g;
  }

  TF1 *MakeLine(const FitResult &r, double xlo, double xhi, const std::string &name)
  {
    TF1 *f = new TF1(name.c_str(), "[0] + [1]*x", xlo, xhi);
    f->SetParameters(r.a, r.b);
    return f;
  }

  // --------------------------------------------------------------------------
  // BuildSlopeOnlyBand
  //
  // Pin the line at the slope/intercept decorrelation pivot x0 = -Cab/Vb (where
  // the two parameters are independent) and let ONLY the slope vary:
  //
  //     band(x) = yhat(x) +/- sqrt(k) * sigma_b * |x - x0|
  //
  // Width is zero at x0 and grows linearly outward. This displays the slope
  // (tilt) uncertainty with the intercept/normalisation uncertainty — into which
  // any fully-correlated scale systematic collects — REMOVED from the band.
  // k = 1.0 is the 1-sigma slope band (Delta-chi2 = 1 for the single parameter).
  // --------------------------------------------------------------------------
  TGraph *BuildSlopeOnlyBand(const FitResult &r, double xlo, double xhi,
                             int color, float alpha, double k = 1.0, int npts = 100)
  {
    const double x0 = -r.Cab / r.Vb;
    std::vector<double> xu, yu, xl, yl;
    for (int i = 0; i <= npts; ++i)
    {
      const double x = xlo + (xhi - xlo) * i / npts;
      const double yhat = r.a + r.b * x;
      const double half = std::sqrt(k) * r.sb * std::fabs(x - x0);
      xu.push_back(x); yu.push_back(yhat + half);
      xl.push_back(x); yl.push_back(yhat - half);
    }
    std::vector<double> px, py;
    px.insert(px.end(), xu.begin(), xu.end());
    py.insert(py.end(), yu.begin(), yu.end());
    for (int i = (int)xl.size() - 1; i >= 0; --i) { px.push_back(xl[i]); py.push_back(yl[i]); }
    px.push_back(xu.front()); py.push_back(yu.front());

    TGraph *g = new TGraph((int)px.size(), px.data(), py.data());
    g->SetFillColorAlpha(color, alpha);
    g->SetLineWidth(0);
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
void figure6_plots_mc(
    const char *dijet_input        = "output_mbd_correct_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root",
    const char *efrac_input        = "output_mbd_correct_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root",
    const char *dijet_herwig_input = "output_mbd_correct_herwig_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root",
    const char *efrac_herwig_input = "output_mbd_correct_herwig_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root",
    const char *output_name        = "h_jet_pt_bin_unc_figure_mctoys.pdf",
    double      ptlo               = 21.,
    double      pthi               = 63.5,
    int         ntoys              = 1000,
    double      band_k             = 2.3,    // 1.0 = pointwise 1-sigma fit band;
                                             // 2.30 = joint 68% region for 2 params
    int         slope_only_band    = 0)      // 1 = slope-only band (intercept/scale
                                             //     uncertainty removed, pinned at pivot)
{
  gROOT->LoadMacro("sPhenixStyle.C");
  gROOT->ProcessLine("SetsPhenixStyle()");

  // ==========================================================================
  //  Yield-weighted <pT> per bin (the "fold with the spectrum" step) is taken
  //  from ProjectionX of the 2D histograms below — see ProjectionXbar. The X
  //  axis of those 2D histograms is the jet pT, so projecting gives dN/dpT and
  //  the per-bin yield-weighted mean is the correct abscissa for the linear fit.
  // ==========================================================================

  TFile *f_dijet        = TFile::Open(dijet_input);
  TFile *f_efrac        = TFile::Open(efrac_input);
  TFile *f_dijet_herwig = TFile::Open(dijet_herwig_input);
  TFile *f_efrac_herwig = TFile::Open(efrac_herwig_input);
  if (!f_dijet || f_dijet->IsZombie() || !f_efrac || f_efrac->IsZombie() ||
      !f_dijet_herwig || f_dijet_herwig->IsZombie() ||
      !f_efrac_herwig || f_efrac_herwig->IsZombie())
  {
    std::cerr << "Could not open input files" << std::endl;
    return;
  }

  // -------------------- load all systematic variation profiles --------------
  const std::vector<std::string> syst = GetSystList();
  std::vector<TH1D *> h_dijet, h_efrac;
  h_dijet.reserve(syst.size());
  h_efrac.reserve(syst.size());

  for (size_t i = 0; i < syst.size() - 1; ++i)
  {
    TH1D *hd = BuildProfileHistogram(f_dijet, syst[i], Form("h_dijet_%zu", i));
    TH1D *he = BuildProfileHistogram(f_efrac, syst[i], Form("h_efrac_%zu", i));
    if (!hd || !he) { std::cerr << "Failed at syst index " << i << std::endl; return; }
    h_dijet.push_back(hd);
    h_efrac.push_back(he);
  }
  TH1D *hd_hw = BuildProfileHistogram(f_dijet_herwig, syst.back(), Form("h_dijet_herwig_%zu", syst.size()-1));
  TH1D *he_hw = BuildProfileHistogram(f_efrac_herwig, syst.back(), Form("h_efrac_herwig_%zu", syst.size()-1));
  if (!hd_hw || !he_hw) { std::cerr << "Failed at Herwig syst index" << std::endl; return; }
  h_dijet.push_back(hd_hw);
  h_efrac.push_back(he_hw);

  // -------------------- truth profiles --------------------------------------
  TH1D *h_truth_pythia = BuildProfileHistogramFromKey(f_dijet,        "h_truth_calib_dijet", "h_truth_pythia");
  TH1D *h_truth_herwig = BuildProfileHistogramFromKey(f_dijet_herwig, "h_truth_calib_dijet", "h_truth_herwig");
  if (!h_truth_pythia || !h_truth_herwig) { std::cerr << "Failed to load truth histograms" << std::endl; return; }

  // -------------------- normalisation ---------------------------------------
  const double norm = 3.0 / (2.2 * 2. * M_PI);
  for (size_t i = 0; i < h_dijet.size(); ++i) { h_dijet[i]->Scale(norm); h_efrac[i]->Scale(norm); }
  h_truth_pythia->Scale(norm);
  h_truth_herwig->Scale(norm);

  TH1D *h_dijet_nom = h_dijet[0];
  TH1D *h_efrac_nom = h_efrac[0];

  // -------------------- resolve per-bin <pT> from ProjectionX ---------------
  // Analysis pT bin edges come from the nominal profile axis. <pT> for the data
  // fits is taken from the unfolded 2D (h_unfold_<nominal>); the truth fits use
  // their own truth 2D projections so each slope sits at its proper abscissa.
  const int nb = h_dijet_nom->GetNbinsX();
  std::vector<double> ptedges(nb + 1);
  for (int j = 0; j < nb; ++j) ptedges[j] = h_dijet_nom->GetXaxis()->GetBinLowEdge(j + 1);
  ptedges[nb] = h_dijet_nom->GetXaxis()->GetBinUpEdge(nb);

  const std::string nom_key = "h_unfold_" + syst[0];
  std::vector<double> xbar_dijet = ProjectionXbar(f_dijet, nom_key, ptedges, "dijet_unfold");
  std::vector<double> xbar_efrac = ProjectionXbar(f_efrac, nom_key, ptedges, "efrac_unfold");
  std::vector<double> xbar_tpy   = ProjectionXbar(f_dijet,        "h_truth_calib_dijet", ptedges, "truth_pythia");
  std::vector<double> xbar_thw   = ProjectionXbar(f_dijet_herwig, "h_truth_calib_dijet", ptedges, "truth_herwig");

  std::cout << "<pT> dijet :";
  for (double v : xbar_dijet) std::cout << " " << v;
  std::cout << "\n<pT> efrac :";
  for (double v : xbar_efrac) std::cout << " " << v;
  std::cout << std::endl;

  // -------------------- MC fits ---------------------------------------------
  const std::vector<std::vector<int>> sources = GetSystSources();
  TRandom3 rng(12345);

  FitResult rd = FitWithMCSystematics(h_dijet_nom, h_dijet, sources, xbar_dijet, ptlo, pthi, ntoys, rng, "dijet");
  FitResult re = FitWithMCSystematics(h_efrac_nom, h_efrac, sources, xbar_efrac, ptlo, pthi, ntoys, rng, "efrac");

  std::cout << "Dijet : a = " << rd.a << " +/- " << rd.sa
            << "  b = " << rd.b << " +/- " << rd.sb
            << "  rho(a,b) = " << rd.rho
            << "  chi2/ndf = " << rd.chi2 << "/" << rd.ndf << std::endl;
  std::cout << "Efrac : a = " << re.a << " +/- " << re.sa
            << "  b = " << re.b << " +/- " << re.sb
            << "  rho(a,b) = " << re.rho
            << "  chi2/ndf = " << re.chi2 << "/" << re.ndf << std::endl;

  // central best-fit lines
  TF1 *f_dijet_nom = MakeLine(rd, ptlo, pthi, "f_dijet_nom");
  TF1 *f_efrac_nom = MakeLine(re, ptlo, pthi, "f_efrac_nom");
  f_dijet_nom->SetLineColor(kBlack);    f_dijet_nom->SetLineWidth(2); f_dijet_nom->SetLineStyle(1);
  f_efrac_nom->SetLineColor(kAzure + 2); f_efrac_nom->SetLineWidth(2); f_efrac_nom->SetLineStyle(1);

  // Uncertainty bands from the TOTAL (stat+syst) covariance.
  //   slope_only_band = 0: full band, width = sqrt(band_k * Var(yhat(x)))
  //       band_k = 1.0  -> pointwise 1-sigma band (comparable to data error bars)
  //       band_k = 2.30 -> envelope of the joint 68% CL region for 2 parameters
  //   slope_only_band = 1: slope-only band pinned at the decorrelation pivot,
  //       showing the tilt uncertainty with the normalisation/scale removed
  //       (use band_k = 1.0 for the 1-sigma slope band)
  TGraph *g_dijet_band, *g_efrac_band;
  if (slope_only_band)
  {
    g_dijet_band = BuildSlopeOnlyBand(rd, ptlo, pthi, kGray + 2, 0.55, band_k);
    g_efrac_band = BuildSlopeOnlyBand(re, ptlo, pthi, kAzure - 9, 0.55, band_k);
  }
  else
  {
    g_dijet_band = BuildContourBand(rd, ptlo, pthi, kGray + 2, 0.55, band_k);
    g_efrac_band = BuildContourBand(re, ptlo, pthi, kAzure - 9, 0.55, band_k);
  }

  // -------------------- truth fit lines (central only) ----------------------
  // (truth uses the same SolveWLS for a consistent slope definition)
  std::vector<double> yt_py(nb), st_py(nb), yt_hw(nb), st_hw(nb);
  for (int j = 0; j < nb; ++j)
  {
    yt_py[j] = h_truth_pythia->GetBinContent(j + 1); st_py[j] = h_truth_pythia->GetBinError(j + 1);
    yt_hw[j] = h_truth_herwig->GetBinContent(j + 1); st_hw[j] = h_truth_herwig->GetBinError(j + 1);
  }
  WLS tpy = SolveWLS(xbar_tpy, yt_py, st_py, ptlo, pthi);
  WLS thw = SolveWLS(xbar_thw, yt_hw, st_hw, ptlo, pthi);

  TF1 *f_truth_pythia = new TF1("f_truth_pythia", "[0] + [1]*x", ptlo, pthi);
  f_truth_pythia->SetParameters(tpy.a, tpy.b);
  f_truth_pythia->SetLineColor(kRed + 1); f_truth_pythia->SetLineWidth(2); f_truth_pythia->SetLineStyle(7);
  TF1 *f_truth_herwig = new TF1("f_truth_herwig", "[0] + [1]*x", ptlo, pthi);
  f_truth_herwig->SetParameters(thw.a, thw.b);
  f_truth_herwig->SetLineColor(kSpring - 6); f_truth_herwig->SetLineWidth(2); f_truth_herwig->SetLineStyle(9);

  h_truth_pythia->SetStats(0); h_truth_pythia->SetLineWidth(0);
  h_truth_pythia->SetMarkerColor(kRed + 1); h_truth_pythia->SetMarkerStyle(21);
  h_truth_herwig->SetStats(0); h_truth_herwig->SetLineWidth(0);
  h_truth_herwig->SetMarkerColor(kSpring - 6); h_truth_herwig->SetMarkerStyle(22);

  // -------------------- canvas ----------------------------------------------
  TCanvas *c = new TCanvas("c_slope_unc", "", 700, 700);
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.12);

  h_dijet_nom->SetStats(0);
  h_dijet_nom->SetLineColor(kBlack); h_dijet_nom->SetLineWidth(2);
  h_dijet_nom->SetMarkerColor(kBlack); h_dijet_nom->SetMarkerStyle(20);
  h_dijet_nom->GetXaxis()->SetTitle("p_{T,lead} [GeV]");
  h_dijet_nom->GetXaxis()->SetTitleSize(25); h_dijet_nom->GetXaxis()->SetTitleFont(43);
  h_dijet_nom->GetXaxis()->SetLabelFont(43); h_dijet_nom->GetXaxis()->SetLabelSize(25);
  h_dijet_nom->GetXaxis()->SetRangeUser(ptlo, pthi);
  h_dijet_nom->GetYaxis()->SetRangeUser(0.3, 0.8);
  h_dijet_nom->GetYaxis()->SetTitle("<#SigmaE_{T}/#delta#eta#delta#phi> [GeV]");
  h_dijet_nom->GetYaxis()->SetTitleSize(25); h_dijet_nom->GetYaxis()->SetTitleFont(43);
  h_dijet_nom->GetYaxis()->SetTitleOffset(1.5);
  h_dijet_nom->GetYaxis()->SetLabelFont(43); h_dijet_nom->GetYaxis()->SetLabelSize(25);

  h_efrac_nom->SetStats(0);
  h_efrac_nom->SetLineColor(kAzure + 2); h_efrac_nom->SetLineWidth(2);
  h_efrac_nom->SetMarkerColor(kAzure + 2); h_efrac_nom->SetMarkerStyle(33); h_efrac_nom->SetMarkerSize(1.4);

  TGraphErrors *g_dijet_pts = HistToGraph(h_dijet_nom);
  TGraphErrors *g_efrac_pts = HistToGraph(h_efrac_nom);

  h_dijet_nom->SetLineWidth(0);

  // -------------------- draw ------------------------------------------------
  // bands behind, then truth, then best-fit lines, then DATA POINTS (stat bars) on top
  h_dijet_nom->Draw("E");
  g_dijet_band->Draw("f same");
  g_efrac_band->Draw("f same");

  f_truth_pythia->Draw("same");
  f_truth_herwig->Draw("same");
  h_truth_pythia->Draw("P same");
  h_truth_herwig->Draw("P same");

  f_dijet_nom->Draw("same");
  f_efrac_nom->Draw("same");

  g_dijet_pts->Draw("P same");   // data points carry STATISTICAL errors only
  g_efrac_pts->Draw("P same");

  // -------------------- legend ----------------------------------------------
  f_truth_pythia->SetMarkerStyle(21); f_truth_pythia->SetMarkerColor(kRed + 1);
  f_truth_herwig->SetMarkerStyle(22); f_truth_herwig->SetMarkerColor(kSpring - 6);

  f_dijet_nom->SetMarkerSize(1.2); f_efrac_nom->SetMarkerSize(1.4);
  f_dijet_nom->SetMarkerStyle(20); f_dijet_nom->SetMarkerColor(kBlack);
  f_efrac_nom->SetMarkerStyle(33); f_efrac_nom->SetMarkerColor(kAzure + 2);


    std::cout << "Dijet : a = " << rd.a << " +/- " << rd.sa
            << "  b = " << rd.b << " +/- " << rd.sb
            << "  rho(a,b) = " << rd.rho
            << "  chi2/ndf = " << rd.chi2 << "/" << rd.ndf << std::endl;
  std::cout << "Efrac : a = " << re.a << " +/- " << re.sa
            << "  b = " << re.b << " +/- " << re.sb
            << "  rho(a,b) = " << re.rho
            << "  chi2/ndf = " << re.chi2 << "/" << re.ndf << std::endl;

  TLegend *leg = new TLegend(0.15, 0.58, 0.4, 0.92);
  leg->SetTextFont(43); leg->SetTextSize(19);
  leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry("", "#bf{#it{sPHENIX}} Internal", "");
  leg->AddEntry("", "200 GeV p+p  anti-k_{t} R=0.4  |#eta_{jet}|<0.7", "");
  leg->AddEntry(f_dijet_nom,
                Form("Ex. Dijet: slope = %.4f #kern[-0.15]{#pm %.4f}, y-inter = %.4f #kern[-0.15]{#pm %.4f}", rd.b, rd.sb, rd.a, rd.sa), "lp");
  leg->AddEntry(h_efrac_nom,
                Form("Inc. Jet: slope = %.4f #kern[-0.15]{#pm %.4f}, y-inter = %.4f #kern[-0.15]{#pm %.4f}", re.b, re.sb, re.a, re.sa), "lp");
  leg->AddEntry(f_truth_pythia, Form("PYTHIA truth: slope = %.4f", tpy.b), "lp");
  leg->AddEntry(f_truth_herwig, Form("HERWIG truth: slope = %.4f", thw.b), "lp");
  leg->AddEntry("", (slope_only_band
                       ? "Linear fit, slope uncertainty only (scale removed)"
                       : Form("Linear fit: %s band",
                              (band_k < 1.5 ? "1#sigma" : "68% CI 2-par"))), "");
  leg->AddEntry(g_dijet_band, "Excluded Dijet", "f");
  leg->AddEntry(g_efrac_band, "Inclusive Jet", "f");
  leg->Draw();

  c->Update();
  c->SaveAs(output_name);
  std::cout << "Saved: " << output_name << std::endl;
}