// plot_comparison.C
// Compares 0mrad vs 1.5mrad for:
//   1) ProfileX of h_ue_pt_transverse_record (from analysis_data_run28_output)
//   2) ProfileX of h_calibjet_pt_dijet_eff   (from analysis_data_run28_output, rebinned)
//   3) ProfileX of h_unfold_calib_dijet_reweight_trim_10_3 (from sphenix_primary_run28_output_files, rebinned)

#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <string>
#include "unfold_Def.h"

// ─── Helper: rebin a TH2D to variable bins, return ownership to caller ───────
TH2D* RebinToVariable(TH2D* hu, const std::string& newname,
                       int nptbins, double* ptbins,
                       int netbins, double* etbins)
{
    TH2D* hv = new TH2D(newname.c_str(), hu->GetTitle(),
                        nptbins, ptbins, netbins, etbins);
    hv->Sumw2();
    for (int ix = 1; ix <= hu->GetNbinsX(); ix++) {
        for (int iy = 1; iy <= hu->GetNbinsY(); iy++) {
            hv->SetBinContent(ix, iy, hu->GetBinContent(ix, iy));
            hv->SetBinError  (ix, iy, hu->GetBinError  (ix, iy));
        }
    }
    return hv;
}

// ─── Helper: style a profile ─────────────────────────────────────────────────
void StyleProfile(TH1* h, Color_t col, Style_t marker, const char* xtitle, const char* ytitle)
{
    h->SetMarkerColor(col);
    h->SetLineColor(col);
    h->SetMarkerStyle(marker);
    h->SetMarkerSize(1.1);
    h->GetXaxis()->SetTitle(xtitle);
    h->GetYaxis()->SetTitle(ytitle);
    h->GetYaxis()->SetTitleOffset(1.4);
}

// ─── Helper: build ratio histogram (h1 / h0), propagate errors ───────────────
TH1D* MakeRatio(TH1D* h1, TH1D* h0, const std::string& name)
{
    // build fresh histogram (no cloning!)
    const TArrayD* bins = h1->GetXaxis()->GetXbins();
    TH1D* hr = nullptr;

    if (bins->GetSize() > 0) {
        hr = new TH1D(name.c_str(), "", h1->GetNbinsX(), bins->GetArray());
    } else {
        hr = new TH1D(name.c_str(), "", h1->GetNbinsX(),
                      h1->GetXaxis()->GetXmin(),
                      h1->GetXaxis()->GetXmax());
    }

    hr->SetDirectory(0);
    hr->Sumw2();

    int nbins = std::min(h1->GetNbinsX(), h0->GetNbinsX());

    for (int i = 1; i <= nbins; i++) {
        double num  = h1->GetBinContent(i);
        double den  = h0->GetBinContent(i);
        double enum_ = h1->GetBinError(i);
        double eden  = h0->GetBinError(i);

        if (den == 0.0) continue;

        double r = num / den;
        double err = (num != 0.0)
            ? r * std::sqrt((enum_/num)*(enum_/num) + (eden/den)*(eden/den))
            : 0.0;

        hr->SetBinContent(i, r);
        hr->SetBinError(i, err);
    }

    return hr;
}

// ─── Draw comparison + ratio canvas ──────────────────────────────────────────
// h0 = 0 mrad (denominator / blue), h1 = 1.5 mrad (numerator / red)
void DrawComparisonWithRatio(TH1D* h0, TH1D* h1,
                              const char* title,
                              const char* label0, const char* label1,
                              const char* xtitle,
                              const char* ytitle,
                              const char* outname)
{
    TCanvas* c = new TCanvas(outname, title, 800, 900);

    // ── Upper pad (comparison) ────────────────────────────────────────────────
    TPad* pTop = new TPad(Form("%s_top", outname), "", 0.0, 0.32, 1.0, 1.0);
    pTop->SetBottomMargin(0.03);
    pTop->SetLeftMargin(0.14);
    pTop->SetTopMargin(0.10);
    pTop->Draw();
    pTop->cd();

    double ymax = std::max(h0->GetMaximum(), h1->GetMaximum()) * 1.3;
    double ymin = std::min({h0->GetMinimum(), h1->GetMinimum(), 0.0}) * 1.1;
    h0->SetMinimum(ymin);
    h0->SetMaximum(ymax);
    h0->SetTitle(title);

    // suppress x-axis labels/title on top pad
    h0->GetXaxis()->SetLabelSize(0);
    h0->GetXaxis()->SetTitleSize(0);
    h0->GetYaxis()->SetTitle(ytitle);
    h0->GetYaxis()->SetTitleSize(0.055);
    h0->GetYaxis()->SetLabelSize(0.050);
    h0->GetYaxis()->SetTitleOffset(1.2);

    h1->GetXaxis()->SetLabelSize(0);
    h1->GetXaxis()->SetTitleSize(0);

    if (!strcmp(title,"h_ue_pt_transverse_record")) {
        h0->GetXaxis()->SetRangeUser(14,63);
    } else {
        h0->GetXaxis()->SetRangeUser(21,63);
    }

    h0->Draw("E1");
    h1->Draw("E1 SAME");

    TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.050);
    leg->AddEntry(h0, label0, "lp");
    leg->AddEntry(h1, label1, "lp");
    leg->Draw();

    // ── Lower pad (ratio) ─────────────────────────────────────────────────────
    c->cd();
    TPad* pBot = new TPad(Form("%s_bot", outname), "", 0.0, 0.0, 1.0, 0.32);
    pBot->SetTopMargin(0.03);
    pBot->SetBottomMargin(0.30);
    pBot->SetLeftMargin(0.14);
    pBot->Draw();
    pBot->cd();

    TH1D* hr = MakeRatio(h1, h0, Form("ratio_%s", outname));

    hr->SetMarkerColor(kBlack);
    hr->SetLineColor(kBlack);
    hr->SetMarkerStyle(20);
    hr->SetMarkerSize(1.0);
    hr->SetTitle("");

    // auto-range the ratio axis with some padding
    double rmax = 1.5, rmin = 0.5;
    for (int i = 1; i <= hr->GetNbinsX(); i++) {
        double v = hr->GetBinContent(i);
        double e = hr->GetBinError(i);
        if (v == 0.0) continue;
        if (v + e > rmax) rmax = v + e;
        if (v - e < rmin) rmin = v - e;
    }
    rmax = rmax * 1.15;
    rmin = (rmin > 0) ? rmin * 0.85 : rmin * 1.15;

    hr->SetMinimum(rmin);
    hr->SetMaximum(rmax);

    hr->GetXaxis()->SetTitle(xtitle);
    hr->GetXaxis()->SetTitleSize(0.13);
    hr->GetXaxis()->SetLabelSize(0.11);
    hr->GetXaxis()->SetTitleOffset(1.0);
    hr->GetYaxis()->SetTitle(Form("%s / %s", label1, label0));
    hr->GetYaxis()->SetTitleSize(0.10);
    hr->GetYaxis()->SetLabelSize(0.10);
    hr->GetYaxis()->SetTitleOffset(0.55);
    hr->GetYaxis()->SetNdivisions(505);
    hr->GetYaxis()->SetRangeUser(0.8,1.2);
    if (!strcmp(title,"h_ue_pt_transverse_record")) {
        hr->GetXaxis()->SetRangeUser(14,63);
    } else {
        hr->GetXaxis()->SetRangeUser(21,63);
    }

    hr->Draw("E1");

    // dashed reference line at ratio = 1
    double xlo = hr->GetXaxis()->GetBinLowEdge(1);
    double xhi = hr->GetXaxis()->GetBinUpEdge(hr->GetNbinsX());
    TLine* line = new TLine(21, 1.0, 63, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray+2);
    line->SetLineWidth(2);
    line->Draw();

    c->SaveAs(Form("%s.png", outname));
    c->SaveAs(Form("%s.pdf", outname));
    std::cout << "[saved] " << outname << ".png / .pdf\n";

    delete c;
}

// ─── Main macro ──────────────────────────────────────────────────────────────
void plot_0mrad_1_5mrad_comparison()
{
    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    // Input files
    const char* f_ana_0   = "analysis_data_run28_output/output_0mrad_pu_correct_data_dijet_bkg_cut.root";
    const char* f_ana_15  = "analysis_data_run28_output/output_1_5mrad_pu_correct_data_dijet_bkg_cut.root";
    const char* f_unf_0   = "sphenix_primary_run28_output_files/output_0mrad_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_run28_iter_3_1000toys.root";
    const char* f_unf_15  = "sphenix_primary_run28_output_files/output_1_5mrad_unfolded_data_8calibetbin_dijet_bkg_cut_calib_dijet_run28_iter_3_1000toys.root";

    TFile* fana0  = TFile::Open(f_ana_0,  "READ");
    TFile* fana15 = TFile::Open(f_ana_15, "READ");
    TFile* funf0  = TFile::Open(f_unf_0,  "READ");
    TFile* funf15 = TFile::Open(f_unf_15, "READ");

    if (!fana0  || fana0->IsZombie())  { std::cerr << "Cannot open " << f_ana_0  << "\n"; return; }
    if (!fana15 || fana15->IsZombie()) { std::cerr << "Cannot open " << f_ana_15 << "\n"; return; }
    if (!funf0  || funf0->IsZombie())  { std::cerr << "Cannot open " << f_unf_0  << "\n"; return; }
    if (!funf15 || funf15->IsZombie()) { std::cerr << "Cannot open " << f_unf_15 << "\n"; return; }

    // ── 1) h_ue_pt_transverse_record (no rebinning) ─────────────────────────
    {
        TH2D* h2_0  = (TH2D*)fana0 ->Get("h_ue_pt_transverse_record");
        TH2D* h2_15 = (TH2D*)fana15->Get("h_ue_pt_transverse_record");
        if (!h2_0 || !h2_15) {
            std::cerr << "Missing h_ue_pt_transverse_record in one of the analysis files\n";
        } else {
            TH1D* p0  = (TH1D*)h2_0 ->ProfileX("prof_ue_0mrad");
            TH1D* p15 = (TH1D*)h2_15->ProfileX("prof_ue_1_5mrad");
            StyleProfile(p0,  kBlue+1,  20, "p_{T} [GeV/c]", "#LT UE p_{T}^{transverse} #GT");
            StyleProfile(p15, kRed+1,   21, "p_{T} [GeV/c]", "#LT UE p_{T}^{transverse} #GT");
            DrawComparisonWithRatio(p0, p15,
                           "h_ue_pt_transverse_record",
                           "0 mrad", "1.5 mrad",
                           "p_{T}^{uncalib} [GeV]",
                           "#LT E_{T}^{reco} #GT [GeV]",
                           "sphenix_primary_run28_output_files/compare_0mrad_1_5mrad_h_ue_pt_transverse_record");
        }
    }

    // ── 2) h_calibjet_pt_dijet_eff (rebin to calib variable bins) ───────────
    {
        TH2D* h2_0  = (TH2D*)fana0 ->Get("h_calibjet_pt_dijet_eff");
        TH2D* h2_15 = (TH2D*)fana15->Get("h_calibjet_pt_dijet_eff");
        if (!h2_0 || !h2_15) {
            std::cerr << "Missing h_calibjet_pt_dijet_eff in one of the analysis files\n";
        } else {
            TH2D* hv0  = RebinToVariable(h2_0,  "h_var_calib_eff_0mrad",  calibnpt, calibptbins, calibnet, calibetbins);
            TH2D* hv15 = RebinToVariable(h2_15, "h_var_calib_eff_1_5mrad", calibnpt, calibptbins, calibnet, calibetbins);
            TH1D* p0   = (TH1D*)hv0 ->ProfileX("prof_calibeff_0mrad");
            TH1D* p15  = (TH1D*)hv15->ProfileX("prof_calibeff_1_5mrad");
            StyleProfile(p0,  kBlue+1, 20, "p_{T}^{calib} [GeV/c]", "#LT E_{T}^{reco} #GT");
            StyleProfile(p15, kRed+1,  21, "p_{T}^{calib} [GeV/c]", "#LT E_{T}^{reco} #GT");
            DrawComparisonWithRatio(p0, p15,
                           "h_calibjet_pt_dijet_eff",
                           "0 mrad", "1.5 mrad",
                            "p_{T}^{calib} [GeV]",
                            "#LT E_{T}^{reco} #GT [GeV]",
                           "sphenix_primary_run28_output_files/compare_0mrad_1_5mrad_h_calibjet_pt_dijet_eff");
        }
    }

    // ── 3) h_unfold_calib_dijet_reweight_trim_10_3 (rebin to truth variable bins) ──
    {
        const char* hname = "h_unfold_calib_dijet_reweight_trim_10_3";
        TH2D* h2_0  = (TH2D*)funf0 ->Get(hname);
        TH2D* h2_15 = (TH2D*)funf15->Get(hname);
        if (!h2_0 || !h2_15) {
            std::cerr << "Missing " << hname << " in one of the unfolded files\n";
        } else {
            TH2D* hv0  = RebinToVariable(h2_0,  "h_var_unfold_0mrad",  truthnpt, truthptbins, truthnet, truthetbins);
            TH2D* hv15 = RebinToVariable(h2_15, "h_var_unfold_1_5mrad", truthnpt, truthptbins, truthnet, truthetbins);
            TH1D* p0   = (TH1D*)hv0 ->ProfileX("prof_unfold_0mrad");
            TH1D* p15  = (TH1D*)hv15->ProfileX("prof_unfold_1_5mrad");
            StyleProfile(p0,  kBlue+1, 20, "p_{T}^{truth} [GeV/c]", "#LT E_{T}^{truth} #GT");
            StyleProfile(p15, kRed+1,  21, "p_{T}^{truth} [GeV/c]", "#LT E_{T}^{truth} #GT");
            DrawComparisonWithRatio(p0, p15,
                           "h_unfold_calib_dijet_reweight_trim_10_3",
                           "0 mrad", "1.5 mrad",
                           "p_{T}^{truth} [GeV]",
                           "#LT E_{T}^{truth} #GT [GeV]",
                           "sphenix_primary_run28_output_files/compare_0mrad_1_5mrad_h_unfold_calib_dijet_reweight_trim_10_3");
        }
    }

    fana0->Close();
    fana15->Close();
    funf0->Close();
    funf15->Close();

    std::cout << "\nDone. Output: compare_*.png / compare_*.pdf\n";
}