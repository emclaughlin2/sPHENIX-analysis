#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TColor.h>
#include <TGraphAsymmErrors.h>
#include <TROOT.h>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iterator>
#include <unfold_Def.h>

void draw_profile(std::vector<TH1D*> truth, std::vector<TH1D*> unfold, bool jet, std::vector<std::string> leg_tags, int max_iter, const char* output_name = nullptr)
{
    TCanvas* canvas = new TCanvas("canvas", "", 600, 800);

    //std::vector<std::string> colors = {"kBlack","kRed","kBlue","kGreen+2","kBlue-3","kMagneta+1"};

    TPad* pad1 = new TPad("pad1", "", 0, 0.5, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();

    truth[0]->SetStats(0);
    truth[0]->SetLineColor(max_iter + 2);
    truth[0]->SetMarkerColor(max_iter + 2);
    truth[0]->GetXaxis()->SetLabelSize(0);
    truth[1]->SetStats(0);
    truth[1]->SetLineColor(max_iter + 4);
    truth[1]->SetMarkerColor(max_iter + 4);
    truth[0]->GetXaxis()->SetRangeUser(17, 63);
    truth[1]->GetXaxis()->SetRangeUser(17, 63);

    for (int i = 0; i < max_iter; ++i) {
        unfold[i]->SetStats(0);
        if (i >= 4) {
            unfold[i]->SetLineColor(i + 2); // Use TColor::GetColor if needed
            unfold[i]->SetMarkerColor(i + 2);
        } else {
            unfold[i]->SetLineColor(i + 1); // Use TColor::GetColor if needed
            unfold[i]->SetMarkerColor(i + 1);
        }
        unfold[i]->GetXaxis()->SetLabelSize(0);
        unfold[i]->GetXaxis()->SetRangeUser(17, 63);
        unfold[i]->GetYaxis()->SetRangeUser(1,7);
    }


    unfold[0]->GetYaxis()->SetTitle("<#SigmaE_{T}> [GeV]");
    unfold[0]->GetXaxis()->SetTitle("p_{T,lead} [GeV]");
    unfold[0]->Draw();
    for (int i = 1; i < max_iter; ++i) {
        unfold[i]->Draw("same");
    }
    for (int i = 0; i < truth.size(); i++) {
        truth[i]->Draw("same");
    }

    TLegend* leg = new TLegend(0.17, 0.72, 0.6, 0.92);
    leg->SetNColumns(2);
    for (int i = 0; i < max_iter; ++i) {
        leg->AddEntry(unfold[i], leg_tags[i].c_str(), "lp");
    }
    leg->AddEntry(truth[0], "Pythia8 Truth", "lp");
    leg->AddEntry(truth[1], "Herwig Truth", "lp");
    leg->SetTextSize(0.04);
    leg->Draw();

    canvas->cd();
    TPad* pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.5);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();
    
    std::vector<TH1D*> ratios;
    /*
    TH1D* r0 = (TH1D*)truth[0]->Clone("truth_ratio0");
    r0->Divide(unfold[0]);
    TH1D* r1 = (TH1D*)truth[1]->Clone("truth_ratio1");
    r1->Divide(unfold[0]);
    for (int i = 1; i < r0->GetNbinsX() + 1; i++) {
        r0->SetBinError(i,0);
        r1->SetBinError(i,0);
    }
    */
    for (int i = 1; i < max_iter; ++i) {
        TH1D* ratio = (TH1D*)unfold[i]->Clone(Form("ratio%d", i));
        ratio->Divide(unfold[0]);
        for (int j = 1; j < ratio->GetNbinsX() + 1; j++) {
            ratio->SetBinError(j,0);
        }
        ratios.push_back(ratio);
    }
    std::cout << ratios.size() << std::endl;
    for (int i = 0; i < ratios.size(); i++) {
        for (int j = 2; j < ratios[i]->GetNbinsX(); j++) {
            std::cout << ratios[i]->GetBinContent(j) << " ";
        }
        std::cout << std::endl;
    }

    ratios[0]->GetYaxis()->SetTitle("Var/Nominal Ratio");
    ratios[0]->GetYaxis()->SetNdivisions(510);
    ratios[0]->GetYaxis()->SetRangeUser(0.5, 1.5);
    ratios[0]->GetYaxis()->SetTitleSize(25);
    ratios[0]->GetYaxis()->SetTitleFont(43);
    ratios[0]->GetYaxis()->SetTitleOffset(1.5);
    ratios[0]->GetYaxis()->SetLabelFont(43);
    ratios[0]->GetYaxis()->SetLabelSize(25);
    ratios[0]->GetXaxis()->SetTitle("p_{T} [GeV]");
    ratios[0]->GetXaxis()->SetTitleSize(25);
    ratios[0]->GetXaxis()->SetTitleFont(43);
    ratios[0]->GetXaxis()->SetTitleOffset(0);
    ratios[0]->GetXaxis()->SetLabelFont(43);
    ratios[0]->GetXaxis()->SetLabelSize(25);
    ratios[0]->GetXaxis()->SetRangeUser(17, 63);

    std::cout << ratios.size() << std::endl;
    for (int i = 0; i < ratios.size(); ++i) {
        ratios[i]->GetXaxis()->SetRangeUser(17, 63);
        ratios[i]->GetYaxis()->SetRangeUser(0.5,1.5);
        std::cout << i << " " << leg_tags[i+1] << " ";
        for (int j = 2; j < ratios[i]->GetNbinsX(); j++) {
            std::cout << ratios[i]->GetBinContent(j) << " ";
        }
        std::cout << std::endl;
        if (i == 0) ratios[i]->Draw("hist");
        else ratios[i]->Draw("hist, same");
    }

    TLine* line0 = new TLine(17, 1.0, 63, 1.0);
    TLine* line1 = new TLine(17, 0.95, 63, 0.95);
    TLine* line2 = new TLine(17, 1.05, 63, 1.05);
    line0->SetLineStyle(1);
    line1->SetLineStyle(2);
    line2->SetLineStyle(2);
    line0->Draw("same");
    line1->Draw("same");
    line2->Draw("same");
    
    canvas->Update();
    canvas->Draw();
    if (output_name) canvas->SaveAs(output_name);
}

void draw_test_profile(std::vector<TH1D*> truth, std::vector<TH1D*> unfold, bool jet, std::vector<std::string> leg_tags, int max_iter, const char* output_name = nullptr)
{
    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);

    //std::vector<std::string> colors = {"kBlack","kRed","kBlue","kGreen+2","kBlue-3","kMagneta+1"};

    truth[0]->SetStats(0);
    truth[0]->SetLineColor(max_iter + 2);
    truth[0]->SetMarkerColor(max_iter + 2);
    truth[1]->SetStats(0);
    truth[1]->SetLineColor(max_iter + 4);
    truth[1]->SetMarkerColor(max_iter + 4);
    truth[0]->GetXaxis()->SetRangeUser(17, 63);
    truth[1]->GetXaxis()->SetRangeUser(17, 63);

    for (int i = 0; i < max_iter; ++i) {
        unfold[i]->SetStats(0);
        if (i >= 4) {
            unfold[i]->SetLineColor(i + 2); // Use TColor::GetColor if needed
            unfold[i]->SetMarkerColor(i + 2);
        } else {
            unfold[i]->SetLineColor(i + 1); // Use TColor::GetColor if needed
            unfold[i]->SetMarkerColor(i + 1);
        }
        unfold[i]->GetXaxis()->SetRangeUser(17, 63);
        unfold[i]->GetYaxis()->SetRangeUser(1,7);
    }


    unfold[0]->GetYaxis()->SetTitle("<#SigmaE_{T}> [GeV]");
    unfold[0]->GetXaxis()->SetTitle("p_{T,lead} [GeV]");
    unfold[0]->Draw();
    for (int i = 1; i < max_iter; ++i) {
        unfold[i]->Draw("same");
    }
    for (int i = 0; i < truth.size(); i++) {
        truth[i]->Draw("same");
    }

    TLegend* leg = new TLegend(0.17, 0.72, 0.92, 0.92);
    //leg->SetNColumns(2);
    for (int i = 0; i < max_iter; ++i) {
        leg->AddEntry(unfold[i], leg_tags[i].c_str(), "lp");
    }
    leg->AddEntry(truth[0], "Pythia8 Truth (PU correction)", "lp");
    leg->AddEntry(truth[1], "Pythia8 Truth (no PU correction)", "lp");
    leg->SetTextSize(0.04);
    leg->Draw();

    canvas->Draw();
    if (output_name) canvas->SaveAs(output_name);
}

// Draw nominal + total syst band vs Pythia8 truth, with data/Pythia8 ratio
// unfold[0] = nominal, unfold[1:] = systematic variations; truth[0] = Pythia8
void draw_result_with_syst(std::vector<TH1D*> truth, std::vector<TH1D*> unfold, const char* output_name = nullptr)
{
    int nbins = unfold[0]->GetNbinsX();

    // Compute asymmetric total systematic uncertainty via quad sum of up/down shifts
    std::vector<double> x(nbins), y(nbins), exl(nbins), exh(nbins), eyl(nbins), eyh(nbins);
    for (int b = 0; b < nbins; b++) {
        x[b]   = unfold[0]->GetBinCenter(b + 1);
        y[b]   = unfold[0]->GetBinContent(b + 1);
        exl[b] = unfold[0]->GetBinWidth(b + 1) / 2.0;
        exh[b] = exl[b];
        double sum_up2 = 0, sum_dn2 = 0;
        for (int i = 1; i < (int)unfold.size(); i++) {
            double delta = unfold[i]->GetBinContent(b + 1) - y[b];
            if (delta >= 0) sum_up2 += delta * delta;
            else            sum_dn2 += delta * delta;
        }
        eyh[b] = std::sqrt(sum_up2);
        eyl[b] = std::sqrt(sum_dn2);
    }

    TGraphAsymmErrors* g_syst = new TGraphAsymmErrors(nbins, x.data(), y.data(),
                                                       exl.data(), exh.data(),
                                                       eyl.data(), eyh.data());
    g_syst->SetFillColorAlpha(kAzure - 9, 0.6);
    g_syst->SetFillStyle(1001);
    g_syst->SetLineWidth(0);

    // --- Canvas ---
    TCanvas* canvas = new TCanvas("canvas_result", "", 600, 800);

    TPad* pad1 = new TPad("pad1_result", "", 0, 0.5, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();

    unfold[0]->SetStats(0);
    unfold[0]->SetLineColor(kBlack);
    unfold[0]->SetMarkerColor(kBlack);
    unfold[0]->SetMarkerStyle(20);
    unfold[0]->GetXaxis()->SetLabelSize(0);
    unfold[0]->GetXaxis()->SetRangeUser(17, 63);
    unfold[0]->GetYaxis()->SetRangeUser(1, 7);
    unfold[0]->GetYaxis()->SetTitle("<#SigmaE_{T}> [GeV]");

    truth[0]->SetStats(0);
    truth[0]->SetLineColor(kRed + 1);
    truth[0]->SetLineWidth(2);
    truth[0]->SetMarkerColor(kRed + 1);
    truth[0]->GetXaxis()->SetRangeUser(17, 63);

    unfold[0]->Draw("E1");
    g_syst->Draw("2 same");
    unfold[0]->Draw("E1 same");
    truth[0]->Draw("hist same");

    TLegend* leg = new TLegend(0.17, 0.65, 0.78, 0.92);
    leg->AddEntry("","#bf{#it{sPHENIX}} Internal","");
    leg->AddEntry("","200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
    leg->AddEntry("","Exclusive dijet","");
    leg->AddEntry(unfold[0], "Unfolded data (stat. unc.)", "lp");
    leg->AddEntry(g_syst, "Total syst. uncertainty", "f");
    leg->AddEntry(truth[0], "Pythia8 Truth", "l");
    leg->SetTextSize(0.04);
    leg->Draw();

    // --- Ratio pad ---
    canvas->cd();
    TPad* pad2 = new TPad("pad2_result", "", 0, 0.0, 1, 0.5);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();

    // Nominal / Pythia8 ratio with stat errors from data only
    TH1D* ratio_nominal = (TH1D*)unfold[0]->Clone("ratio_nominal_result");
    for (int b = 1; b <= nbins; b++) {
        double pyt  = truth[0]->GetBinContent(b);
        double nom  = unfold[0]->GetBinContent(b);
        double stat = unfold[0]->GetBinError(b);
        if (pyt > 0) {
            ratio_nominal->SetBinContent(b, nom / pyt);
            ratio_nominal->SetBinError(b, stat / pyt);
        } else {
            ratio_nominal->SetBinContent(b, 0);
            ratio_nominal->SetBinError(b, 0);
        }
    }

    // Syst band for ratio: same absolute syst / pythia8
    std::vector<double> ry(nbins), reyl(nbins), reyh(nbins);
    for (int b = 0; b < nbins; b++) {
        double pyt = truth[0]->GetBinContent(b + 1);
        ry[b]   = (pyt > 0) ? y[b]   / pyt : 0;
        reyl[b] = (pyt > 0) ? eyl[b] / pyt : 0;
        reyh[b] = (pyt > 0) ? eyh[b] / pyt : 0;
    }
    TGraphAsymmErrors* g_ratio_syst = new TGraphAsymmErrors(nbins, x.data(), ry.data(),
                                                             exl.data(), exh.data(),
                                                             reyl.data(), reyh.data());
    g_ratio_syst->SetFillColorAlpha(kAzure - 9, 0.6);
    g_ratio_syst->SetFillStyle(1001);
    g_ratio_syst->SetLineWidth(0);

    ratio_nominal->SetStats(0);
    ratio_nominal->GetYaxis()->SetTitle("Data / Pythia8");
    ratio_nominal->GetYaxis()->SetNdivisions(510);
    ratio_nominal->GetYaxis()->SetRangeUser(0.5, 1.5);
    ratio_nominal->GetYaxis()->SetTitleSize(25);
    ratio_nominal->GetYaxis()->SetTitleFont(43);
    ratio_nominal->GetYaxis()->SetTitleOffset(1.5);
    ratio_nominal->GetYaxis()->SetLabelFont(43);
    ratio_nominal->GetYaxis()->SetLabelSize(25);
    ratio_nominal->GetXaxis()->SetTitle("p_{T,lead} [GeV]");
    ratio_nominal->GetXaxis()->SetTitleSize(25);
    ratio_nominal->GetXaxis()->SetTitleFont(43);
    ratio_nominal->GetXaxis()->SetTitleOffset(4);
    ratio_nominal->GetXaxis()->SetLabelFont(43);
    ratio_nominal->GetXaxis()->SetLabelSize(25);
    ratio_nominal->GetXaxis()->SetRangeUser(17, 63);

    ratio_nominal->Draw("E1");
    g_ratio_syst->Draw("2 same");
    ratio_nominal->Draw("E1 same");

    TLine* line0 = new TLine(17, 1.0, 63, 1.0);
    TLine* line1 = new TLine(17, 0.95, 63, 0.95);
    TLine* line2 = new TLine(17, 1.05, 63, 1.05);
    line0->SetLineStyle(1);
    line1->SetLineStyle(2);
    line2->SetLineStyle(2);
    line0->Draw("same");
    line1->Draw("same");
    line2->Draw("same");

    canvas->Update();
    canvas->Draw();
    if (output_name) canvas->SaveAs(output_name);
}

void plot_result() {
    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    std::vector<std::string> syst = {"calib_dijet_reweight_trim_10","calib_dijet_jesdown_reweight_trim_10","calib_dijet_jesup_reweight_trim_10","calib_dijet_jerdown_reweight_trim_10","calib_dijet_jerup_reweight_trim_10","calib_dijet_trim_10"};
    std::vector<std::string> syst_tags = {"Nominal", "JES Down", "JES Up", "JER Down", "JER Up", "No reweight", "Iter 3", "Iter 7"};

    std::vector<std::string> truth_syst = {"pythia","herwig"};

    TFile* fpileup = TFile::Open("run21_output_files/output_unfolded_data_pu_correct_run21_iter_3_1000toys.root");
    TFile* f = TFile::Open("run28_output_files/output_unfolded_data_dijet_bkg_cut_run28_iter_3_1000toys.root");
    TFile* fherwig = TFile::Open("run21_output_files/output_unfolded_data_herwig_calib_dijet_run21_iter_3_1000toys.root");

    std::vector<TH2D*> h_uni_truth_2D;
    std::vector<TH2D*> h_truth_2D;
    std::vector<TProfile*> h_truth_prof;
    std::vector<TH1D*> h_truth;
    std::vector<TH2D*> h_unfold_2D;
    std::vector<TProfile*> h_unfold_prof;
    std::vector<TH1D*> h_unfold;
    TH2D* h[2];

    // --- Truth histograms (uniform -> variable binning) ---
    for (int i = 0; i < 2; i++) {
        if (i == 0) {
            h_uni_truth_2D.push_back((TH2D*)f->Get("h_truth_calib_dijet"));
            h[i] = new TH2D("h_var_truth_calib_dijet","",truthnpt, truthptbins, truthnet, truthetbins);
        }
        if (i == 1) {
            h_uni_truth_2D.push_back((TH2D*)fherwig->Get("h_truth_calib_dijet"));
            h[i] = new TH2D("h_var_truth_herwig_calib_dijet","",truthnpt, truthptbins, truthnet, truthetbins);
        }
        for (int ix = 1; ix <= h_uni_truth_2D[i]->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= h_uni_truth_2D[i]->GetNbinsY(); iy++) {
                h[i]->SetBinContent(ix,iy,h_uni_truth_2D[i]->GetBinContent(ix,iy));
                h[i]->SetBinError(ix,iy,h_uni_truth_2D[i]->GetBinError(ix,iy));
            }
        }
        h_truth_2D.push_back(h[i]);
    }

    for (int i = 0; i < 2; i++) {
        h_truth_prof.push_back(h_truth_2D[i]->ProfileX(("truth_prof_"+truth_syst[i]).c_str(),1,h_truth_2D[i]->GetNbinsY()-1));
        int tnbins = h_truth_prof[i]->GetNbinsX();
        const TAxis* txaxis = h_truth_prof[i]->GetXaxis();
        std::vector<double> tedges(tnbins + 1);
        for (int b = 0; b <= tnbins; ++b) { tedges[b] = txaxis->GetBinLowEdge(b + 1); }
        tedges[tnbins] = txaxis->GetBinUpEdge(tnbins);
        TH1D* truth_hist = new TH1D(("truth_hist_"+truth_syst[i]).c_str(), "", tnbins, tedges.data());
        for (int b = 1; b <= tnbins; ++b) {
            truth_hist->SetBinContent(b, h_truth_prof[i]->GetBinContent(b));
            truth_hist->SetBinError(b, h_truth_prof[i]->GetBinError(b));
        }
        h_truth.push_back(truth_hist);
    }

    // --- Unfold histograms ---

    // All syst from nominal file f (trim_10_5): uniform -> variable binning
    for (int i = 0; i < (int)syst.size(); i++) {
        TH2D* hu = (TH2D*)f->Get(("h_unfold_"+syst[i]+"_5").c_str());
        TH2D* hv = new TH2D(("h_var_unfold_"+syst[i]+"_5").c_str(),"",truthnpt, truthptbins, truthnet, truthetbins);
        for (int ix = 1; ix <= hu->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= hu->GetNbinsY(); iy++) {
                hv->SetBinContent(ix,iy,hu->GetBinContent(ix,iy));
                hv->SetBinError(ix,iy,hu->GetBinError(ix,iy));
            }
        }
        h_unfold_2D.push_back(hv);
    }

    // f: iter 3 and iter 7 of nominal reweight_trim_10 (uniform -> variable binning)
    {
        TH2D* hu = (TH2D*)f->Get("h_unfold_calib_dijet_reweight_trim_10_3");
        TH2D* hv = new TH2D("h_var_unfold_calib_dijet_reweight_trim_10_3","",truthnpt, truthptbins, truthnet, truthetbins);
        for (int ix = 1; ix <= hu->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= hu->GetNbinsY(); iy++) {
                hv->SetBinContent(ix,iy,hu->GetBinContent(ix,iy));
                hv->SetBinError(ix,iy,hu->GetBinError(ix,iy));
            }
        }
        h_unfold_2D.push_back(hv);
    }
    {
        TH2D* hu = (TH2D*)f->Get("h_unfold_calib_dijet_reweight_trim_10_7");
        TH2D* hv = new TH2D("h_var_unfold_calib_dijet_reweight_trim_10_7","",truthnpt, truthptbins, truthnet, truthetbins);
        for (int ix = 1; ix <= hu->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= hu->GetNbinsY(); iy++) {
                hv->SetBinContent(ix,iy,hu->GetBinContent(ix,iy));
                hv->SetBinError(ix,iy,hu->GetBinError(ix,iy));
            }
        }
        h_unfold_2D.push_back(hv);
    }

    // fpileup: h_unfold_calib_dijet_reweight_trim10_4 (uniform -> variable binning)
    {
        TH2D* hu = (TH2D*)fpileup->Get("h_unfold_calib_dijet_reweight_trim_10_4");
        TH2D* hv = new TH2D("h_var_unfold_pu_trim_10_4","",truthnpt, truthptbins, truthnet, truthetbins);
        for (int ix = 1; ix <= hu->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= hu->GetNbinsY(); iy++) {
                hv->SetBinContent(ix,iy,hu->GetBinContent(ix,iy));
                hv->SetBinError(ix,iy,hu->GetBinError(ix,iy));
            }
        }
        //h_unfold_2D.push_back(hv);
    }

    // fherwig: h_unfold_calib_dijet_reweight_trim_5_8 (already variable binning)
    {
        TH2D* hv = (TH2D*)fherwig->Get("h_unfold_calib_dijet_reweight_trim_5_8");
        //h_unfold_2D.push_back(hv);
    }

    // Profile and extract 1D for all unfold histograms
    auto get_unfold_prof_name = [&](int i) -> std::string {
        if (i < (int)syst.size()) return "unfold_prof_"+syst[i]+"_5";
        if (i == (int)syst.size())     return "unfold_prof_calib_dijet_reweight_trim_10_3";
        if (i == (int)syst.size() + 1) return "unfold_prof_calib_dijet_reweight_trim_10_7";
        return "unfold_prof_unknown_" + std::to_string(i);
    };
    auto get_unfold_hist_name = [&](int i) -> std::string {
        if (i < (int)syst.size()) return "unfold_hist_"+syst[i]+"_5";
        if (i == (int)syst.size())     return "unfold_hist_calib_dijet_reweight_trim_10_3";
        if (i == (int)syst.size() + 1) return "unfold_hist_calib_dijet_reweight_trim_10_7";
        return "unfold_hist_unknown_" + std::to_string(i);
    };

    for (int i = 0; i < (int)h_unfold_2D.size(); i++) {
        h_unfold_prof.push_back(h_unfold_2D[i]->ProfileX(get_unfold_prof_name(i).c_str(),1,h_truth_2D[0]->GetNbinsY()-1));
        int nbins = h_unfold_prof[i]->GetNbinsX();
        const TAxis* xaxis = h_unfold_prof[i]->GetXaxis();
        std::vector<double> edges(nbins + 1);
        for (int b = 0; b <= nbins; ++b) { edges[b] = xaxis->GetBinLowEdge(b + 1); }
        edges[nbins] = xaxis->GetBinUpEdge(nbins);
        TH1D* hist = new TH1D(get_unfold_hist_name(i).c_str(), "", nbins, edges.data());
        for (int b = 1; b <= nbins; ++b) {
            hist->SetBinContent(b, h_unfold_prof[i]->GetBinContent(b));
            hist->SetBinError(b, h_unfold_prof[i]->GetBinError(b));
        }
        h_unfold.push_back(hist);
    }

    draw_profile(h_truth, h_unfold, true, syst_tags, (int)h_unfold.size(), "h_profile_run28_dijet_bkg_cut_figure.png");

    draw_result_with_syst(h_truth, h_unfold, "h_result_run28_dijet_bkg_cut_syst_band_figure.png");

}
