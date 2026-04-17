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
#include <iomanip>

int rgb[25][3] = {{0,0,0}, {230, 25, 75}, {60, 180, 75}, {255, 225, 25}, {0, 130, 200},
                 {245, 130, 48}, {145, 30, 180}, {70, 240, 240}, {240, 50, 230}, {210, 245, 60},
                {250, 190, 212}, {0, 128, 128}, {220, 190, 255}, {170, 110, 40}, {128, 128, 128}, 
                {128, 0, 0}, {128, 128, 0}, {255, 215, 180}, {0, 0, 128}, {34, 139, 34}, 
                {59, 213, 224},{173, 117, 219},{209, 62, 109},{108, 137, 204},{189, 10, 67}};
int colors[25];
for (int i = 0; i < 25; i++) colors[i] = TColor::GetColor(rgb[i][0],rgb[i][1],rgb[i][2]);

void draw_profile(std::vector<TH1D*> truth, std::vector<TH1D*> unfold, bool jet, std::vector<std::string> leg_tags, int max_iter, int dijet, const char* output_name = nullptr)
{
    TCanvas* canvas = new TCanvas("canvas", "", 600, 800);

    //std::vector<std::string> colors = {"kBlack","kRed","kBlue","kGreen+2","kBlue-3","kMagneta+1"};

    TPad* pad1 = new TPad("pad1", "", 0, 0.5, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();

    for (int i = max_iter - 1; i >= 0; --i) {
        unfold[i]->SetStats(0); // Use TColor::GetColor if needed
        unfold[i]->SetLineColor(colors[i]); // Use TColor::GetColor if needed
        unfold[i]->SetMarkerColor(colors[i]);
        unfold[i]->GetXaxis()->SetLabelSize(0);
        unfold[i]->GetXaxis()->SetRangeUser(21, 63);
        unfold[i]->GetYaxis()->SetRangeUser(0.0, 1.2);
        unfold[i]->GetYaxis()->SetTitle("<#SigmaE_{T}/#delta#eta#delta#phi> [GeV]");
        unfold[i]->GetXaxis()->SetTitle("p_{T,lead} [GeV]");
    }

    unfold[1]->Draw();
    for (int i = 2; i < max_iter; ++i) {
        unfold[i]->Draw("same");
    } 
    unfold[0]->Draw("same");

    TLegend* leg = new TLegend(0.17, 0.5, 0.92, 0.92);
    leg->AddEntry("","#bf{#it{sPHENIX}} Internal","");
    leg->AddEntry("","200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
    if (dijet) { leg->AddEntry("","Exclusive dijet",""); }
    else { leg->AddEntry("","Inclusive jet",""); }
    leg->SetNColumns(2);
    for (int i = 0; i < max_iter; ++i) {
        leg->AddEntry(unfold[i], leg_tags[i].c_str(), "lp");
    }
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
    ratios[0]->GetXaxis()->SetTitle("p_{T,lead} [GeV]");
    ratios[0]->GetXaxis()->SetTitleSize(25);
    ratios[0]->GetXaxis()->SetTitleFont(43);
    ratios[0]->GetXaxis()->SetTitleOffset(0);
    ratios[0]->GetXaxis()->SetLabelFont(43);
    ratios[0]->GetXaxis()->SetLabelSize(25);
    ratios[0]->GetXaxis()->SetRangeUser(21, 63);

    std::cout << ratios.size() << std::endl;
    for (int i = 0; i < ratios.size(); ++i) {
        ratios[i]->GetXaxis()->SetRangeUser(21, 63);
        ratios[i]->GetYaxis()->SetRangeUser(0.5,1.5);
        std::cout << i << " " << leg_tags[i+1] << " ";
        float avgbin = 0;
        int nbins = 0;
        for (int j = 2; j < ratios[i]->GetNbinsX(); j++) {
            avgbin += ratios[i]->GetBinContent(j) - 1.0; 
            nbins++;
            //std::cout << ratios[i]->GetBinContent(j) << " ";
        }
        std::cout << std::setprecision(3) << (avgbin*100.0)/nbins << "% " << std::endl;
        if (i == 0) ratios[i]->Draw("hist");
        else ratios[i]->Draw("hist, same");
    }

    TLine* line0 = new TLine(21, 1.0, 63, 1.0);
    TLine* line1 = new TLine(21, 0.95, 63, 0.95);
    TLine* line2 = new TLine(21, 1.05, 63, 1.05);
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
    truth[0]->SetLineColor(max_iter + 1);
    truth[0]->SetMarkerColor(max_iter + 1);
    truth[1]->SetStats(0);
    truth[1]->SetLineColor(max_iter + 2);
    truth[1]->SetMarkerColor(max_iter + 2);
    truth[0]->GetXaxis()->SetRangeUser(21, 63);
    truth[1]->GetXaxis()->SetRangeUser(21, 63);

    for (int i = 0; i < max_iter; ++i) {
        unfold[i]->SetStats(0);
        unfold[i]->SetLineColor(colors[i]); // Use TColor::GetColor if needed
        unfold[i]->SetMarkerColor(colors[i]);
        unfold[i]->GetXaxis()->SetRangeUser(21, 63);
        unfold[i]->GetYaxis()->SetRangeUser(0.0, 1.2);
    }


    unfold[0]->GetYaxis()->SetTitle("<#SigmaE_{T}/#delta#eta#delta#phi> [GeV]");
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
void draw_result_with_syst(std::vector<TH1D*> truth, std::vector<TH1D*> unfold, int dijet, const char* output_name = nullptr)
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

    g_syst->Print();

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
    unfold[0]->GetXaxis()->SetRangeUser(21, 63);
    unfold[0]->GetYaxis()->SetRangeUser(0.0, 1.2);
    unfold[0]->GetYaxis()->SetTitle("<#SigmaE_{T}/#delta#eta#delta#phi> [GeV]");

    truth[0]->SetStats(0);
    truth[0]->SetLineColor(kRed + 1);
    truth[0]->SetLineWidth(2);
    truth[0]->SetMarkerColor(kRed + 1);
    truth[0]->GetXaxis()->SetRangeUser(21, 63);

    unfold[0]->Draw("E1");
    g_syst->Draw("2 same");
    unfold[0]->Draw("E1 same");
    truth[0]->Draw("hist same");

    TLegend* leg = new TLegend(0.17, 0.65, 0.78, 0.92);
    leg->AddEntry("","#bf{#it{sPHENIX}} Internal","");
    leg->AddEntry("","200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
    if (dijet) { leg->AddEntry("","Exclusive dijet",""); }
    else { leg->AddEntry("","Inclusive jet",""); }
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
    //ratio_nominal->GetXaxis()->SetTitleOffset(4);
    ratio_nominal->GetXaxis()->SetLabelFont(43);
    ratio_nominal->GetXaxis()->SetLabelSize(25);
    ratio_nominal->GetXaxis()->SetRangeUser(21, 63);

    ratio_nominal->Draw("E1");
    g_ratio_syst->Draw("2 same");
    ratio_nominal->Draw("E1 same");

    TLine* line0 = new TLine(21, 1.0, 63, 1.0);
    TLine* line1 = new TLine(21, 0.95, 63, 0.95);
    TLine* line2 = new TLine(21, 1.05, 63, 1.05);
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

void plot_result(int dijet = 1) {
    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    std::vector<std::string> syst = {
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
    "calib_dijet_reweight_trim_10_2_var_down_etEffCorrected" 
    };
    std::vector<std::string> syst_tags = {"Nominal", "JES Down", "JES Up", "JER Down", "JER Up", "No reweight", "Jet Trigger Eff down", "Jet Trigger Eff up", "PU corr", "Timing eff",
     "EMCal scale up", "EMCal scale down", "IHCal scale up", "IHCal scale down", "OHCal scale up", "OHCal scale down", "Had resp up", "Had resp down", "Cluster Phi Res",
     "Cluster E res","Noise down","Noise up","OHCal MC/Data var","MBD eff up", "MBD eff down", "Herwig"}; 

    std::vector<std::string> truth_syst = {"pythia","herwig"};

    TFile* f = nullptr;
    if (dijet) {
        f = TFile::Open("sphenix_primary_run28_output_files/output_mbd_correct_unfolded_data_8calibetbin_dijet_bkg_cut_run28_iter_3_1000toys.root");
    } else {
        f = TFile::Open("sphenix_primary_run28_output_files/output_mbd_correct_unfolded_data_8calibetbin_efrac_bkg_cut_run28_iter_3_1000toys.root");
    }
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
            h_uni_truth_2D.push_back((TH2D*)f->Get("h_truth_calib_dijet"));
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
        std::cout << i << std::endl;
        TH2D* hu = (TH2D*)f->Get(("h_unfold_"+syst[i]).c_str());
        TH2D* hv = new TH2D(("h_var_unfold_"+syst[i]+"_2").c_str(),"",truthnpt, truthptbins, truthnet, truthetbins);
        for (int ix = 1; ix <= hu->GetNbinsX(); ix++) {
            for (int iy = 1; iy <= hu->GetNbinsY(); iy++) {
                hv->SetBinContent(ix,iy,hu->GetBinContent(ix,iy));
                hv->SetBinError(ix,iy,hu->GetBinError(ix,iy));
            }
        }
        h_unfold_2D.push_back(hv);
    }

    // fherwig: h_unfold_calib_dijet_reweight_trim_5_8 (already variable binning)
    {
        TH2D* hv = (TH2D*)fherwig->Get("h_unfold_calib_dijet_reweight_trim_5_8");
        //h_unfold_2D.push_back(hv);
    }

    // Profile and extract 1D for all unfold histograms
    auto get_unfold_prof_name = [&](int i) -> std::string {
        if (i < (int)syst.size()) return "unfold_prof_"+syst[i]+"_5";
        return "unfold_prof_herwig_trim_5_8";
    };
    auto get_unfold_hist_name = [&](int i) -> std::string {
        if (i < (int)syst.size()) return "unfold_hist_"+syst[i]+"_5";
        return "unfold_hist_herwig_trim_5_8";
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
    for (int i = 0; i < (int)h_unfold.size(); i++) {
        h_unfold[i]->Scale(3.0 / (2.2*2*M_PI));
    }
    for (int i = 0; i < (int)h_truth.size(); i++) {
        h_truth[i]->Scale(3.0 / (2.2*2*M_PI));
    }

    if (dijet) {
        draw_profile(h_truth, h_unfold, true, syst_tags, (int)h_unfold.size(), dijet, "sphenix_primary_run28_output_files/h_profile_run28_w_calo_syst_8calibetbin_dijet_bkg_cut_figure.png");
        draw_result_with_syst(h_truth, h_unfold, dijet, "sphenix_primary_run28_output_files/h_result_run28_w_calo_syst_8calibetbin_dijet_bkg_cut_syst_band_figure.png");
    } else {
        draw_profile(h_truth, h_unfold, true, syst_tags, (int)h_unfold.size(), dijet, "sphenix_primary_run28_output_files/h_profile_run28_w_calo_syst_8calibetbin_efrac_bkg_cut_figure.png");
        draw_result_with_syst(h_truth, h_unfold, dijet, "sphenix_primary_run28_output_files/h_result_run28_w_calo_syst_8calibetbin_efrac_bkg_cut_syst_band_figure.png");
    }
}
