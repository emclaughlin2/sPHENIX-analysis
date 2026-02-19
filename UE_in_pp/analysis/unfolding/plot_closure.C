#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TColor.h>
#include <TROOT.h>
#include <vector>
#include <string>
#include <sstream>
#include <unfold_Def.h>

int colors[] = {kBlue, kGreen+2, kMagenta, kOrange, kCyan+1,
                    kViolet, kSpring+8, kAzure-4, kPink+8,
                    kYellow+2, kTeal+3, kGray+2, kBlue+3,
                    kOrange+7, kGreen+4, kMagenta+2, kRed+3,
                    kAzure+6, kBlack, kRed+1, kBlue+1, kGreen+1, 
                    kMagenta+1, kOrange+1, 
                    kCyan+2, kViolet+1, kPink+1, kYellow+3};

void draw_unfolded_spectra_one_unfold(TH1D* truth, TH1D* measure, TH1D* unfold, bool jet, std::string leg_tags, int max_iter, std::pair<double, double> x_range, std::string trim, const char* output_name = nullptr)
{
    TCanvas* canvas = new TCanvas("canvas", "", 600, 800);

    TPad* pad1 = new TPad("pad1", "", 0, 0.5, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    //if (jet) 
    pad1->SetLogy(1);
    //else pad1->SetLogx(1);

    truth->SetStats(0);
    measure->SetStats(0);
    truth->SetLineColor(2);
    truth->SetMarkerColor(2);
    measure->SetLineColor(1);
    measure->SetMarkerColor(1);
    truth->GetXaxis()->SetLabelSize(0);
    measure->GetXaxis()->SetLabelSize(0);

    unfold->SetStats(0);
    unfold->SetLineColor(colors[0]); // Use TColor::GetColor if needed
    unfold->SetMarkerColor(colors[0]);
    unfold->GetXaxis()->SetLabelSize(0);

    vector<double> bins;
    vector<double> bin_size;
    for (int i = 1; i <= truth->GetNbinsX() + 1; i++) {
        bins.push_back(truth->GetBinLowEdge(i));
    }
    for (int i = 1; i < bins.size(); i++) {
        bin_size.push_back(bins[i] - bins[i-1]);
    }

    for (int i = 1; i < truth->GetNbinsX() + 1; i++) {
        truth->SetBinContent(i,truth->GetBinContent(i)/bin_size[i-1]);
        truth->SetBinError(i,truth->GetBinError(i)/bin_size[i-1]);
    }

    bins.clear();
    bin_size.clear();
    for (int i = 1; i <= measure->GetNbinsX() + 1; i++) {
        bins.push_back(measure->GetBinLowEdge(i));
    }
    for (int i = 1; i < bins.size(); i++) {
        bin_size.push_back(bins[i] - bins[i-1]);
    }

    for (int i = 1; i < measure->GetNbinsX() + 1; i++) {
        measure->SetBinContent(i,measure->GetBinContent(i)/bin_size[i-1]);
        measure->SetBinError(i,measure->GetBinError(i)/bin_size[i-1]);
    }

    bins.clear();
    bin_size.clear();
    for (int i = 1; i <= unfold->GetNbinsX() + 1; i++) {
        bins.push_back(unfold->GetBinLowEdge(i));
    }
    for (int i = 1; i < bins.size(); i++) {
        bin_size.push_back(bins[i] - bins[i-1]);
    }

    for (int i = 1; i < unfold->GetNbinsX() + 1; i++) {
        unfold->SetBinContent(i,unfold->GetBinContent(i)/bin_size[i-1]);
        unfold->SetBinError(i,unfold->GetBinError(i)/bin_size[i-1]);
    }

    truth->Scale(1.0/truth->Integral());
    measure->Scale(1.0/measure->Integral());
    unfold->Scale(1.0/unfold->Integral());

    truth->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
    measure->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
    unfold->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
    if (jet) truth->GetYaxis()->SetRangeUser(0.0000001,10);
    if (!jet) truth->GetYaxis()->SetRangeUser(0.0005,0.7);

    truth->Draw();
    measure->Draw("same");
    unfold->Draw("same");

    TLegend* leg = new TLegend(0.55, 0.65, 0.9, 0.9);
    leg->AddEntry("","#bf{#it{sPHENIX}} Simulation Internal","");
    leg->AddEntry("",trim.c_str(),"");
    leg->AddEntry("","Pythia8 200 GeV p+p","");
    leg->AddEntry(truth, "Truth", "lp");
    leg->AddEntry(measure, "Measured", "lp");
    leg->AddEntry(unfold, leg_tags.c_str(), "lp");
    leg->SetTextSize(0.04);
    leg->Draw();

    canvas->cd();
    TPad* pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.5);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();
    //if (!jet) pad2->SetLogx(1);

    TH1D* ratio = (TH1D*)unfold->Clone("ratio");
    ratio->Divide(truth);
    ratio->GetYaxis()->SetTitle("Unfolded/Truth Ratio");
    ratio->GetYaxis()->SetNdivisions(510);
    ratio->GetYaxis()->SetRangeUser(0.9, 1.1);
    ratio->GetYaxis()->SetTitleSize(25);
    ratio->GetYaxis()->SetTitleFont(43);
    ratio->GetYaxis()->SetTitleOffset(1.5);
    ratio->GetYaxis()->SetLabelFont(43);
    ratio->GetYaxis()->SetLabelSize(25);
    ratio->GetXaxis()->SetTitle(jet ? "p_{T} [GeV]" : "#SigmaE_{T} [GeV]");
    ratio->GetXaxis()->SetTitleSize(25);
    ratio->GetXaxis()->SetTitleFont(43);
    ratio->GetXaxis()->SetTitleOffset(0);
    ratio->GetXaxis()->SetLabelFont(43);
    ratio->GetXaxis()->SetLabelSize(25);
    ratio->GetXaxis()->SetRangeUser(x_range.first, x_range.second);

    ratio->Draw("ep");

    TLine* line0 = new TLine(x_range.first, 1.0, x_range.second, 1.0);
    TLine* line1 = new TLine(x_range.first, 0.95, x_range.second, 0.95);
    TLine* line2 = new TLine(x_range.first, 1.05, x_range.second, 1.05);
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

void draw_unfolded_spectra(TH1D* truth, TH1D* measure, std::vector<TH1D*> unfold, bool jet, std::vector<std::string> leg_tags, int max_iter, std::pair<double, double> x_range, std::string trim, const char* output_name = nullptr)
{
    TCanvas* canvas = new TCanvas("canvas", "", 600, 800);

    TPad* pad1 = new TPad("pad1", "", 0, 0.5, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    //if (jet) 
    pad1->SetLogy(1);
    //else pad1->SetLogx(1);

    truth->SetStats(0);
    measure->SetStats(0);
    truth->SetLineColor(2);
    truth->SetMarkerColor(2);
    measure->SetLineColor(1);
    measure->SetMarkerColor(1);
    truth->GetXaxis()->SetLabelSize(0);
    measure->GetXaxis()->SetLabelSize(0);

    vector<double> bins;
    vector<double> bin_size;
    for (int i = 1; i <= truth->GetNbinsX() + 1; i++) {
        bins.push_back(truth->GetBinLowEdge(i));
    }
    for (int i = 1; i < bins.size(); i++) {
        bin_size.push_back(bins[i] - bins[i-1]);
    }

    for (int i = 1; i < truth->GetNbinsX() + 1; i++) {
        truth->SetBinContent(i,truth->GetBinContent(i)/bin_size[i-1]);
        truth->SetBinError(i,truth->GetBinError(i)/bin_size[i-1]);
    }

    bins.clear();
    bin_size.clear();
    for (int i = 1; i <= measure->GetNbinsX() + 1; i++) {
        bins.push_back(measure->GetBinLowEdge(i));
    }
    for (int i = 1; i < bins.size(); i++) {
        bin_size.push_back(bins[i] - bins[i-1]);
    }

    for (int i = 1; i < measure->GetNbinsX() + 1; i++) {
        measure->SetBinContent(i,measure->GetBinContent(i)/bin_size[i-1]);
        measure->SetBinError(i,measure->GetBinError(i)/bin_size[i-1]);
    }

    bins.clear();
    bin_size.clear();
    for (int i = 1; i <= unfold[0]->GetNbinsX() + 1; i++) {
        bins.push_back(unfold[0]->GetBinLowEdge(i));
    }
    for (int i = 1; i < bins.size(); i++) {
        bin_size.push_back(bins[i] - bins[i-1]);
    }

    for (int u = 0; u < max_iter; u++) {
        for (int i = 1; i < unfold[u]->GetNbinsX() + 1; i++) {
            unfold[u]->SetBinContent(i,unfold[u]->GetBinContent(i)/bin_size[i-1]);
            unfold[u]->SetBinError(i,unfold[u]->GetBinError(i)/bin_size[i-1]);
        }
    }

    truth->Scale(1.0/truth->Integral());
    measure->Scale(1.0/measure->Integral());
    for (int i = 0; i < max_iter; ++i) {
        unfold[i]->Scale(1.0/unfold[i]->Integral());
    }

    for (int i = 0; i < max_iter; ++i) {
        unfold[i]->SetStats(0);
        unfold[i]->SetLineColor(colors[i]); // Use TColor::GetColor if needed
        unfold[i]->SetMarkerColor(colors[i]);
        unfold[i]->GetXaxis()->SetLabelSize(0);
    }

    truth->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
    measure->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
    for (int i = 0; i < max_iter; i++) {
        unfold[i]->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
    }
    if (jet) truth->GetYaxis()->SetRangeUser(0.0000001,10);
    if (!jet) truth->GetYaxis()->SetRangeUser(0.0005,0.7);

    truth->Draw();
    measure->Draw("same");
    for (int i = 0; i < max_iter; ++i) {
        unfold[i]->Draw("same");
    }

    TLegend* leg = new TLegend(0.55, 0.5, 0.9, 0.9);
    leg->AddEntry("","#bf{#it{sPHENIX}} Simulation Internal","");
    leg->AddEntry("",trim.c_str(),"");
    leg->AddEntry("","Pythia8 200 GeV p+p","");
    leg->AddEntry(truth, "Truth", "lp");
    leg->AddEntry(measure, "Measured", "lp");
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
    //if (!jet) pad2->SetLogx(1);

    std::vector<TH1D*> ratios;
    for (int i = 0; i < max_iter; ++i) {
        TH1D* ratio = (TH1D*)unfold[i]->Clone(Form("ratio%d", i));
        ratio->Divide(truth);
        ratios.push_back(ratio);
    }

    TH1D* r0 = ratios[0];
    r0->GetYaxis()->SetTitle("Unfolded/Truth Ratio");
    r0->GetYaxis()->SetNdivisions(510);
    r0->GetYaxis()->SetRangeUser(0.9, 1.1);
    r0->GetYaxis()->SetTitleSize(25);
    r0->GetYaxis()->SetTitleFont(43);
    r0->GetYaxis()->SetTitleOffset(1.5);
    r0->GetYaxis()->SetLabelFont(43);
    r0->GetYaxis()->SetLabelSize(25);
    r0->GetXaxis()->SetTitle(jet ? "p_{T} [GeV]" : "#SigmaE_{T} [GeV]");
    r0->GetXaxis()->SetTitleSize(25);
    r0->GetXaxis()->SetTitleFont(43);
    r0->GetXaxis()->SetTitleOffset(0);
    r0->GetXaxis()->SetLabelFont(43);
    r0->GetXaxis()->SetLabelSize(25);
    r0->GetXaxis()->SetRangeUser(x_range.first, x_range.second);

    r0->Draw("ep");
    for (int i = 1; i < max_iter; ++i) {
        ratios[i]->Draw("ep same");
    }

    TLine* line0 = new TLine(x_range.first, 1.0, x_range.second, 1.0);
    TLine* line1 = new TLine(x_range.first, 0.95, x_range.second, 0.95);
    TLine* line2 = new TLine(x_range.first, 1.05, x_range.second, 1.05);
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

void plot_closure(const char* closurefile = "output_closure_sim_run28_iter_2_1000toys.root") {
    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    std::vector<std::string> full_leg_tags = { "Unfold, iter 1"};
    //std::vector<std::string> half_leg_tags = { "Unfold, iter 1", "Unfold, iter 2", "Unfold, iter 3", "Unfold, iter 4", "Unfold, iter 5", "Unfold, iter 6", 
    //"Unfold, iter 7", "Unfold, iter 8", "Unfold, iter 9", "Unfold, iter 10"};
    std::vector<std::string> half_leg_tags = { "Unfold, iter 1", "Unfold, iter 2", "Unfold, iter 3", "Unfold, iter 4", "Unfold, iter 5", "Unfold, iter 6"};
    std::vector<std::string> syst = {"calib_dijet","calib_dijet_jesdown","calib_dijet_jesup","calib_dijet_jerdown","calib_dijet_jerup","calib_dijet_half1","calib_dijet_half2"};;
    std::vector<std::string> syst1D = {"jetpt_respmatrix","caloet_respmatrix"};
    std::vector<std::string> trim = {"","_trim_5","_trim_10","_reweight","_reweight_trim_5","_reweight_trim_10"};

    TFile *f = new TFile(closurefile, "READ");
    string bkg_cut;
    if (strstr(closurefile, "dijet")) { bkg_cut = "_dijet"; }
    if (strstr(closurefile, "none")) { bkg_cut = "_none"; }

    /*
    TH2D* h_truth_2D[7][6];
    TH2D* h_measure_2D[7][6];
    TH2D* h_full_unfold_2D[7][6];
    TH1D* hj_truth[7][6];
    TH1D* hc_truth[7][6];
    TH1D* hj_measure[7][6];
    TH1D* hc_measure[7][6];
    TH1D* hj_full_unfold[7][6];
    TH1D* hc_full_unfold[7][6];

    TH2D* h_half_unfold_2D[6][10];
    TH1D* hj_half_unfold[6][10];
    TH1D* hc_half_unfold[6][10];
    */

    std::vector<std::vector<TH2D*>> h_uni_truth_2D(7, std::vector<TH2D*>(6, nullptr));
    std::vector<std::vector<TH2D*>> h_uni_measure_2D(7, std::vector<TH2D*>(6, nullptr));
    std::vector<std::vector<TH2D*>> h_uni_full_unfold_2D(7, std::vector<TH2D*>(6, nullptr));

    std::vector<std::vector<TH2D*>> h_truth_2D(7, std::vector<TH2D*>(6, nullptr));
    std::vector<std::vector<TH2D*>> h_measure_2D(7, std::vector<TH2D*>(6, nullptr));
    std::vector<std::vector<TH2D*>> h_full_unfold_2D(7, std::vector<TH2D*>(6, nullptr));

    std::vector<std::vector<TH1D*>> hj_truth(7, std::vector<TH1D*>(6, nullptr));
    std::vector<std::vector<TH1D*>> hc_truth(7, std::vector<TH1D*>(6, nullptr));
    std::vector<std::vector<TH1D*>> hj_measure(7, std::vector<TH1D*>(6, nullptr));
    std::vector<std::vector<TH1D*>> hc_measure(7, std::vector<TH1D*>(6, nullptr));
    std::vector<std::vector<TH1D*>> hj_full_unfold(7, std::vector<TH1D*>(6, nullptr));
    std::vector<std::vector<TH1D*>> hc_full_unfold(7, std::vector<TH1D*>(6, nullptr));

    std::vector<std::vector<TH2D*>> h_uni_half_unfold_2D(6, std::vector<TH2D*>(10, nullptr));
    std::vector<std::vector<TH2D*>> h_half_unfold_2D(6, std::vector<TH2D*>(10, nullptr));
    std::vector<std::vector<TH1D*>> hj_half_unfold(6, std::vector<TH1D*>(10, nullptr));
    std::vector<std::vector<TH1D*>> hc_half_unfold(6, std::vector<TH1D*>(10, nullptr));

    std::vector<std::vector<TH1D*>> h_truth_1D(2, std::vector<TH1D*>(6, nullptr));
    std::vector<std::vector<TH1D*>> h_measure_1D(2, std::vector<TH1D*>(6, nullptr));
    std::vector<std::vector<TH1D*>> h_full_unfold_1D(2, std::vector<TH1D*>(6, nullptr));

    
    for (int i = 0; i < syst.size(); i++) {
        for (int j = 0; j < trim.size(); j++) {
            h_uni_truth_2D[i][j] = (TH2D*)f->Get(("h_truth_"+syst[i]+trim[j]).c_str());
            h_uni_measure_2D[i][j] = (TH2D*)f->Get(("h_measure_"+syst[i]+trim[j]).c_str());
            h_truth_2D[i][j] = new TH2D(("h_var_truth_"+syst[i]+trim[j]+bkg_cut).c_str(),"",truthnpt, truthptbins, truthnet, truthetbins);
            h_measure_2D[i][j] = new TH2D(("h_var_measure_"+syst[i]+trim[j]+bkg_cut).c_str(),"",calibnpt, calibptbins, calibnet, calibetbins);
            for (int ix = 1; ix < h_uni_truth_2D[i][j]->GetNbinsX() + 1; ix++) {
                for (int iy = 1; iy < h_uni_truth_2D[i][j]->GetNbinsY() + 1; iy++) {
                    h_truth_2D[i][j]->SetBinContent(ix,iy,h_uni_truth_2D[i][j]->GetBinContent(ix,iy));
                    h_truth_2D[i][j]->SetBinError(ix,iy,h_uni_truth_2D[i][j]->GetBinError(ix,iy));
                }
            }
            for (int ix = 1; ix < h_uni_measure_2D[i][j]->GetNbinsX() + 1; ix++) {
                for (int iy = 1; iy < h_uni_measure_2D[i][j]->GetNbinsY() + 1; iy++) {
                    h_measure_2D[i][j]->SetBinContent(ix,iy,h_uni_measure_2D[i][j]->GetBinContent(ix,iy));
                    h_measure_2D[i][j]->SetBinError(ix,iy,h_uni_measure_2D[i][j]->GetBinError(ix,iy));
                }
            }
        }
    }
    for (int i = 0; i < syst.size(); i++) {
        for (int j = 0; j < trim.size(); j++) {
            hj_truth[i][j] = h_truth_2D[i][j]->ProjectionX(("hj_truth_"+syst[i]+trim[j]+bkg_cut).c_str());
            hc_truth[i][j] = h_truth_2D[i][j]->ProjectionY(("hc_truth_"+syst[i]+trim[j]+bkg_cut).c_str());
            hj_measure[i][j] = h_measure_2D[i][j]->ProjectionX(("hj_measure_"+syst[i]+trim[j]+bkg_cut).c_str());
            hc_measure[i][j] = h_measure_2D[i][j]->ProjectionY(("hc_measure_"+syst[i]+trim[j]+bkg_cut).c_str());
        }
    }
    for (int i = 0; i < syst.size(); i++) {
        for (int j = 0; j < trim.size(); j++) {
            h_uni_full_unfold_2D[i][j] = (TH2D*)f->Get(("h_unfold_full_"+syst[i]+trim[j]+"_1").c_str());
            h_full_unfold_2D[i][j] = new TH2D(("h_var_unfold_full_"+syst[i]+trim[j]+bkg_cut).c_str(),"",truthnpt, truthptbins, truthnet, truthetbins);
            for (int ix = 1; ix < h_uni_full_unfold_2D[i][j]->GetNbinsX() + 1; ix++) {
                for (int iy = 1; iy < h_uni_full_unfold_2D[i][j]->GetNbinsY() + 1; iy++) {
                    h_full_unfold_2D[i][j]->SetBinContent(ix,iy,h_uni_full_unfold_2D[i][j]->GetBinContent(ix,iy));
                    h_full_unfold_2D[i][j]->SetBinError(ix,iy,h_uni_full_unfold_2D[i][j]->GetBinError(ix,iy));
                }
            }
        }
    }
    for (int i = 0; i < syst.size(); i++) {
        for (int j = 0; j < trim.size(); j++) {
            hj_full_unfold[i][j] = h_full_unfold_2D[i][j]->ProjectionX(("hj_full_unfold_"+syst[i]+trim[j]+bkg_cut).c_str());
            hc_full_unfold[i][j] = h_full_unfold_2D[i][j]->ProjectionY(("hc_full_unfold_"+syst[i]+trim[j]+bkg_cut).c_str());
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < trim.size(); j++) {
            h_truth_1D[i][j] = (TH1D*)f->Get(("h_truth_"+syst1D[i]+trim[j]).c_str());
            h_measure_1D[i][j] = (TH1D*)f->Get(("h_measure_"+syst1D[i]+trim[j]).c_str());
            h_full_unfold_1D[i][j] = (TH1D*)f->Get(("h_unfold_full_"+syst1D[i]+trim[j]+"_1").c_str());
        }
    }

    for (int i = 0; i < trim.size(); i++) {
        for (int j = 0; j < 6; j++) {
            h_uni_half_unfold_2D[i][j] = (TH2D*)f->Get(("h_unfold_half"+trim[i]+"_"+to_string(j+1)).c_str());
            h_half_unfold_2D[i][j] = new TH2D(("h_var_unfold_half_"+trim[i]+bkg_cut+"_"+to_string(j+1)).c_str(),"",truthnpt, truthptbins, truthnet, truthetbins);
            for (int ix = 1; ix < h_uni_half_unfold_2D[i][j]->GetNbinsX() + 1; ix++) {
                for (int iy = 1; iy < h_uni_half_unfold_2D[i][j]->GetNbinsY() + 1; iy++) {
                    h_half_unfold_2D[i][j]->SetBinContent(ix,iy,h_uni_half_unfold_2D[i][j]->GetBinContent(ix,iy));
                    h_half_unfold_2D[i][j]->SetBinError(ix,iy,h_uni_half_unfold_2D[i][j]->GetBinError(ix,iy));
                }
            }
        }
    }
    for (int i = 0; i < trim.size(); i++) {
        for (int j = 0; j < 6; j++) {
            hj_half_unfold[i][j] = h_half_unfold_2D[i][j]->ProjectionX(("hj_half_unfold"+trim[i]+bkg_cut+to_string(j)).c_str());
            hc_half_unfold[i][j] = h_half_unfold_2D[i][j]->ProjectionY(("hc_half_unfold"+trim[i]+bkg_cut+to_string(j)).c_str());
        }
    }
    /*
    for (int i = 0; i < syst.size(); i++) {
        for (int j = 0; j < trim.size(); j++) {
            string jet_outfile = "sim_plots_run28/h_jet_spectrum_" + syst[i] + trim[j]+bkg_cut + "_full_closure_iter_2_1000toys_extended_reco_bin.png";
            string et_outfile = "sim_plots_run28/h_et_spectrum_" + syst[i] + trim[j]+bkg_cut + "_full_closure_iter_2_1000toys_extended_reco_bin.png";
            draw_unfolded_spectra_one_unfold(hj_truth[i][j], hj_measure[i][j], hj_full_unfold[i][j], true, full_leg_tags[0], 1, std::make_pair(17.0, 82.0), trim[i]+bkg_cut, jet_outfile.c_str());
            draw_unfolded_spectra_one_unfold(hc_truth[i][j], hc_measure[i][j], hc_full_unfold[i][j], false, full_leg_tags[0], 1, std::make_pair(0, 15.0), trim[i]+bkg_cut, et_outfile.c_str());
        }
    }
    */
    for (int i = 0; i < trim.size(); i++) {
        string jet_outfile = "sim_plots_run28/h_jet_spectrum" + trim[i]+bkg_cut + "_half_closure_iter_2_1000toys_extended_reco_bin_feb.png";
        string et_outfile = "sim_plots_run28/h_et_spectrum" + trim[i]+bkg_cut + "_half_closure_iter_2_1000toys_extended_reco_bin_feb.png";
        draw_unfolded_spectra(hj_truth[5][i], hj_measure[6][i], hj_half_unfold[i], true, half_leg_tags, 6, std::make_pair(17.0, 82.0), trim[i]+bkg_cut, jet_outfile.c_str());
        draw_unfolded_spectra(hc_truth[5][i], hc_measure[6][i], hc_half_unfold[i], false, half_leg_tags, 6, std::make_pair(0, 15.0), trim[i]+bkg_cut, et_outfile.c_str());
    }
    
    /*
    for (int i = 0; i < 6; i++) {
        string jet_outfile = "sim_plots_run21/h_jet_spectrum_1D" + trim[i]+bkg_cut + "_full_closure_iter_3_1000toys.png";
        string et_outfile = "sim_plots_run21/h_et_spectrum_1D" + trim[i]+bkg_cut + "_full_closure_iter_3_1000toys.png";
        draw_unfolded_spectra_one_unfold(h_truth_1D[0][i], h_measure_1D[0][i], h_full_unfold_1D[0][i], true, full_leg_tags[0], 1, std::make_pair(0.0,1.0), jet_outfile.c_str());
        draw_unfolded_spectra_one_unfold(h_truth_1D[1][i], h_measure_1D[1][i], h_full_unfold_1D[1][i], false, full_leg_tags[0], 1, std::make_pair(0.0,1.0), et_outfile.c_str());
    }
    */
    f->Close();
}
