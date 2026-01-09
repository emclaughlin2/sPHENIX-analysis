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

void draw_unfolded_spectra(TH1D* truth, TH1D* measure, std::vector<TH1D*> unfold, bool jet, std::vector<std::string> leg_tags, int max_iter, std::pair<double, double> x_range, const char* output_name = nullptr)
{
    TCanvas* canvas = new TCanvas("canvas", "", 600, 800);

    TPad* pad1 = new TPad("pad1", "", 0, 0.5, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    if (jet) pad1->SetLogy(1);
    else pad1->SetLogx(1);

    truth->SetStats(0);
    measure->SetStats(0);
    truth->SetLineColor(2);
    truth->SetMarkerColor(2);
    measure->SetLineColor(1);
    measure->SetMarkerColor(1);
    truth->GetXaxis()->SetLabelSize(0);
    measure->GetXaxis()->SetLabelSize(0);
    truth->GetXaxis()->SetRangeUser(x_range.first, x_range.second);

    for (int i = 0; i < max_iter; ++i) {
        unfold[i]->SetStats(0);
        unfold[i]->SetLineColor(i + 3); // Use TColor::GetColor if needed
        unfold[i]->SetMarkerColor(i + 3);
        unfold[i]->GetXaxis()->SetLabelSize(0);
    }

    truth->Draw();
    measure->Draw("same");
    for (int i = 0; i < max_iter; ++i) {
        unfold[i]->Draw("same");
    }

    TLegend* leg = new TLegend(0.77, 0.65, 0.9, 0.9);
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
    if (!jet) pad2->SetLogx(1);

    std::vector<TH1D*> ratios;
    for (int i = 0; i < max_iter; ++i) {
        TH1D* ratio = (TH1D*)unfold[i]->Clone(Form("ratio%d", i));
        ratio->Divide(truth);
        ratios.push_back(ratio);
    }

    TH1D* r0 = ratios[0];
    r0->GetYaxis()->SetTitle("Unfolded/Truth Ratio");
    r0->GetYaxis()->SetNdivisions(510);
    r0->GetYaxis()->SetRangeUser(0.5, 1.5);
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

void closure_plots() {
    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    std::vector<std::string> full_leg_tags = { "Unfold, iter 1", "Unfold, iter 2", "Unfold, iter 1", "Unfold, iter 2" };

    std::vector<std::string> half_leg_tags = { "Unfold, iter 1", "Unfold, iter 2", "Unfold, iter 3", "Unfold, iter 4", "Unfold, iter 5", "Unfold, iter 6", 
    "Unfold, iter 7", "Unfold, iter 8", "Unfold, iter 9", "Unfold, iter 10", "Unfold, iter 11", "Unfold, iter 12", "Unfold, iter 13", 
    "Unfold, iter 14", "Unfold, iter 15", "Unfold, iter 16", "Unfold, iter 17", "Unfold, iter 18", "Unfold, iter 19", "Unfold, iter 20"
    };

    TFile* f = TFile::Open("output_closure_sim_4bin_newetbin_dijet_1000toys.root");

    TH2D* h_full_truth_2D = (TH2D*)f->Get("h_truth_calib_dijet");
    TH2D* h_full_measure_2D = (TH2D*)f->Get("h_measure_calib_dijet");
    TH2D* h_half_truth_2D = (TH2D*)f->Get("h_truth_calib_dijet_half1");
    TH2D* h_half_measure_2D = (TH2D*)f->Get("h_measure_calib_dijet_half2");
    TH1D* hj_full_truth = h_full_truth_2D->ProjectionX("hj_full_truth");
    TH1D* hj_full_measure = h_full_measure_2D->ProjectionX("hj_full_measure");
    TH1D* hc_full_truth = h_full_truth_2D->ProjectionY("hc_full_truth");
    TH1D* hc_full_measure = h_full_measure_2D->ProjectionY("hc_full_measure");
    TH1D* hj_half_truth = h_half_truth_2D->ProjectionX("hj_half_truth");
    TH1D* hj_half_measure = h_half_measure_2D->ProjectionX("hj_half_measure");
    TH1D* hc_half_truth = h_half_truth_2D->ProjectionY("hc_half_truth");
    TH1D* hc_half_measure = h_half_measure_2D->ProjectionY("hc_half_measure");

    std::vector<TH1D*> hj_full_unfold; std::vector<TH1D*> hc_full_unfold;
    std::vector<TH1D*> hj_full_unfold_trim_5; std::vector<TH1D*> hc_full_unfold_trim_5;
    std::vector<TH1D*> hj_full_unfold_trim_10; std::vector<TH1D*> hc_full_unfold_trim_10;
    for (int i = 0; i < 4; ++i) {
        std::ostringstream name, name_trim5, name_trim10;
        if (i % 2 == 0) {
            name << "h_unfold_calib_full_" << (i / 2 + 1);
            name_trim5 << "h_unfold_calib_trim_5_full_" << (i / 2 + 1);
            name_trim10 << "h_unfold_calib_trim_10_full_" << (i / 2 + 1);
        } else {
            name << "h_unfold_calib_full_extra_" << (i / 2 + 1);
            name_trim5 << "h_unfold_calib_trim_5_full_extra_" << (i / 2 + 1);
            name_trim10 << "h_unfold_calib_trim_10_full_extra_" << (i / 2 + 1);
        }

        TH2D* h2d = (TH2D*)f->Get(name.str().c_str());
        if (!h2d) { std::cerr << "Could not find histogram " << name.str() << std::endl; continue; }
        hj_full_unfold.push_back(h2d->ProjectionX(("projx_" + name.str()).c_str()));
        hc_full_unfold.push_back(h2d->ProjectionY(("projy_" + name.str()).c_str()));

        TH2D* h2d_trim5 = (TH2D*)f->Get(name_trim5.str().c_str());
        if (!h2d_trim5) { std::cerr << "Could not find histogram " << name_trim5.str() << std::endl; continue; }
        hj_full_unfold_trim_5.push_back(h2d_trim5->ProjectionX(("projx_" + name_trim5.str()).c_str()));
        hc_full_unfold_trim_5.push_back(h2d_trim5->ProjectionY(("projy_" + name_trim5.str()).c_str()));

        TH2D* h2d_trim10 = (TH2D*)f->Get(name_trim10.str().c_str());
        if (!h2d_trim10) { std::cerr << "Could not find histogram " << name_trim10.str() << std::endl; continue; }
        hj_full_unfold_trim_10.push_back(h2d_trim10->ProjectionX(("projx_" + name_trim10.str()).c_str()));
        hc_full_unfold_trim_10.push_back(h2d_trim10->ProjectionY(("projy_" + name_trim10.str()).c_str()));
    }

    std::vector<TH1D*> hj_unfold; std::vector<TH1D*> hc_unfold;
    std::vector<TH1D*> hj_unfold_trim_5; std::vector<TH1D*> hc_unfold_trim_5;
    std::vector<TH1D*> hj_unfold_trim_10; std::vector<TH1D*> hc_unfold_trim_10;
    for (int i = 0; i < 20; ++i) {
        std::ostringstream name, name_trim5, name_trim10;
        name << "h_unfold_calib_half_" << (i + 1);
        name_trim5 << "h_unfold_calib_trim_5_half_" << (i + 1);
        name_trim10 << "h_unfold_calib_trim_10_half_" << (i + 1);

        TH2D* h2d = (TH2D*)f->Get(name.str().c_str());
        if (!h2d) { std::cerr << "Could not find histogram " << name.str() << std::endl; continue; }
        hj_unfold.push_back(h2d->ProjectionX(("projx_" + name.str()).c_str()));
        hc_unfold.push_back(h2d->ProjectionY(("projy_" + name.str()).c_str()));

        TH2D* h2d_trim5 = (TH2D*)f->Get(name_trim5.str().c_str());
        if (!h2d_trim5) { std::cerr << "Could not find histogram " << name_trim5.str() << std::endl; continue; }
        hj_unfold_trim_5.push_back(h2d_trim5->ProjectionX(("projx_" + name_trim5.str()).c_str()));
        hc_unfold_trim_5.push_back(h2d_trim5->ProjectionY(("projy_" + name_trim5.str()).c_str()));

        TH2D* h2d_trim10 = (TH2D*)f->Get(name_trim10.str().c_str());
        if (!h2d_trim10) { std::cerr << "Could not find histogram " << name_trim10.str() << std::endl; continue; }
        hj_unfold_trim_10.push_back(h2d_trim10->ProjectionX(("projx_" + name_trim10.str()).c_str()));
        hc_unfold_trim_10.push_back(h2d_trim10->ProjectionY(("projy_" + name_trim10.str()).c_str()));
    }

    draw_unfolded_spectra(hj_full_truth, hj_full_measure, hj_full_unfold, true, full_leg_tags, 4, std::make_pair(14.0, 82.0), "figure/h_jet_spectrum_full_closure_1000toys.png");
    draw_unfolded_spectra(hc_full_truth, hc_full_measure, hc_full_unfold, false, full_leg_tags, 4, std::make_pair(0.1,35), "figure/h_et_spectrum_full_closure_1000toys.png");
    draw_unfolded_spectra(hj_full_truth, hj_full_measure, hj_full_unfold_trim_5, true, full_leg_tags, 4, std::make_pair(14.0, 82.0), "figure/h_jet_spectrum_full_trim_5_closure_1000toys.png");
    draw_unfolded_spectra(hc_full_truth, hc_full_measure, hc_full_unfold_trim_5, false, full_leg_tags, 4, std::make_pair(0.1,35), "figure/h_et_spectrum_full_trim_5_closure_1000toys.png");
    draw_unfolded_spectra(hj_full_truth, hj_full_measure, hj_full_unfold_trim_10, true, full_leg_tags, 4, std::make_pair(14.0, 82.0), "figure/h_jet_spectrum_full_trim_10_closure_1000toys.png");
    draw_unfolded_spectra(hc_full_truth, hc_full_measure, hc_full_unfold_trim_10, false, full_leg_tags, 4, std::make_pair(0.1,35), "figure/h_et_spectrum_full_trim_10_closure_1000toys.png");

    draw_unfolded_spectra(hj_half_truth, hj_half_measure, hj_unfold, true, half_leg_tags, 6, std::make_pair(14.0, 82.0), "figure/h_jet_spectrum_half_closure_1000toys.png");
    draw_unfolded_spectra(hc_half_truth, hc_half_measure, hc_unfold, false, half_leg_tags, 6, std::make_pair(0.1, 35), "figure/h_et_spectrum_half_closure_1000toys.png");
    draw_unfolded_spectra(hj_half_truth, hj_half_measure, hj_unfold_trim_5, true, half_leg_tags, 6, std::make_pair(14.0, 82.0), "figure/h_jet_spectrum_half_trim_5_closure_1000toys.png");
    draw_unfolded_spectra(hc_half_truth, hc_half_measure, hc_unfold_trim_5, false, half_leg_tags, 6, std::make_pair(0.1, 35), "figure/h_et_spectrum_half_trim_5_closure_1000toys.png");
    draw_unfolded_spectra(hj_half_truth, hj_half_measure, hj_unfold_trim_10, true, half_leg_tags, 6, std::make_pair(14.0, 82.0), "figure/h_jet_spectrum_half_trim_10_closure_1000toys.png");
    draw_unfolded_spectra(hc_half_truth, hc_half_measure, hc_unfold_trim_10, false, half_leg_tags, 6, std::make_pair(0.1, 35), "figure/h_et_spectrum_half_trim_10_closure_1000toys.png");


    f->Close();
}
