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

void draw_unfolded_spectra(TH1D* truth, TH1D* measure, std::vector<TH1D*> unfold, bool jet, std::vector<std::string> leg_tags, int max_iter, std::pair<double, double> x_range, const char* output_name = nullptr, bool first = true)
{
    TH1D* truth_scale = dynamic_cast<TH1D*>(truth->Clone("truth_scale"));
    TCanvas* canvas = new TCanvas("canvas", "", 600, 800);

    TPad* pad1 = new TPad("pad1", "", 0, 0.5, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    if (jet) pad1->SetLogy(1);
    else pad1->SetLogx(1);

    vector<double> bins;
    vector<double> bin_size;
    for (int i = 1; i <= truth->GetNbinsX() + 1; i++) {
        bins.push_back(truth->GetBinLowEdge(i));
    }
    for (int i = 1; i < bins.size(); i++) {
        bin_size.push_back(bins[i] - bins[i-1]);
    }

    for (int i = 1; i < truth->GetNbinsX() + 1; i++) {
        if (first) truth->SetBinContent(i,truth->GetBinContent(i)/bin_size[i-1]);
        if (first) truth->SetBinError(i,truth->GetBinError(i)/bin_size[i-1]);
        for (int u = 0; u < max_iter; u++) {
            unfold[u]->SetBinContent(i,unfold[u]->GetBinContent(i)/bin_size[i-1]);
            unfold[u]->SetBinError(i,unfold[u]->GetBinError(i)/bin_size[i-1]);
        }
    }

    if (first) {
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
    }

    if (first) {
        truth->Scale(1.0/truth->GetEntries());
        measure->Scale(1.0/measure->GetEntries());
    }
    for (int i = 0; i < max_iter; i++) {
        unfold[i]->Scale(1.0/unfold[i]->GetEntries());
    }

    truth->SetStats(0);
    measure->SetStats(0);
    truth->SetLineColor(2);
    truth->SetMarkerColor(2);
    measure->SetLineColor(1);
    measure->SetMarkerColor(1);
    truth->GetXaxis()->SetLabelSize(0);
    measure->GetXaxis()->SetLabelSize(0);
    truth->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
    //if (!jet) truth->GetYaxis()->SetRangeUser(0,30000);

    for (int i = 0; i < max_iter; ++i) {
        unfold[i]->SetStats(0);
        unfold[i]->SetLineColor(i + 3); // Use TColor::GetColor if needed
        unfold[i]->SetMarkerColor(i + 3);
        unfold[i]->GetXaxis()->SetLabelSize(0);
    }

    measure->Draw();
    truth->Draw("same");
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
    r0->GetYaxis()->SetRangeUser(0, 2);
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

void unfolding_plots() {
    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    std::vector<std::string> leg_tags = { "Unfold, iter 1", "Unfold, iter 2", "Unfold, iter 3", "Unfold, iter 4", "Unfold, iter 5", "Unfold, iter 6", 
    "Unfold, iter 7", "Unfold, iter 8", "Unfold, iter 9", "Unfold, iter 10", "Unfold, iter 11", "Unfold, iter 12", "Unfold, iter 13", 
    "Unfold, iter 14", "Unfold, iter 15", "Unfold, iter 16", "Unfold, iter 17", "Unfold, iter 18", "Unfold, iter 19", "Unfold, iter 20"
    };

    TFile* f = TFile::Open("output_unfolded_data_4bin_newetbin_dijet_1000toys.root");
    TH2D* h_truth_2D = (TH2D*)f->Get("h_truth_calib_dijet");
    TH2D* h_measure_2D = (TH2D*)f->Get("h_calibjet_pt_dijet_eff");
    TH1D* hj_truth = h_truth_2D->ProjectionX("hj_truth");
    TH1D* hj_measure = h_measure_2D->ProjectionX("hj_measure");
    TH1D* hc_truth = h_truth_2D->ProjectionY("hc_truth");
    TH1D* hc_measure = h_measure_2D->ProjectionY("hc_measure");

    //hj_truth->Scale(hj_measure->Integral()/hj_truth->Integral());
    //hc_truth->Scale(hc_measure->Integral()/hc_truth->Integral());

    std::vector<TH2D*> h_unfold, h_unfold_trim_5, h_unfold_trim_10;
    std::vector<TH1D*> hj_unfold; std::vector<TH1D*> hc_unfold;
    std::vector<TH1D*> hj_unfold_trim_5; std::vector<TH1D*> hc_unfold_trim_5;
    std::vector<TH1D*> hj_unfold_trim_10; std::vector<TH1D*> hc_unfold_trim_10;
    for (int i = 0; i < 20; ++i) {
        std::ostringstream name, name_trim5, name_trim10;
        name << "h_unfold_calib_all_" << (i + 1);
        name_trim5 << "h_unfold_calib_trim_5_" << (i + 1);
        name_trim10 << "h_unfold_calib_trim_10_" << (i + 1);

        TH2D* h2d = (TH2D*)f->Get(name.str().c_str());
        if (!h2d) { std::cerr << "Could not find histogram " << name.str() << std::endl; continue; }
        h_unfold.push_back((TH2D*)h2d->Clone(Form("h_unfold_%d",i)));
        hj_unfold.push_back(h2d->ProjectionX(("projx_" + name.str()).c_str()));
        hc_unfold.push_back(h2d->ProjectionY(("projy_" + name.str()).c_str()));

        TH2D* h2d_trim5 = (TH2D*)f->Get(name_trim5.str().c_str());
        if (!h2d_trim5) { std::cerr << "Could not find histogram " << name_trim5.str() << std::endl; continue; }
        h_unfold_trim_5.push_back((TH2D*)h2d_trim5->Clone(Form("h_unfold_trim_5_%d",i)));
        hj_unfold_trim_5.push_back(h2d_trim5->ProjectionX(("projx_" + name_trim5.str()).c_str()));
        hc_unfold_trim_5.push_back(h2d_trim5->ProjectionY(("projy_" + name_trim5.str()).c_str()));

        TH2D* h2d_trim10 = (TH2D*)f->Get(name_trim10.str().c_str());
        if (!h2d_trim10) { std::cerr << "Could not find histogram " << name_trim10.str() << std::endl; continue; }
        h_unfold_trim_10.push_back((TH2D*)h2d_trim10->Clone(Form("h_unfold_trim_10_%d",i)));
        hj_unfold_trim_10.push_back(h2d_trim10->ProjectionX(("projx_" + name_trim10.str()).c_str()));
        hc_unfold_trim_10.push_back(h2d_trim10->ProjectionY(("projy_" + name_trim10.str()).c_str()));
    }
    
    draw_unfolded_spectra(hj_truth, hj_measure, hj_unfold, true, leg_tags, 6, std::make_pair(14.0, 82.0), "figure/h_unfolded_jet_spectrum_1000toys.png",true);
    draw_unfolded_spectra(hc_truth, hc_measure, hc_unfold, false, leg_tags, 6, std::make_pair(0.1, 35), "figure/h_unfolded_et_spectrum_1000toys.png",true);
    draw_unfolded_spectra(hj_truth, hj_measure, hj_unfold_trim_5, true, leg_tags, 6, std::make_pair(14.0, 82.0), "figure/h_unfolded_jet_spectrum_trim_5_1000toys.png",false);
    draw_unfolded_spectra(hc_truth, hc_measure, hc_unfold_trim_5, false, leg_tags, 6, std::make_pair(0.1, 35), "figure/h_unfolded_et_spectrum_trim_5_1000toys.png",false);
    draw_unfolded_spectra(hj_truth, hj_measure, hj_unfold_trim_10, true, leg_tags, 6, std::make_pair(14.0, 82.0), "figure/h_unfolded_jet_spectrum_trim_10_1000toys.png",false);
    draw_unfolded_spectra(hc_truth, hc_measure, hc_unfold_trim_10, false, leg_tags, 6, std::make_pair(0.1, 35), "figure/h_unfolded_et_spectrum_trim_10_1000toys.png",false);
    
    
    // For iteration optmization need to find: 
    // Sum of the bin error of unfolded distribution 
    // Variation of the sum of bin content between current unfolding distribution and last unfolding distribution 

    std::cout << h_unfold.size() << std::endl;

    std::vector<double> iteration;
    std::vector<double> sigma_val, sigma_error, sigma_val_diff;
    std::vector<double> sigma_val_trim_5, sigma_error_trim_5, sigma_val_diff_trim_5;
    std::vector<double> sigma_val_trim_10, sigma_error_trim_10, sigma_val_diff_trim_10;
    for (int it = 0; it < 20; it++) {
        if (it > 0) iteration.push_back(it+1);
        double error = 0;
        double val = 0;
        double error_5 = 0;
        double val_5 = 0;
        double error_10 = 0;
        double val_10 = 0;
        for (int i = 2; i < h_unfold[it]->GetNbinsX(); i++) {
            for (int j = 1; j < h_unfold[it]->GetNbinsY(); j++) {
                val += h_unfold[it]->GetBinContent(i,j);
                val_5 += h_unfold_trim_5[it]->GetBinContent(i,j);
                val_10 += h_unfold_trim_10[it]->GetBinContent(i,j);
                if (it > 0) {
                    error += h_unfold[it]->GetBinError(i,j);
                    error_5 += h_unfold_trim_5[it]->GetBinError(i,j);
                    error_10 += h_unfold_trim_10[it]->GetBinError(i,j);
                }
            }
        }
        sigma_val.push_back(val);
        sigma_val_trim_5.push_back(val_5);
        sigma_val_trim_10.push_back(val_10);
        if (it > 0) {
            sigma_error.push_back(error/val);
            sigma_error_trim_5.push_back(error_5/val_5);
            sigma_error_trim_10.push_back(error_10/val_10);
        }
    }

    for (int it = 1; it < 20; it++) {
        sigma_val_diff.push_back(fabs(sigma_val[it] - sigma_val[it-1])/sigma_val[it]);
        sigma_val_diff_trim_5.push_back(fabs(sigma_val_trim_5[it] - sigma_val_trim_5[it-1])/sigma_val[it]);
        sigma_val_diff_trim_10.push_back(fabs(sigma_val_trim_10[it] - sigma_val_trim_10[it-1])/sigma_val[it]);
    }

    sigma_val.erase(sigma_val.begin());
    sigma_val_trim_5.erase(sigma_val_trim_5.begin());
    sigma_val_trim_10.erase(sigma_val_trim_10.begin());

    std::cout << "it \t val \t error \t total" << std::endl;
    for (int it = 0; it < 19; it++) {
        std::cout << iteration[it] << "\t" << sigma_val_diff[it] << "\t" << sigma_error[it] << "\t" << sigma_val_diff[it] + sigma_error[it] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "it \t val \t error \t total" << std::endl;
    for (int it = 0; it < 19; it++) {
        std::cout << iteration[it] << "\t" << sigma_val_diff_trim_5[it] << "\t" << sigma_error_trim_5[it] << "\t" << sigma_val_diff_trim_5[it] + sigma_error_trim_5[it] << std::endl;
    }
    std::cout << std::endl;

    std::cout << "it \t val \t error \t total" << std::endl;
    for (int it = 0; it < 19; it++) {
        std::cout << iteration[it] << "\t" << sigma_val_diff_trim_10[it] << "\t" << sigma_error_trim_10[it] << "\t" << sigma_val_diff_trim_10[it] + sigma_error_trim_10[it] << std::endl;
    }

    int n = iteration.size();
    TGraph* g_val_diff = new TGraph(n, &iteration[0], &sigma_val_diff[0]);
    TGraph* g_val_diff_trim_5 = new TGraph(n, &iteration[0], &sigma_val_diff_trim_5[0]);
    TGraph* g_val_diff_trim_10 = new TGraph(n, &iteration[0], &sigma_val_diff_trim_10[0]);
    TGraph* g_error = new TGraph(n, &iteration[0], &sigma_error[0]);
    TGraph* g_error_trim_5 = new TGraph(n, &iteration[0], &sigma_error_trim_5[0]);
    TGraph* g_error_trim_10 = new TGraph(n, &iteration[0], &sigma_error_trim_10[0]);

    // Sum in quadrature
    std::vector<double> total_error(n), total_error_trim_5(n), total_error_trim_10(n);
    for (int i = 0; i < n; ++i) {
        total_error[i] = std::sqrt(std::pow(sigma_val_diff[i],2) + std::pow(sigma_error[i],2));
        total_error_trim_5[i] = std::sqrt(std::pow(sigma_val_diff_trim_5[i],2) + std::pow(sigma_error_trim_5[i],2));
        total_error_trim_10[i] = std::sqrt(std::pow(sigma_val_diff_trim_10[i],2) + std::pow(sigma_error_trim_10[i],2));
    }
    TGraph* g_total = new TGraph(n, &iteration[0], &total_error[0]);
    TGraph* g_total_trim_5 = new TGraph(n, &iteration[0], &total_error_trim_5[0]);
    TGraph* g_total_trim_10 = new TGraph(n, &iteration[0], &total_error_trim_10[0]);

    // Styling
    g_val_diff->SetLineColor(kRed); g_val_diff->SetMarkerColor(kRed);
    g_val_diff_trim_5->SetLineColor(kRed); g_val_diff_trim_5->SetMarkerColor(kRed);
    g_val_diff_trim_10->SetLineColor(kRed); g_val_diff_trim_10->SetMarkerColor(kRed);
    g_error->SetLineColor(kBlue); g_error->SetMarkerColor(kBlue);
    g_error_trim_5->SetLineColor(kBlue); g_error_trim_5->SetMarkerColor(kBlue);
    g_error_trim_10->SetLineColor(kBlue); g_error_trim_10->SetMarkerColor(kBlue);
    g_total->SetLineColor(kBlack); g_total->SetMarkerColor(kBlack);
    g_total_trim_5->SetLineColor(kBlack); g_total_trim_5->SetMarkerColor(kBlack);
    g_total_trim_10->SetLineColor(kBlack); g_total_trim_10->SetMarkerColor(kBlack);

    // Draw
    TCanvas* c1 = new TCanvas("c1", "Val Diff and Errors", 800, 600);
    g_val_diff->GetYaxis()->SetRangeUser(0,0.07);
    g_val_diff->GetXaxis()->SetTitle("Iteration");
    g_val_diff->GetYaxis()->SetTitle("#sigma");
    g_val_diff->Draw("ALP");
    g_error->Draw("LP SAME");
    g_total->Draw("LP SAME");

    TLegend* leg1 = new TLegend(0.2, 0.65, 0.5, 0.92);
    leg1->SetTextSize(0.035);
    leg1->AddEntry("","#it{#bf{sPHENIX}} Internal");
    leg1->AddEntry("","Untrimmed resp matrix","");
    leg1->AddEntry(g_val_diff, "#Sigma#delta_{it}", "lp");
    leg1->AddEntry(g_error, "#Sigma(#sigma_{stat}^{2} #oplus #sigma_{unfold}^{2})", "lp");
    leg1->AddEntry(g_total, "#Sigma#delta_{it} #oplus #Sigma(#sigma_{stat}^{2} #oplus #sigma_{unfold}^{2})", "lp");
    leg1->Draw();

    c1->SaveAs("figure/h_unfolding_iterations.png");

    TCanvas* c2 = new TCanvas("c2", "Val Diff and Errors", 800, 600);
    g_val_diff_trim_5->GetYaxis()->SetRangeUser(0,0.07);
    g_val_diff_trim_5->GetXaxis()->SetTitle("Iteration");
    g_val_diff_trim_5->GetYaxis()->SetTitle("#sigma");
    g_val_diff_trim_5->Draw("ALP");
    g_error_trim_5->Draw("LP SAME");
    g_total_trim_5->Draw("LP SAME");

    TLegend* leg2 = new TLegend(0.2, 0.65, 0.5, 0.92);
    leg2->SetTextSize(0.035);
    leg2->AddEntry("","#it{#bf{sPHENIX}} Internal");
    leg2->AddEntry("","Trim < 5 entries resp matrix","");
    leg2->AddEntry(g_val_diff_trim_5, "#Sigma#delta_{it}", "lp");
    leg2->AddEntry(g_error_trim_5, "#Sigma#sqrt{#sigma_{stat}^{2} + #sigma_{unfold}^{2}}", "lp");
    leg2->AddEntry(g_total_trim_5, "#Sigma#delta_{it} #oplus #Sigma(#sigma_{stat}^{2} #oplus #sigma_{unfold}^{2})", "lp");
    leg2->Draw();

    c2->SaveAs("figure/h_unfolding_trim_5_iterations.png");


    TCanvas* c3 = new TCanvas("c3", "Val Diff and Errors", 800, 600);
    g_val_diff_trim_10->GetYaxis()->SetRangeUser(0,0.07);
    g_val_diff_trim_10->GetXaxis()->SetTitle("Iteration");
    g_val_diff_trim_10->GetYaxis()->SetTitle("#sigma");
    g_val_diff_trim_10->Draw("ALP");
    g_error_trim_10->Draw("LP SAME");
    g_total_trim_10->Draw("LP SAME");

    TLegend* leg3 = new TLegend(0.2, 0.65, 0.5, 0.92);
    leg3->SetTextSize(0.035);
    leg3->AddEntry("","#it{#bf{sPHENIX}} Internal");
    leg3->AddEntry("","Trim < 10 entries resp matrix","");
    leg3->AddEntry(g_val_diff_trim_10, "#Sigma#delta_{it}", "lp");
    leg3->AddEntry(g_error_trim_10, "#Sigma#sqrt{#sigma_{stat}^{2} + #sigma_{unfold}^{2}}", "lp");
    leg3->AddEntry(g_total_trim_10, "#Sigma#delta_{it} #oplus #Sigma(#sigma_{stat}^{2} #oplus #sigma_{unfold}^{2})", "lp");
    leg3->Draw();

    c3->SaveAs("figure/h_unfolding_trim_10_iterations.png");

    g_total->SetLineColor(kBlue); g_total->SetMarkerColor(kBlue);
    g_total_trim_5->SetLineColor(kRed); g_total_trim_5->SetMarkerColor(kRed);
    g_total_trim_10->SetLineColor(kGreen+2); g_total_trim_10->SetMarkerColor(kGreen+2);    

    TCanvas* c = new TCanvas("c", "Val Diff and Errors", 800, 600);
    g_total->GetYaxis()->SetRangeUser(0,0.07);
    g_total->GetXaxis()->SetTitle("Iteration");
    g_total->GetYaxis()->SetTitle("#sigma");
    g_total->Draw("ALP");
    g_total_trim_5->Draw("LP SAME");
    g_total_trim_10->Draw("LP SAME");

    TLegend* leg = new TLegend(0.2, 0.7, 0.5, 0.88);
    leg->SetTextSize(0.035);
    leg->AddEntry("","#it{#bf{sPHENIX}} Internal");
    leg->AddEntry(g_total, "Untrimmed resp matrix", "lp");
    leg->AddEntry(g_total_trim_5, "Trim < 5 entries resp matrix", "lp");
    leg->AddEntry(g_total_trim_10, "Trim < 10 entries resp matrix", "lp");
    leg->Draw();

    c->SaveAs("figure/h_unfolding_compare_trim_iterations.png");











}
