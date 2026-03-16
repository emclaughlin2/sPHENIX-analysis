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
#include <iterator>

void normalize_hist(TH1D* h) {
    const int nX = h->GetNbinsX();
    double integral = 0;

    // Temporary storage
    std::vector<double>contents(nX) ;
    std::vector<double>errors(nX);

    // First pass: normalize by bin area and compute total
    for (int ix = 1; ix <= nX; ++ix) {
        double wX = h->GetXaxis()->GetBinWidth(ix);
        double value = h->GetBinContent(ix);
        double error = h->GetBinError(ix);
        double norm_value = value / wX;
        double norm_error = error / wX;
        contents[ix-1] = norm_value;
        errors[ix-1] = norm_error;
        integral += value;
    }

    // Second pass: divide by total integral if desired
    for (int ix = 1; ix <= nX; ++ix) {
        double norm_value = contents[ix-1] / integral;
        double norm_error = errors[ix-1] / integral;
        h->SetBinContent(ix, norm_value);
        h->SetBinError(ix, norm_error);
    }
}

void draw_unfolded_spectra(TH1D* truth, TH1D* unfold, bool jet, std::pair<double, double> x_range, const char* output_name = nullptr)
{ 
    TCanvas* canvas = new TCanvas("canvas", "", 600, 800);
    TPad* pad1 = new TPad("pad1", "", 0, 0.5, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->Draw();
    pad1->cd();
    pad1->SetLogy(1);

    truth->SetStats(0);
    unfold->SetStats(0);
    truth->SetLineColor(2);
    truth->SetMarkerColor(2);
    unfold->SetLineColor(1);
    unfold->SetMarkerColor(1);
    truth->GetXaxis()->SetLabelSize(0);
    unfold->GetXaxis()->SetLabelSize(0);
    unfold->GetXaxis()->SetRangeUser(x_range.first, x_range.second);
    unfold->Draw();
    truth->Draw("same");

    TLegend* leg = new TLegend(0.77, 0.65, 0.9, 0.75);
    leg->AddEntry(truth, "2D Unfold", "lp");
    leg->AddEntry(unfold, "1D Unfold", "lp");
    leg->SetTextSize(0.04);
    leg->Draw();

    canvas->cd();
    TPad* pad2 = new TPad("pad2", "", 0, 0.0, 1, 0.5);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();

    TH1D* ratio = (TH1D*)unfold->Clone("ratio");
    ratio->Divide(truth);
    ratio->GetYaxis()->SetTitle("1D/2D Ratio");
    ratio->GetYaxis()->SetNdivisions(510);
    ratio->GetYaxis()->SetRangeUser(0.8, 1.2);
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

void draw_iteration_graph(std::vector<TGraph*> graphs, std::vector<std::string> leg_tags, std::string syst_tag, std::string syst, string fname = "") {
    TCanvas* c1 = new TCanvas("c1", "Val Diff and Errors", 800, 600);
    graphs[0]->GetYaxis()->SetRangeUser(0,0.2);
    graphs[0]->GetXaxis()->SetTitle("Iteration");
    graphs[0]->GetYaxis()->SetTitle("#sigma");
    graphs[0]->Draw("ALP");
    for (int i = 1; i < graphs.size(); i++) {
        graphs[i]->Draw("LP SAME");
    }

    double y_min = 0.92-0.035*(graphs.size()+2);
    TLegend* leg1 = new TLegend(0.2, y_min, 0.5, 0.92);
    leg1->SetTextSize(0.035);
    leg1->AddEntry("","#it{#bf{sPHENIX}} Internal");
    if (syst_tag != "") leg1->AddEntry("",syst_tag.c_str(),"");
    for (int i = 0; i < graphs.size(); i++) {
        leg1->AddEntry(graphs[i], leg_tags[i].c_str(), "lp");
    }
    leg1->Draw();

    string outfile = fname + "_h_1D_unfolding_iterations_all_respmatrices.png";
    if (syst != "") outfile = fname + "_h_1D_unfolding_iterations_" + syst + "_.png";

    c1->SaveAs(outfile.c_str());
}

void plot_1D_unfold(const char* unfoldfile = "run28_output_files/output_unfolded_data_pu_correct_calib_dijet_run28_iter_3_1000toys.root") {
    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    std::vector<std::string> syst = {"jetpt_respmatrix","caloet_respmatrix"};
    TFile* f_data = TFile::Open(unfoldfile);

    TH1D* h_truth1D[2];
    TH1D* h_measure1D[2];
    TH1D* h_unfold1D[2][20];
    h_measure1D[0] = (TH1D*)f_data->Get("h_jetpt");
    h_measure1D[1] = (TH1D*)f_data->Get("h_caloet");
    for (int j = 0; j < 2; j++) { 
        h_truth1D[j] = (TH1D*)f_data->Get(("h_truth_"+syst[j]+"reweight_trim_10").c_str()); 
        for (int n = 0; n < 20; n++) {
            h_unfold1D[j][n] = (TH1D*)f_data->Get(("h_unfold_"+syst[j]+"_reweight_trim_10_"+to_string(n+1)).c_str());
        }
    }

    // For iteration optmization need to find: 
    // Sum of the bin error of unfolded distribution 
    // Variation of the sum of bin content between current unfolding distribution and last unfolding distribution 
    for (int it = 0; it < 20; it++) {
        for (int i = 0; i < syst.size(); i++) {
            normalize_hist(h_unfold1D[i][it]);
        }
    }

    std::vector<double> iteration;
    std::vector<std::vector<double>> sigma_error(2, std::vector<double>(19, 0.0));
    std::vector<std::vector<double>> sigma_val_diff(2, std::vector<double>(19, 0.0));
    std::vector<std::vector<double>> total_error(2, std::vector<double>(19, 0.0));
    
    // absolute difference and error calcuation 
    for (int it = 1; it < 20; it++) {
        iteration.push_back(it);
        for (int i = 0; i < syst.size(); i++) {
            int xstart = 1;
            if (i == 0) { xstart = 2; }
            for (int x = xstart; x < h_unfold1D[i][it]->GetNbinsX(); x++) {
                sigma_error[i][it-1] += pow(h_unfold1D[i][it]->GetBinError(x),2);
                sigma_val_diff[i][it-1] += pow(h_unfold1D[i][it]->GetBinContent(x) - h_unfold1D[i][it-1]->GetBinContent(x),2);
            }
        }
    }
    
    /*
    // relative difference and error calculation 
    for (int it = 1; it < 20; it++) {
        iteration.push_back(it);
        for (int i = 0; i < syst.size(); i++) {
            int xstart = 1;
            if (i == 0) { xstart = 2; }
            for (int x = xstart; x < h_unfold1D[i][it]->GetNbinsX(); x++) {
                sigma_error[i][it-1] += pow(h_unfold1D[i][it]->GetBinError(x)/h_unfold1D[i][it]->GetBinContent(x),2);
                sigma_val_diff[i][it-1] += pow((h_unfold1D[i][it]->GetBinContent(x) - h_unfold1D[i][it-1]->GetBinContent(x))/h_unfold1D[i][it]->GetBinContent(x),2);
            }
        }
    }
    */
    for (int i = 0; i < syst.size(); i++) {
        for (int it = 0; it < sigma_val_diff[i].size(); it++) {
            total_error[i][it] = std::sqrt(sigma_val_diff[i][it] + sigma_error[i][it]);
        }
    }

    for (int i = 0; i < syst.size(); i++) {
        for (int it = 0; it < sigma_val_diff[i].size(); it++) {
            sigma_val_diff[i][it] = sqrt(sigma_val_diff[i][it]);
            sigma_error[i][it] = sqrt(sigma_error[i][it]);
        }
    }    

    for (int i = 0; i < syst.size(); i++) {
        std::cout << "it \t val \t error \t total" << std::endl;
        for (int it = 0; it < 19; it++) {
            std::cout << iteration[it] << "\t" << sigma_val_diff[i][it] << "\t" << sigma_error[i][it] << "\t" << total_error[i][it] << std::endl;
        }
        std::cout << std::endl;
    }

    int n = iteration.size();
    TGraph* g_val_diff[2]; TGraph* g_error[2]; TGraph* g_total[2];
    for (int i = 0; i < syst.size(); i++) {
        g_val_diff[i] = new TGraph(n, &iteration[0], &sigma_val_diff[i][0]);
        g_error[i] = new TGraph(n, &iteration[0], &sigma_error[i][0]);
        g_total[i] = new TGraph(n, &iteration[0], &total_error[i][0]);
        g_val_diff[i]->SetLineColor(kRed); g_val_diff[i]->SetMarkerColor(kRed);
        g_error[i]->SetLineColor(kBlue); g_error[i]->SetMarkerColor(kBlue);
        g_total[i]->SetLineColor(kBlack); g_total[i]->SetMarkerColor(kBlack);
    }

    std::vector<std::string> iter_leg_tags = {"#Sigma#delta_{it}", "#Sigma(#sigma_{stat}^{2} #oplus #sigma_{unfold}^{2})", "#Sigma#delta_{it} #oplus #Sigma(#sigma_{stat}^{2} #oplus #sigma_{unfold}^{2})"};
    std::vector<std::string> trim_iter_leg_tags = {"Jet pT trim < 10 entries reweighted resp matrix","Calo ET trim < 10 entries reweighted resp matrix"};

    for (int i = 0; i < syst.size(); i++) {
        draw_iteration_graph({g_val_diff[i], g_error[i], g_total[i]}, iter_leg_tags, trim_iter_leg_tags[i], syst[i], unfoldfile);
    }
    std::string total_syst_tag = "#Sigma#delta_{it} #oplus #Sigma(#sigma_{stat}^{2} #oplus #sigma_{unfold}^{2})";
    std::string total_syst = "";
    for (int i = 0; i < syst.size(); i++) {
        if (i < 4) { g_total[i]->SetLineColor(i+1); g_total[i]->SetMarkerColor(i+1); }
        else { g_total[i]->SetLineColor(i+2); g_total[i]->SetMarkerColor(i+2); }
    }
    std::vector<TGraph*> graphs(std::begin(g_total), std::end(g_total));
    draw_iteration_graph(graphs, trim_iter_leg_tags, total_syst_tag, total_syst, unfoldfile);
   
    //string jet_outfile = "uniform_bin_run21_figure/h_unfolded_jet_spectrum_1D_vs_2D_iter_3_1000toys.png";
    //string et_outfile = "uniform_bin_run21_figure/h_unfolded_et_spectrum_1D_vs_2D_iter_3_1000toys.png";
    //draw_unfolded_spectra(hj_unfold, hj_unfold_1D, true, std::make_pair(0.0, 1.0), jet_outfile.c_str());
    //draw_unfolded_spectra(hc_unfold, hc_unfold_1D, false, std::make_pair(0.0, 1.0), et_outfile.c_str());

}
