#include <TFile.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <TGraphErrors.h>
#include <TVirtualFitter.h>

int color[4] = {1,2,4,6};

int make_jet_bkg_eff_fits() {

    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    float xMin = 22; float xMax = 62; int nFinePoints = 100;

    // --- Load TEfficiency ---
    TFile* f2 = TFile::Open("analysis_sim_run28_output/output_sim_efrac_bkg_cut.root"); // replace with your file
    TH2D* h2[2]; TH1D* h[2];
    std::vector<std::string> hist_name = {"h_total_measure_tight","h_pass_cut_measure_tight"};
    std::vector<std::string> leg_tags = {"No bkg cut","E fraction cut"};
    for (int i = 0; i < 2; i++) {
        h2[i] = (TH2D*) f2->Get(hist_name[i].c_str());
    }
    TH1D* totalHisto = (TH1D*) h2[0]->ProjectionX("totalHisto");
    TH1D* eff_hist = (TH1D*) h2[1]->ProjectionX("passedHisto");
    h[0] = (TH1D*) h2[0]->ProjectionX("h0");
    h[1] = (TH1D*) h2[1]->ProjectionX("h1");
    eff_hist->Divide(eff_hist, totalHisto, 1.0, 1.0, "B");

    // --- Define fit function ---
    TF1* fit_func = new TF1("fit_func","[0]*x*x + [1]*x + [2]",xMin, xMax);
    eff_hist->Fit(fit_func,"REMQN", "", xMin, xMax); 

    // --- Create graph for fine points (for confidence band) ---
    TGraphErrors *grFineConf = new TGraphErrors(nFinePoints);
    grFineConf->SetName("grFineConf");

    for (int i = 0; i < nFinePoints; ++i) {
        double x = xMin + i * (xMax - xMin) / (nFinePoints - 1);
        grFineConf->SetPoint(i, x, fit_func->Eval(x)); 
    }

    // --- Compute 1σ confidence intervals ---
    TVirtualFitter *fitter = TVirtualFitter::GetFitter();
    if (!fitter) {
        std::cerr << "Error: TVirtualFitter not available!" << std::endl;
        return 0;
    }

    fitter->GetConfidenceIntervals(grFineConf, 0.683);

    // We'll use a smooth interpolation (spline) for the ±1σ curves
    TGraph* gr_plus1  = new TGraph(nFinePoints);
    TGraph* gr_minus1 = new TGraph(nFinePoints);

    for (int i = 0; i < nFinePoints; ++i) {
        double x, y;
        grFineConf->GetPoint(i, x, y);
        double ey = grFineConf->GetErrorY(i);
        gr_plus1->SetPoint(i, x, y + ey);
        gr_minus1->SetPoint(i, x, y - ey);
    }

    // Optionally, fit splines to them (for smooth TF1s)
    TSpline3* spline_plus1  = new TSpline3("fit_plus1_spline",  gr_plus1);
    TSpline3* spline_minus1 = new TSpline3("fit_minus1_spline", gr_minus1);

    // --- Save to ROOT file ---
    TFile outFile("tefficiency_fits.root","RECREATE");
    fit_func->Write("fit_nominal");        // nominal curve
    spline_plus1->Write("fit_plus1");      // upper 1σ
    spline_minus1->Write("fit_minus1");    // lower 1σ
    outFile.Close();
    std::cout << "Saved smooth ±1σ TF1s to tefficiency_fits.root" << std::endl;

    // --- Optional: draw ---
    TCanvas* c = new TCanvas("c","TEfficiency fit with smooth ±1σ",800,600);
    eff_hist->GetYaxis()->SetRangeUser(0.4,1.0);
    //eff_hist->GetXaxis()->SetRangeUser(22,62);
    eff_hist->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
    eff_hist->GetYaxis()->SetTitle("Efficiency");
    eff_hist->Draw();
    fit_func->SetLineColor(kRed);
    fit_func->SetLineWidth(2);
    fit_func->Draw("L SAME");

    spline_plus1->SetLineColor(kGreen+2);
    spline_plus1->SetLineWidth(2);
    spline_plus1->Draw("L SAME");

    spline_minus1->SetLineColor(kBlue+2);
    spline_minus1->SetLineWidth(2);
    spline_minus1->Draw("L SAME");
    TLegend* legend = new TLegend(0.55,0.2,0.85,0.4);
    legend->AddEntry(eff_hist, "Efrac Cut Efficiency", "p");
    legend->AddEntry(fit_func, "Fit", "l");
    legend->AddEntry(spline_plus1, "+1 #sigma", "l");
    legend->AddEntry(spline_minus1, "-1 #sigma", "l");
    legend->Draw();

    c->Update();
    c->SaveAs("tefficiency_fit.png");

        // --- Optional: draw ---
    TCanvas* c2 = new TCanvas("c2","TEfficiency fit with smooth ±1σ",800,600);
    c2->SetLogy(1);
    TLegend *leg = new TLegend(.5,.7,.9,.9);
    for (int i = 0; i < 2; i++) {
        h[i]->GetXaxis()->SetRangeUser(22,62);
        h[i]->SetLineColor(color[i]);
        h[i]->SetMarkerColor(color[i]);
        leg->AddEntry(h[i],leg_tags[i].c_str(),"pl");
        if (i == 0) {
            h[i]->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
            h[i]->GetYaxis()->SetTitle("N_{events}");
            h[i]->Draw();
        } else {
            h[i]->Draw("same");
        }
    }
    leg->Draw("same");
    c2->SaveAs("all_jet_events.png");

    return 0;
}
