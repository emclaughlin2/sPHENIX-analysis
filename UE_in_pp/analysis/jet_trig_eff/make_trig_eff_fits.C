#include <TFile.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <vector>
#include <algorithm>
#include <iostream>

int color[4] = {1,2,4,6};

int make_trig_eff_fits() {

    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    // --- Load TEfficiency ---
    TFile* f2 = TFile::Open("analysis_output/output.root"); // replace with your file
    TH1D* h_jets[4];
    std::vector<std::string> hist_name = {"h_leadjet_zcut_mb","h_leadjet_zcut_jet8","h_leadjet_zcut_jet10","h_leadjet_zcut_jet12"};
    std::vector<std::string> leg_tags = {"MB (trig12)","Jet8 (trig33)","Jet10 (trig 34)","Jet12 (trig35)"};
    for (int i = 0; i < 4; i++) {
        h_jets[i] = (TH1D*) f2->Get(hist_name[i].c_str());
    }
    TEfficiency* jet10 = (TEfficiency*) f2->Get("eff_h_leadjet_zcut_jet10");
    TH1D* passedHisto = (TH1D*)jet10->GetPassedHistogram();
    TH1D* totalHisto = (TH1D*)jet10->GetTotalHistogram();
    TH1D* eff_hist = new TH1D(); // Or other appropriate TH1 type
    eff_hist->Divide(passedHisto, totalHisto, 1.0, 1.0, "B");
    if (!jet10) { std::cerr << "TEfficiency not found!" << std::endl; return 1; }

    // --- Define fit function ---
    TF1* fit_func = new TF1("fit_func","[0] + [1]*exp([2]*(x-[3]))",10,35);
    fit_func->SetParameters(0.965, -0.52, -0.38, 10.87);

    // --- Binomial likelihood fit ---
    jet10->Fit(fit_func,"RLS"); // likelihood fit

    // --- Extract fit parameters and errors ---
    int npar = fit_func->GetNpar();
    std::vector<double> p(npar), ep(npar);
    for (int i=0;i<npar;i++) { 
        p[i] = fit_func->GetParameter(i);
        ep[i] = fit_func->GetParError(i);
    }

    // --- Toy MC for ±1σ bands ---
    const int nToys = 1000;
    const int nPoints = 500;
    double xMin = 10, xMax = 35;

    TGraph* gr_plus1   = new TGraph(nPoints);
    TGraph* gr_minus1  = new TGraph(nPoints);

    TRandom3 rnd(0);

    for (int i=0;i<nPoints;i++) {
        double x = xMin + i*(xMax-xMin)/(nPoints-1);

        // Toy evaluation
        std::vector<double> toy_vals;
        for (int t=0;t<nToys;t++) {
            for (int j=0;j<npar;j++) 
                fit_func->SetParameter(j, rnd.Gaus(p[j], ep[j]));
            toy_vals.push_back(fit_func->Eval(x));
        }

        std::sort(toy_vals.begin(), toy_vals.end());
        double y_low  = toy_vals[int(0.16*toy_vals.size())];
        double y_high = toy_vals[int(0.84*toy_vals.size())];

        gr_plus1->SetPoint(i, x, y_high);
        gr_minus1->SetPoint(i, x, y_low);
    }

    // --- Wrap TGraphs with TSpline3 for smooth ±1σ curves ---
    TSpline3* spline_plus1   = new TSpline3("spline_plus1",   gr_plus1);
    TSpline3* spline_minus1  = new TSpline3("spline_minus1",  gr_minus1);

    // --- Styling ---
    fit_func->SetLineColor(kBlack);       // nominal TF1
    spline_plus1->SetLineColor(kRed);
    spline_minus1->SetLineColor(kBlue);

    // --- Save to ROOT file ---
    TFile outFile("tefficiency_fits.root","RECREATE");
    fit_func->Write("fit_nominal");        // nominal curve
    spline_plus1->Write("fit_plus1");      // upper 1σ
    spline_minus1->Write("fit_minus1");    // lower 1σ
    outFile.Close();
    std::cout << "Saved smooth ±1σ TF1s to tefficiency_fits.root" << std::endl;

    // --- Optional: draw ---
    TCanvas* c = new TCanvas("c","TEfficiency fit with smooth ±1σ",800,600);
    jet10->Draw("AP");
    c->Update();
    jet10->GetPaintedGraph()->GetXaxis()->SetRangeUser(0,40);
    c->Update();
    jet10->GetPaintedGraph()->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
    c->Update();
    jet10->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency");
    c->Update();
    fit_func->Draw("L SAME");
    spline_plus1->Draw("L SAME");
    spline_minus1->Draw("L SAME");
    c->Update();
    c->SaveAs("tefficiency_fit.png");

        // --- Optional: draw ---
    TCanvas* c2 = new TCanvas("c2","TEfficiency fit with smooth ±1σ",800,600);
    c2->SetLogy(1);
    TLegend *leg = new TLegend(.5,.7,.9,.9);
    for (int i = 0; i < 4; i++) {
        h_jets[i]->GetXaxis()->SetRangeUser(0,40);
        h_jets[i]->SetLineColor(color[i]);
        h_jets[i]->SetMarkerColor(color[i]);
        leg->AddEntry(h_jets[i],leg_tags[i].c_str(),"pl");
        if (i == 0) {
            h_jets[i]->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
            h_jets[i]->GetYaxis()->SetTitle("N_{events}");
            h_jets[i]->Draw();
        } else {
            h_jets[i]->Draw("same");
        }
    }
    leg->Draw("same");
    c2->SaveAs("all_jet_events.png");

    return 0;
}
