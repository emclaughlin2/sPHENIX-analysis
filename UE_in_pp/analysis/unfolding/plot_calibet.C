//double calibptbins[] = {21, 26, 32.5, 40.5, 63.5}; // 52.0
//double truthptbins[] = {17, 21, 26, 32.5, 40.5, 63.5, 82}; // 72
//double calibetbins[] = {-1.08, -0.1, 0.0, 0.1, 1.08, 1.97, 3.05, 4.68, 6.2, 15.0}; // calib bins set 1
//double truthetbins[] = {0.0, 0.5, 1.08, 1.97, 3.05, 4.68, 6.2, 15.0, 35.0};
//int calibnet = sizeof(calibetbins) / sizeof(calibetbins[0]) - 1;
//int truthnet = sizeof(truthetbins) / sizeof(truthetbins[0]) - 1;
//int calibnpt = sizeof(calibptbins) / sizeof(calibptbins[0]) - 1;
//int truthnpt = sizeof(truthptbins) / sizeof(truthptbins[0]) - 1;
#include "unfold_Def.h"

void plot_calibet(bool dijet = 1) 
{
    string infile;
    if (dijet) {
        infile = "analysis_data_run28_output/output_pu_correct_data_dijet_bkg_cut.root";
    } else {
        infile = "analysis_data_run28_output/output_pu_correct_data_efrac_bkg_cut.root";
    }

    // Open file
    TFile* f = TFile::Open(infile.c_str(), "READ");
    if (!f || f->IsZombie()) {
        Error("plot_calibjets", "Could not open file %s", infile.c_str());
        return;
    }

    // Input histogram names
    std::vector<std::string> names = {
        "h_calibjet_pt_dijet_eff",
        "h_calibjet_pt_dijet_pu_correct_et",
        "h_calibjet_pt_pu_up_et",
        "h_calibjet_pt_pu_down_et"
    };

    // Labels
    std::vector<std::string> labels = {
        "No PU correction",
        "PU correction",
        "Dijet PU correction",
        "EFrac PU correction"
    };

    std::vector<TH1D*> projY_hists;
    std::vector<TLine*> mean_lines;
    std::vector<float> means;

    for (size_t i = 0; i < names.size(); i++) {
        TH2D* h_uni = (TH2D*)f->Get(names[i].c_str());
        if (!h_uni) {
            Warning("plot_calibjets_manual","Histogram %s not found!",names[i].c_str());
            continue;
        }

        // Build variable-binning histogram (copy content manually)
        int nx = h_uni->GetNbinsX();
        double xmin = h_uni->GetXaxis()->GetXmin();
        double xmax = h_uni->GetXaxis()->GetXmax();

        TH2D* h_var = new TH2D((names[i] + "_var").c_str(), "",
                               nx, xmin, xmax, calibnet, calibetbins);

        for (int ix = 1; ix <= nx; ix++) {
            for (int iy = 1; iy <= h_uni->GetNbinsY(); iy++) {
                double val = h_uni->GetBinContent(ix,iy);
                double err = h_uni->GetBinError(ix,iy);
                h_var->SetBinContent(ix,iy,val);
                h_var->SetBinError(ix,iy,err);
            }
        }

        // Projection on Y (variable-binned axis)
        TH1D* py = h_var->ProjectionY((names[i] + "_py").c_str());
        double mean = py->GetMean();
        double ymax = py->GetBinContent(py->FindBin(py->GetMean()));
        means.push_back(mean);

        // Scale by bin width
        for (int b = 1; b <= py->GetNbinsX(); b++) {
            double c = py->GetBinContent(b);
            double e = py->GetBinError(b);
            double w = py->GetBinWidth(b);
            if (w > 0) {
                py->SetBinContent(b, c / w);
                py->SetBinError(b, e / w);
            }
        }

        // Normalize to unit integral
        double integral = py->Integral("width");
        if (integral > 0) py->Scale(1.0 / integral);

        py->SetLineWidth(2);
        py->SetLineColor(i+1);
        py->SetMarkerColor(i+1);
        //py->SetTitle(labels[i].c_str());

        projY_hists.push_back(py);

        // Compute mean and draw vertical line
        TLine* line = new TLine(mean, 0, mean, 1);
        line->SetLineColor(i+1);
        line->SetLineStyle(2); // dashed
        line->SetLineWidth(2);
        mean_lines.push_back(line);
    }

    // Plot
    TCanvas* c = new TCanvas("c","PU corrections",600,500);
    gStyle->SetOptStat(0);

    TLegend* leg = new TLegend(0.3,0.75,0.88,0.88);
    if (dijet) {
        leg->AddEntry("","Exclusive Dijet","");
    } else {
        leg->AddEntry("","Inclusive Jet","");
    }
    bool first = true;
    for (size_t i = 0; i < projY_hists.size(); i++) {
        if (!projY_hists[i]) continue;
        if (first) {
            projY_hists[i]->GetXaxis()->SetTitle("#SigmaE_{T} [GeV]");
            projY_hists[i]->GetYaxis()->SetTitle("1/N dN/d#SigmaE_{T} [GeV^{-1}]");
            projY_hists[i]->GetYaxis()->SetRangeUser(0,1);
            //projY_hists[i]->GetXaxis()->SetRangeUser(-1.08,6.2);
            projY_hists[i]->SetMarkerStyle(20);
            projY_hists[i]->Draw();
            first = false;
        } else {
            if (i == 3 || i == 2) {continue;}
            projY_hists[i]->SetMarkerStyle(20);
            projY_hists[i]->Draw("SAME"); 
        }
        mean_lines[i]->Draw("SAME");
        leg->AddEntry(projY_hists[i], (labels[i]+" mean = "+to_string(means[i])).c_str(), "lpe");
    }
    leg->Draw();

    if (dijet) {
        c->SaveAs("sphenix_primary_run28_output_files/calibet_neg_1GeV_et_cut_dijet_et_binning_comparison.png");
    } else {
        c->SaveAs("sphenix_primary_run28_output_files/calibet_neg_1GeV_et_cut_efrac_et_binning_comparison.png");
    }
}
