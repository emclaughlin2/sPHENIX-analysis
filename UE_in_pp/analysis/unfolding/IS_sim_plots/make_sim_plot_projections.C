#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <vector>
#include <iostream>

int make_sim_plot_projections(string infile = "analysis_sim_run21_output/output_dijet_sim_iter_1.root", string outname = "pythia_truth_ue_pt_projections") {
    // Open file and get your TH2D
    TFile* f = TFile::Open(infile.c_str());
    TH2D* h2 = (TH2D*) f->Get("h_ue_pt_truth_transverse_record");
    if(!h2) { std::cerr << "TH2D not found!" << std::endl; return 1; }

    int nX = h2->GetNbinsX();
    int nY = h2->GetNbinsY();

    int rgb[20][3] = {
        {230, 25, 75}, {60, 180, 75}, {255, 225, 25}, {0, 130, 200},
        {245, 130, 48}, {145, 30, 180}, {70, 240, 240}, {240, 50, 230},
        {210, 245, 60}, {250, 190, 212}, {0, 128, 128}, {220, 190, 255},
        {170, 110, 40}, {128, 128, 128}, {128, 0, 0}, {0, 0, 0},
        {128, 128, 0}, {255, 215, 180}, {0, 0, 128}, {34, 139, 34}
    };

    std::vector<int> colors;
    for(int i=0;i<20;i++) {
        colors.push_back(TColor::GetColor(rgb[i][0], rgb[i][1], rgb[i][2]));
    }

    // --- Create vectors for projections ---
    std::vector<TH1D*> projections;
    std::vector<TH1D*> projections_norm;

    std::vector<std::string> leg_names = {"p_{T}=17-21GeV","p_{T}=21-26GeV","p_{T}=26-32.5GeV","p_{T}=32.5-40.5GeV","p_{T}=40.5-52GeV","p_{T}=52-63.5GeV","p_{T}=63.5-82GeV"};

    for(int ix = 1; ix < nX; ix++) {
        // Regular Y projection
        TH1D* h_proj = h2->ProjectionY(Form("proj_y_bin%d", ix), ix, ix, "e");
        for(int iy = 1; iy <= nY; iy++) {
            double bw = h_proj->GetXaxis()->GetBinWidth(iy);
            double c = h_proj->GetBinContent(iy);
            double e = h_proj->GetBinError(iy);
            h_proj->SetBinContent(iy, c / bw);
            h_proj->SetBinError(iy, e / bw);
        }
        projections.push_back(h_proj);

        // Self-normalized projection
        TH1D* h_norm = (TH1D*) h_proj->Clone(Form("proj_y_bin%d_norm", ix));
        double integral = h_norm->Integral("width"); // include bin widths in integral
        if(integral > 0) h_norm->Scale(1.0 / integral);
        projections_norm.push_back(h_norm);
    }

    // --- Canvas 1: regular projections ---
    TCanvas* c1 = new TCanvas("c1","Y projections",800,600);
    TLegend* leg1 = new TLegend(0.7,0.1,0.9,0.9);
    leg1->SetBorderSize(1); leg1->SetFillColor(0); leg1->SetTextSize(0.03);
    leg1->AddEntry("","Herwig7 200 GeV p+p","");

    for(size_t i=0;i<projections.size();i++) {
        TH1D* h = projections[i];
        h->SetLineColor(colors[i % colors.size()]);
        h->SetLineWidth(2);
        h->SetStats(0);
        h->GetXaxis()->SetTitle("#SigmaE_{T} [GeV]");
        if(i==0) h->Draw("HIST");
        else h->Draw("HIST SAME");
        leg1->AddEntry(h, leg_names[i].c_str(), "l");
    }
    leg1->Draw();
    c1->SaveAs((outname+".png").c_str());

    // --- Canvas 2: self-normalized projections ---
    TCanvas* c2 = new TCanvas("c2","Y projections normalized",800,600);
    TLegend* leg2 = new TLegend(0.7,0.1,0.9,0.9);
    leg2->SetBorderSize(1); leg2->SetFillColor(0); leg2->SetTextSize(0.03);
    leg2->AddEntry("","Herwig7 200 GeV p+p","");

    for(size_t i=0;i<projections_norm.size();i++) {
        TH1D* h = projections_norm[i];
        h->SetMarkerColor(colors[i % colors.size()]);
        h->SetLineColor(colors[i % colors.size()]);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1.2);
        h->SetStats(0);
        h->GetXaxis()->SetTitle("#SigmaE_{T} [GeV]");
        h->GetYaxis()->SetRangeUser(0,1);
        if(i==0) h->Draw("ep");
        else h->Draw("ep SAME");
        leg2->AddEntry(h, leg_names[i].c_str(), "lpe");
    }
    leg2->Draw();
    c2->SaveAs((outname+"_normalized.png").c_str());

    // --- Save histograms to ROOT file ---
    TFile* fout = TFile::Open((outname+".root").c_str(),"RECREATE");
    for(auto h : projections) h->Write();
    for(auto h : projections_norm) h->Write();
    fout->Close();

    std::cout << "Saved " << projections.size() << " projections and " 
              << projections_norm.size() << " normalized projections." << std::endl;

    return 0;
}
