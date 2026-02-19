#include <TFile.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <string>

int compare_jet10_efficiencies() {

    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    // --- Open input files ---
    TFile* f_ppg09_0mrad = TFile::Open("ppg09_trig10_zvtx_lt_60_output.root");
    TFile* f_ppg09_1_5mrad = TFile::Open("ppg09_trig12_zvtx_lt_10_output.root");
    TFile* f_ana509_0mrad = TFile::Open("ana509_analysis_output/0mrad_new_dijet_output.root");
    TFile* f_ana509_1_5mrad = TFile::Open("ana509_analysis_output/1.5mrad_new_dijet_output.root");
    
    if (!f_ppg09_0mrad) {
        std::cerr << "Error: Could not open ppg09_trig10_zvtx_lt_60_output.root" << std::endl;
        return 1;
    }
    if (!f_ppg09_1_5mrad) {
        std::cerr << "Error: Could not open ppg09_trig12_zvtx_lt_10_output.root" << std::endl;
        return 1;
    }
    if (!f_ana509_0mrad) {
        std::cerr << "Error: Could not open 0mrad_new_dijet_output.root" << std::endl;
        return 1;
    }
    if (!f_ana509_1_5mrad) {
        std::cerr << "Error: Could not open 1.5mrad_new_dijet_output.root" << std::endl;
        return 1;
    }

    // --- Define efficiency names and labels ---
    std::vector<std::string> eff_names = {
        "eff_04_pt_trig18",           // from ppg09_trig10_zvtx_lt_60_output.root
        "eff_04_pt_trig34",           // from ppg09_trig12_zvtx_lt_10_output.root
        "eff_h_leadjet_jet10",        // from ana509_analysis_output/0mrad_new_dijet_output.root
        "eff_h_leadjet_zcut_jet10"    // from ana509_analysis_output/1.5mrad_new_dijet_output.root
    };
    std::vector<std::string> eff_labels = {
        "PPG09: Trig18/Trig10, offline |zvtx| < 60 cm",
        "PPG09: Trig34/Trig12, offline |zvtx| < 10 cm",
        "PPG10: Trig18/Trig10, offline |zvtx| < 60 cm",
        "PPG10: Trig34/Trig12, offline |zvtx| < 10 cm"
    };
    std::vector<TFile*> eff_files = {
        f_ppg09_0mrad,
        f_ppg09_1_5mrad,
        f_ana509_0mrad,
        f_ana509_1_5mrad
    };

    // --- Load efficiencies ---
    TEfficiency* efficiencies[4] = {nullptr};
    int color_array[4] = {kBlack, kRed, kBlue, kGreen+2};
    
    for (size_t i = 0; i < eff_names.size(); i++) {
        efficiencies[i] = (TEfficiency*) eff_files[i]->Get(eff_names[i].c_str());
        if (!efficiencies[i]) {
            std::cerr << "Error: TEfficiency " << eff_names[i] << " not found!" << std::endl;
            return 1;
        }
        std::cout << "Loaded " << eff_names[i] << " from file " << i << std::endl;
    }

    // --- Create comparison canvas ---
    TCanvas* c = new TCanvas("c", "Efficiency Comparison: Jet10 Efficiencies", 800, 600);
    TLegend *leg = new TLegend(0.15, 0.72, 0.5, 0.92);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);

    // Create blank histogram to set up axes
    TH1D* h_frame = new TH1D("h_frame_compare", "", 100, 0, 60);
    h_frame->SetMinimum(0);
    h_frame->SetMaximum(1.5);
    h_frame->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
    h_frame->GetYaxis()->SetTitle("Efficiency");
    h_frame->Draw();

    // Draw all efficiencies
    for (int i = 0; i < 4; i++) {
        // Verify efficiency has data
        TH1D* passed = (TH1D*)efficiencies[i]->GetPassedHistogram();
        TH1D* total = (TH1D*)efficiencies[i]->GetTotalHistogram();
        if (passed && total) {
            std::cout << "Efficiency " << i << " (" << eff_labels[i] << "): passed entries=" 
                      << passed->GetEntries() << ", total entries=" << total->GetEntries() << std::endl;
        }
        
        // Style the efficiency before drawing
        efficiencies[i]->SetLineColor(color_array[i]);
        efficiencies[i]->SetMarkerColor(color_array[i]);
        efficiencies[i]->SetMarkerStyle(20 + i);
        efficiencies[i]->SetMarkerSize(1.5);
        efficiencies[i]->SetLineWidth(2);
        
        // Draw all with "P SAME" since frame histogram sets up the axes
        efficiencies[i]->Draw("P SAME");
        c->Update(); // Update after each draw to ensure it's rendered
        leg->AddEntry(efficiencies[i], eff_labels[i].c_str(), "lep");
        std::cout << "Drew efficiency " << i << ": " << eff_labels[i] << " (color=" << color_array[i] << ")" << std::endl;
    }
    
    leg->Draw("same");
    c->Update();
    
    // Save the plot
    c->SaveAs("jet10_efficiency_comparison.png");
    std::cout << "Saved efficiency comparison to jet10_efficiency_comparison.png" << std::endl;

    // Clean up
    f_ppg09_0mrad->Close();
    f_ppg09_1_5mrad->Close();
    f_ana509_0mrad->Close();
    f_ana509_1_5mrad->Close();

    return 0;
}
