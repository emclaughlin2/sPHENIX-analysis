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

namespace {
    TEfficiency* RebinEfficiencyBy(TEfficiency* efficiency, int rebin_factor) {
        if (!efficiency) {
            return nullptr;
        }

        const TH1* passed = efficiency->GetPassedHistogram();
        const TH1* total = efficiency->GetTotalHistogram();
        if (!passed || !total) {
            return efficiency;
        }

        auto* passed_rebinned = static_cast<TH1*>(passed->Clone(
            Form("%s_passed_rebinned_%d", efficiency->GetName(), rebin_factor)));
        auto* total_rebinned = static_cast<TH1*>(total->Clone(
            Form("%s_total_rebinned_%d", efficiency->GetName(), rebin_factor)));

        passed_rebinned->Rebin(rebin_factor);
        total_rebinned->Rebin(rebin_factor);

        if (!TEfficiency::CheckConsistency(*passed_rebinned, *total_rebinned)) {
            std::cerr << "Warning: Rebinned histograms are inconsistent for "
                      << efficiency->GetName() << ". Using original efficiency." << std::endl;
            return efficiency;
        }

        auto* rebinned_efficiency = new TEfficiency(*passed_rebinned, *total_rebinned);
        rebinned_efficiency->SetName(Form("%s_rebinned_%d", efficiency->GetName(), rebin_factor));
        rebinned_efficiency->SetTitle(efficiency->GetTitle());
        return rebinned_efficiency;
    }
}

int compare_jet10_time_efficiencies() {

    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    // --- Open input files ---
    TFile* f_0mrad_early = TFile::Open("ana509_analysis_output/0mrad_early_dijet_output.root");
    TFile* f_0mrad_mid = TFile::Open("ana509_analysis_output/0mrad_mid_dijet_output.root");
    TFile* f_0mrad_1_5mrad = TFile::Open("ana509_analysis_output/trig18_trig34_dijet_output.root");
    TFile* f_1_5mrad = TFile::Open("ana509_analysis_output/1.5mrad_only_dijet_output.root");
    
    if (!f_0mrad_early) {
        std::cerr << "Error: Could not open ppg09_trig10_zvtx_lt_60_output.root" << std::endl;
        return 1;
    }
    if (!f_0mrad_mid) {
        std::cerr << "Error: Could not open ppg09_trig12_zvtx_lt_10_output.root" << std::endl;
        return 1;
    }
    if (!f_0mrad_1_5mrad) {
        std::cerr << "Error: Could not open 0mrad_new_dijet_output.root" << std::endl;
        return 1;
    }
    if (!f_1_5mrad) {
        std::cerr << "Error: Could not open 1.5mrad_new_dijet_output.root" << std::endl;
        return 1;
    }

    // --- Define efficiency names and labels ---
    std::vector<std::string> eff_names = {
        "eff_h_leadjet_jet10",
        "eff_h_leadjet_jet10",
        //"eff_h_leadjet_jet10",
        //"eff_h_leadjet_zcut_jet10",
        "eff_h_leadjet_zcut_jet10"};

    std::vector<std::string> eff_labels = {
        "0mrad early: Trig18/Trig10, offline |zvtx| < 60 cm",
        "0mrad mid: Trig18/Trig10, offline |zvtx| < 60 cm",
        //"0mrad 1.5mrad: Trig18/Trig10, offline |zvtx| < 60 cm",
        //"0mrad 1.5mrad: Trig34/Trig12, offline |zvtx| < 10 cm",
        "1.5mrad: Trig34/Trig12, offline |zvtx| < 10 cm"
    };
    std::vector<TFile*> eff_files = {
        f_0mrad_early,
        f_0mrad_mid,
        //f_0mrad_1_5mrad,
        //f_0mrad_1_5mrad,
        f_1_5mrad
    };

    // --- Load efficiencies ---
    TEfficiency* efficiencies[5] = {nullptr};
    int color_array[5] = {kBlack, kRed, kBlue, kGreen+2, kMagenta};
    
    for (size_t i = 0; i < eff_names.size(); i++) {
        efficiencies[i] = (TEfficiency*) eff_files[i]->Get(eff_names[i].c_str());
        if (!efficiencies[i]) {
            std::cerr << "Error: TEfficiency " << eff_names[i] << " not found!" << std::endl;
            return 1;
        }
        //if (eff_files[i] == f_0mrad_1_5mrad
        //    && (eff_names[i] == "eff_h_leadjet_jet10"
        //        || eff_names[i] == "eff_h_leadjet_zcut_jet10")) {
        //    efficiencies[i] = RebinEfficiencyBy(efficiencies[i], 2);
        //}
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
    for (int i = 0; i < 3; i++) {
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
    c->SaveAs("jet10_vs_time_efficiency_comparison_clean.png");
    std::cout << "Saved efficiency comparison to jet10_vs_time_efficiency_comparison.png" << std::endl;

    // Clean up
    f_0mrad_early->Close();
    f_0mrad_mid->Close();
    f_0mrad_1_5mrad->Close();
    f_1_5mrad->Close();

    return 0;
}
