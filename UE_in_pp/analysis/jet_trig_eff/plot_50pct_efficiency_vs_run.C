#include <TFile.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

// Forward declarations
double Find50PercentPoint(TEfficiency* eff, double* error_low = nullptr, double* error_high = nullptr);
void PlotEfficiency(TEfficiency* eff, int run, int trig);
TEfficiency* RebinEfficiencyBy(TEfficiency* efficiency, int rebin_factor);

int plot_50pct_efficiency_vs_run() {

    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    // --- Read run numbers from list files ---
    std::vector<int> trig18_runs;
    std::vector<int> trig34_runs;
    
    // Read trig18_runs.list
    std::ifstream file18("trig18_runs.list");
    if (!file18.is_open()) {
        std::cerr << "Error: Could not open trig18_runs.list" << std::endl;
        return 1;
    }
    std::string line;
    while (std::getline(file18, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip empty lines and comments
        std::istringstream iss(line);
        int run;
        if (iss >> run) {
            trig18_runs.push_back(run);
        }
    }
    file18.close();
    std::cout << "Read " << trig18_runs.size() << " run numbers from trig18_runs.list" << std::endl;
    
    // Sort trig18 runs
    std::sort(trig18_runs.begin(), trig18_runs.end());
    
    // Read trig34_runs.list
    std::ifstream file34("trig34_runs.list");
    if (!file34.is_open()) {
        std::cerr << "Error: Could not open trig34_runs.list" << std::endl;
        return 1;
    }
    while (std::getline(file34, line)) {
        if (line.empty() || line[0] == '#') continue; // Skip empty lines and comments
        std::istringstream iss(line);
        int run;
        if (iss >> run) {
            trig34_runs.push_back(run);
        }
    }
    file34.close();
    std::cout << "Read " << trig34_runs.size() << " run numbers from trig34_runs.list" << std::endl;
    
    // Sort trig34 runs
    std::sort(trig34_runs.begin(), trig34_runs.end());

    // Batch size for combining runs
    const int batch_size = 100;

    // --- Process trig18 runs in batches of 10 ---
    std::vector<int> trig18_run_numbers;
    std::vector<double> trig18_50pct_pt;
    std::vector<double> trig18_50pct_pt_err;
    std::vector<double> trig18_passed_events;
    std::vector<double> trig18_total_events;
    for (size_t batch_start = 0; batch_start < trig18_runs.size(); batch_start += batch_size) {
        size_t batch_end = TMath::Min(batch_start + batch_size, trig18_runs.size());
        int group_run = trig18_runs[batch_start]; // Use first run as group identifier
        
        // Create combined histograms
        TH1F* h_passed_combined = nullptr;
        TH1F* h_total_combined = nullptr;
        int files_loaded = 0;
        
        for (size_t i = batch_start; i < batch_end; i++) {
            int run = trig18_runs[i];
            std::string filename = Form("test_analysis_output/trig10_zvtx_lt_60_output_%d.root", run);
            
            TFile* f = TFile::Open(filename.c_str());
            if (!f || f->IsZombie()) {
                std::cerr << "Warning: Could not open " << filename << std::endl;
                continue;
            }
            
            TH1F* h_passed = (TH1F*) f->Get("h_leadingjet04pt_all_18");
            TH1F* h_total = (TH1F*) f->Get("h_leadingjet04pt_all_10");
            
            if (!h_passed || !h_total) {
                std::cerr << "Warning: Could not find histograms in " << filename << std::endl;
                f->Close();
                continue;
            }
            
            // Clone first histogram, add subsequent ones
            if (files_loaded == 0) {
                h_passed_combined = (TH1F*) h_passed->Clone(Form("h_passed_combined_%d", group_run));
                h_total_combined = (TH1F*) h_total->Clone(Form("h_total_combined_%d", group_run));
                if (!h_passed_combined || !h_total_combined) {
                    std::cerr << "Error: Failed to clone histograms for run " << run << std::endl;
                    f->Close();
                    break;
                }
                // Reset directory to avoid ownership issues
                h_passed_combined->SetDirectory(0);
                h_total_combined->SetDirectory(0);
            } else {
                if (h_passed_combined && h_total_combined) {
                    h_passed_combined->Add(h_passed);
                    h_total_combined->Add(h_total);
                }
            }
            files_loaded++;
            f->Close();
        }
        
        if (files_loaded == 0 || !h_passed_combined || !h_total_combined) {
            std::cerr << "Warning: No valid files in batch starting at run " << group_run << std::endl;
            if (h_passed_combined) delete h_passed_combined;
            if (h_total_combined) delete h_total_combined;
            continue;
        }
        
        std::cout << "Combined " << files_loaded << " files for trig18 batch starting at run " << group_run << std::endl;
        std::cout << "  Combined passed entries: " << h_passed_combined->GetEntries() << std::endl;
        std::cout << "  Combined total entries: " << h_total_combined->GetEntries() << std::endl;
        
        // Create efficiency
        if (!TEfficiency::CheckConsistency(*h_passed_combined, *h_total_combined)) {
            std::cerr << "Warning: Combined histograms not consistent for batch starting at run " << group_run << std::endl;
            delete h_passed_combined;
            delete h_total_combined;
            continue;
        }
        
        TEfficiency* eff = nullptr;
        try {
            eff = new TEfficiency(*h_passed_combined, *h_total_combined);
            if (!eff) {
                std::cerr << "Error: Failed to create TEfficiency for batch starting at run " << group_run << std::endl;
                delete h_passed_combined;
                delete h_total_combined;
                continue;
            }
        } catch (...) {
            std::cerr << "Exception: Failed to create TEfficiency for batch starting at run " << group_run << std::endl;
            delete h_passed_combined;
            delete h_total_combined;
            continue;
        }
        
        // Rebin efficiency before finding 50% point
        TEfficiency* eff_rebinned = RebinEfficiencyBy(eff, 2);
        TEfficiency* eff_to_use = (eff_rebinned) ? eff_rebinned : eff;
        
        // Find 50% efficiency point with errors
        double pt_50_err_low = 0, pt_50_err_high = 0;
        double pt_50 = Find50PercentPoint(eff_to_use, &pt_50_err_low, &pt_50_err_high);
        
        // Store event counts
        double passed_events = h_passed_combined->GetEntries();
        double total_events = h_total_combined->GetEntries();
        
        if (pt_50 > 0) {
            trig18_run_numbers.push_back(group_run);
            trig18_50pct_pt.push_back(pt_50);
            // Use average of low and high errors, or the larger one for symmetric error
            double avg_err = (pt_50_err_low + pt_50_err_high) / 2.0;
            trig18_50pct_pt_err.push_back(avg_err);
            trig18_passed_events.push_back(passed_events);
            trig18_total_events.push_back(total_events);
            std::cout << "Run group " << group_run << " (trig18, " << files_loaded << " runs): 50% efficiency at pT = " << pt_50 << " GeV" << std::endl;
            PlotEfficiency(eff_to_use, group_run, 18);
        } else {
            std::cerr << "Warning: Could not find 50% point for run group " << group_run << " (trig18)" << std::endl;
            // Plot and save the efficiency
            PlotEfficiency(eff_to_use, group_run, 18);
        }
        
        // Clean up efficiencies
        if (eff_rebinned && eff_rebinned != eff) {
            delete eff_rebinned;
        }
        delete eff;
        delete h_passed_combined;
        delete h_total_combined;
    }
    
    // --- Process trig34 runs in batches of 10 ---
    std::vector<int> trig34_run_numbers;
    std::vector<double> trig34_50pct_pt;
    std::vector<double> trig34_50pct_pt_err;
    std::vector<double> trig34_passed_events;
    std::vector<double> trig34_total_events;
    
    for (size_t batch_start = 0; batch_start < trig34_runs.size(); batch_start += batch_size) {
        size_t batch_end = TMath::Min(batch_start + batch_size, trig34_runs.size());
        int group_run = trig34_runs[batch_start]; // Use first run as group identifier
        
        // Create combined histograms
        TH1F* h_passed_combined = nullptr;
        TH1F* h_total_combined = nullptr;
        int files_loaded = 0;
        
        for (size_t i = batch_start; i < batch_end; i++) {
            int run = trig34_runs[i];
            std::string filename = Form("test_analysis_output/trig12_zvtx_lt_10_output_%d.root", run);
            
            TFile* f = TFile::Open(filename.c_str());
            if (!f || f->IsZombie()) {
                std::cerr << "Warning: Could not open " << filename << std::endl;
                continue;
            }
            
            TH1F* h_passed = (TH1F*) f->Get("h_leadingjet04pt_all_34");
            TH1F* h_total = (TH1F*) f->Get("h_leadingjet04pt_all_12");
            
            if (!h_passed || !h_total) {
                std::cerr << "Warning: Could not find histograms in " << filename << std::endl;
                f->Close();
                continue;
            }
            
            // Clone first histogram, add subsequent ones
            if (files_loaded == 0) {
                h_passed_combined = (TH1F*) h_passed->Clone(Form("h_passed_combined_%d", group_run));
                h_total_combined = (TH1F*) h_total->Clone(Form("h_total_combined_%d", group_run));
                if (!h_passed_combined || !h_total_combined) {
                    std::cerr << "Error: Failed to clone histograms for run " << run << std::endl;
                    f->Close();
                    break;
                }
                // Reset directory to avoid ownership issues
                h_passed_combined->SetDirectory(0);
                h_total_combined->SetDirectory(0);
            } else {
                if (h_passed_combined && h_total_combined) {
                    h_passed_combined->Add(h_passed);
                    h_total_combined->Add(h_total);
                }
            }
            files_loaded++;
            f->Close();
        }
        
        if (files_loaded == 0 || !h_passed_combined || !h_total_combined) {
            std::cerr << "Warning: No valid files in batch starting at run " << group_run << std::endl;
            if (h_passed_combined) delete h_passed_combined;
            if (h_total_combined) delete h_total_combined;
            continue;
        }
        
        std::cout << "Combined " << files_loaded << " files for trig34 batch starting at run " << group_run << std::endl;
        std::cout << "  Combined passed entries: " << h_passed_combined->GetEntries() << std::endl;
        std::cout << "  Combined total entries: " << h_total_combined->GetEntries() << std::endl;
        
        // Create efficiency
        if (!TEfficiency::CheckConsistency(*h_passed_combined, *h_total_combined)) {
            std::cerr << "Warning: Combined histograms not consistent for batch starting at run " << group_run << std::endl;
            delete h_passed_combined;
            delete h_total_combined;
            continue;
        }
        
        TEfficiency* eff = nullptr;
        try {
            eff = new TEfficiency(*h_passed_combined, *h_total_combined);
            if (!eff) {
                std::cerr << "Error: Failed to create TEfficiency for batch starting at run " << group_run << std::endl;
                delete h_passed_combined;
                delete h_total_combined;
                continue;
            }
        } catch (...) {
            std::cerr << "Exception: Failed to create TEfficiency for batch starting at run " << group_run << std::endl;
            delete h_passed_combined;
            delete h_total_combined;
            continue;
        }
        
        // Rebin efficiency before finding 50% point
        TEfficiency* eff_rebinned = RebinEfficiencyBy(eff, 4);
        TEfficiency* eff_to_use = (eff_rebinned) ? eff_rebinned : eff;
        
        // Find 50% efficiency point with errors
        double pt_50_err_low = 0, pt_50_err_high = 0;
        double pt_50 = Find50PercentPoint(eff_to_use, &pt_50_err_low, &pt_50_err_high);
        
        // Store event counts
        double passed_events = h_passed_combined->GetEntries();
        double total_events = h_total_combined->GetEntries();
        
        if (pt_50 > 0) {
            trig34_run_numbers.push_back(group_run);
            trig34_50pct_pt.push_back(pt_50);
            // Use average of low and high errors, or the larger one for symmetric error
            double avg_err = (pt_50_err_low + pt_50_err_high) / 2.0;
            trig34_50pct_pt_err.push_back(avg_err);
            trig34_passed_events.push_back(passed_events);
            trig34_total_events.push_back(total_events);
            std::cout << "Run group " << group_run << " (trig34, " << files_loaded << " runs): 50% efficiency at pT = " << pt_50 << " GeV" << std::endl;
            PlotEfficiency(eff_to_use, group_run, 34);
        } else {
            std::cerr << "Warning: Could not find 50% point for run group " << group_run << " (trig34)" << std::endl;
            // Plot and save the efficiency
            PlotEfficiency(eff_to_use, group_run, 34);
        }
        
        // Clean up efficiencies
        if (eff_rebinned && eff_rebinned != eff) {
            delete eff_rebinned;
        }
        delete eff;
        delete h_passed_combined;
        delete h_total_combined;
    }
    
    // --- Create plots ---
    TCanvas* c = new TCanvas("c", "50% Efficiency pT vs Run Number", 800, 600);
    
    // Create graphs with error bars
    TGraphErrors* gr_trig18 = nullptr;
    TGraphErrors* gr_trig34 = nullptr;
    
    if (trig18_run_numbers.size() > 0) {
        gr_trig18 = new TGraphErrors();
        gr_trig18->SetName("gr_trig18");
        gr_trig18->SetTitle("50% Efficiency p_{T} vs Run Number");
        gr_trig18->GetXaxis()->SetTitle("Run Number");
        gr_trig18->GetYaxis()->SetTitle("p_{T} at 50% Efficiency [GeV]");
        gr_trig18->SetMarkerColor(kRed);
        gr_trig18->SetLineColor(kRed);
        gr_trig18->SetMarkerStyle(20);
        gr_trig18->SetMarkerSize(1.2);
        gr_trig18->SetLineWidth(2);
        
        for (size_t i = 0; i < trig18_run_numbers.size(); i++) {
            gr_trig18->SetPoint(i, trig18_run_numbers[i], trig18_50pct_pt[i]);
            gr_trig18->SetPointError(i, 0, trig18_50pct_pt_err[i]); // x error = 0, y error from efficiency
        }
    }
    
    if (trig34_run_numbers.size() > 0) {
        gr_trig34 = new TGraphErrors();
        gr_trig34->SetName("gr_trig34");
        gr_trig34->SetTitle("50% Efficiency p_{T} vs Run Number");
        gr_trig34->GetXaxis()->SetTitle("Run Number");
        gr_trig34->GetYaxis()->SetTitle("p_{T} at 50% Efficiency [GeV]");
        gr_trig34->SetMarkerColor(kBlue);
        gr_trig34->SetLineColor(kBlue);
        gr_trig34->SetMarkerStyle(21);
        gr_trig34->SetMarkerSize(1.2);
        gr_trig34->SetLineWidth(2);
        
        for (size_t i = 0; i < trig34_run_numbers.size(); i++) {
            gr_trig34->SetPoint(i, trig34_run_numbers[i], trig34_50pct_pt[i]);
            gr_trig34->SetPointError(i, 0, trig34_50pct_pt_err[i]); // x error = 0, y error from efficiency
        }
    }
    
    // Create frame histogram to set up axes
    TH1D* h_frame = new TH1D("h_frame_run", "", 100, 47000.0, 54000.0);
    
    // Determine y-axis range from data
    double y_min = 0;
    double y_max = 0;
    if (trig18_50pct_pt.size() > 0) {
        y_min = *std::min_element(trig18_50pct_pt.begin(), trig18_50pct_pt.end());
        y_max = *std::max_element(trig18_50pct_pt.begin(), trig18_50pct_pt.end());
    }
    if (trig34_50pct_pt.size() > 0) {
        double y_min_34 = *std::min_element(trig34_50pct_pt.begin(), trig34_50pct_pt.end());
        double y_max_34 = *std::max_element(trig34_50pct_pt.begin(), trig34_50pct_pt.end());
        if (trig18_50pct_pt.size() == 0) {
            y_min = y_min_34;
            y_max = y_max_34;
        } else {
            y_min = TMath::Min(y_min, y_min_34);
            y_max = TMath::Max(y_max, y_max_34);
        }
    }
    
    // Add some padding to y-axis
    double y_range = y_max - y_min;
    if (y_range > 0) {
        y_min = y_min - 0.1 * y_range;
        y_max = y_max + 0.1 * y_range;
    } else {
        y_min = 0;
        y_max = 20; // Default range if no data
    }
    
    h_frame->SetMinimum(y_min);
    h_frame->SetMaximum(y_max);
    h_frame->GetXaxis()->SetTitle("Run Number");
    h_frame->GetXaxis()->SetRangeUser(47000, 54000);
    h_frame->GetYaxis()->SetTitle("p_{T} at 50% Efficiency [GeV]");
    h_frame->Draw();
    
    // Draw graphs
    TLegend* leg = new TLegend(0.15, 0.7, 0.4, 0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    
    bool has_data = false;
    
    if (gr_trig18 && gr_trig18->GetN() > 0) {
        gr_trig18->Draw("P SAME");
        leg->AddEntry(gr_trig18, "Trig18/Trig10", "lep");
        has_data = true;
    }
    
    if (gr_trig34 && gr_trig34->GetN() > 0) {
        gr_trig34->Draw("P SAME");
        leg->AddEntry(gr_trig34, "Trig34/Trig12", "lep");
        has_data = true;
    }
    
    if (has_data) {
        leg->Draw("same");
        c->Update();
        c->SaveAs("50pct_efficiency_vs_run.png");
        std::cout << "Saved plot to 50pct_efficiency_vs_run.png" << std::endl;
    } else {
        std::cerr << "Error: No data points to plot!" << std::endl;
        return 1;
    }
    
    // --- Create plot for passed events vs run number ---
    TCanvas* c_passed = new TCanvas("c_passed", "Passed Events vs Run Number", 800, 600);
    TGraph* gr_passed_18 = nullptr;
    TGraph* gr_passed_34 = nullptr;
    
    // Create frame for passed events plot
    double max_passed = 0;
    if (trig18_passed_events.size() > 0) {
        double max_18 = *std::max_element(trig18_passed_events.begin(), trig18_passed_events.end());
        max_passed = max_18;
    }
    if (trig34_passed_events.size() > 0) {
        double max_34 = *std::max_element(trig34_passed_events.begin(), trig34_passed_events.end());
        max_passed = TMath::Max(max_passed, max_34);
    }
    max_passed = max_passed * 1.1; // Add 10% padding
    
    TH1D* h_frame_passed = new TH1D("h_frame_passed", "", 100, 47000, 54000);
    h_frame_passed->SetMinimum(0);
    h_frame_passed->SetMaximum(max_passed);
    h_frame_passed->GetXaxis()->SetTitle("Run Number");
    h_frame_passed->GetXaxis()->SetRangeUser(47000, 54000);
    h_frame_passed->GetYaxis()->SetTitle("Number of Passed Events");
    h_frame_passed->Draw();
    
    TLegend* leg_passed = new TLegend(0.15, 0.7, 0.4, 0.9);
    leg_passed->SetFillStyle(0);
    leg_passed->SetBorderSize(0);
    
    if (trig18_passed_events.size() > 0) {
        gr_passed_18 = new TGraph();
        for (size_t i = 0; i < trig18_run_numbers.size(); i++) {
            gr_passed_18->SetPoint(i, trig18_run_numbers[i], trig18_passed_events[i]);
        }
        gr_passed_18->SetMarkerColor(kRed);
        gr_passed_18->SetLineColor(kRed);
        gr_passed_18->SetMarkerStyle(20);
        gr_passed_18->SetMarkerSize(1.2);
        gr_passed_18->SetLineWidth(2);
        gr_passed_18->Draw("P SAME");
        leg_passed->AddEntry(gr_passed_18, "Trig18/Trig10", "lep");
    }
    
    if (trig34_passed_events.size() > 0) {
        gr_passed_34 = new TGraph();
        for (size_t i = 0; i < trig34_run_numbers.size(); i++) {
            gr_passed_34->SetPoint(i, trig34_run_numbers[i], trig34_passed_events[i]);
        }
        gr_passed_34->SetMarkerColor(kBlue);
        gr_passed_34->SetLineColor(kBlue);
        gr_passed_34->SetMarkerStyle(21);
        gr_passed_34->SetMarkerSize(1.2);
        gr_passed_34->SetLineWidth(2);
        gr_passed_34->Draw("P SAME");
        leg_passed->AddEntry(gr_passed_34, "Trig34/Trig12", "lep");
    }
    
    leg_passed->Draw("same");
    c_passed->Update();
    c_passed->SaveAs("passed_events_vs_run.png");
    std::cout << "Saved plot to passed_events_vs_run.png" << std::endl;
    
    // --- Create plot for total events vs run number ---
    TCanvas* c_total = new TCanvas("c_total", "Total Events vs Run Number", 800, 600);
    TGraph* gr_total_18 = nullptr;
    TGraph* gr_total_34 = nullptr;
    
    // Create frame for total events plot
    double max_total = 0;
    if (trig18_total_events.size() > 0) {
        double max_18 = *std::max_element(trig18_total_events.begin(), trig18_total_events.end());
        max_total = max_18;
    }
    if (trig34_total_events.size() > 0) {
        double max_34 = *std::max_element(trig34_total_events.begin(), trig34_total_events.end());
        max_total = TMath::Max(max_total, max_34);
    }
    max_total = max_total * 1.1; // Add 10% padding
    
    TH1D* h_frame_total = new TH1D("h_frame_total", "", 100, 47000, 54000);
    h_frame_total->SetMinimum(0);
    h_frame_total->SetMaximum(max_total);
    h_frame_total->GetXaxis()->SetTitle("Run Number");
    h_frame_total->GetXaxis()->SetRangeUser(47000, 54000);
    h_frame_total->GetYaxis()->SetTitle("Number of Total Events");
    h_frame_total->Draw();
    
    TLegend* leg_total = new TLegend(0.15, 0.7, 0.4, 0.9);
    leg_total->SetFillStyle(0);
    leg_total->SetBorderSize(0);
    
    if (trig18_total_events.size() > 0) {
        gr_total_18 = new TGraph();
        for (size_t i = 0; i < trig18_run_numbers.size(); i++) {
            gr_total_18->SetPoint(i, trig18_run_numbers[i], trig18_total_events[i]);
        }
        gr_total_18->SetMarkerColor(kRed);
        gr_total_18->SetLineColor(kRed);
        gr_total_18->SetMarkerStyle(20);
        gr_total_18->SetMarkerSize(1.2);
        gr_total_18->SetLineWidth(2);
        gr_total_18->Draw("P SAME");
        leg_total->AddEntry(gr_total_18, "Trig18/Trig10", "lep");
    }
    
    if (trig34_total_events.size() > 0) {
        gr_total_34 = new TGraph();
        for (size_t i = 0; i < trig34_run_numbers.size(); i++) {
            gr_total_34->SetPoint(i, trig34_run_numbers[i], trig34_total_events[i]);
        }
        gr_total_34->SetMarkerColor(kBlue);
        gr_total_34->SetLineColor(kBlue);
        gr_total_34->SetMarkerStyle(21);
        gr_total_34->SetMarkerSize(1.2);
        gr_total_34->SetLineWidth(2);
        gr_total_34->Draw("P SAME");
        leg_total->AddEntry(gr_total_34, "Trig34/Trig12", "lep");
    }
    
    leg_total->Draw("same");
    c_total->Update();
    c_total->SaveAs("total_events_vs_run.png");
    std::cout << "Saved plot to total_events_vs_run.png" << std::endl;
    
    return 0;
}

// Helper function to find 50% efficiency point with error estimation
double Find50PercentPoint(TEfficiency* eff, double* error_low, double* error_high) {
    if (!eff) {
        std::cerr << "Find50PercentPoint: eff is null" << std::endl;
        return -1;
    }
    
    // Initialize error pointers if provided
    if (error_low) *error_low = 0;
    if (error_high) *error_high = 0;
    
    TGraphAsymmErrors* gr = nullptr;
    try {
        gr = eff->CreateGraph();
    } catch (...) {
        std::cerr << "Find50PercentPoint: Exception creating graph" << std::endl;
        return -1;
    }
    
    if (!gr) {
        std::cerr << "Find50PercentPoint: Failed to create graph" << std::endl;
        return -1;
    }
    
    if (gr->GetN() == 0) {
        std::cerr << "Find50PercentPoint: Graph has no points" << std::endl;
        delete gr;
        return -1;
    }
    
    double target_eff = 0.5;
    double pt_50 = -1;
    int crossing_index = -1;
    
    // Find the point where efficiency crosses 0.5
    for (int i = 0; i < gr->GetN() - 1; i++) {
        double x1, y1, x2, y2;
        gr->GetPoint(i, x1, y1);
        gr->GetPoint(i + 1, x2, y2);
        
        // Check if 0.5 is between y1 and y2
        if ((y1 <= target_eff && y2 >= target_eff) || (y1 >= target_eff && y2 <= target_eff)) {
            // Linear interpolation
            if (TMath::Abs(y2 - y1) > 1e-6) {
                pt_50 = x1 + (target_eff - y1) * (x2 - x1) / (y2 - y1);
            } else {
                pt_50 = (x1 + x2) / 2.0;
            }
            crossing_index = i;
            break;
        }
    }
    
    // If we didn't find a crossing, check if any point is exactly at 0.5
    if (pt_50 < 0) {
        for (int i = 0; i < gr->GetN(); i++) {
            double x, y;
            gr->GetPoint(i, x, y);
            if (TMath::Abs(y - target_eff) < 0.01) {
                pt_50 = x;
                crossing_index = i;
                break;
            }
        }
    }
    
    // Estimate error in x-value using slope method
    if (pt_50 > 0 && crossing_index >= 0 && (error_low || error_high)) {
        // Get efficiency errors at the crossing point
        double err_low = gr->GetErrorYlow(crossing_index);
        double err_high = gr->GetErrorYhigh(crossing_index);
        
        // Get the points around the crossing
        double x1, y1, x2, y2;
        if (crossing_index < gr->GetN() - 1) {
            gr->GetPoint(crossing_index, x1, y1);
            gr->GetPoint(crossing_index + 1, x2, y2);
        } else {
            gr->GetPoint(crossing_index - 1, x1, y1);
            gr->GetPoint(crossing_index, x2, y2);
        }
        
        // Calculate the slope of the efficiency curve: dx/dy = (x2-x1)/(y2-y1)
        double slope = 0;
        if (TMath::Abs(y2 - y1) > 1e-6) {
            slope = (x2 - x1) / (y2 - y1);
        } else {
            // If slope is undefined, use bin width as estimate
            slope = TMath::Abs(x2 - x1) / 0.1; // Assume 0.1 efficiency change over bin
        }
        
        // Error in x is approximately: error_x = |slope| * error_y
        // The absolute value ensures positive errors
        double err_x_low = TMath::Abs(slope * err_low);
        double err_x_high = TMath::Abs(slope * err_high);
        
        // Set the errors
        if (error_low) *error_low = err_x_low;
        if (error_high) *error_high = err_x_high;
        
        // Minimum error: use half bin width if error is too small
        double bin_width = TMath::Abs(x2 - x1);
        if (bin_width < 0.1) bin_width = 0.5; // Default if bin width is too small
        if (error_low && *error_low < bin_width * 0.1) *error_low = bin_width * 0.5;
        if (error_high && *error_high < bin_width * 0.1) *error_high = bin_width * 0.5;
    }
    
    delete gr;
    return pt_50;
}

// Helper function to rebin efficiency
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
        delete passed_rebinned;
        delete total_rebinned;
        return efficiency;
    }

    auto* rebinned_efficiency = new TEfficiency(*passed_rebinned, *total_rebinned);
    rebinned_efficiency->SetName(Form("%s_rebinned_%d", efficiency->GetName(), rebin_factor));
    rebinned_efficiency->SetTitle(efficiency->GetTitle());
    
    // Clean up temporary histograms
    delete passed_rebinned;
    delete total_rebinned;
    
    return rebinned_efficiency;
}

// Helper function to plot efficiency when 50% point not found
void PlotEfficiency(TEfficiency* eff, int run, int trig) {
    if (!eff) return;
    
    // Efficiency is already rebinned when passed to this function
    TCanvas* c_eff = new TCanvas(Form("c_eff_run%d_trig%d", run, trig), 
                                  Form("Efficiency Run %d Trig%d", run, trig), 
                                  800, 600);
    
    // Create blank histogram to set up axes
    TH1D* h_frame = new TH1D(Form("h_frame_run%d_trig%d", run, trig), "", 100, 0, 60);
    h_frame->SetMinimum(0);
    h_frame->SetMaximum(1.1);
    h_frame->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
    h_frame->GetYaxis()->SetTitle("Efficiency");
    h_frame->Draw();
    
    // Style and draw efficiency
    eff->SetLineColor(kBlue);
    eff->SetMarkerColor(kBlue);
    eff->SetMarkerStyle(20);
    eff->SetMarkerSize(1.2);
    eff->SetLineWidth(2);
    eff->Draw("P SAME");
    
    // Add text label
    TText* t_label = new TText();
    t_label->SetNDC();
    t_label->SetTextFont(42);
    t_label->SetTextSize(0.04);
    t_label->SetTextAlign(13);
    t_label->DrawText(0.15, 0.95, Form("Run %d, Trig%d (No 50%% point found)", run, trig));
    
    c_eff->Update();
    
    // Save the plot
    std::string filename = Form("test_no_50pct/100batch/run%d_trig%d_efficiency.png", run, trig);
    c_eff->SaveAs(filename.c_str());
    std::cout << "Saved efficiency plot to " << filename << std::endl;
    
    delete c_eff;
    delete h_frame;
    delete t_label;
}
