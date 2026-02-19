#include <TFile.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TText.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <string>

int compare_ppg09_mrad_efficiencies() {

    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    // --- Open input files ---
    TFile* f_0mrad = TFile::Open("ppg09_trig10_zvtx_lt_60_output.root");
    TFile* f_1_5mrad = TFile::Open("ppg09_trig12_zvtx_lt_10_output.root");
    
    if (!f_0mrad) {
        std::cerr << "Error: Could not open 0mrad_new_output.root" << std::endl;
        return 1;
    }
    if (!f_1_5mrad) {
        std::cerr << "Error: Could not open 1.5mrad_new_output.root" << std::endl;
        return 1;
    }

    // --- Define efficiency names and labels ---
    // From 0mrad_new_output.root
    std::vector<std::string> eff_names_0mrad = {
        "eff_04_pt_trig18",
        "eff_04_pt_trig22",
    };
    std::vector<std::string> eff_labels_0mrad = {
        "Trig18/Trig10, offline |zvtx| < 60 cm",
        "Trig22/Trig10, offline |zvtx| < 60 cm"};

    // From 1.5mrad_new_output.root
    std::vector<std::string> eff_names_1_5mrad = {
        "eff_04_pt_trig34",
        "eff_04_pt_trig22"
    };
    std::vector<std::string> eff_labels_1_5mrad = {
        "Trig34/Trig12, offline |zvtx| < 10 cm",
        "Trig22/Trig12, offline |zvtx| < 10 cm"
    };

    // --- Load efficiencies ---
    TEfficiency* efficiencies[4] = {nullptr};
    std::vector<std::string> all_labels;
    int color_array[6] = {kBlack, kRed, kBlue, kGreen+2, kMagenta+2, kOrange+7}; // Colors for the 5 efficiencies
    
    // Load from 0mrad file
    for (size_t i = 0; i < eff_names_0mrad.size(); i++) {
        efficiencies[i] = (TEfficiency*) f_0mrad->Get(eff_names_0mrad[i].c_str());
        if (!efficiencies[i]) {
            std::cerr << "Error: TEfficiency " << eff_names_0mrad[i] << " not found in 0mrad file!" << std::endl;
            return 1;
        }
        std::cout << "Loaded " << eff_names_0mrad[i] << " from 0mrad file" << std::endl;
        all_labels.push_back(eff_labels_0mrad[i]);
    }
    
    // Load from 1.5mrad file
    for (size_t i = 0; i < eff_names_1_5mrad.size(); i++) {
        efficiencies[2 + i] = (TEfficiency*) f_1_5mrad->Get(eff_names_1_5mrad[i].c_str());
        if (!efficiencies[2 + i]) {
            std::cerr << "Error: TEfficiency " << eff_names_1_5mrad[i] << " not found in 1.5mrad file!" << std::endl;
            return 1;
        }
        std::cout << "Loaded " << eff_names_1_5mrad[i] << " from 1.5mrad file" << std::endl;
        all_labels.push_back(eff_labels_1_5mrad[i]);
    }

    // --- Create comparison canvas ---
    TCanvas* c = new TCanvas("c", "Efficiency Comparison: 0mrad vs 1.5mrad", 800, 600);
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

    // Check if efficiencies 0 and 1 overlap
    if (efficiencies[0] && efficiencies[1]) {
        TGraphAsymmErrors* gr0 = efficiencies[0]->CreateGraph();
        TGraphAsymmErrors* gr1 = efficiencies[1]->CreateGraph();
        if (gr0 && gr1) {
            int n0 = gr0->GetN();
            int n1 = gr1->GetN();
            std::cout << "\n=== Checking overlap between efficiency 0 and 1 ===" << std::endl;
            std::cout << "Efficiency 0 has " << n0 << " points" << std::endl;
            std::cout << "Efficiency 1 has " << n1 << " points" << std::endl;
            
            int overlap_count = 0;
            double max_diff = 0;
            for (int i = 0; i < n0 && i < n1; i++) {
                double x0, y0, x1, y1;
                gr0->GetPoint(i, x0, y0);
                gr1->GetPoint(i, x1, y1);
                double diff = TMath::Abs(y0 - y1);
                if (diff < 0.001) { // Very close values
                    overlap_count++;
                }
                if (diff > max_diff) max_diff = diff;
                if (i < 10) { // Print first 10 points
                    std::cout << "Point " << i << ": x0=" << x0 << " y0=" << y0 
                              << ", x1=" << x1 << " y1=" << y1 
                              << ", diff=" << diff << std::endl;
                }
            }
            std::cout << "Overlapping points (diff < 0.001): " << overlap_count << " out of " 
                      << TMath::Min(n0, n1) << std::endl;
            std::cout << "Maximum difference: " << max_diff << std::endl;
            std::cout << "================================================\n" << std::endl;
        }
    }

    // Draw all efficiencies
    for (int i = 0; i < 4; i++) {
        // Verify efficiency has data
        TH1D* passed = (TH1D*)efficiencies[i]->GetPassedHistogram();
        TH1D* total = (TH1D*)efficiencies[i]->GetTotalHistogram();
        if (passed && total) {
            std::cout << "Efficiency " << i << " (" << all_labels[i] << "): passed entries=" 
                      << passed->GetEntries() << ", total entries=" << total->GetEntries() << std::endl;
        }
        
        // Style the efficiency before drawing
        efficiencies[i]->SetLineColor(color_array[i]);
        efficiencies[i]->SetMarkerColor(color_array[i]);
        efficiencies[i]->SetMarkerStyle(20 + i);
        if (i == 1) { efficiencies[i]->SetMarkerStyle(24); }
        if (i == 4) { efficiencies[i]->SetMarkerStyle(20); }
        if (i == 5) { efficiencies[i]->SetMarkerStyle(24); }
        efficiencies[i]->SetMarkerSize(1.5); // Slightly larger markers
        efficiencies[i]->SetLineWidth(2);
        
        // Draw all with "P SAME" since frame histogram sets up the axes
        efficiencies[i]->Draw("P SAME");
        c->Update(); // Update after each draw to ensure it's rendered
        leg->AddEntry(efficiencies[i], all_labels[i].c_str(), "lep");
        std::cout << "Drew efficiency " << i << ": " << all_labels[i] << " (color=" << color_array[i] << ")" << std::endl;
    }
    
    leg->Draw("same");
    c->Update();
    
    // Save the plot
    c->SaveAs("ppg09_wzcut_mrad_efficiency_comparison.png");
    std::cout << "Saved efficiency comparison to mrad_efficiency_comparison.png" << std::endl;

    // --- Ratio Plot ---
    // Load numerator and denominator efficiencies
    TEfficiency* eff_num = (TEfficiency*) f_1_5mrad->Get("eff_04_pt_trig34");
    TEfficiency* eff_denom = (TEfficiency*) f_0mrad->Get("eff_04_pt_trig18");
    
    if (!eff_num || !eff_denom) {
        std::cerr << "Warning: Could not load efficiencies for ratio plot! Skipping." << std::endl;
    } else {
        std::cout << "Loaded numerator: eff_04_pt_trig34" << std::endl;
        std::cout << "Loaded denominator: eff_04_pt_trig18" << std::endl;
        
        // Create graphs
        TGraphAsymmErrors* gr_num = eff_num->CreateGraph();
        TGraphAsymmErrors* gr_denom = eff_denom->CreateGraph();
        
        if (!gr_num || !gr_denom) {
            std::cerr << "Warning: Could not create graphs for ratio plot! Skipping." << std::endl;
        } else {
            // Create ratio graph
            TGraphAsymmErrors* gr_ratio = new TGraphAsymmErrors();
            gr_ratio->SetName("ratio_graph_trig34_trig18");
            
            int n_num = gr_num->GetN();
            int n_denom = gr_denom->GetN();
            int n_points = TMath::Min(n_num, n_denom);
            int point_count = 0;
            
            for (int j = 0; j < n_points; j++) {
                double x_num, y_num, x_denom, y_denom;
                gr_num->GetPoint(j, x_num, y_num);
                gr_denom->GetPoint(j, x_denom, y_denom);
                
                // Only add point if both have valid values and x values match (within tolerance)
                if (TMath::Abs(x_num - x_denom) < 0.1 && y_denom > 0 && y_num >= 0) {
                    double ratio = y_num / y_denom;
                    
                    // Calculate errors using error propagation
                    double err_num_low = gr_num->GetErrorYlow(j);
                    double err_num_high = gr_num->GetErrorYhigh(j);
                    double err_denom_low = gr_denom->GetErrorYlow(j);
                    double err_denom_high = gr_denom->GetErrorYhigh(j);
                    
                    // Relative errors
                    double rel_err_num_low = (y_num > 0) ? err_num_low / y_num : 0;
                    double rel_err_num_high = (y_num > 0) ? err_num_high / y_num : 0;
                    double rel_err_denom_low = err_denom_low / y_denom;
                    double rel_err_denom_high = err_denom_high / y_denom;
                    
                    // Combined relative error (assuming uncorrelated)
                    double rel_err_low = TMath::Sqrt(rel_err_num_low*rel_err_num_low + rel_err_denom_low*rel_err_denom_low);
                    double rel_err_high = TMath::Sqrt(rel_err_num_high*rel_err_num_high + rel_err_denom_high*rel_err_denom_high);
                    
                    double err_low = ratio * rel_err_low;
                    double err_high = ratio * rel_err_high;
                    
                    gr_ratio->SetPoint(point_count, x_num, ratio);
                    gr_ratio->SetPointError(point_count, 
                                           gr_num->GetErrorXlow(j), gr_num->GetErrorXhigh(j),
                                           err_low, err_high);
                    point_count++;
                }
            }
            
            if (point_count > 0) {
                std::cout << "Created ratio graph with " << point_count << " points" << std::endl;
                
                // Create ratio plot canvas
                TCanvas* c_ratio = new TCanvas("c_ratio", "Efficiency Ratio: Trig34/Trig18", 800, 600);
                TLegend *leg_ratio = new TLegend(0.15, 0.72, 0.5, 0.92);
                leg_ratio->SetFillStyle(0);
                leg_ratio->SetBorderSize(0);
                leg_ratio->SetTextSize(0.035);
                
                // Create blank histogram to set up axes
                TH1D* h_frame_ratio = new TH1D("h_frame_ratio", "", 100, 0, 60);
                h_frame_ratio->SetMinimum(0.0);
                h_frame_ratio->SetMaximum(2.0);
                h_frame_ratio->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
                h_frame_ratio->GetYaxis()->SetTitle("Ratio");
                h_frame_ratio->Draw();
                
                // Style and draw ratio graph
                gr_ratio->SetLineColor(kRed);
                gr_ratio->SetMarkerColor(kRed);
                gr_ratio->SetMarkerStyle(20);
                gr_ratio->SetMarkerSize(1.5);
                gr_ratio->SetLineWidth(2);
                gr_ratio->Draw("P SAME");
                leg_ratio->AddEntry(gr_ratio, "Trig34/Trig12 / Trig18/Trig10", "lep");
                
                // Draw reference line at y=1
                TLine* line_one = new TLine(0, 1, 60, 1);
                line_one->SetLineStyle(2);
                line_one->SetLineColor(kGray+1);
                line_one->Draw("SAME");
                
                leg_ratio->Draw("same");
                c_ratio->Update();
                
                // Save the plot
                c_ratio->SaveAs("ppg09_wzcut_mrad_efficiency_ratio.png");
                std::cout << "Saved ratio plot to ppg09_wzcut_mrad_efficiency_ratio.png" << std::endl;
            } else {
                std::cerr << "Warning: Ratio graph has no valid points" << std::endl;
            }
        }
    }

    // Clean up
    f_0mrad->Close();
    f_1_5mrad->Close();

    return 0;
}
