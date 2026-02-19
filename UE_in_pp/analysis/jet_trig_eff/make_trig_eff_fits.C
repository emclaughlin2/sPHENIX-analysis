#include <TFile.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TSpline.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TText.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TMath.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

int color[4] = {1,2,4,6};
int eff_color[5] = {1,2,4,6,8}; // Colors for efficiency comparison
int quad_color[4] = {1,2,4,6}; // Colors for quadrant comparison

int make_trig_eff_fits() {

    gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    gROOT->ProcessLine("SetsPhenixStyle()");

    // --- Load TEfficiency ---
    TFile* f2 = TFile::Open("ana509_analysis_output/1.5mrad_only_dijet_output.root");  // replace with your file
    TH1D* h_jets[4];
    std::vector<std::string> hist_name = {"h_leadjet_zcut_mb","h_leadjet_zcut_jet8","h_leadjet_zcut_jet10","h_leadjet_zcut_jet12"};
    std::vector<std::string> leg_tags = {"MB (trig12)","Jet8 (trig33)","Jet10 (trig 34)","Jet12 (trig35)"};
    //std::vector<std::string> hist_name = {"h_leadjet_mb","h_leadjet_jet8","h_leadjet_jet10","h_leadjet_jet12"};
    //std::vector<std::string> leg_tags = {"MB (trig10)","Jet8 (trig17)","Jet10 (trig 18)","Jet12 (trig19)"};
    for (int i = 0; i < 4; i++) {
        h_jets[i] = (TH1D*) f2->Get(hist_name[i].c_str());
    }
    TEfficiency* jet10 = (TEfficiency*) f2->Get("eff_h_leadjet_zcut_jet10");
    //TEfficiency* jet10 = (TEfficiency*) f2->Get("eff_h_leadjet_jet10");
    TH1D* passedHisto = (TH1D*)jet10->GetPassedHistogram();
    TH1D* totalHisto = (TH1D*)jet10->GetTotalHistogram();
    TH1D* eff_hist = new TH1D(); // Or other appropriate TH1 type
    eff_hist->Divide(passedHisto, totalHisto, 1.0, 1.0, "B");
    if (!jet10) { std::cerr << "TEfficiency not found!" << std::endl; return 1; }

    // --- Define fit function ---
    //TF1* fit_func = new TF1("fit_func","[0] + [1]*exp([2]*(x-[3]))",10,35);
    //fit_func->SetParameters(0.965,  -0.52, -0.38, 10.87);

    TF1* fit_func0 = new TF1(
        "fit_func",
        "[0] + [1]/pow(1+exp(-[3]*(x-[2])),[4])",
        6, 30
    );
    fit_func0->SetLineColor(6);
    fit_func0->SetParameters(
        0.0022,   // [0] baseline
        0.98,   // [1] plateau height
        10.87,   // [2] turn-on location
        0.85,    // [3] slope
        1.0    // [4] sharpness
    );

    jet10->Fit(fit_func0,"RLS"); // likelihood fit

    TF1* fit_func = new TF1(
        "fit_func",
        "[0] + [1]/pow(1+exp(-[3]*(x-[2])),[4]) + [5]*pow(x,[6])",
        6, 30
    );
    fit_func->SetLineColor(3);
    fit_func->SetParameters(
        0.0022,   // [0] baseline
        0.98,   // [1] plateau height
        10.87,   // [2] turn-on location
        0.85,    // [3] slope
        1.0,     // [4] sharpness
        0.002,
        0.8
    );
    fit_func->FixParameter(0, fit_func0->GetParameter(0));
    fit_func->FixParameter(1, fit_func0->GetParameter(1));
    fit_func->FixParameter(2, fit_func0->GetParameter(2));
    fit_func->FixParameter(3, fit_func0->GetParameter(3));

    //TF1* fit_func = new TF1("fit_func", "[0] + [1]/pow(1+exp(-[3]*(x-[2])),[4]) + [5]*pow(x,[6])", 6, 100);
    //fit_func->SetParameters(0.0022, 0.965, 7.36, 0.467, 5.801, 0.0019, 0.8);

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
    double xMin = 6, xMax = 30;

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
    TFile outFile("1.5mrad_tefficiency_fits.root","RECREATE");
    fit_func->Write("fit_nominal");        // nominal curve
    spline_plus1->Write("fit_plus1");      // upper 1σ
    spline_minus1->Write("fit_minus1");    // lower 1σ
    outFile.Close();
    std::cout << "Saved smooth ±1σ TF1s to tefficiency_fits.root" << std::endl;

    // --- Optional: draw ---
    TCanvas* c = new TCanvas("c","TEfficiency fit with smooth ±1σ",800,600);
    jet10->Draw("AP");
    c->Update();
    jet10->GetPaintedGraph()->GetXaxis()->SetRangeUser(0,60);
    c->Update();
    jet10->GetPaintedGraph()->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
    c->Update();
    jet10->GetPaintedGraph()->GetYaxis()->SetTitle("Efficiency");
    c->Update();
    fit_func->Draw("L SAME");
    spline_plus1->Draw("L SAME");
    spline_minus1->Draw("L SAME");
    c->Update();
    c->SaveAs("1.5mrad_tefficiency_fit.png");

        // --- Optional: draw ---
    TCanvas* c2 = new TCanvas("c2","TEfficiency fit with smooth ±1σ",800,600);
    c2->SetLogy(1);
    TLegend *leg = new TLegend(.5,.7,.9,.9);
    for (int i = 0; i < 4; i++) {
        h_jets[i]->GetXaxis()->SetRangeUser(0,60);
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
    c2->SaveAs("1.5mrad_all_jet_events.png");

    // --- Efficiency Comparison ---
    std::vector<std::string> eff_names = {
        "eff_h_leadjet_jet10",
        "eff_h_leadjet_trig22",
        "eff_h_leadjet_trig22_wzcut",
        "eff_h_leadjet_zcut_jet10",
        "eff_h_leadjet_zcut_trig22"
    };
    std::vector<std::string> eff_labels = {
        "leadjet_jet10",
        "leadjet_trig22",
        "leadjet_trig22_wzcut",
        "leadjet_zcut_jet10",
        "leadjet_zcut_trig22"
    };

    TEfficiency* efficiencies[5] = {nullptr};
    for (int i = 0; i < 5; i++) {
        efficiencies[i] = (TEfficiency*) f2->Get(eff_names[i].c_str());
        if (!efficiencies[i]) {
            std::cerr << "TEfficiency " << eff_names[i] << " not found!" << std::endl;
            return 1;
        }
    }

    // Create comparison canvas
    TCanvas* c3 = new TCanvas("c3","Efficiency Comparison",800,600);
    TLegend *leg_eff = new TLegend(0.15,0.15,0.45,0.4);
    leg_eff->SetFillStyle(0);
    leg_eff->SetBorderSize(0);

    // Create blank histogram to set up axes
    TH1D* h_frame = new TH1D("h_frame_eff", "", 100, 0, 60);
    h_frame->SetMinimum(0);
    h_frame->SetMaximum(1.1);
    h_frame->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
    h_frame->GetYaxis()->SetTitle("Efficiency");
    h_frame->Draw();

    for (int i = 0; i < 5; i++) {
        // Get the painted graph and style it
        efficiencies[i]->SetLineColor(eff_color[i]);
        efficiencies[i]->SetMarkerColor(eff_color[i]);
        efficiencies[i]->SetMarkerStyle(20+i);
        efficiencies[i]->SetMarkerSize(1.2);
        efficiencies[i]->Draw("P SAME");
        leg_eff->AddEntry(efficiencies[i], eff_labels[i].c_str(), "lep");
    }
    leg_eff->Draw("same");
    c3->Update();
    c3->SaveAs("1.5mrad_efficiency_comparison.png");
    std::cout << "Saved efficiency comparison to efficiency_comparison.png" << std::endl;

    // --- Quadrant Efficiency Comparisons ---
    // Define efficiency base names (without quadrant suffix) and their labels
    // Base names are like "eff_h_leadjet_zcut_jet12", quadrants are _0, _1, _2, _3
    std::vector<std::string> quad_eff_bases = {
        "eff_h_leadjet_jet8",
        "eff_h_leadjet_jet10",
        "eff_h_leadjet_jet12",
        "eff_h_leadjet_zcut_jet8",
        "eff_h_leadjet_zcut_jet10",
        "eff_h_leadjet_zcut_jet12"
    };
    std::vector<std::string> quad_eff_labels = {
        "leadjet_jet8",
        "leadjet_jet10",
        "leadjet_jet12",
        "leadjet_zcut_jet8",
        "leadjet_zcut_jet10",
        "leadjet_zcut_jet12"
    };
    std::vector<std::string> quad_labels = {"Quadrant 0", "Quadrant 1", "Quadrant 2", "Quadrant 3"};

    // Create comparison plots for each efficiency type
    for (size_t eff_idx = 0; eff_idx < quad_eff_bases.size(); eff_idx++) {
        TEfficiency* quad_effs[4] = {nullptr};
        bool all_found = true;
        
        // Load all 4 quadrants for this efficiency type
        // Quadrants are _0, _1, _2, _3
        for (int quad = 0; quad < 4; quad++) {
            std::string eff_name = quad_eff_bases[eff_idx] + "_" + std::to_string(quad+1);
            quad_effs[quad] = (TEfficiency*) f2->Get(eff_name.c_str());
            
            if (!quad_effs[quad]) {
                std::cerr << "Warning: TEfficiency " << eff_name << " not found!" << std::endl;
                all_found = false;
            }
        }
        
        if (!all_found) {
            std::cerr << "Skipping quadrant comparison for " << quad_eff_bases[eff_idx] << std::endl;
            continue;
        }

        // Create canvas for this efficiency type
        TCanvas* c_quad = new TCanvas(Form("c_quad_%zu", eff_idx), 
                                      Form("Quadrant Comparison: %s", quad_eff_labels[eff_idx].c_str()),
                                      800, 600);
        TLegend *leg_quad = new TLegend(0.15, 0.15, 0.45, 0.4);
        leg_quad->SetFillStyle(0);
        leg_quad->SetBorderSize(0);

        // Create blank histogram to set up axes
        TH1D* h_frame_quad = new TH1D(Form("h_frame_quad_%zu", eff_idx), "", 100, 0, 60);
        h_frame_quad->SetMinimum(0);
        h_frame_quad->SetMaximum(1.1);
        h_frame_quad->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
        h_frame_quad->GetYaxis()->SetTitle("Efficiency");
        h_frame_quad->Draw();

        for (int quad = 0; quad < 4; quad++) {
            // Style the efficiency
            quad_effs[quad]->SetLineColor(quad_color[quad]);
            quad_effs[quad]->SetMarkerColor(quad_color[quad]);
            quad_effs[quad]->SetMarkerStyle(20 + quad);
            quad_effs[quad]->SetMarkerSize(1.2);
            quad_effs[quad]->Draw("P SAME");
            leg_quad->AddEntry(quad_effs[quad], quad_labels[quad].c_str(), "lep");
        }
        leg_quad->Draw("same");
        
        // Add text label showing the base efficiency name
        TText* t_quad_label = new TText();
        t_quad_label->SetNDC();
        t_quad_label->SetTextFont(42);
        t_quad_label->SetTextSize(0.04);
        t_quad_label->SetTextAlign(13); // left-top alignment
        t_quad_label->DrawText(0.17, 0.95, quad_eff_labels[eff_idx].c_str());
        
        c_quad->Update();
        
        // Save the plot
        std::string out_name = "1.5mrad_quadrant_comparison_" + quad_eff_labels[eff_idx] + ".png";
        c_quad->SaveAs(out_name.c_str());
        std::cout << "Saved quadrant comparison to " << out_name << std::endl;
        
        // --- Create ratio plot for quadrants (divide by quadrant 0) ---
        if (quad_effs[0]) {
            TEfficiency* eff_denom = quad_effs[0]; // Use quadrant 0 as denominator
            
            // Create ratio graphs for quadrants 1, 2, 3
            std::vector<TGraphAsymmErrors*> ratio_graphs;
            std::vector<std::string> ratio_labels;
            TGraphAsymmErrors* gr_denom = eff_denom->CreateGraph();
            
            for (int quad = 1; quad < 4; quad++) {
                if (!quad_effs[quad]) continue;
                
                TGraphAsymmErrors* gr_num = quad_effs[quad]->CreateGraph();
                if (!gr_num || !gr_denom) {
                    std::cerr << "Warning: Could not create graphs for quadrant " << quad << " ratio" << std::endl;
                    continue;
                }
                
                // Create ratio graph
                TGraphAsymmErrors* gr_ratio = new TGraphAsymmErrors();
                gr_ratio->SetName(Form("ratio_graph_quad%d", quad));
                
                int n_points = TMath::Min(gr_num->GetN(), gr_denom->GetN());
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
                    ratio_graphs.push_back(gr_ratio);
                    ratio_labels.push_back(Form("Quadrant %d / Quadrant 0", quad));
                    std::cout << "Created ratio graph for quadrant " << quad << " with " << point_count << " points" << std::endl;
                } else {
                    std::cerr << "Warning: Ratio graph for quadrant " << quad << " has no valid points" << std::endl;
                }
            }
            
            if (ratio_graphs.size() > 0) {
                // Create ratio plot canvas
                TCanvas* c_ratio = new TCanvas(Form("c_ratio_%zu", eff_idx), 
                                              Form("Quadrant Ratio: %s", quad_eff_labels[eff_idx].c_str()),
                                              800, 600);
                TLegend *leg_ratio = new TLegend(0.15, 0.72, 0.5, 0.92);
                leg_ratio->SetFillStyle(0);
                leg_ratio->SetBorderSize(0);
                leg_ratio->SetTextSize(0.035);
                
                // Create blank histogram to set up axes
                TH1D* h_frame_ratio = new TH1D(Form("h_frame_ratio_%zu", eff_idx), "", 100, 0, 60);
                h_frame_ratio->SetMinimum(0.0);
                h_frame_ratio->SetMaximum(2.0);
                h_frame_ratio->GetXaxis()->SetTitle("Leading Jet p_{T} [GeV]");
                h_frame_ratio->GetYaxis()->SetTitle("Ratio");
                h_frame_ratio->Draw();
                
                // Draw ratio graphs
                for (size_t i = 0; i < ratio_graphs.size(); i++) {
                    ratio_graphs[i]->SetLineColor(quad_color[i+1]);
                    ratio_graphs[i]->SetMarkerColor(quad_color[i+1]);
                    ratio_graphs[i]->SetMarkerStyle(20 + i);
                    ratio_graphs[i]->SetMarkerSize(1.5);
                    ratio_graphs[i]->SetLineWidth(2);
                    ratio_graphs[i]->Draw("P SAME");
                    leg_ratio->AddEntry(ratio_graphs[i], ratio_labels[i].c_str(), "lep");
                }
                
                // Draw reference line at y=1
                TLine* line_one = new TLine(0, 1, 60, 1);
                line_one->SetLineStyle(2);
                line_one->SetLineColor(kGray+1);
                line_one->Draw("SAME");
                
                // Add text label showing the base efficiency name
                TText* t_ratio_label = new TText();
                t_ratio_label->SetNDC();
                t_ratio_label->SetTextFont(42);
                t_ratio_label->SetTextSize(0.04);
                t_ratio_label->SetTextAlign(13);
                t_ratio_label->DrawText(0.17, 0.95, quad_eff_labels[eff_idx].c_str());
                
                leg_ratio->Draw("same");
                c_ratio->Update();
                
                // Save the plot
                std::string out_name_ratio = "1.5mrad_quadrant_ratio_" + quad_eff_labels[eff_idx] + ".png";
                c_ratio->SaveAs(out_name_ratio.c_str());
                std::cout << "Saved quadrant ratio plot to " << out_name_ratio << std::endl;
            }
        }
    }

    return 0;
}
