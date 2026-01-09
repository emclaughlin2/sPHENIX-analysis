#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <vector>
#include <iostream>

// Colors for cluster outlines
int colors[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange, kCyan+1,
                    kViolet, kSpring+8, kAzure-4, kPink+8,
                    kYellow+2, kTeal+3, kGray+2, kBlue+3,
                    kOrange+7, kGreen+4, kMagenta+2, kRed+3,
                    kAzure+6, kBlack, kRed+1, kBlue+1, kGreen+1, 
                    kMagenta+1, kOrange+1, 
                    kCyan+2, kViolet+1, kPink+1, kYellow+3};

/// Draw the boxy outline of a cluster histogram
void draw_cluster_outline(TH2D* h, int color) {
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();

    for (int ix = 1; ix <= nx; ix++) {
        for (int iy = 1; iy <= ny; iy++) {
            if (h->GetBinContent(ix, iy) > 0) {
                double x1 = h->GetXaxis()->GetBinLowEdge(ix);
                double x2 = h->GetXaxis()->GetBinUpEdge(ix);
                double y1 = h->GetYaxis()->GetBinLowEdge(iy);
                double y2 = h->GetYaxis()->GetBinUpEdge(iy);

                // Only draw edges that border empty space
                if (h->GetBinContent(ix-1, iy) <= 0) {
                    TLine* l = new TLine(x1, y1, x1, y2);
                    l->SetLineColor(color); l->SetLineWidth(2); l->Draw("SAME");
                }
                if (h->GetBinContent(ix+1, iy) <= 0) {
                    TLine* l = new TLine(x2, y1, x2, y2);
                    l->SetLineColor(color); l->SetLineWidth(2); l->Draw("SAME");
                }
                if (h->GetBinContent(ix, iy-1) <= 0) {
                    TLine* l = new TLine(x1, y1, x2, y1);
                    l->SetLineColor(color); l->SetLineWidth(2); l->Draw("SAME");
                }
                if (h->GetBinContent(ix, iy+1) <= 0) {
                    TLine* l = new TLine(x1, y2, x2, y2);
                    l->SetLineColor(color); l->SetLineWidth(2); l->Draw("SAME");
                }
            }
        }
    }
}

void draw_topo_event_display(const char* filename = "event_display.root")
{
    //gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
    //gROOT->ProcessLine("SetsPhenixStyle()");
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLightTemperature);   // nice colormap for the background

    // Open file
    TFile* infile = TFile::Open(filename, "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return;
    }

    float leading_jet_phi = 2.18291;
    float delta_phi = (2.0*M_PI)/64.0;
    float emcal_delta_phi = (2.0*M_PI)/256.0;
    float trans_start1 = (leading_jet_phi - 2.0*M_PI/3.0) / delta_phi;
    float trans_end1 = (leading_jet_phi - M_PI/3.0) / delta_phi;
    float trans_start2 = (leading_jet_phi + M_PI/3.0) / delta_phi;
    float trans_end2 = (leading_jet_phi + 2.0*M_PI/3.0) / delta_phi;
    float emcal_trans_start1 = (leading_jet_phi - 2.0*M_PI/3.0) / emcal_delta_phi;
    float emcal_trans_end1 = (leading_jet_phi - M_PI/3.0) / emcal_delta_phi;
    float emcal_trans_start2 = (leading_jet_phi + M_PI/3.0) / emcal_delta_phi;
    float emcal_trans_end2 = (leading_jet_phi + 2.0*M_PI/3.0) / emcal_delta_phi;

    // Load total event display histogram
    TH2D* h_total_emcal = (TH2D*)infile->Get("h_total_topo_emcal_event_display");

    // Load cluster histograms
    std::vector<TH2D*> emcal_clusters;
    for (int i = 0; i <= 29; i++) {
        TString hname = Form("h_topo_emcal_event_display_%d", i);
        TH2D* h = (TH2D*)infile->Get(hname);
        if (h) emcal_clusters.push_back(h);
    }

    // Canvas
    TCanvas* emcal_canvas = new TCanvas("emcal_canvas","Topo Event Display",600,1000);
    TLegend* emcal_legend = new TLegend(.05,.9,.55,.99);
    emcal_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Preliminary","");
    emcal_legend->AddEntry("","Pythia8 200 GeV p+p Topoclusters","");
    emcal_legend->SetFillStyle(0);
    emcal_legend->SetBorderSize(0);
    emcal_legend->SetTextSize(0.04);
    emcal_canvas->SetLogz(1);
    emcal_canvas->SetLeftMargin(0.12);
    emcal_canvas->SetRightMargin(0.2);
    // Draw base heatmap
    h_total_emcal->GetXaxis()->SetTitle("EMCal i#eta");
    h_total_emcal->GetYaxis()->SetTitle("EMCal i#phi");
    h_total_emcal->GetZaxis()->SetTitle("Tower E [GeV]");
    h_total_emcal->Draw("COLZ");
    gPad->Update();
    // Access the palette axis
    TPaletteAxis* emcal_palette = (TPaletteAxis*)h_total_emcal->GetListOfFunctions()->FindObject("palette");
    if (emcal_palette) {
        // Move and resize: Xmin, Ymin, Xmax, Ymax in NDC (normalized device coords)
        emcal_palette->SetX1NDC(0.82);
        emcal_palette->SetX2NDC(0.9);
    }

    // Draw transparent transition boxes
    Double_t alpha = 0.1;
    Int_t lightblue = TColor::GetColorTransparent(kGray, alpha);

    TBox* emcalbox1 = new TBox(0.1, emcal_trans_start1, 95.9, emcal_trans_end1);
    emcalbox1->SetFillColor(lightblue);
    emcalbox1->SetFillStyle(3004);
    emcalbox1->SetLineColor(0);
    emcalbox1->Draw("same");

    TBox* emcalbox2 = new TBox(0.1, emcal_trans_start2, 95.9, emcal_trans_end2);
    emcalbox2->SetFillColor(lightblue);
    emcalbox2->SetFillStyle(3004);
    emcalbox2->SetLineColor(0);
    emcalbox2->Draw("same");
    // Draw each cluster as a boxy outline
    for (size_t ic = 0; ic < emcal_clusters.size(); ic++) {
        draw_cluster_outline(emcal_clusters[ic], colors[ic % 29]);
    }
    emcal_legend->Draw("same");
    emcal_canvas->Update();
    emcal_canvas->SaveAs("sim_plots_run21/emcal_topo_event_display_update.png");
    emcal_canvas->SaveAs("sim_plots_run21/emcal_topo_event_display_update.pdf");

        // Load total event display histogram
    TH2D* h_total_ihcal = (TH2D*)infile->Get("h_total_topo_ihcal_event_display");

    // Load cluster histograms
    std::vector<TH2D*> ihcal_clusters;
    for (int i = 0; i <= 29; i++) {
        TString hname = Form("h_topo_ihcal_event_display_%d", i);
        TH2D* h = (TH2D*)infile->Get(hname);
        if (h) ihcal_clusters.push_back(h);
    }

    // Canvas
    TCanvas* ihcal_canvas = new TCanvas("ihcal_canvas","Topo Event Display",600,1000);
    TLegend* ihcal_legend = new TLegend(.05,.9,.55,.99);
    ihcal_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Preliminary","");
    ihcal_legend->AddEntry("","Pythia8 200 GeV p+p Topoclusters","");
    ihcal_legend->SetFillStyle(0);
    ihcal_legend->SetBorderSize(0);
    ihcal_legend->SetTextSize(0.04);
    ihcal_canvas->SetLogz(1);
    ihcal_canvas->SetRightMargin(0.2);
    // Draw base heatmap
    h_total_ihcal->GetXaxis()->SetTitle("IHCal i#eta");
    h_total_ihcal->GetYaxis()->SetTitle("IHCal i#phi");
    h_total_ihcal->GetZaxis()->SetTitle("Tower E [GeV]");
    h_total_ihcal->Draw("COLZ");
    gPad->Update();
    // Access the palette axis
    TPaletteAxis* ihcal_palette = (TPaletteAxis*)h_total_ihcal->GetListOfFunctions()->FindObject("palette");
    if (ihcal_palette) {
        // Move and resize: Xmin, Ymin, Xmax, Ymax in NDC (normalized device coords)
        ihcal_palette->SetX1NDC(0.82);
        ihcal_palette->SetX2NDC(0.9);
    }
    // Draw each cluster as a boxy outline
    for (size_t ic = 0; ic < ihcal_clusters.size(); ic++) {
        draw_cluster_outline(ihcal_clusters[ic], colors[ic % 29]);
    }
    TBox* ihcalbox1 = new TBox(0, trans_start1, 24, trans_end1);
    ihcalbox1->SetFillColor(lightblue);
    ihcalbox1->SetFillStyle(3004);
    ihcalbox1->SetLineColor(0);
    ihcalbox1->Draw("same");

    TBox* ihcalbox2 = new TBox(0, trans_start2, 24, trans_end2);
    ihcalbox2->SetFillColor(lightblue);
    ihcalbox2->SetFillStyle(3004);
    ihcalbox2->SetLineColor(0);
    ihcalbox2->Draw("same");
    ihcal_legend->Draw("same");
    ihcal_canvas->Update();
    ihcal_canvas->SaveAs("sim_plots_run21/ihcal_topo_event_display_update.png");
    ihcal_canvas->SaveAs("sim_plots_run21/ihcal_topo_event_display_update.pdf");

    // Load total event display histogram
    TH2D* h_total_ohcal = (TH2D*)infile->Get("h_total_topo_ohcal_event_display");
    // Load cluster histograms
    std::vector<TH2D*> ohcal_clusters;
    for (int i = 0; i <= 29; i++) {
        TString hname = Form("h_topo_ohcal_event_display_%d", i);
        TH2D* h = (TH2D*)infile->Get(hname);
        if (h) ohcal_clusters.push_back(h);
    }
    // Canvas
    TCanvas* ohcal_canvas = new TCanvas("ohcal_canvas","Topo Event Display",600,1000);
    TLegend* ohcal_legend = new TLegend(.05,.9,.55,.99);
    ohcal_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Preliminary","");
    ohcal_legend->AddEntry("","Pythia8 200 GeV p+p Topoclusters","");
    ohcal_legend->SetFillStyle(0);
    ohcal_legend->SetBorderSize(0);
    ohcal_legend->SetTextSize(0.04);
    ohcal_canvas->SetLogz(1);
    ohcal_canvas->SetRightMargin(0.2);
    // Draw base heatmap
    h_total_ohcal->GetXaxis()->SetTitle("OHCal i#eta");
    h_total_ohcal->GetYaxis()->SetTitle("OHCal i#phi");
    h_total_ohcal->GetZaxis()->SetTitle("Tower E [GeV]");
    h_total_ohcal->Draw("COLZ");
    gPad->Update();
    // Access the palette axis
    TPaletteAxis* ohcal_palette = (TPaletteAxis*)h_total_ohcal->GetListOfFunctions()->FindObject("palette");
    if (ohcal_palette) {
        // Move and resize: Xmin, Ymin, Xmax, Ymax in NDC (normalized device coords)
        ohcal_palette->SetX1NDC(0.82);
        ohcal_palette->SetX2NDC(0.9);
    }
    TBox* ohcalbox1 = new TBox(0, trans_start1, 24, trans_end1);
    ohcalbox1->SetFillColor(lightblue);
    ohcalbox1->SetFillStyle(3004);
    ohcalbox1->SetLineColor(0);
    ohcalbox1->Draw("same");

    TBox* ohcalbox2 = new TBox(0, trans_start2, 24, trans_end2);
    ohcalbox2->SetFillColor(lightblue);
    ohcalbox2->SetFillStyle(3004);
    ohcalbox2->SetLineColor(0);
    ohcalbox2->Draw("same");
    // Draw each cluster as a boxy outline
    for (size_t ic = 0; ic < ohcal_clusters.size(); ic++) {
        draw_cluster_outline(ohcal_clusters[ic], colors[ic % 29]);
    }
    ohcal_legend->Draw("same");
    ohcal_canvas->Update();
    ohcal_canvas->SaveAs("sim_plots_run21/ohcal_topo_event_display_update.png");
    ohcal_canvas->SaveAs("sim_plots_run21/ohcal_topo_event_display_update.pdf");

        // Load total event display histogram
    TH2D* h_total = (TH2D*)infile->Get("h_total_topo_event_display");

    // Load cluster histograms
    std::vector<TH2D*> clusters;
    for (int i = 0; i <= 29; i++) {
        TString hname = Form("h_topo_event_display_%d", i);
        TH2D* h = (TH2D*)infile->Get(hname);
        if (h) clusters.push_back(h);
    }

    // Canvas
    TCanvas* c1 = new TCanvas("c1","Topo Event Display",600,1000);
    TLegend* legend = new TLegend(.05,.9,.55,.99);
    legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Preliminary","");
    legend->AddEntry("","Pythia8 200 GeV p+p Topoclusters","");
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.04);
    c1->SetLogz(1);
    c1->SetRightMargin(0.2);

    // Draw base heatmap
    h_total->GetXaxis()->SetTitle("Calo i#eta");
    h_total->GetYaxis()->SetTitle("Calo i#phi");
    h_total->GetZaxis()->SetTitle("Tower E [GeV]");
    h_total->Draw("COLZ");

    gPad->Update();

    // Access the palette axis
    TPaletteAxis* palette = (TPaletteAxis*)h_total->GetListOfFunctions()->FindObject("palette");
    if (palette) {
        // Move and resize: Xmin, Ymin, Xmax, Ymax in NDC (normalized device coords)
        palette->SetX1NDC(0.82);
        palette->SetX2NDC(0.9);
    }
    TBox* box1 = new TBox(0, trans_start1, 24, trans_end1);
    box1->SetFillColor(lightblue);
    box1->SetFillStyle(3004);
    box1->SetLineColor(0);
    box1->Draw("same");

    TBox* box2 = new TBox(0, trans_start2, 24, trans_end2);
    box2->SetFillColor(lightblue);
    box2->SetFillStyle(3004);
    box2->SetLineColor(0);
    box2->Draw("same");

    // Draw each cluster as a boxy outline
    for (size_t ic = 0; ic < clusters.size(); ic++) {
        draw_cluster_outline(clusters[ic], colors[ic % 29]);
    }
    legend->Draw("same");
    c1->Update();
    c1->SaveAs("sim_plots_run21/topo_event_display_update.png");
    c1->SaveAs("sim_plots_run21/topo_event_display_update.pdf");

    infile->Close();
}
