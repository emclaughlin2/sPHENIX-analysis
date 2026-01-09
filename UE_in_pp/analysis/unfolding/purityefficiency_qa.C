#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>
#include "unfold_Def.h"

void purityefficiency_qa(string infilename = "analysis_sim_run21_output/output_efrac_bkg_cut_sim_iter_1.root", string outfilename = "output_purityefficiency_qa_efrac_bkg_cut.root", string bkg_cut_type = "efrac") {
  // Read Files
  TFile *f_in = new TFile(infilename.c_str(), "READ");
  if (!f_in) {
    std::cout << "Error: cannot open output_sim.root" << std::endl;
    return;
  }

  // Load and apply sPHENIX style (if available)
  {
    const char* sphenixStylePath = "/sphenix/u/egm2153/spring_2023/sPhenixStyle.C";
    // Attempt to load and apply; if file not present ROOT will print an error but code continues
    gROOT->ProcessLine(Form(".L %s+", sphenixStylePath));
    // If the function exists, call it
    // Note: SetsPhenixStyle() is defined in the style file
    gROOT->ProcessLine("SetsPhenixStyle();");
  }



  // Read in QA histograms
  TH2D* h_truth_qa = (TH2D*)f_in->Get("h_truth_qa");
  TH2D* h_measure_qa = (TH2D*)f_in->Get("h_measure_qa");
  TH2D* h_purity_match = (TH2D*)f_in->Get("h_purity_match");
  TH2D* h_purity_truth_pt_cut = (TH2D*)f_in->Get("h_purity_truth_pt_cut");
  TH2D* h_efficiency_match = (TH2D*)f_in->Get("h_efficiency_match");
  TH2D* h_efficiency_reco_pt_cut = (TH2D*)f_in->Get("h_efficiency_reco_pt_cut");
  TH2D* h_efficiency_reco_bkg_cut = (TH2D*)f_in->Get("h_efficiency_reco_bkg_cut");
  // Additional efficiency variations
  TH2D* h_efficiency_match_reco_bkg_cut = (TH2D*)f_in->Get("h_efficiency_match_reco_bkg_cut");
  TH2D* h_efficiency_reco_pt_cut_reco_bkg_cut = (TH2D*)f_in->Get("h_efficiency_reco_pt_cut_reco_bkg_cut");
  TH2D* h_efficiency_match_reco_pt_cut_reco_bkg_cut = (TH2D*)f_in->Get("h_efficiency_match_reco_pt_cut_reco_bkg_cut");
  TH2D* h_efficiency_match_reco_pt_cut = (TH2D*)f_in->Get("h_efficiency_match_reco_pt_cut");

  // --- Create 1D projections first (X = jet pT, Y = calo ET) ---
  if (!h_measure_qa) { std::cout << "Warning: h_measure_qa not found\n"; }
  if (!h_truth_qa)   { std::cout << "Warning: h_truth_qa not found\n"; }
  if (!h_purity_match) { std::cout << "Warning: h_purity_match not found\n"; }
  if (!h_purity_truth_pt_cut) { std::cout << "Warning: h_purity_truth_pt_cut not found\n"; }
  if (!h_efficiency_match) { std::cout << "Warning: h_efficiency_match not found\n"; }
  if (!h_efficiency_reco_pt_cut) { std::cout << "Warning: h_efficiency_reco_pt_cut not found\n"; }
  if (!h_efficiency_reco_bkg_cut) { std::cout << "Warning: h_efficiency_reco_bkg_cut not found\n"; }

  // Projections of measure and truth (denominators)
  TH1D* h_measure_x = h_measure_qa ? h_measure_qa->ProjectionX("h_measure_qa_projX") : nullptr;
  TH1D* h_measure_y = h_measure_qa ? h_measure_qa->ProjectionY("h_measure_qa_projY") : nullptr;
  TH1D* h_truth_x = h_truth_qa ? h_truth_qa->ProjectionX("h_truth_qa_projX") : nullptr;
  TH1D* h_truth_y = h_truth_qa ? h_truth_qa->ProjectionY("h_truth_qa_projY") : nullptr;

  // Projections of purity 2D histograms (numerators)
  TH1D* purity_match_x = h_purity_match ? h_purity_match->ProjectionX("purity_match_x") : nullptr;
  TH1D* purity_truth_pt_cut_x = h_purity_truth_pt_cut ? h_purity_truth_pt_cut->ProjectionX("purity_truth_pt_cut_x") : nullptr;
  TH1D* purity_match_y = h_purity_match ? h_purity_match->ProjectionY("purity_match_y") : nullptr;
  TH1D* purity_truth_pt_cut_y = h_purity_truth_pt_cut ? h_purity_truth_pt_cut->ProjectionY("purity_truth_pt_cut_y") : nullptr;

  // Projections of efficiency 2D histograms (numerators)
  TH1D* efficiency_match_x = h_efficiency_match ? h_efficiency_match->ProjectionX("efficiency_match_x") : nullptr;
  TH1D* efficiency_reco_pt_cut_x = h_efficiency_reco_pt_cut ? h_efficiency_reco_pt_cut->ProjectionX("efficiency_reco_pt_cut_x") : nullptr;
  TH1D* efficiency_reco_bkg_cut_x = h_efficiency_reco_bkg_cut ? h_efficiency_reco_bkg_cut->ProjectionX("efficiency_reco_bkg_cut_x") : nullptr;
  // Additional efficiency projections
  TH1D* efficiency_match_reco_bkg_cut_x = h_efficiency_match_reco_bkg_cut ? h_efficiency_match_reco_bkg_cut->ProjectionX("efficiency_match_reco_bkg_cut_x") : nullptr;
  TH1D* efficiency_reco_pt_cut_reco_bkg_cut_x = h_efficiency_reco_pt_cut_reco_bkg_cut ? h_efficiency_reco_pt_cut_reco_bkg_cut->ProjectionX("efficiency_reco_pt_cut_reco_bkg_cut_x") : nullptr;
  TH1D* efficiency_match_reco_pt_cut_reco_bkg_cut_x = h_efficiency_match_reco_pt_cut_reco_bkg_cut ? h_efficiency_match_reco_pt_cut_reco_bkg_cut->ProjectionX("efficiency_match_reco_pt_cut_reco_bkg_cut_x") : nullptr;
  TH1D* efficiency_match_reco_pt_cut_x = h_efficiency_match_reco_pt_cut ? h_efficiency_match_reco_pt_cut->ProjectionX("efficiency_match_reco_pt_cut_x") : nullptr;
  TH1D* efficiency_match_y = h_efficiency_match ? h_efficiency_match->ProjectionY("efficiency_match_y") : nullptr;
  TH1D* efficiency_reco_pt_cut_y = h_efficiency_reco_pt_cut ? h_efficiency_reco_pt_cut->ProjectionY("efficiency_reco_pt_cut_y") : nullptr;
  TH1D* efficiency_reco_bkg_cut_y = h_efficiency_reco_bkg_cut ? h_efficiency_reco_bkg_cut->ProjectionY("efficiency_reco_bkg_cut_y") : nullptr;
  TH1D* efficiency_match_reco_bkg_cut_y = h_efficiency_match_reco_bkg_cut ? h_efficiency_match_reco_bkg_cut->ProjectionY("efficiency_match_reco_bkg_cut_y") : nullptr;
  TH1D* efficiency_reco_pt_cut_reco_bkg_cut_y = h_efficiency_reco_pt_cut_reco_bkg_cut ? h_efficiency_reco_pt_cut_reco_bkg_cut->ProjectionY("efficiency_reco_pt_cut_reco_bkg_cut_y") : nullptr;
  TH1D* efficiency_match_reco_pt_cut_reco_bkg_cut_y = h_efficiency_match_reco_pt_cut_reco_bkg_cut ? h_efficiency_match_reco_pt_cut_reco_bkg_cut->ProjectionY("efficiency_match_reco_pt_cut_reco_bkg_cut_y") : nullptr;
  TH1D* efficiency_match_reco_pt_cut_y = h_efficiency_match_reco_pt_cut ? h_efficiency_match_reco_pt_cut->ProjectionY("efficiency_match_reco_pt_cut_y") : nullptr;

  // --- Now divide the 1D histograms to form purity/efficiency vs X and Y ---
  TH1D* purity_match_x_div = nullptr;
  TH1D* purity_truth_pt_cut_x_div = nullptr;
  TH1D* purity_match_y_div = nullptr;
  TH1D* purity_truth_pt_cut_y_div = nullptr;

  if (purity_match_x && h_measure_x) { purity_match_x_div = (TH1D*)purity_match_x->Clone("purity_match_x_div"); purity_match_x_div->Divide(h_measure_x); }
  if (purity_truth_pt_cut_x && h_measure_x) { purity_truth_pt_cut_x_div = (TH1D*)purity_truth_pt_cut_x->Clone("purity_truth_pt_cut_x_div"); purity_truth_pt_cut_x_div->Divide(h_measure_x); }
  if (purity_match_y && h_measure_y) { purity_match_y_div = (TH1D*)purity_match_y->Clone("purity_match_y_div"); purity_match_y_div->Divide(h_measure_y); }
  if (purity_truth_pt_cut_y && h_measure_y) { purity_truth_pt_cut_y_div = (TH1D*)purity_truth_pt_cut_y->Clone("purity_truth_pt_cut_y_div"); purity_truth_pt_cut_y_div->Divide(h_measure_y); }

  TH1D* efficiency_match_x_div = nullptr;
  TH1D* efficiency_reco_pt_cut_x_div = nullptr;
  TH1D* efficiency_reco_bkg_cut_x_div = nullptr;
  TH1D* efficiency_match_reco_bkg_cut_x_div = nullptr;
  TH1D* efficiency_reco_pt_cut_reco_bkg_cut_x_div = nullptr;
  TH1D* efficiency_match_reco_pt_cut_reco_bkg_cut_x_div = nullptr;
  TH1D* efficiency_match_reco_pt_cut_x_div = nullptr;
  TH1D* efficiency_match_y_div = nullptr;
  TH1D* efficiency_reco_pt_cut_y_div = nullptr;
  TH1D* efficiency_reco_bkg_cut_y_div = nullptr;
  TH1D* efficiency_match_reco_bkg_cut_y_div = nullptr;
  TH1D* efficiency_reco_pt_cut_reco_bkg_cut_y_div = nullptr;
  TH1D* efficiency_match_reco_pt_cut_reco_bkg_cut_y_div = nullptr;
  TH1D* efficiency_match_reco_pt_cut_y_div = nullptr;

  if (efficiency_match_x && h_truth_x) { efficiency_match_x_div = (TH1D*)efficiency_match_x->Clone("efficiency_match_x_div"); efficiency_match_x_div->Divide(h_truth_x); }
  if (efficiency_reco_pt_cut_x && h_truth_x) { efficiency_reco_pt_cut_x_div = (TH1D*)efficiency_reco_pt_cut_x->Clone("efficiency_reco_pt_cut_x_div"); efficiency_reco_pt_cut_x_div->Divide(h_truth_x); }
  if (efficiency_reco_bkg_cut_x && h_truth_x) { efficiency_reco_bkg_cut_x_div = (TH1D*)efficiency_reco_bkg_cut_x->Clone("efficiency_reco_bkg_cut_x_div"); efficiency_reco_bkg_cut_x_div->Divide(h_truth_x); }
  if (efficiency_match_reco_bkg_cut_x && h_truth_x) { efficiency_match_reco_bkg_cut_x_div = (TH1D*)efficiency_match_reco_bkg_cut_x->Clone("efficiency_match_reco_bkg_cut_x_div"); efficiency_match_reco_bkg_cut_x_div->Divide(h_truth_x); }
  if (efficiency_reco_pt_cut_reco_bkg_cut_x && h_truth_x) { efficiency_reco_pt_cut_reco_bkg_cut_x_div = (TH1D*)efficiency_reco_pt_cut_reco_bkg_cut_x->Clone("efficiency_reco_pt_cut_reco_bkg_cut_x_div"); efficiency_reco_pt_cut_reco_bkg_cut_x_div->Divide(h_truth_x); }
  if (efficiency_match_reco_pt_cut_reco_bkg_cut_x && h_truth_x) { efficiency_match_reco_pt_cut_reco_bkg_cut_x_div = (TH1D*)efficiency_match_reco_pt_cut_reco_bkg_cut_x->Clone("efficiency_match_reco_pt_cut_reco_bkg_cut_x_div"); efficiency_match_reco_pt_cut_reco_bkg_cut_x_div->Divide(h_truth_x); }
  if (efficiency_match_reco_pt_cut_x && h_truth_x) { efficiency_match_reco_pt_cut_x_div = (TH1D*)efficiency_match_reco_pt_cut_x->Clone("efficiency_match_reco_pt_cut_x_div"); efficiency_match_reco_pt_cut_x_div->Divide(h_truth_x); }
  if (efficiency_match_y && h_truth_y) { efficiency_match_y_div = (TH1D*)efficiency_match_y->Clone("efficiency_match_y_div"); efficiency_match_y_div->Divide(h_truth_y); }
  if (efficiency_reco_pt_cut_y && h_truth_y) { efficiency_reco_pt_cut_y_div = (TH1D*)efficiency_reco_pt_cut_y->Clone("efficiency_reco_pt_cut_y_div"); efficiency_reco_pt_cut_y_div->Divide(h_truth_y); }
  if (efficiency_reco_bkg_cut_y && h_truth_y) { efficiency_reco_bkg_cut_y_div = (TH1D*)efficiency_reco_bkg_cut_y->Clone("efficiency_reco_bkg_cut_y_div"); efficiency_reco_bkg_cut_y_div->Divide(h_truth_y); }
  if (efficiency_match_reco_bkg_cut_y && h_truth_y) { efficiency_match_reco_bkg_cut_y_div = (TH1D*)efficiency_match_reco_bkg_cut_y->Clone("efficiency_match_reco_bkg_cut_y_div"); efficiency_match_reco_bkg_cut_y_div->Divide(h_truth_y); }
  if (efficiency_reco_pt_cut_reco_bkg_cut_y && h_truth_y) { efficiency_reco_pt_cut_reco_bkg_cut_y_div = (TH1D*)efficiency_reco_pt_cut_reco_bkg_cut_y->Clone("efficiency_reco_pt_cut_reco_bkg_cut_y_div"); efficiency_reco_pt_cut_reco_bkg_cut_y_div->Divide(h_truth_y); }
  if (efficiency_match_reco_pt_cut_reco_bkg_cut_y && h_truth_y) { efficiency_match_reco_pt_cut_reco_bkg_cut_y_div = (TH1D*)efficiency_match_reco_pt_cut_reco_bkg_cut_y->Clone("efficiency_match_reco_pt_cut_reco_bkg_cut_y_div"); efficiency_match_reco_pt_cut_reco_bkg_cut_y_div->Divide(h_truth_y); }
  if (efficiency_match_reco_pt_cut_y && h_truth_y) { efficiency_match_reco_pt_cut_y_div = (TH1D*)efficiency_match_reco_pt_cut_y->Clone("efficiency_match_reco_pt_cut_y_div"); efficiency_match_reco_pt_cut_y_div->Divide(h_truth_y); }

  // --- Draw the four 1D plots using plain ROOT calls (no helper) ---
  // Purity vs jet pT (X)
  TCanvas* c_purity_pt = new TCanvas("c_purity_pt","Purity vs pT",800,600);
  if (purity_match_x_div) {
    purity_match_x_div->SetLineColor(kRed); purity_match_x_div->SetMarkerStyle(20); purity_match_x_div->SetMarkerColor(kRed);
    purity_match_x_div->SetTitle("");
    purity_match_x_div->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
    purity_match_x_div->GetYaxis()->SetTitle("Purity");
    purity_match_x_div->GetYaxis()->SetRangeUser(0.0,1.5);
    purity_match_x_div->Draw("E1");
  }
  if (purity_truth_pt_cut_x_div) { purity_truth_pt_cut_x_div->SetLineColor(kBlue); purity_truth_pt_cut_x_div->SetMarkerStyle(21); purity_truth_pt_cut_x_div->SetMarkerColor(kBlue); purity_truth_pt_cut_x_div->Draw("E1 SAME"); }
  TLegend* leg_purity_pt = new TLegend(0.6,0.8,0.88,0.88);
  leg_purity_pt->SetTextSize(0.04);
  if (purity_match_x_div) leg_purity_pt->AddEntry(purity_match_x_div,"purity_match","lep");
  if (purity_truth_pt_cut_x_div) leg_purity_pt->AddEntry(purity_truth_pt_cut_x_div,"purity_truth_pt_cut","lep");
  leg_purity_pt->Draw();
  c_purity_pt->SaveAs(("test_figure/purity_vs_pt_"+bkg_cut_type+".pdf").c_str());

  // Purity vs calo ET (Y)
  TCanvas* c_purity_et = new TCanvas("c_purity_et","Purity vs calo ET",800,600);
  if (purity_match_y_div) {
    purity_match_y_div->SetLineColor(kRed); purity_match_y_div->SetMarkerStyle(20); purity_match_y_div->SetMarkerColor(kRed);
    purity_match_y_div->SetTitle("");
    purity_match_y_div->GetXaxis()->SetTitle("E_{T}^{calo} [GeV]");
    purity_match_y_div->GetYaxis()->SetTitle("Purity");
    purity_match_y_div->GetYaxis()->SetRangeUser(0.0,1.5);
    purity_match_y_div->Draw("E1");
  }
  if (purity_truth_pt_cut_y_div) { purity_truth_pt_cut_y_div->SetLineColor(kBlue); purity_truth_pt_cut_y_div->SetMarkerStyle(21); purity_truth_pt_cut_y_div->SetMarkerColor(kBlue); purity_truth_pt_cut_y_div->Draw("E1 SAME"); }
  TLegend* leg_purity_et = new TLegend(0.6,0.8,0.88,0.88);
  leg_purity_et->SetTextSize(0.04);
  if (purity_match_y_div) leg_purity_et->AddEntry(purity_match_y_div,"purity_match","lep");
  if (purity_truth_pt_cut_y_div) leg_purity_et->AddEntry(purity_truth_pt_cut_y_div,"purity_truth_pt_cut","lep");
  leg_purity_et->Draw();
  c_purity_et->SaveAs(("test_figure/purity_vs_et_"+bkg_cut_type+".pdf").c_str());

  // Efficiency vs jet pT (X)
  TCanvas* c_eff_pt = new TCanvas("c_eff_pt","Efficiency vs pT",800,600);
  if (efficiency_match_x_div) {
    efficiency_match_x_div->SetLineColor(kRed); efficiency_match_x_div->SetMarkerStyle(20); efficiency_match_x_div->SetMarkerColor(kRed);
    efficiency_match_x_div->SetTitle("");
    efficiency_match_x_div->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
    efficiency_match_x_div->GetYaxis()->SetTitle("Efficiency");
    efficiency_match_x_div->GetYaxis()->SetRangeUser(0.0,1.5);
    efficiency_match_x_div->Draw("E1");
  }
  if (efficiency_reco_pt_cut_x_div) { efficiency_reco_pt_cut_x_div->SetLineColor(kBlue); efficiency_reco_pt_cut_x_div->SetMarkerStyle(21); efficiency_reco_pt_cut_x_div->SetMarkerColor(kBlue); efficiency_reco_pt_cut_x_div->Draw("E1 SAME"); }
  if (efficiency_reco_bkg_cut_x_div) { efficiency_reco_bkg_cut_x_div->SetLineColor(kGreen+2); efficiency_reco_bkg_cut_x_div->SetMarkerStyle(22); efficiency_reco_bkg_cut_x_div->SetMarkerColor(kGreen+2); efficiency_reco_bkg_cut_x_div->Draw("E1 SAME"); }
  if (efficiency_match_reco_bkg_cut_x_div) { efficiency_match_reco_bkg_cut_x_div->SetLineColor(kMagenta); efficiency_match_reco_bkg_cut_x_div->SetMarkerStyle(23); efficiency_match_reco_bkg_cut_x_div->SetMarkerColor(kMagenta); efficiency_match_reco_bkg_cut_x_div->Draw("E1 SAME"); }
  if (efficiency_reco_pt_cut_reco_bkg_cut_x_div) { efficiency_reco_pt_cut_reco_bkg_cut_x_div->SetLineColor(kCyan); efficiency_reco_pt_cut_reco_bkg_cut_x_div->SetMarkerStyle(24); efficiency_reco_pt_cut_reco_bkg_cut_x_div->SetMarkerColor(kCyan); efficiency_reco_pt_cut_reco_bkg_cut_x_div->Draw("E1 SAME"); }
  if (efficiency_match_reco_pt_cut_reco_bkg_cut_x_div) { efficiency_match_reco_pt_cut_reco_bkg_cut_x_div->SetLineColor(kBlack); efficiency_match_reco_pt_cut_reco_bkg_cut_x_div->SetMarkerStyle(20); efficiency_match_reco_pt_cut_reco_bkg_cut_x_div->SetMarkerColor(kBlack); efficiency_match_reco_pt_cut_reco_bkg_cut_x_div->Draw("E1 SAME"); }
  if (efficiency_match_reco_pt_cut_x_div) { efficiency_match_reco_pt_cut_x_div->SetLineColor(kViolet); efficiency_match_reco_pt_cut_x_div->SetMarkerStyle(26); efficiency_match_reco_pt_cut_x_div->SetMarkerColor(kViolet); efficiency_match_reco_pt_cut_x_div->Draw("E1 SAME"); }
  TLegend* leg_eff_pt = new TLegend(0.55,0.75,0.92,0.92);
  leg_eff_pt->SetTextSize(0.035);
  if (efficiency_match_x_div) leg_eff_pt->AddEntry(efficiency_match_x_div,"efficiency_match","lep");
  if (efficiency_reco_pt_cut_x_div) leg_eff_pt->AddEntry(efficiency_reco_pt_cut_x_div,"efficiency_reco_pt_cut","lep");
  if (efficiency_reco_bkg_cut_x_div) leg_eff_pt->AddEntry(efficiency_reco_bkg_cut_x_div,"efficiency_reco_bkg_cut","lep");
  if (efficiency_match_reco_bkg_cut_x_div) leg_eff_pt->AddEntry(efficiency_match_reco_bkg_cut_x_div,"efficiency_match_reco_bkg_cut","lep");
  if (efficiency_reco_pt_cut_reco_bkg_cut_x_div) leg_eff_pt->AddEntry(efficiency_reco_pt_cut_reco_bkg_cut_x_div,"efficiency_reco_pt_cut_reco_bkg_cut","lep");
  if (efficiency_match_reco_pt_cut_reco_bkg_cut_x_div) leg_eff_pt->AddEntry(efficiency_match_reco_pt_cut_reco_bkg_cut_x_div,"efficiency_match_reco_pt_cut_reco_bkg_cut","lep");
  if (efficiency_match_reco_pt_cut_x_div) leg_eff_pt->AddEntry(efficiency_match_reco_pt_cut_x_div,"efficiency_match_reco_pt_cut","lep");
  leg_eff_pt->Draw();
  c_eff_pt->SaveAs(("test_figure/efficiency_vs_pt_"+bkg_cut_type+".pdf").c_str());

  // Efficiency vs calo ET (Y)
  TCanvas* c_eff_et = new TCanvas("c_eff_et","Efficiency vs calo ET",800,600);
  if (efficiency_match_y_div) {
    efficiency_match_y_div->SetLineColor(kRed); efficiency_match_y_div->SetMarkerStyle(20); efficiency_match_y_div->SetMarkerColor(kRed);
    efficiency_match_y_div->SetTitle("");
    efficiency_match_y_div->GetXaxis()->SetTitle("E_{T}^{calo} [GeV]");
    efficiency_match_y_div->GetYaxis()->SetTitle("Efficiency");
    efficiency_match_y_div->GetYaxis()->SetRangeUser(0.0,1.5);
    efficiency_match_y_div->Draw("E1");
  }
  if (efficiency_reco_pt_cut_y_div) { efficiency_reco_pt_cut_y_div->SetLineColor(kBlue); efficiency_reco_pt_cut_y_div->SetMarkerStyle(21); efficiency_reco_pt_cut_y_div->SetMarkerColor(kBlue); efficiency_reco_pt_cut_y_div->Draw("E1 SAME"); }
  if (efficiency_reco_bkg_cut_y_div) { efficiency_reco_bkg_cut_y_div->SetLineColor(kGreen+2); efficiency_reco_bkg_cut_y_div->SetMarkerStyle(22); efficiency_reco_bkg_cut_y_div->SetMarkerColor(kGreen+2); efficiency_reco_bkg_cut_y_div->Draw("E1 SAME"); }
  if (efficiency_match_reco_bkg_cut_y_div) { efficiency_match_reco_bkg_cut_y_div->SetLineColor(kMagenta); efficiency_match_reco_bkg_cut_y_div->SetMarkerStyle(23); efficiency_match_reco_bkg_cut_y_div->SetMarkerColor(kMagenta); efficiency_match_reco_bkg_cut_y_div->Draw("E1 SAME"); }
  if (efficiency_reco_pt_cut_reco_bkg_cut_y_div) { efficiency_reco_pt_cut_reco_bkg_cut_y_div->SetLineColor(kCyan); efficiency_reco_pt_cut_reco_bkg_cut_y_div->SetMarkerStyle(24); efficiency_reco_pt_cut_reco_bkg_cut_y_div->SetMarkerColor(kCyan); efficiency_reco_pt_cut_reco_bkg_cut_y_div->Draw("E1 SAME"); }
  if (efficiency_match_reco_pt_cut_reco_bkg_cut_y_div) { efficiency_match_reco_pt_cut_reco_bkg_cut_y_div->SetLineColor(kOrange+7); efficiency_match_reco_pt_cut_reco_bkg_cut_y_div->SetMarkerStyle(25); efficiency_match_reco_pt_cut_reco_bkg_cut_y_div->SetMarkerColor(kOrange+7); efficiency_match_reco_pt_cut_reco_bkg_cut_y_div->Draw("E1 SAME"); }
  if (efficiency_match_reco_pt_cut_y_div) { efficiency_match_reco_pt_cut_y_div->SetLineColor(kViolet); efficiency_match_reco_pt_cut_y_div->SetMarkerStyle(26); efficiency_match_reco_pt_cut_y_div->SetMarkerColor(kViolet); efficiency_match_reco_pt_cut_y_div->Draw("E1 SAME"); }
  TLegend* leg_eff_et = new TLegend(0.55,0.75,0.92,0.92);
  leg_eff_et->SetTextSize(0.035);
  if (efficiency_match_y_div) leg_eff_et->AddEntry(efficiency_match_y_div,"match","lep");
  if (efficiency_reco_pt_cut_y_div) leg_eff_et->AddEntry(efficiency_reco_pt_cut_y_div,"reco_pt_cut","lep");
  if (efficiency_reco_bkg_cut_y_div) leg_eff_et->AddEntry(efficiency_reco_bkg_cut_y_div,"reco_bkg_cut","lep");
  if (efficiency_match_reco_bkg_cut_y_div) leg_eff_et->AddEntry(efficiency_match_reco_bkg_cut_y_div,"match && reco_bkg_cut","lep");
  if (efficiency_reco_pt_cut_reco_bkg_cut_y_div) leg_eff_et->AddEntry(efficiency_reco_pt_cut_reco_bkg_cut_y_div,"reco_pt_cut && reco_bkg_cut","lep");
  if (efficiency_match_reco_pt_cut_reco_bkg_cut_y_div) leg_eff_et->AddEntry(efficiency_match_reco_pt_cut_reco_bkg_cut_y_div,"match && reco_pt_cut && reco_bkg_cut","lep");
  if (efficiency_match_reco_pt_cut_y_div) leg_eff_et->AddEntry(efficiency_match_reco_pt_cut_y_div,"match && reco_pt_cut","lep");
  leg_eff_et->Draw();
  c_eff_et->SaveAs(("test_figure/efficiency_vs_et_"+bkg_cut_type+".pdf").c_str());

  // Write 1D histograms to output file
  TFile* f_out = new TFile(outfilename.c_str(), "RECREATE");
  if (f_out && f_out->IsOpen()) {
    if (purity_match_x_div) purity_match_x_div->Write();
    if (purity_truth_pt_cut_x_div) purity_truth_pt_cut_x_div->Write();
    if (purity_match_y_div) purity_match_y_div->Write();
    if (purity_truth_pt_cut_y_div) purity_truth_pt_cut_y_div->Write();
    if (efficiency_match_x_div) efficiency_match_x_div->Write();
    if (efficiency_reco_pt_cut_x_div) efficiency_reco_pt_cut_x_div->Write();
    if (efficiency_reco_bkg_cut_x_div) efficiency_reco_bkg_cut_x_div->Write();
    if (efficiency_match_reco_bkg_cut_x_div) efficiency_match_reco_bkg_cut_x_div->Write();
    if (efficiency_reco_pt_cut_reco_bkg_cut_x_div) efficiency_reco_pt_cut_reco_bkg_cut_x_div->Write();
    if (efficiency_match_reco_pt_cut_reco_bkg_cut_x_div) efficiency_match_reco_pt_cut_reco_bkg_cut_x_div->Write();
    if (efficiency_match_reco_pt_cut_x_div) efficiency_match_reco_pt_cut_x_div->Write();
    if (efficiency_match_y_div) efficiency_match_y_div->Write();
    if (efficiency_reco_pt_cut_y_div) efficiency_reco_pt_cut_y_div->Write();
    if (efficiency_reco_bkg_cut_y_div) efficiency_reco_bkg_cut_y_div->Write();
    if (efficiency_match_reco_bkg_cut_y_div) efficiency_match_reco_bkg_cut_y_div->Write();
    if (efficiency_reco_pt_cut_reco_bkg_cut_y_div) efficiency_reco_pt_cut_reco_bkg_cut_y_div->Write();
    if (efficiency_match_reco_pt_cut_reco_bkg_cut_y_div) efficiency_match_reco_pt_cut_reco_bkg_cut_y_div->Write();
    if (efficiency_match_reco_pt_cut_y_div) efficiency_match_reco_pt_cut_y_div->Write();
    f_out->Close();
  } else {
    std::cout << "Warning: could not open output file " << outfilename << std::endl;
  }

}