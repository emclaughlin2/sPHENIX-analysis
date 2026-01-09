#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include <TMatrixD.h>
#include "unfold_Def.h"
//#include "/sphenix/user/hanpuj/CaloDataAna24_skimmed/src/draw_template.C" 

int colors[] = {kRed, kBlue, kGreen+2, kMagenta, kOrange, kCyan+1,
                    kViolet, kSpring+8, kAzure-4, kPink+8,
                    kYellow+2, kTeal+3, kGray+2, kBlue+3,
                    kOrange+7, kGreen+4, kMagenta+2, kRed+3,
                    kAzure+6, kBlack, kRed+1, kBlue+1, kGreen+1, 
                    kMagenta+1, kOrange+1, 
                    kCyan+2, kViolet+1, kPink+1, kYellow+3};

void make_simplots() {
  gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
  gROOT->ProcessLine("SetsPhenixStyle()");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  //********** General Set up **********//
  const float PI = TMath::Pi();
  const float jet_radius = 0.4;

  std::vector<TH1F*> h_input;
  std::vector<int> color;
  std::vector<int> markerstyle;
  std::vector<std::string> text;
  std::vector<std::string> legend;
 
  //********** Files **********//
  TFile *f_sim = new TFile("analysis_sim_run21_output/sim_plots_output_dijet_sim_iter_2.root", "READ");
  TFile *f_sim_towers = new TFile("analysis_sim_run21_output/sim_plots_towers_output_dijet_sim_iter_1.root", "READ");
  TFile *f_sim_herwig = new TFile("analysis_sim_run21_output/sim_plots_output_dijet_sim_iter_1_herwig.root", "READ");

  //TFile *f_sim = new TFile("analysis_sim_run21_output/sim_plots_output_dijet_sim_iter_1_finebins.root", "READ");
  //TFile *f_sim_towers = new TFile("analysis_sim_run21_output/sim_plots_towers_output_dijet_sim_iter_1.root", "READ");
  //TFile *f_sim_herwig = new TFile("analysis_sim_run21_output/output_dijet_sim_iter_1_herwig.root", "READ");

  // histograms for ET distribution plots 
  TH1D* h_et[4];
  TH1D* h_et_bigbin[4];
  h_et_bigbin[0] = (TH1D*)f_sim->Get("h_et_transverse_record"); 
  h_et_bigbin[1] = (TH1D*)f_sim_towers->Get("h_et_transverse_record");
  h_et_bigbin[2] = (TH1D*)f_sim->Get("h_et_truth_transverse_record"); 
  h_et_bigbin[3] = (TH1D*)f_sim->Get("h_reco_thres_et_truth_transverse_record");
  h_et[0] = (TH1D*)f_sim->Get("h_nw_et_transverse_record"); 
  h_et[1] = (TH1D*)f_sim_towers->Get("h_nw_et_transverse_record");
  h_et[2] = (TH1D*)f_sim->Get("h_nw_et_truth_transverse_record"); 
  h_et[3] = (TH1D*)f_sim->Get("h_reco_thres_nw_et_truth_transverse_record");

  for (int i = 0; i < 4; i++) { h_et[i]->Rebin(16); }
  vector<double> bins; vector<double> bin_size;
  for (int i = 1; i <= h_et_bigbin[0]->GetNbinsX() + 1; i++) { bins.push_back(h_et_bigbin[0]->GetBinLowEdge(i)); }
  for (int i = 1; i < bins.size(); i++) { bin_size.push_back(bins[i] - bins[i-1]); }
  for (int b = 0; b < 2; b++) {
    for (int i = 1; i < h_et_bigbin[b]->GetNbinsX() + 1; i++) {
      h_et_bigbin[b]->SetBinContent(i,h_et_bigbin[b]->GetBinContent(i)/bin_size[i-1]);
      h_et_bigbin[b]->SetBinError(i,h_et_bigbin[b]->GetBinError(i)/bin_size[i-1]);
    }
  }
  bins.clear(); bin_size.clear();
  for (int i = 1; i <= h_et_bigbin[2]->GetNbinsX() + 1; i++) { bins.push_back(h_et_bigbin[2]->GetBinLowEdge(i)); }
  for (int i = 1; i < bins.size(); i++) { bin_size.push_back(bins[i] - bins[i-1]); }
  for (int b = 2; b < 4; b++) {
    for (int i = 1; i < h_et_bigbin[b]->GetNbinsX() + 1; i++) {
      h_et_bigbin[b]->SetBinContent(i,h_et_bigbin[b]->GetBinContent(i)/bin_size[i-1]);
      h_et_bigbin[b]->SetBinError(i,h_et_bigbin[b]->GetBinError(i)/bin_size[i-1]);
    }
  }
  
  for (int i = 0; i < 4; i++) {
    h_et[i]->Scale(1.0/h_et[0]->GetBinWidth(1));
  }

  float h_et_mean[4];
  for (int i = 0; i < 4; i++) {
    h_et_mean[i] = h_et[i]->GetMean();
  }
  
  std::vector<std::string> et_dist_leg_tag = {"<#SigmaE_{T}^{topo}> = ","<#SigmaE_{T}^{towers}> = ","<#SigmaE_{T}^{truth}> = ","<#SigmaE_{T}^{truth,topo reco}> = "};
  TCanvas* et_dist_canvas = new TCanvas("et_dist_canvas","et_dist_canvas", 800, 600);
  et_dist_canvas->SetLogy(1);
  TLegend* et_dist_legend = new TLegend(.54,.65,.9,.92);
  et_dist_legend->SetTextSize(0.038);
  et_dist_legend->AddEntry("","#bf{Transverse Region}","");
  for (int i = 0; i < 4; i++) {
    h_et[i]->SetLineColor(colors[i]);
    h_et[i]->SetMarkerColor(colors[i]);
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2) << h_et_mean[i];
    std::string str_value = oss.str();
    et_dist_legend->AddEntry(h_et[i],(et_dist_leg_tag[i]+str_value+" GeV").c_str(),"le");
    if (i == 0) {
      h_et[i]->GetXaxis()->SetRangeUser(-5,20);
      h_et[i]->GetYaxis()->SetRangeUser(0.1,200000);
      h_et[i]->SetYTitle("dN/d#SigmaE_{T} [GeV^{-1}]");
      h_et[i]->Draw("hist");
    } else {
      h_et[i]->Draw("hist,same");
    }
  }
  et_dist_legend->Draw("same");
  TLatex latex;
  latex.SetTextSize(0.038);
  latex.SetTextAlign(13);  //align at top
  latex.DrawLatex(-4.2,100000,"#bf{#it{sPHENIX}} Simulation Preliminary");
  TLatex latex2;
  latex2.SetTextSize(0.038);
  latex2.SetTextAlign(13);  //align at top
  latex2.DrawLatex(-4.2,50000,"Pythia8 200 GeV p+p");
  TLine* line[4];
  for (int i = 0; i < 4; i++) {
    line[i] = new TLine(h_et_mean[i],0.1,h_et_mean[i],h_et[i]->GetBinContent(h_et[i]->FindBin(h_et_mean[i])));
    line[i]->SetLineColor(colors[i]);
    line[i]->SetLineWidth(2);
    line[i]->SetLineStyle(2);
    line[i]->Draw("same");
  }
  et_dist_canvas->SaveAs("sim_plots_run21/et_distribution.png");
  et_dist_canvas->SaveAs("sim_plots_run21/et_distribution.pdf");

/*
  TCanvas* et_dist_bigbin_canvas = new TCanvas("et_dist_bigbin_canvas","et_dist_bigbin_canvas", 700, 600);
  TLegend* et_dist_bigbin_legend = new TLegend(.4,.5,.9,.9);
  et_dist_bigbin_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Preliminary","");
  et_dist_bigbin_legend->AddEntry("","Pythia8 200 GeV p+p","");
  et_dist_bigbin_legend->SetFillStyle(0);
  et_dist_bigbin_legend->SetTextSize(0.04);
  for (int i = 0; i < 4; i++) {
    h_et_bigbin[i]->SetLineColor(i+1);
    h_et_bigbin[i]->SetMarkerColor(i+1);
    et_dist_bigbin_legend->AddEntry(h_et_bigbin[i],et_dist_leg_tag[i].c_str(),"lep");
  }
  h_et_bigbin[0]->Draw("hist");
  h_et_bigbin[1]->Draw("hist,same");
  h_et_bigbin[2]->Draw("hist,same");
  h_et_bigbin[3]->Draw("hist,same");
  et_dist_bigbin_legend->Draw("same");
  et_dist_bigbin_canvas->SaveAs("sim_plots_run21/et_distribution_big_bins.png");
  */
  /*
  // histograms for truth UE as a function of jet pT 
  TH2D* h_ue_pt_transverse_record = (TH2D*)f_sim->Get("h_ue_pt_transverse_record"); 
  TH2D* h_truth_ue_pt_transverse_record = (TH2D*)f_sim->Get("h_ue_pt_truth_transverse_record"); 
  TH2D* h_topo_thres_truth_ue_pt_transverse_record = (TH2D*)f_sim->Get("h_reco_thres_ue_pt_truth_transverse_record");
  TH2D* h_herwig_ue_pt_transverse_record = (TH2D*)f_sim_herwig->Get("h_ue_pt_truth_transverse_record");
  TH1D* h_pythia_ue_pt = (TH1D*)h_truth_ue_pt_transverse_record->ProfileX("h_pythia_ue_pt",1,16);
  TH1D* h_herwig_ue_pt = (TH1D*)h_herwig_ue_pt_transverse_record->ProfileX("h_herwig_ue_pt",1,16);

  TH1D* h_pythia = (TH1D*)h_pythia_ue_pt->Clone("h_pythia");
  TH1D* h_herwig = (TH1D*)h_herwig_ue_pt->Clone("h_herwig");
  TH2D* h_pythia_2D = (TH2D*)h_truth_ue_pt_transverse_record->Clone("h_pythia_2D");
  TH2D* h_herwig_2D = (TH2D*)h_herwig_ue_pt_transverse_record->Clone("h_herwig_2D");

  h_pythia_ue_pt->Scale(1.0/(secteta*sectphi));
  h_herwig_ue_pt->Scale(1.0/(secteta*sectphi));

  TCanvas* truth_ue_vs_pt_canvas = new TCanvas("truth_ue_vs_pt_canvas","truth_ue_vs_pt_canvas", 750, 700);
  TLegend* truth_ue_vs_pt_legend = new TLegend(.2,.2,.5,.5);
  truth_ue_vs_pt_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Preliminary","");
  truth_ue_vs_pt_legend->AddEntry("","anti-k_{t}#it{R}=0.4, |#eta_{jet}| < 0.7, p_{T}^{lead jet} > 17 GeV","");
  truth_ue_vs_pt_legend->AddEntry("","p_{T}^{sublead jet} > 0.3 p_{T}^{lead jet}, |#Delta#phi| > 3#pi/4","");
  truth_ue_vs_pt_legend->AddEntry(h_pythia_ue_pt,"Pythia8 200 GeV p+p","lep");
  truth_ue_vs_pt_legend->AddEntry(h_herwig_ue_pt,"Herwig7 200 GeV p+p","lep");
  truth_ue_vs_pt_legend->SetTextSize(0.038);
  h_pythia_ue_pt->SetLineColor(2);
  h_pythia_ue_pt->SetMarkerColor(2);
  h_herwig_ue_pt->SetLineColor(4);
  h_herwig_ue_pt->SetMarkerColor(4);
  h_pythia_ue_pt->GetYaxis()->SetRangeUser(0,0.8);
  h_pythia_ue_pt->SetYTitle("Truth <#SigmaE_{T}/#delta#eta#delta#phi> [GeV]");
  h_pythia_ue_pt->SetXTitle("p_{T}^{truth lead jet} [GeV]");
  h_pythia_ue_pt->GetXaxis()->SetRangeUser(17,62);
  h_pythia_ue_pt->Draw();
  h_herwig_ue_pt->Draw("same");
  truth_ue_vs_pt_legend->Draw("same");
  truth_ue_vs_pt_canvas->SaveAs("sim_plots_run21/truth_ue_vs_pt_finebins.png");
  truth_ue_vs_pt_canvas->SaveAs("sim_plots_run21/truth_ue_vs_pt_finebins.pdf");

  TCanvas* raw_truth_pythia_canvas = new TCanvas("raw_truth_pythia_canvas","raw_truth_pythia_canvas", 800, 600);
  raw_truth_pythia_canvas->SetTopMargin(0.18);
  raw_truth_pythia_canvas->SetLogz(1);
  TLegend* raw_truth_pythia_legend = new TLegend(.05,.82,.55,.99);
  raw_truth_pythia_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Internal","");
  raw_truth_pythia_legend->AddEntry("","Pythia8 200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
  raw_truth_pythia_legend->AddEntry("","p_{T}^{sublead jet} > 0.3 p_{T}^{lead jet}, |#Delta#phi| > 3#pi/4","");
  raw_truth_pythia_legend->SetFillStyle(0);
  raw_truth_pythia_legend->SetTextSize(0.035);
  h_pythia_2D->SetStats(0);
  h_pythia_2D->GetZaxis()->SetRangeUser(0.001,10000);
  h_pythia_2D->GetYaxis()->SetRangeUser(0,20);
  h_pythia_2D->GetXaxis()->SetRangeUser(17,62);
  h_pythia_2D->SetYTitle("Truth #SigmaE_{T} [GeV]");
  h_pythia_2D->SetXTitle("p_{T}^{truth lead jet} [GeV]");
  h_pythia_2D->Draw("colz");
  h_pythia->SetLineColor(4);
  h_pythia->SetMarkerColor(4);
  h_pythia->SetLineWidth(2);
  h_pythia->Draw("hist,same");
  raw_truth_pythia_canvas->SetRightMargin(0.2);
  raw_truth_pythia_canvas->Update();
  TPaletteAxis *raw_truth_pythia_palette = (TPaletteAxis*)h_pythia_2D->GetListOfFunctions()->FindObject("palette");
  raw_truth_pythia_palette->SetX1NDC(0.82);
  raw_truth_pythia_palette->SetX2NDC(0.9);
  raw_truth_pythia_legend->Draw("same");
  raw_truth_pythia_canvas->Update();
  raw_truth_pythia_canvas->SaveAs("sim_plots_run21/unscaled_truth_pythia_ue_pt_2D.png");
  raw_truth_pythia_canvas->SaveAs("sim_plots_run21/unscaled_truth_pythia_ue_pt_2D.pdf");

  TCanvas* raw_truth_herwig_canvas = new TCanvas("raw_truth_herwig_canvas","raw_truth_herwig_canvas", 800, 600);
  raw_truth_herwig_canvas->SetTopMargin(0.18);
  raw_truth_herwig_canvas->SetLogz(1);
  TLegend* raw_truth_herwig_legend = new TLegend(.05,.82,.55,.99);
  raw_truth_herwig_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Internal","");
  raw_truth_herwig_legend->AddEntry("","Herwig7 200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
  raw_truth_herwig_legend->AddEntry("","p_{T}^{sublead jet} > 0.3 p_{T}^{lead jet}, |#Delta#phi| > 3#pi/4","");
  raw_truth_herwig_legend->SetFillStyle(0);
  raw_truth_herwig_legend->SetTextSize(0.035);
  h_herwig_2D->SetStats(0);
  h_herwig_2D->GetZaxis()->SetRangeUser(0.001,10000);
  h_herwig_2D->GetYaxis()->SetRangeUser(0,20);
  h_herwig_2D->GetXaxis()->SetRangeUser(17,62);
  h_herwig_2D->SetYTitle("Truth #SigmaE_{T} [GeV]");
  h_herwig_2D->SetXTitle("p_{T}^{truth lead jet} [GeV]");
  h_herwig_2D->Draw("colz");
  h_herwig->SetLineColor(4);
  h_herwig->SetMarkerColor(4);
  h_herwig->SetLineWidth(2);
  h_herwig->Draw("hist,same");
  raw_truth_herwig_canvas->SetRightMargin(0.2);
  raw_truth_herwig_canvas->Update();
  TPaletteAxis *raw_truth_herwig_palette = (TPaletteAxis*)h_herwig_2D->GetListOfFunctions()->FindObject("palette");
  raw_truth_herwig_palette->SetX1NDC(0.82);
  raw_truth_herwig_palette->SetX2NDC(0.9);
  raw_truth_herwig_legend->Draw("same");
  raw_truth_herwig_canvas->Update();
  raw_truth_herwig_canvas->SaveAs("sim_plots_run21/unscaled_truth_herwig_ue_pt_2D.png");
  raw_truth_herwig_canvas->SaveAs("sim_plots_run21/unscaled_truth_herwig_ue_pt_2D.pdf");

  for (int i = 1; i < h_truth_ue_pt_transverse_record->GetNbinsX() + 1; i++) {
    float x_tot = 0.0;
    float hx_tot = 0.0;
    for (int j = 1; j < h_truth_ue_pt_transverse_record->GetNbinsY() + 1; j++) {
      x_tot += h_truth_ue_pt_transverse_record->GetBinContent(i,j);
      hx_tot += h_herwig_ue_pt_transverse_record->GetBinContent(i,j);
    }
    for (int j = 1; j < h_truth_ue_pt_transverse_record->GetNbinsY() + 1; j++) {
      h_truth_ue_pt_transverse_record->SetBinContent(i,j,h_truth_ue_pt_transverse_record->GetBinContent(i,j)/x_tot);
      h_herwig_ue_pt_transverse_record->SetBinContent(i,j,h_herwig_ue_pt_transverse_record->GetBinContent(i,j)/hx_tot);
      h_truth_ue_pt_transverse_record->SetBinError(i,j,h_truth_ue_pt_transverse_record->GetBinError(i,j)/x_tot);
      h_herwig_ue_pt_transverse_record->SetBinError(i,j,h_herwig_ue_pt_transverse_record->GetBinError(i,j)/hx_tot);
    }
  }

  for (int i = 1; i < h_truth_ue_pt_transverse_record->GetNbinsX() + 1; i++) {
    for (int j = 1; j < h_truth_ue_pt_transverse_record->GetNbinsY() + 1; j++) {
      float ybin = h_truth_ue_pt_transverse_record->GetYaxis()->GetBinWidth(j);
      float hybin = h_herwig_ue_pt_transverse_record->GetYaxis()->GetBinWidth(j);
      h_truth_ue_pt_transverse_record->SetBinContent(i,j,h_truth_ue_pt_transverse_record->GetBinContent(i,j)/ybin);
      h_herwig_ue_pt_transverse_record->SetBinContent(i,j,h_herwig_ue_pt_transverse_record->GetBinContent(i,j)/hybin);
      h_truth_ue_pt_transverse_record->SetBinError(i,j,h_truth_ue_pt_transverse_record->GetBinError(i,j)/ybin);
      h_herwig_ue_pt_transverse_record->SetBinError(i,j,h_herwig_ue_pt_transverse_record->GetBinError(i,j)/hybin);
    }
  }

  TCanvas* truth_pythia_canvas = new TCanvas("truth_pythia_canvas","truth_pythia_canvas", 800, 600);
  truth_pythia_canvas->SetTopMargin(0.18);
  truth_pythia_canvas->SetLogz(1);
  TLegend* truth_pythia_legend = new TLegend(.05,.82,.55,.99);
  truth_pythia_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Internal","");
  truth_pythia_legend->AddEntry("","Pythia8 200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
  truth_pythia_legend->AddEntry("","p_{T}^{sublead jet} > 0.3 p_{T}^{lead jet}, |#Delta#phi| > 3#pi/4","");
  truth_pythia_legend->SetFillStyle(0);
  truth_pythia_legend->SetTextSize(0.035);
  h_truth_ue_pt_transverse_record->SetStats(0);
  h_truth_ue_pt_transverse_record->GetZaxis()->SetRangeUser(0.00005,2);
  h_truth_ue_pt_transverse_record->GetYaxis()->SetRangeUser(0,20);
  h_truth_ue_pt_transverse_record->GetXaxis()->SetRangeUser(17,62);
  h_truth_ue_pt_transverse_record->SetYTitle("Truth #SigmaE_{T} [GeV]");
  h_truth_ue_pt_transverse_record->SetXTitle("p_{T}^{truth lead jet} [GeV]");
  h_truth_ue_pt_transverse_record->Draw("colz");
  h_pythia->SetLineColor(4);
  h_pythia->SetMarkerColor(4);
  h_pythia->SetLineWidth(2);
  h_pythia->Draw("hist,same");
  truth_pythia_canvas->SetRightMargin(0.2);
  truth_pythia_canvas->Update();
  TPaletteAxis *truth_pythia_palette = (TPaletteAxis*)h_truth_ue_pt_transverse_record->GetListOfFunctions()->FindObject("palette");
  truth_pythia_palette->SetX1NDC(0.82);
  truth_pythia_palette->SetX2NDC(0.9);
  truth_pythia_legend->Draw("same");
  truth_pythia_canvas->Update();
  truth_pythia_canvas->SaveAs("sim_plots_run21/truth_pythia_ue_pt_2D.png");
  truth_pythia_canvas->SaveAs("sim_plots_run21/truth_pythia_ue_pt_2D.pdf");

  TCanvas* truth_herwig_canvas = new TCanvas("truth_herwig_canvas","truth_herwig_canvas", 800, 600);
  truth_herwig_canvas->SetTopMargin(0.18);
  truth_herwig_canvas->SetLogz(1);
  TLegend* truth_herwig_legend = new TLegend(.05,.82,.55,.99);
  truth_herwig_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Internal","");
  truth_herwig_legend->AddEntry("","Herwig7 200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
  truth_herwig_legend->AddEntry("","p_{T}^{sublead jet} > 0.3 p_{T}^{lead jet}, |#Delta#phi| > 3#pi/4","");
  truth_herwig_legend->SetFillStyle(0);
  truth_herwig_legend->SetTextSize(0.035);
  h_herwig_ue_pt_transverse_record->SetStats(0);
  h_herwig_ue_pt_transverse_record->GetZaxis()->SetRangeUser(0.00005,2);
  h_herwig_ue_pt_transverse_record->GetYaxis()->SetRangeUser(0,20);
  h_herwig_ue_pt_transverse_record->GetXaxis()->SetRangeUser(17,62);
  h_herwig_ue_pt_transverse_record->SetYTitle("Truth #SigmaE_{T} [GeV]");
  h_herwig_ue_pt_transverse_record->SetXTitle("p_{T}^{truth lead jet} [GeV]");
  h_herwig_ue_pt_transverse_record->Draw("colz");
  h_herwig->SetLineColor(4);
  h_herwig->SetMarkerColor(4);
  h_herwig->SetLineWidth(2);
  h_herwig->Draw("hist,same");
  truth_herwig_canvas->SetRightMargin(0.2);
  truth_herwig_canvas->Update();
  TPaletteAxis *truth_herwig_palette = (TPaletteAxis*)h_herwig_ue_pt_transverse_record->GetListOfFunctions()->FindObject("palette");
  truth_herwig_palette->SetX1NDC(0.82);
  truth_herwig_palette->SetX2NDC(0.9);
  truth_herwig_legend->Draw("same");
  truth_herwig_canvas->Update();
  truth_herwig_canvas->SaveAs("sim_plots_run21/truth_herwig_ue_pt_2D.png");
  truth_herwig_canvas->SaveAs("sim_plots_run21/truth_herwig_ue_pt_2D.pdf");
*/
  // histograms for unfolding plots
  RooUnfoldResponse* h[4];
  h[0] = (RooUnfoldResponse*)f_sim->Get("h_jetpt_respmatrix");
  h[1] = (RooUnfoldResponse*)f_sim->Get("h_caloet_respmatrix");
  h[2] = (RooUnfoldResponse*)f_sim->Get("h_jetpt_respmatrix_counts");
  h[3] = (RooUnfoldResponse*)f_sim->Get("h_caloet_respmatrix_counts");
  TH2D* h_jet_pt_respmatrix_calib_dijet = (TH2D*)h[0]->Hresponse(); 
  TH2D* h_calo_et_respmatrix_calib_dijet = (TH2D*)h[1]->Hresponse(); 
  TH2D* h_nw_jet_pt_respmatrix_calib_dijet = (TH2D*)h[2]->Hresponse(); 
  TH2D* h_nw_calo_et_respmatrix_calib_dijet = (TH2D*)h[3]->Hresponse(); 
  TH2D* h_truth_calib_dijet = (TH2D*)f_sim->Get("h_truth_calib_dijet"); 
  TH2D* h_measure_calib_dijet = (TH2D*)f_sim->Get("h_measure_calib_dijet"); 
  RooUnfoldResponse* h_respmatrix_calib_dijet = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_trim_10");  
  RooUnfoldResponse* h_counts_calib_dijet = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_counts"); 
  const TMatrixD& count_matrix = h_counts_calib_dijet->Mresponse(false);
  const TMatrixD& resp_matrix = h_respmatrix_calib_dijet->Mresponse(false);  // use normalized=true if needed

  // histograms for jet QA
  TH1D* h_lead_spectra_record = (TH1D*)f_sim->Get("h_lead_spectra_record"); 
  TH1D* h_sub_spectra_record = (TH1D*)f_sim->Get("h_sub_spectra_record"); 
  TH1D* h_truth_lead_spectra_record = (TH1D*)f_sim->Get("h_lead_truth_spectra_record"); 
  TH1D* h_truth_sub_spectra_record = (TH1D*)f_sim->Get("h_sub_truth_spectra_record"); 
  
  int nxr = h_measure_calib_dijet->GetNbinsX(); // reco x bins
  int nyr = h_measure_calib_dijet->GetNbinsY(); // reco y bins
  int nxt = h_truth_calib_dijet->GetNbinsX(); // truth x bins
  int nyt = h_truth_calib_dijet->GetNbinsY(); // truth y bins

  TH2D* hist2d = new TH2D("response2D", "Response Matrix;Reco Bin;Truth Bin", nxr*nyr, 0, nxr*nyr, nxt*nyt, 0, nxt*nyt);
  TH2D* counts_hist2d = new TH2D("counts_response2D", "Response Matrix;Reco Bin;Truth Bin", nxr*nyr, 0, nxr*nyr, nxt*nyt, 0, nxt*nyt);
  TH2D* jet_pt_respmatrix = new TH2D("jet_pt_respmatrix",";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", calibnpt, calibptbins, truthnpt, truthptbins);
  TH2D* calo_et_respmatrix = new TH2D("calo_et_respmatrix",";#SigmaE_{T}^{Reco} [GeV];#SigmaE_{T}^{Truth} [GeV]", calibnet, calibetbins, truthnet, truthetbins);
  TH2D* counts_jet_pt_respmatrix = new TH2D("counts_jet_pt_respmatrix",";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", calibnpt, calibptbins, truthnpt, truthptbins);
  TH2D* counts_calo_et_respmatrix = new TH2D("counts_calo_et_respmatrix",";#SigmaE_{T}^{Reco} [GeV];#SigmaE_{T}^{Truth} [GeV]", calibnet, calibetbins, truthnet, truthetbins);

  for (int ixr = 0; ixr < nxr; ++ixr) {
    for (int iyr = 0; iyr < nyr; ++iyr) {
      int recoBin = iyr * nxr + ixr;
      for (int ixt = 0; ixt < nxt; ++ixt) {
        for (int iyt = 0; iyt < nyt; ++iyt) {
          int truthBin = iyt * nxt + ixt;
          hist2d->SetBinContent(recoBin + 1, truthBin + 1, resp_matrix(recoBin, truthBin));
          counts_hist2d->SetBinContent(recoBin + 1, truthBin + 1, count_matrix(recoBin, truthBin));
        }
      }
    }
  }

  // Project pt_reco vs pt_truth
  for (int ixr = 0; ixr < nxr; ++ixr) {
    for (int ixt = 0; ixt < nxt; ++ixt) {
      double sum = 0.0;
      double rw_sum = 0.0;
      for (int iyr = 0; iyr < nyr; ++iyr) {
        for (int iyt = 0; iyt < nyt; ++iyt) {
          int row = iyr * nxr + ixr;  // flatten (pt_reco, et_reco)
          int col = iyt * nxt + ixt;  // flatten (pt_truth, et_truth)
          sum += resp_matrix(row, col);
          rw_sum += count_matrix(row, col);
        }
      }
      jet_pt_respmatrix->SetBinContent(ixr + 1, ixt + 1, sum);
      counts_jet_pt_respmatrix->SetBinContent(ixr + 1, ixt + 1, rw_sum);
    }
  }

  // Project et_reco vs et_truth
  for (int iyr = 0; iyr < nyr; ++iyr) {
    for (int iyt = 0; iyt < nyt; ++iyt) {
      double sum = 0.0;
      double rw_sum = 0.0;
      for (int ixr = 0; ixr < nxr; ++ixr) {
        for (int ixt = 0; ixt < nxt; ++ixt) {
          int row = iyr * nxr + ixr;
          int col = iyt * nxt + ixt;
          sum += resp_matrix(row, col);
          rw_sum += count_matrix(row, col);
        }
      }
      calo_et_respmatrix->SetBinContent(iyr + 1, iyt + 1, sum);
      counts_calo_et_respmatrix->SetBinContent(iyr + 1, iyt + 1, rw_sum);
    }
  }

  
  TCanvas* resp_canvas = new TCanvas("resp_canvas","resp_canvas", 700, 650);
  resp_canvas->SetTopMargin(0.18);
  TLegend* resp_legend = new TLegend(.05,.82,.55,.99);
  resp_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Preliminary Pythia8 200 GeV p+p","");
  resp_legend->AddEntry("","anti-k_{t}#it{R}=0.4, |#eta_{jet}| < 0.7, p_{T}^{truth lead jet} > 17 GeV","");
  resp_legend->AddEntry("","p_{T}^{reco lead jet} > 21 GeV, p_{T}^{sublead jet} > 0.3 p_{T}^{lead jet}, |#Delta#phi| > 3#pi/4","");
  resp_legend->SetFillStyle(0);
  resp_legend->SetTextSize(0.035);
  resp_canvas->SetLogz(1);
  resp_canvas->SetRightMargin(0.2);
  hist2d->GetZaxis()->SetRangeUser(0.000001,800);
  hist2d->Draw("colz");
  resp_canvas->Update();
  TPaletteAxis *resp_palette = (TPaletteAxis*)hist2d->GetListOfFunctions()->FindObject("palette");
  resp_palette->SetX1NDC(0.82);
  resp_palette->SetX2NDC(0.9);
  resp_legend->Draw("same");
  resp_canvas->Update();
  resp_canvas->SaveAs("sim_plots_run21/raw_respmatrix.png");
  resp_canvas->SaveAs("sim_plots_run21/raw_respmatrix.pdf");

  TCanvas* counts_resp_canvas = new TCanvas("counts_resp_canvas","counts_resp_canvas", 700, 650);
  counts_resp_canvas->SetTopMargin(0.18);
  TLegend* counts_resp_legend = new TLegend(.05,.82,.55,.99);
  counts_resp_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Internal","");
  counts_resp_legend->AddEntry("","Pythia8 200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
  counts_resp_legend->AddEntry("","p_{T}^{sublead jet} > 0.3 p_{T}^{lead jet}, |#Delta#phi| > 3#pi/4","");
  counts_resp_legend->SetFillStyle(0);
  counts_resp_legend->SetTextSize(0.035);
  counts_resp_canvas->SetLogz(1);
  counts_resp_canvas->SetRightMargin(0.2);
  counts_hist2d->GetZaxis()->SetRangeUser(0,10000);
  counts_hist2d->Draw("colz");
  counts_resp_canvas->Update();
  TPaletteAxis *counts_resp_palette = (TPaletteAxis*)counts_hist2d->GetListOfFunctions()->FindObject("palette");
  counts_resp_palette->SetX1NDC(0.82);
  counts_resp_palette->SetX2NDC(0.9);
  counts_resp_legend->Draw("same");
  counts_resp_canvas->Update();
  counts_resp_canvas->SaveAs("sim_plots_run21/counts_respmatrix.png");
  counts_resp_canvas->SaveAs("sim_plots_run21/counts_respmatrix.pdf");

  TCanvas* jet_pt_resp_canvas = new TCanvas("jet_pt_resp_canvas","resp_canvas", 700, 650);
  jet_pt_resp_canvas->SetTopMargin(0.18);
  TLegend* jet_pt_resp_legend = new TLegend(.05,.82,.55,.99);
  jet_pt_resp_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Preliminary","");
  jet_pt_resp_legend->AddEntry("","Pythia8 200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
  jet_pt_resp_legend->AddEntry("","p_{T}^{sublead jet} > 0.3 p_{T}^{lead jet}, |#Delta#phi| > 3#pi/4","");
  jet_pt_resp_legend->SetFillStyle(0);
  jet_pt_resp_legend->SetTextSize(0.035);
  jet_pt_resp_canvas->SetLogz(1);
  jet_pt_resp_canvas->SetRightMargin(0.2);
  jet_pt_respmatrix->GetZaxis()->SetRangeUser(0.0001,3000);
  jet_pt_respmatrix->Draw("colz");
  jet_pt_resp_canvas->Update();
  TPaletteAxis *jet_pt_resp_palette = (TPaletteAxis*)jet_pt_respmatrix->GetListOfFunctions()->FindObject("palette");
  jet_pt_resp_palette->SetX1NDC(0.82);
  jet_pt_resp_palette->SetX2NDC(0.9);
  jet_pt_resp_legend->Draw("same");
  jet_pt_resp_canvas->Update();
  jet_pt_resp_canvas->SaveAs("sim_plots_run21/jet_pt_respmatrix.png");
  jet_pt_resp_canvas->SaveAs("sim_plots_run21/jet_pt_respmatrix.pdf");

  TCanvas* calo_et_resp_canvas = new TCanvas("calo_et_resp_canvas","resp_canvas", 700, 650);
  calo_et_resp_canvas->SetTopMargin(0.18);
  TLegend* calo_et_resp_legend = new TLegend(.05,.82,.55,.99);
  calo_et_resp_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Preliminary","");
  calo_et_resp_legend->AddEntry("","Pythia8 200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
  calo_et_resp_legend->AddEntry("","p_{T}^{sublead jet} > 0.3 p_{T}^{lead jet}, |#Delta#phi| > 3#pi/4","");
  calo_et_resp_legend->SetFillStyle(0);
  calo_et_resp_legend->SetTextSize(0.035);
  calo_et_resp_canvas->SetLogz(1);
  calo_et_resp_canvas->SetRightMargin(0.2);
  calo_et_respmatrix->GetZaxis()->SetRangeUser(0.00001,500);
  calo_et_respmatrix->Draw("colz");
  calo_et_resp_canvas->Update();
  TPaletteAxis *calo_et_resp_palette = (TPaletteAxis*)calo_et_respmatrix->GetListOfFunctions()->FindObject("palette");
  calo_et_resp_palette->SetX1NDC(0.82);
  calo_et_resp_palette->SetX2NDC(0.9);
  calo_et_resp_legend->Draw("same");
  calo_et_resp_canvas->Update();
  calo_et_resp_canvas->SaveAs("sim_plots_run21/calo_et_respmatrix.png");
  calo_et_resp_canvas->SaveAs("sim_plots_run21/calo_et_respmatrix.pdf");
  
  TCanvas* counts_jet_pt_resp_canvas = new TCanvas("counts_jet_pt_resp_canvas","counts_resp_canvas", 700, 650);
  counts_jet_pt_resp_canvas->SetTopMargin(0.18);
  TLegend* counts_jet_pt_resp_legend = new TLegend(.05,.82,.55,.99);
  counts_jet_pt_resp_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Internal","");
  counts_jet_pt_resp_legend->AddEntry("","Pythia8 200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
  counts_jet_pt_resp_legend->AddEntry("","p_{T}^{sublead jet} > 0.3 p_{T}^{lead jet}, |#Delta#phi| > 3#pi/4","");
  counts_jet_pt_resp_legend->SetFillStyle(0);
  counts_jet_pt_resp_legend->SetTextSize(0.035);
  counts_jet_pt_resp_canvas->SetLogz(1);
  counts_jet_pt_resp_canvas->SetRightMargin(0.2);
  counts_jet_pt_respmatrix->GetZaxis()->SetRangeUser(0,100000);
  counts_jet_pt_respmatrix->Draw("colz");
  counts_jet_pt_resp_canvas->Update();
  TPaletteAxis *counts_jet_pt_resp_palette = (TPaletteAxis*)counts_jet_pt_respmatrix->GetListOfFunctions()->FindObject("palette");
  counts_jet_pt_resp_palette->SetX1NDC(0.82);
  counts_jet_pt_resp_palette->SetX2NDC(0.9);
  counts_jet_pt_resp_legend->Draw("same");
  counts_jet_pt_resp_canvas->Update();
  counts_jet_pt_resp_canvas->SaveAs("sim_plots_run21/counts_jet_pt_respmatrix.png");
  counts_jet_pt_resp_canvas->SaveAs("sim_plots_run21/counts_jet_pt_respmatrix.pdf");

  TCanvas* counts_calo_et_resp_canvas = new TCanvas("counts_calo_et_resp_canvas","counts_resp_canvas", 700, 650);
  counts_calo_et_resp_canvas->SetTopMargin(0.18);
  TLegend* counts_calo_et_resp_legend = new TLegend(.05,.82,.55,.99);
  counts_calo_et_resp_legend->AddEntry("","#bf{#it{sPHENIX}} Simulation Internal","");
  counts_calo_et_resp_legend->AddEntry("","Pythia8 200 GeV p+p anti-k_{t}#it{R}=0.4 |#eta_{jet}| < 0.7","");
  counts_calo_et_resp_legend->AddEntry("","p_{T}^{sublead jet} > 0.3 p_{T}^{lead jet}, |#Delta#phi| > 3#pi/4","");
  counts_calo_et_resp_legend->SetFillStyle(0);
  counts_calo_et_resp_legend->SetTextSize(0.035);
  counts_calo_et_resp_canvas->SetLogz(1);
  counts_calo_et_resp_canvas->SetRightMargin(0.2);
  counts_calo_et_respmatrix->GetZaxis()->SetRangeUser(0,60000);
  counts_calo_et_respmatrix->Draw("colz");
  counts_calo_et_resp_canvas->Update();
  TPaletteAxis *counts_calo_et_resp_palette = (TPaletteAxis*)counts_calo_et_respmatrix->GetListOfFunctions()->FindObject("palette");
  counts_calo_et_resp_palette->SetX1NDC(0.82);
  counts_calo_et_resp_palette->SetX2NDC(0.9);
  counts_calo_et_resp_legend->Draw("same");
  counts_calo_et_resp_canvas->Update();
  counts_calo_et_resp_canvas->SaveAs("sim_plots_run21/counts_calo_et_respmatrix.png");
  counts_calo_et_resp_canvas->SaveAs("sim_plots_run21/counts_calo_et_respmatrix.pdf");
  
  /*
  TCanvas* jet_spectra_canvas = new TCanvas("jet_spectra_canvas", "jet_spectra_canvas", 800, 600);
  jet_spectra_canvas->cd();
  jet_spectra_canvas->SetLogy(1);
  h_truth_lead_spectra_record->Rebin(20);
  h_lead_spectra_record->Rebin(20);
  h_truth_sub_spectra_record->Rebin(20);
  h_sub_spectra_record->Rebin(20);
  h_truth_lead_spectra_record->SetYTitle("N_{jet}");
  h_truth_lead_spectra_record->SetXTitle("p_{T} [GeV]");
  h_truth_lead_spectra_record->GetXaxis()->SetRangeUser(0,90);
  h_truth_lead_spectra_record->GetYaxis()->SetRangeUser(0.000001,80000);
  h_sub_spectra_record->SetMarkerStyle(24);
  h_truth_sub_spectra_record->SetMarkerStyle(24);
  h_truth_lead_spectra_record->SetMarkerColor(2);
  h_truth_lead_spectra_record->SetLineColor(2);
  h_truth_sub_spectra_record->SetMarkerColor(2);
  h_truth_sub_spectra_record->SetLineColor(2);
  h_truth_lead_spectra_record->Draw("");
  h_lead_spectra_record->Draw("same");
  h_truth_sub_spectra_record->Draw("same");
  h_sub_spectra_record->Draw("same");
  TLegend* jet_spectra_leg = new TLegend(.5,.72,.9,.92);
  jet_spectra_leg->AddEntry(h_lead_spectra_record, "Reco Leading Jet", "pel");
  jet_spectra_leg->AddEntry(h_truth_lead_spectra_record, "Truth Leading Jet", "pel");
  jet_spectra_leg->AddEntry(h_sub_spectra_record, "Reco Subleading Jet", "pel");
  jet_spectra_leg->AddEntry(h_truth_sub_spectra_record, "Truth Subleading Jet", "pel");
  jet_spectra_leg->Draw("same");
  jet_spectra_canvas->SaveAs("sim_plots_run21/jet_spectra_QA.png");
  
  TCanvas* nw_jet_pt_resp_canvas = new TCanvas("nw_jet_pt_resp_canvas", "nw_jet_pt_resp_canvas", 700, 600);
  nw_jet_pt_resp_canvas->SetLogz(1);
  nw_jet_pt_resp_canvas->SetRightMargin(0.2);
  h_nw_jet_pt_respmatrix_calib_dijet->GetZaxis()->SetRangeUser(1,80000);
  h_nw_jet_pt_respmatrix_calib_dijet->Draw("colz");
  nw_jet_pt_resp_canvas->Update();
  TPaletteAxis *nw_jet_pt_resp_palette = (TPaletteAxis*)h_nw_jet_pt_respmatrix_calib_dijet->GetListOfFunctions()->FindObject("palette");
  nw_jet_pt_resp_palette->SetX1NDC(0.82);
  nw_jet_pt_resp_palette->SetX2NDC(0.9);
  nw_jet_pt_resp_canvas->Update();
  nw_jet_pt_resp_canvas->SaveAs("sim_plots_run21/nw_jet_pt_respmatrix.png");

  TCanvas* nw_calo_et_resp_canvas = new TCanvas("nw_calo_et_resp_canvas", "nw_calo_et_resp_canvas", 700, 600);
  nw_calo_et_resp_canvas->SetLogz(1);
  nw_calo_et_resp_canvas->SetRightMargin(0.2);
  h_nw_calo_et_respmatrix_calib_dijet->GetZaxis()->SetRangeUser(1,80000);
  h_nw_calo_et_respmatrix_calib_dijet->Draw("colz");
  nw_calo_et_resp_canvas->Update();
  TPaletteAxis *nw_calo_et_resp_palette = (TPaletteAxis*)h_nw_calo_et_respmatrix_calib_dijet->GetListOfFunctions()->FindObject("palette");
  nw_calo_et_resp_palette->SetX1NDC(0.82);
  nw_calo_et_resp_palette->SetX2NDC(0.9);
  nw_calo_et_resp_canvas->Update();
  nw_calo_et_resp_canvas->SaveAs("sim_plots_run21/nw_calo_et_respmatrix.png");
  */
/*
  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  TFile *f_out = new TFile("sim_plots_run21/all_sim_plots.root", "RECREATE");
  f_out->cd();
  et_dist_canvas->Write();
  //et_dist_bigbin_canvas->Write();
  truth_ue_vs_pt_canvas->Write();
  resp_canvas->Write();
  counts_resp_canvas->Write();
  jet_pt_resp_canvas->Write();
  calo_et_resp_canvas->Write();
  counts_jet_pt_resp_canvas->Write();
  counts_calo_et_resp_canvas->Write();
  //jet_spectra_canvas->Write();
  //nw_jet_pt_resp_canvas->Write();
  //nw_calo_et_resp_canvas->Write();
  f_out->Close();
  std::cout << "All done!" << std::endl;
  */
  
}