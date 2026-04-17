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
#include "/sphenix/user/hanpuj/CaloDataAna24_skimmed/src/draw_template.C" 

void respmatrix_QA(bool dijet = 1) {
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
  TFile *f_sim;
  if (dijet) {
    f_sim = new TFile("sphenix_primary_analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_3.root", "READ");
  } else {
    f_sim = new TFile("sphenix_primary_analysis_sim_run28_output/output_none_bkg_cut_sim_iter_3.root", "READ");
  }
  //TFile *f_purityefficiency = new TFile("output_purityefficiency.root", "READ");

  TH1D* h_zvertex = (TH1D*)f_sim->Get("h_zvertex"); 
  
  //TH1D* h_deltaphi_record = (TH1D*)f_sim->Get("h_deltaphi_record"); 
  //TH1D* h_truth_deltaphi_record = (TH1D*)f_sim->Get("h_truth_deltaphi_record"); 
  
  //TH1D* h_xj_record = (TH1D*)f_sim->Get("h_xj_record"); 
  //TH1D* h_truth_xj_record = (TH1D*)f_sim->Get("h_truth_xj_record"); 
  
  TH1D* h_lead_spectra_record = (TH1D*)f_sim->Get("h_lead_spectra_record"); 
  TH1D* h_sub_spectra_record = (TH1D*)f_sim->Get("h_sub_spectra_record"); 
  TH1D* h_truth_lead_spectra_record = (TH1D*)f_sim->Get("h_lead_truth_spectra_record"); 
  TH1D* h_truth_sub_spectra_record = (TH1D*)f_sim->Get("h_sub_truth_spectra_record"); 
  
  TH1D* h_et_transverse_record = (TH1D*)f_sim->Get("h_et_transverse_record"); 
  TH1D* h_nw_et_transverse_record = (TH1D*)f_sim->Get("h_nw_et_transverse_record"); 
  TH1D* h_truth_et_transverse_record = (TH1D*)f_sim->Get("h_et_truth_transverse_record"); 
  TH1D* h_truth_nw_et_transverse_record = (TH1D*)f_sim->Get("h_nw_et_truth_transverse_record"); 
  
  TH2D* h_ue_pt_transverse_record = (TH2D*)f_sim->Get("h_ue_pt_transverse_record"); 
  TH2D* h_truth_ue_pt_transverse_record = (TH2D*)f_sim->Get("h_ue_pt_truth_transverse_record"); 

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
  RooUnfoldResponse* h_respmatrix_calib_dijet = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet");  
  RooUnfoldResponse* h_respmatrix_calib_dijet_trim = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_trim_10");  
  RooUnfoldResponse* h_respmatrix_calib_dijet_rw_trim = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_reweight_trim_10");  
  TH2D* h_fake_calib_dijet = (TH2D*)f_sim->Get("h_fake_calib_dijet"); 
  TH2D* h_miss_calib_dijet = (TH2D*)f_sim->Get("h_miss_calib_dijet"); 
  RooUnfoldResponse* h_counts_calib_dijet = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_counts");  

  /*
  for (int x = 1; x <= h_nw_jet_pt_respmatrix_calib_dijet->GetNbinsX(); x++) {
    for (int y = 1; y <= h_nw_jet_pt_respmatrix_calib_dijet->GetNbinsY(); y++) {
      if (h_nw_jet_pt_respmatrix_calib_dijet->GetBinContent(x,y) > 0 && h_nw_jet_pt_respmatrix_calib_dijet->GetBinContent(x,y) <= 5) {
        std::cout << "Bin with " << h_nw_jet_pt_respmatrix_calib_dijet->GetBinContent(x,y) << " entries: x " << x << " y " << y << std::endl; 
      }
    }
  }
  std::cout << std::endl;
  for (int x = 1; x <= h_nw_calo_et_respmatrix_calib_dijet->GetNbinsX(); x++) {
    for (int y = 1; y <= h_nw_calo_et_respmatrix_calib_dijet->GetNbinsY(); y++) {
      if (h_nw_calo_et_respmatrix_calib_dijet->GetBinContent(x,y) > 0 && h_nw_calo_et_respmatrix_calib_dijet->GetBinContent(x,y) <= 5) {
        std::cout << "Bin with " << h_nw_calo_et_respmatrix_calib_dijet->GetBinContent(x,y) << " entries: x " << x << " y " << y << std::endl; 
      }
    }
  }

  const TMatrixD& counts_matrix = h_counts_calib_dijet->Mresponse(false);  // use normalized=true if needed
  const TMatrixD& resp_matrix = h_respmatrix_calib_dijet->Mresponse(false);  // use normalized=true if needed
  
  int nxr = h_measure_calib_dijet->GetNbinsX(); // reco x bins
  int nyr = h_measure_calib_dijet->GetNbinsY(); // reco y bins
  int nxt = h_truth_calib_dijet->GetNbinsX(); // truth x bins
  int nyt = h_truth_calib_dijet->GetNbinsY(); // truth y bins

  TH2D* hist2d = new TH2D("response2D", "Response Matrix;Reco Bin;Truth Bin", nxr*nyr, 0, nxr*nyr, nxt*nyt, 0, nxt*nyt);

  for (int ixr = 0; ixr < nxr; ++ixr) {
    for (int iyr = 0; iyr < nyr; ++iyr) {
      int recoBin = ixr * nyr + iyr;
      for (int ixt = 0; ixt < nxt; ++ixt) {
        for (int iyt = 0; iyt < nyt; ++iyt) {
          int truthBin = ixt * nyt + iyt;
          if (counts_matrix(recoBin, truthBin) < 20) { continue; }
          double value = resp_matrix(recoBin, truthBin);
          hist2d->SetBinContent(recoBin + 1, truthBin + 1, value);
        }
      }
    }
  }

  */

  const TMatrixD& count_matrix = h_counts_calib_dijet->Mresponse(false);
  const TMatrixD& trim_matrix = h_respmatrix_calib_dijet_trim->Mresponse(false);  // use normalized=true if needed
  const TMatrixD& rw_trim_matrix = h_respmatrix_calib_dijet_rw_trim->Mresponse(false);
  const TMatrixD& resp_matrix = h_respmatrix_calib_dijet->Mresponse(false);  // use normalized=true if needed
  
  int nxr = h_measure_calib_dijet->GetNbinsX(); // reco x bins
  int nyr = h_measure_calib_dijet->GetNbinsY(); // reco y bins
  int nxt = h_truth_calib_dijet->GetNbinsX(); // truth x bins
  int nyt = h_truth_calib_dijet->GetNbinsY(); // truth y bins

  TH2D* hist2d = new TH2D("response2D", "Response Matrix;Reco Bin;Truth Bin", nxr*nyr, 0, nxr*nyr, nxt*nyt, 0, nxt*nyt);
  TH2D* trim_hist2d = new TH2D("trim_response2D", "Trim Response Matrix;Reco Bin;Truth Bin", nxr*nyr, 0, nxr*nyr, nxt*nyt, 0, nxt*nyt);
  TH2D* rw_trim_hist2d = new TH2D("rw_trim_response2D","Reweight Trim Response Matrix;Reco Bin;Truth Bin",nxr*nyr, 0, nxr*nyr, nxt*nyt, 0, nxt*nyt);
  TH2D* trim_jet_pt_respmatrix = new TH2D("trim_jet_pt_respmatrix",";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", calibnpt, calibptbins, truthnpt, truthptbins);
  TH2D* trim_calo_et_respmatrix = new TH2D("trim_calo_et_respmatrix",";#SigmaE_{T}^{Reco} [GeV];#SigmaE_{T}^{Truth} [GeV]", calibnet, calibetbins, truthnet, truthetbins);
  TH2D* rw_trim_jet_pt_respmatrix = new TH2D("rw_trim_jet_pt_respmatrix",";p_{T}^{Calib jet} [GeV];p_{T}^{Truth jet} [GeV]", calibnpt, calibptbins, truthnpt, truthptbins);
  TH2D* rw_trim_calo_et_respmatrix = new TH2D("rw_trim_calo_et_respmatrix",";#SigmaE_{T}^{Reco} [GeV];#SigmaE_{T}^{Truth} [GeV]", calibnet, calibetbins, truthnet, truthetbins);


  for (int ixr = 0; ixr < nxr; ++ixr) {
    for (int iyr = 0; iyr < nyr; ++iyr) {
      int recoBin = iyr * nxr + ixr;
      for (int ixt = 0; ixt < nxt; ++ixt) {
        for (int iyt = 0; iyt < nyt; ++iyt) {
          int truthBin = iyt * nxt + ixt;
          hist2d->SetBinContent(recoBin + 1, truthBin + 1, resp_matrix(recoBin, truthBin));
          trim_hist2d->SetBinContent(recoBin + 1, truthBin + 1, trim_matrix(recoBin, truthBin));
          rw_trim_hist2d->SetBinContent(recoBin + 1, truthBin + 1, rw_trim_matrix(recoBin, truthBin));
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
          sum += trim_matrix(row, col);
          rw_sum += rw_trim_matrix(row, col);
        }
      }
      trim_jet_pt_respmatrix->SetBinContent(ixr + 1, ixt + 1, sum);
      rw_trim_jet_pt_respmatrix->SetBinContent(ixr + 1, ixt + 1, rw_sum);
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
          sum += trim_matrix(row, col);
          rw_sum += rw_trim_matrix(row, col);
        }
      }
      trim_calo_et_respmatrix->SetBinContent(iyr + 1, iyt + 1, sum);
      rw_trim_calo_et_respmatrix->SetBinContent(iyr + 1, iyt + 1, rw_sum);
    }
  }
  /*
  TCanvas* resp_canvas = new TCanvas("resp_canvas","resp_canvas", 700, 650);
  resp_canvas->SetTopMargin(0.18);
  TLegend* resp_legend = new TLegend(.05,.82,.55,.99);
  resp_legend->AddEntry("","#bf{#it{sPHENIX}} Internal Pythia8 200 GeV p+p","");
  resp_legend->AddEntry("","anti-k_{t}#it{R}=0.4, |#eta_{jet}| < 0.7, p_{T}^{truth lead jet} > 17 GeV","");
  if (dijet) resp_legend->AddEntry("","Exclusive Dijet","");
  else { resp_legend->AddEntry("","Inclusive Jet",""); }
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
  if (dijet) resp_canvas->SaveAs("sphenix_primary_run28_output_files/raw_respmatrix_dijet.png");
  else resp_canvas->SaveAs("sphenix_primary_run28_output_files/raw_respmatrix_none.png");

  TCanvas* trim_resp_canvas = new TCanvas("trim_resp_canvas","trim_resp_canvas", 700, 600);
  trim_resp_canvas->SetLogz(1);
  trim_resp_canvas->SetRightMargin(0.2);
  trim_hist2d->GetZaxis()->SetRangeUser(0.000001,500);
  trim_hist2d->Draw("colz");
  trim_resp_canvas->Update();
  TPaletteAxis *trim_resp_palette = (TPaletteAxis*)trim_hist2d->GetListOfFunctions()->FindObject("palette");
  trim_resp_palette->SetX1NDC(0.82);
  trim_resp_palette->SetX2NDC(0.9);
  trim_resp_canvas->Update();
  trim_resp_canvas->SaveAs("sphenix_primary_run28_output_files/respmatrix_QA_herwig_trim_respmatrix.png");

  TCanvas* trim_jet_pt_resp_canvas = new TCanvas("trim_jet_pt_resp_canvas","trim_resp_canvas", 700, 600);
  trim_jet_pt_resp_canvas->SetLogz(1);
  trim_jet_pt_resp_canvas->SetRightMargin(0.2);
  trim_jet_pt_respmatrix->GetZaxis()->SetRangeUser(0.0000001,2000);
  trim_jet_pt_respmatrix->Draw("colz");
  trim_jet_pt_resp_canvas->Update();
  TPaletteAxis *trim_jet_pt_resp_palette = (TPaletteAxis*)trim_jet_pt_respmatrix->GetListOfFunctions()->FindObject("palette");
  trim_jet_pt_resp_palette->SetX1NDC(0.82);
  trim_jet_pt_resp_palette->SetX2NDC(0.9);
  trim_jet_pt_resp_canvas->Update();
  trim_jet_pt_resp_canvas->SaveAs("sphenix_primary_run28_output_files/respmatrix_QA_herwig_trim_jet_pt_respmatrix.png");

  TCanvas* trim_calo_et_resp_canvas = new TCanvas("trim_calo_et_resp_canvas","trim_resp_canvas", 700, 600);
  trim_calo_et_resp_canvas->SetLogz(1);
  trim_calo_et_resp_canvas->SetRightMargin(0.2);
  trim_calo_et_respmatrix->GetZaxis()->SetRangeUser(0.0000001,500);
  trim_calo_et_respmatrix->Draw("colz");
  trim_calo_et_resp_canvas->Update();
  TPaletteAxis *trim_calo_et_resp_palette = (TPaletteAxis*)trim_calo_et_respmatrix->GetListOfFunctions()->FindObject("palette");
  trim_calo_et_resp_palette->SetX1NDC(0.82);
  trim_calo_et_resp_palette->SetX2NDC(0.9);
  trim_calo_et_resp_canvas->Update();
  trim_calo_et_resp_canvas->SaveAs("sphenix_primary_run28_output_files/respmatrix_QA_herwig_trim_calo_et_respmatrix.png");
  */
  TCanvas* rw_trim_resp_canvas = new TCanvas("rw_trim_resp_canvas","rw_trim_resp_canvas", 700, 650);
  rw_trim_resp_canvas->SetTopMargin(0.18);
  TLegend* rw_trim_resp_legend = new TLegend(.05,.82,.55,.99);
  rw_trim_resp_legend->AddEntry("","#bf{#it{sPHENIX}} Internal Pythia8 200 GeV p+p","");
  rw_trim_resp_legend->AddEntry("","anti-k_{t}#it{R}=0.4, |#eta_{jet}| < 0.7, p_{T}^{truth lead jet} > 17 GeV","");
  if (dijet) rw_trim_resp_legend->AddEntry("","Exclusive Dijet","");
  else { rw_trim_resp_legend->AddEntry("","Inclusive Jet",""); }
  rw_trim_resp_legend->SetFillStyle(0);
  rw_trim_resp_legend->SetTextSize(0.035);
  rw_trim_resp_canvas->SetLogz(1);
  rw_trim_resp_canvas->SetRightMargin(0.2);
  rw_trim_hist2d->GetZaxis()->SetRangeUser(0.0000001,800);
  rw_trim_hist2d->Draw("colz");
  rw_trim_resp_canvas->Update();
  TPaletteAxis *rw_trim_resp_palette = (TPaletteAxis*)rw_trim_hist2d->GetListOfFunctions()->FindObject("palette");
  rw_trim_resp_palette->SetX1NDC(0.82);
  rw_trim_resp_palette->SetX2NDC(0.9);
  rw_trim_resp_legend->Draw("same");
  rw_trim_resp_canvas->Update();
  if (dijet) rw_trim_resp_canvas->SaveAs("sphenix_primary_run28_output_files/respmatrix_dijet.png");
  else rw_trim_resp_canvas->SaveAs("sphenix_primary_run28_output_files/respmatrix_none.png");

  TCanvas* rw_trim_jet_pt_resp_canvas = new TCanvas("rw_trim_jet_pt_resp_canvas","rw_trim_jet_pt_resp_canvas", 700, 650);
  rw_trim_jet_pt_resp_canvas->SetTopMargin(0.18);
  TLegend* rw_trim_jet_pt_resp_legend = new TLegend(.05,.82,.55,.99);
  rw_trim_jet_pt_resp_legend->AddEntry("","#bf{#it{sPHENIX}} Internal Pythia8 200 GeV p+p","");
  rw_trim_jet_pt_resp_legend->AddEntry("","anti-k_{t}#it{R}=0.4, |#eta_{jet}| < 0.7, p_{T}^{truth lead jet} > 17 GeV","");
  if (dijet) rw_trim_jet_pt_resp_legend->AddEntry("","Exclusive Dijet","");
  else { rw_trim_jet_pt_resp_legend->AddEntry("","Inclusive Jet",""); }
  rw_trim_jet_pt_resp_legend->SetFillStyle(0);
  rw_trim_jet_pt_resp_legend->SetTextSize(0.035);
  rw_trim_jet_pt_resp_canvas->SetLogz(1);
  rw_trim_jet_pt_resp_canvas->SetRightMargin(0.2);
  rw_trim_jet_pt_respmatrix->GetZaxis()->SetRangeUser(0.0001,100000);
  rw_trim_jet_pt_respmatrix->Draw("colz");
  rw_trim_jet_pt_resp_canvas->Update();
  TPaletteAxis *rw_trim_jet_pt_resp_palette = (TPaletteAxis*)rw_trim_jet_pt_respmatrix->GetListOfFunctions()->FindObject("palette");
  rw_trim_jet_pt_resp_palette->SetX1NDC(0.82);
  rw_trim_jet_pt_resp_palette->SetX2NDC(0.9);
  rw_trim_jet_pt_resp_legend->Draw("same");
  rw_trim_jet_pt_resp_canvas->Update();
  if (dijet) rw_trim_jet_pt_resp_canvas->SaveAs("sphenix_primary_run28_output_files/jet_pt_respmatrix_dijet.png");
  else rw_trim_jet_pt_resp_canvas->SaveAs("sphenix_primary_run28_output_files/jet_pt_respmatrix_none.png");

  TCanvas* rw_trim_calo_et_resp_canvas = new TCanvas("rw_trim_calo_et_resp_canvas","rw_trim_calo_et_resp_canvas", 700, 650);
  rw_trim_calo_et_resp_canvas->SetTopMargin(0.18);
  TLegend* rw_trim_calo_et_resp_legend = new TLegend(.05,.82,.55,.99);
  rw_trim_calo_et_resp_legend->AddEntry("","#bf{#it{sPHENIX}} Internal Pythia8 200 GeV p+p","");
  rw_trim_calo_et_resp_legend->AddEntry("","anti-k_{t}#it{R}=0.4, |#eta_{jet}| < 0.7, p_{T}^{truth lead jet} > 17 GeV","");
  if (dijet) rw_trim_calo_et_resp_legend->AddEntry("","Exclusive Dijet","");
  else { rw_trim_calo_et_resp_legend->AddEntry("","Inclusive Jet",""); }
  rw_trim_calo_et_resp_legend->SetFillStyle(0);
  rw_trim_calo_et_resp_legend->SetTextSize(0.035);
  rw_trim_calo_et_resp_canvas->SetLogz(1);
  rw_trim_calo_et_resp_canvas->SetRightMargin(0.2);
  rw_trim_calo_et_respmatrix->GetZaxis()->SetRangeUser(0.0001,20000);
  rw_trim_calo_et_respmatrix->Draw("colz");
  rw_trim_calo_et_resp_canvas->Update();
  TPaletteAxis *rw_trim_calo_et_resp_palette = (TPaletteAxis*)rw_trim_calo_et_respmatrix->GetListOfFunctions()->FindObject("palette");
  rw_trim_calo_et_resp_palette->SetX1NDC(0.82);
  rw_trim_calo_et_resp_palette->SetX2NDC(0.9);
  rw_trim_calo_et_resp_legend->Draw("same");
  rw_trim_calo_et_resp_canvas->Update();
  if (dijet) rw_trim_calo_et_resp_canvas->SaveAs("sphenix_primary_run28_output_files/calo_et_respmatrix_dijet.png");
  else rw_trim_calo_et_resp_canvas->SaveAs("sphenix_primary_run28_output_files/calo_et_respmatrix_none.png");

  /*
  TCanvas* zvtx_canvas = new TCanvas("zvtx_canvas", "zvtx_canvas", 800, 600);
  zvtx_canvas->cd();
  h_zvertex->Scale(1.0/h_zvertex->Integral());
  h_zvertex->SetYTitle("1/N dN/dv_{z} [cm^{-1}]");
  h_zvertex->SetXTitle("v_{z} [cm]");
  h_zvertex->GetXaxis()->SetRangeUser(-40,40);
  h_zvertex->Draw();
  zvtx_canvas->SaveAs("sphenix_primary_run28_output_files/respmatrix_QA_herwig_zvertex.png");
  
  TCanvas* deltaphi_canvas = new TCanvas("deltaphi_canvas", "deltaphi_canvas", 800, 600);
  deltaphi_canvas->cd();
  h_deltaphi_record->Scale(1.0/h_deltaphi_record->Integral());
  h_truth_deltaphi_record->Scale(1.0/h_truth_deltaphi_record->Integral());
  h_truth_deltaphi_record->SetYTitle("1/N dN/d(#Delta#phi)");
  h_truth_deltaphi_record->SetXTitle("#Delta#phi");
  h_truth_deltaphi_record->GetXaxis()->SetRangeUser(0,4);
  h_truth_deltaphi_record->SetMarkerColor(2);
  h_truth_deltaphi_record->SetLineColor(2);
  h_truth_deltaphi_record->Draw("hist");
  h_deltaphi_record->Draw("hist,same");
  TLegend* deltaphi_leg = new TLegend(.2,.65,.5,.85);
  deltaphi_leg->AddEntry(h_deltaphi_record, "Reco", "l");
  deltaphi_leg->AddEntry(h_truth_deltaphi_record, "Truth", "l");
  deltaphi_leg->Draw("same");
  deltaphi_canvas->SaveAs("uniform_bin_run21_figure/respmatrix_QA_herwig_deltaphi.png");

  TCanvas* xj_canvas = new TCanvas("xj_canvas", "xj_canvas", 800, 600);
  xj_canvas->cd();
  h_xj_record->Scale(1.0/h_xj_record->Integral());
  h_truth_xj_record->Scale(1.0/h_truth_xj_record->Integral());
  h_xj_record->SetYTitle("1/N dN/d(x_{j})");
  h_xj_record->SetXTitle("x_{j}");
  h_xj_record->GetXaxis()->SetRangeUser(0,1);
  h_truth_xj_record->SetMarkerColor(2);
  h_truth_xj_record->SetLineColor(2);
  h_xj_record->Draw("hist");
  h_truth_xj_record->Draw("hist,same");
  TLegend* xj_leg = new TLegend(.2,.65,.5,.85);
  xj_leg->AddEntry(h_xj_record, "Reco", "l");
  xj_leg->AddEntry(h_truth_xj_record, "Truth", "l");
  xj_leg->Draw("same");
  xj_canvas->SaveAs("uniform_bin_run21_figure/respmatrix_QA_herwig_xj.png");
  */
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
  h_truth_lead_spectra_record->GetYaxis()->SetRangeUser(0.000001,50000);
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
  jet_spectra_canvas->SaveAs("uniform_bin_run21_figure/respmatrix_QA_herwig_jet_spectra.png");

  TCanvas* et_spectra_canvas = new TCanvas("et_spectra_canvas", "et_spectra_canvas", 800, 600);
  et_spectra_canvas->SetLogy(1);
  h_truth_et_transverse_record->Rebin(20);
  h_et_transverse_record->Rebin(20);
  h_truth_nw_et_transverse_record->Rebin(20);
  h_nw_et_transverse_record->Rebin(20);
  h_truth_et_transverse_record->Scale(1.0/h_truth_et_transverse_record->Integral());
  h_et_transverse_record->Scale(1.0/h_et_transverse_record->Integral());
  h_truth_nw_et_transverse_record->Scale(1.0/h_truth_nw_et_transverse_record->Integral());
  h_nw_et_transverse_record->Scale(1.0/h_nw_et_transverse_record->Integral());
  h_et_transverse_record->SetYTitle("1/N dN/d#SigmaE_{T} [GeV^{-1}]");
  h_et_transverse_record->SetXTitle("#SigmaE_{T} [GeV]");
  h_et_transverse_record->GetXaxis()->SetRangeUser(0,20);
  h_nw_et_transverse_record->SetMarkerStyle(24);
  h_truth_nw_et_transverse_record->SetMarkerStyle(24);
  h_truth_et_transverse_record->SetMarkerColor(2);
  h_truth_et_transverse_record->SetLineColor(2);
  h_truth_nw_et_transverse_record->SetMarkerColor(2);
  h_truth_nw_et_transverse_record->SetLineColor(2);
  h_et_transverse_record->Draw("");
  h_truth_et_transverse_record->Draw("same");
  h_truth_nw_et_transverse_record->Draw("same");
  h_nw_et_transverse_record->Draw("same");
  TLegend* et_spectra_leg = new TLegend(.5,.72,.9,.92);
  et_spectra_leg->AddEntry(h_et_transverse_record, "Reco", "pel");
  et_spectra_leg->AddEntry(h_truth_et_transverse_record, "Truth", "pel");
  et_spectra_leg->AddEntry(h_nw_et_transverse_record, "Reco (No jet weighting)", "pel");
  et_spectra_leg->AddEntry(h_truth_nw_et_transverse_record, "Truth (No jet weighting)", "pel");
  et_spectra_leg->Draw("same");
  et_spectra_canvas->SaveAs("uniform_bin_run21_figure/respmatrix_QA_herwig_et_spectra.png");

  TCanvas* jet_pt_resp_canvas = new TCanvas("jet_pt_resp_canvas", "jet_pt_resp_canvas", 700, 600);
  jet_pt_resp_canvas->SetLogz(1);
  jet_pt_resp_canvas->SetRightMargin(0.2);
  h_jet_pt_respmatrix_calib_dijet->GetZaxis()->SetRangeUser(0.0000001,2000);
  h_jet_pt_respmatrix_calib_dijet->Draw("colz");
  jet_pt_resp_canvas->Update();
  TPaletteAxis *jet_pt_resp_palette = (TPaletteAxis*)h_jet_pt_respmatrix_calib_dijet->GetListOfFunctions()->FindObject("palette");
  jet_pt_resp_palette->SetX1NDC(0.82);
  jet_pt_resp_palette->SetX2NDC(0.9);
  jet_pt_resp_canvas->Update();
  jet_pt_resp_canvas->SaveAs("uniform_bin_run21_figure/respmatrix_QA_herwig_jet_pt_respmatrix.png");

  TCanvas* jet_pt_resp_canvas_zoom = new TCanvas("jet_pt_resp_canvas_zoom", "jet_pt_resp_canvas_zoom", 700, 600);
  jet_pt_resp_canvas_zoom->SetLogz(1);
  jet_pt_resp_canvas_zoom->SetRightMargin(0.2);
  h_jet_pt_respmatrix_calib_dijet->GetZaxis()->SetRangeUser(0.0000001,2000);
  h_jet_pt_respmatrix_calib_dijet->GetYaxis()->SetRangeUser(18,72);
  h_jet_pt_respmatrix_calib_dijet->Draw("colz");
  jet_pt_resp_canvas_zoom->Update();
  TPaletteAxis *jet_pt_resp_palette_zoom = (TPaletteAxis*)h_jet_pt_respmatrix_calib_dijet->GetListOfFunctions()->FindObject("palette");
  jet_pt_resp_palette_zoom->SetX1NDC(0.82);
  jet_pt_resp_palette_zoom->SetX2NDC(0.9);
  jet_pt_resp_canvas_zoom->Update();
  jet_pt_resp_canvas_zoom->SaveAs("uniform_bin_run21_figure/respmatrix_QA_herwig_jet_pt_respmatrix_zoom.png");

  TCanvas* calo_et_resp_canvas = new TCanvas("calo_et_resp_canvas", "calo_et_resp_canvas", 700, 600);
  calo_et_resp_canvas->SetLogz(1);
  calo_et_resp_canvas->SetRightMargin(0.2);
  h_calo_et_respmatrix_calib_dijet->GetZaxis()->SetRangeUser(0.0000001,500);
  h_calo_et_respmatrix_calib_dijet->Draw("colz");
  calo_et_resp_canvas->Update();
  TPaletteAxis *calo_et_resp_palette = (TPaletteAxis*)h_calo_et_respmatrix_calib_dijet->GetListOfFunctions()->FindObject("palette");
  calo_et_resp_palette->SetX1NDC(0.82);
  calo_et_resp_palette->SetX2NDC(0.9);
  calo_et_resp_canvas->Update();
  calo_et_resp_canvas->SaveAs("uniform_bin_run21_figure/respmatrix_QA_herwig_calo_et_respmatrix.png");
  
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
  nw_jet_pt_resp_canvas->SaveAs("uniform_bin_run21_figure/respmatrix_QA_herwig_nw_jet_pt_respmatrix.png");

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
  nw_calo_et_resp_canvas->SaveAs("uniform_bin_run21_figure/respmatrix_QA_herwig_nw_calo_et_respmatrix.png");

  
  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  TFile *f_out = new TFile("respmatrix_QA_run21_herwig_iter_3.root", "RECREATE");
  f_out->cd();
  resp_canvas->Write();
  trim_resp_canvas->Write();
  trim_jet_pt_resp_canvas->Write();
  trim_calo_et_resp_canvas->Write();
  rw_trim_resp_canvas->Write();
  rw_trim_jet_pt_resp_canvas->Write();
  rw_trim_calo_et_resp_canvas->Write();
  zvtx_canvas->Write();
  deltaphi_canvas->Write();
  xj_canvas->Write();
  jet_spectra_canvas->Write();
  et_spectra_canvas->Write();
  jet_pt_resp_canvas->Write();
  jet_pt_resp_canvas_zoom->Write();
  calo_et_resp_canvas->Write();
  nw_jet_pt_resp_canvas->Write();
  nw_calo_et_resp_canvas->Write();
  f_out->Close();
  std::cout << "All done!" << std::endl;

  */
  
}