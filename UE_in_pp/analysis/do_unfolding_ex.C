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

////////////////////////////////////////// Main Function //////////////////////////////////////////
void do_unfolding() {
  //********** General Set up **********//
  const float PI = TMath::Pi();
  const float jet_radius = 0.4;
 
  //********** Files **********//
  TFile *f_out = new TFile("output_unfolded.root", "RECREATE");

  TFile *f_data = new TFile("output_data.root", "READ");
  TFile *f_in_rm = new TFile("output_reweightmatrix.root", "READ");
  TFile *f_efficiency = new TFile("output_purityefficiency.root", "READ");

  TH1D* h_recojet_pt = (TH1D*)f_data->Get("h_recojet_pt_puritycorr");
  TH1D* h_calibjet_pt = (TH1D*)f_data->Get("h_calibjet_pt_puritycorr");

  TH2D* h_respmatrix_reco_record = (TH2D*)f_in_rm->Get("h_respmatrix_reco_all_record");
  TH2D* h_respmatrix_calib_record = (TH2D*)f_in_rm->Get("h_respmatrix_calib_all_record");
  TH2D* h_respmatrix_reco_reweighted = (TH2D*)f_in_rm->Get("h_respmatrix_reco_all_reweighted");
  TH2D* h_respmatrix_calib_reweighted = (TH2D*)f_in_rm->Get("h_respmatrix_calib_all_reweighted");

  TH1D* h_efficiency_reco = (TH1D*)f_efficiency->Get("h_efficiency_reco_all");
  TH1D* h_efficiency_calib = (TH1D*)f_efficiency->Get("h_efficiency_calib_all");

  TH1D* h_meas_reco = (TH1D*)h_respmatrix_reco_record->ProjectionX("h_meas_reco");
  TH1D* h_truth_reco = (TH1D*)h_respmatrix_reco_record->ProjectionY("h_truth_reco");
  RooUnfoldResponse* resp_reco = new RooUnfoldResponse(h_meas_reco, h_truth_reco, h_respmatrix_reco_record, "resp_reco", "Response matrix for reco jet");
  TH1D* h_meas_calib = (TH1D*)h_respmatrix_calib_record->ProjectionX("h_meas_calib");
  TH1D* h_truth_calib = (TH1D*)h_respmatrix_calib_record->ProjectionY("h_truth_calib");
  RooUnfoldResponse* resp_calib = new RooUnfoldResponse(h_meas_calib, h_truth_calib, h_respmatrix_calib_record, "resp_calib", "Response matrix for calib jet");
  TH1D* h_meas_reco_reweighted = (TH1D*)h_respmatrix_reco_reweighted->ProjectionX("h_meas_reco_reweighted");
  TH1D* h_truth_reco_reweighted = (TH1D*)h_respmatrix_reco_reweighted->ProjectionY("h_truth_reco_reweighted");
  RooUnfoldResponse* resp_reco_reweighted = new RooUnfoldResponse(h_meas_reco_reweighted, h_truth_reco_reweighted, h_respmatrix_reco_reweighted, "resp_reco_reweighted", "Response matrix for reco jet (reweighted)");
  TH1D* h_meas_calib_reweighted = (TH1D*)h_respmatrix_calib_reweighted->ProjectionX("h_meas_calib_reweighted");
  TH1D* h_truth_calib_reweighted = (TH1D*)h_respmatrix_calib_reweighted->ProjectionY("h_truth_calib_reweighted");
  RooUnfoldResponse* resp_calib_reweighted = new RooUnfoldResponse(h_meas_calib_reweighted, h_truth_calib_reweighted, h_respmatrix_calib_reweighted, "resp_calib_reweighted", "Response matrix for calib jet (reweighted)");

  //********** Unfold **********//
  RooUnfoldBayes unfold_reco_1(resp_reco, h_recojet_pt, 1); TH1D* h_unfold_reco_1 = (TH1D*)unfold_reco_1.Hunfold(RooUnfold::kErrors); h_unfold_reco_1->SetName("h_unfold_reco_1");
  RooUnfoldBayes unfold_reco_2(resp_reco, h_recojet_pt, 2); TH1D* h_unfold_reco_2 = (TH1D*)unfold_reco_2.Hunfold(RooUnfold::kErrors); h_unfold_reco_2->SetName("h_unfold_reco_2");
  RooUnfoldBayes unfold_reco_3(resp_reco, h_recojet_pt, 3); TH1D* h_unfold_reco_3 = (TH1D*)unfold_reco_3.Hunfold(RooUnfold::kErrors); h_unfold_reco_3->SetName("h_unfold_reco_3");
  RooUnfoldBayes unfold_reco_4(resp_reco, h_recojet_pt, 4); TH1D* h_unfold_reco_4 = (TH1D*)unfold_reco_4.Hunfold(RooUnfold::kErrors); h_unfold_reco_4->SetName("h_unfold_reco_4");
  RooUnfoldBayes unfold_reco_5(resp_reco, h_recojet_pt, 5); TH1D* h_unfold_reco_5 = (TH1D*)unfold_reco_5.Hunfold(RooUnfold::kErrors); h_unfold_reco_5->SetName("h_unfold_reco_5");
  RooUnfoldBayes unfold_reco_6(resp_reco, h_recojet_pt, 6); TH1D* h_unfold_reco_6 = (TH1D*)unfold_reco_6.Hunfold(RooUnfold::kErrors); h_unfold_reco_6->SetName("h_unfold_reco_6");

  RooUnfoldBayes unfold_calib_1(resp_calib, h_calibjet_pt, 1); TH1D* h_unfold_calib_1 = (TH1D*)unfold_calib_1.Hunfold(RooUnfold::kErrors); h_unfold_calib_1->SetName("h_unfold_calib_1");
  RooUnfoldBayes unfold_calib_2(resp_calib, h_calibjet_pt, 2); TH1D* h_unfold_calib_2 = (TH1D*)unfold_calib_2.Hunfold(RooUnfold::kErrors); h_unfold_calib_2->SetName("h_unfold_calib_2");
  RooUnfoldBayes unfold_calib_3(resp_calib, h_calibjet_pt, 3); TH1D* h_unfold_calib_3 = (TH1D*)unfold_calib_3.Hunfold(RooUnfold::kErrors); h_unfold_calib_3->SetName("h_unfold_calib_3");
  RooUnfoldBayes unfold_calib_4(resp_calib, h_calibjet_pt, 4); TH1D* h_unfold_calib_4 = (TH1D*)unfold_calib_4.Hunfold(RooUnfold::kErrors); h_unfold_calib_4->SetName("h_unfold_calib_4");
  RooUnfoldBayes unfold_calib_5(resp_calib, h_calibjet_pt, 5); TH1D* h_unfold_calib_5 = (TH1D*)unfold_calib_5.Hunfold(RooUnfold::kErrors); h_unfold_calib_5->SetName("h_unfold_calib_5");
  RooUnfoldBayes unfold_calib_6(resp_calib, h_calibjet_pt, 6); TH1D* h_unfold_calib_6 = (TH1D*)unfold_calib_6.Hunfold(RooUnfold::kErrors); h_unfold_calib_6->SetName("h_unfold_calib_6");

  RooUnfoldBayes unfold_reco_reweighted_1(resp_reco_reweighted, h_recojet_pt, 1); TH1D* h_unfold_reco_reweighted_1 = (TH1D*)unfold_reco_reweighted_1.Hunfold(RooUnfold::kErrors); h_unfold_reco_reweighted_1->SetName("h_unfold_reco_reweighted_1");
  RooUnfoldBayes unfold_reco_reweighted_2(resp_reco_reweighted, h_recojet_pt, 2); TH1D* h_unfold_reco_reweighted_2 = (TH1D*)unfold_reco_reweighted_2.Hunfold(RooUnfold::kErrors); h_unfold_reco_reweighted_2->SetName("h_unfold_reco_reweighted_2");
  RooUnfoldBayes unfold_reco_reweighted_3(resp_reco_reweighted, h_recojet_pt, 3); TH1D* h_unfold_reco_reweighted_3 = (TH1D*)unfold_reco_reweighted_3.Hunfold(RooUnfold::kErrors); h_unfold_reco_reweighted_3->SetName("h_unfold_reco_reweighted_3");
  RooUnfoldBayes unfold_reco_reweighted_4(resp_reco_reweighted, h_recojet_pt, 4); TH1D* h_unfold_reco_reweighted_4 = (TH1D*)unfold_reco_reweighted_4.Hunfold(RooUnfold::kErrors); h_unfold_reco_reweighted_4->SetName("h_unfold_reco_reweighted_4");
  RooUnfoldBayes unfold_reco_reweighted_5(resp_reco_reweighted, h_recojet_pt, 5); TH1D* h_unfold_reco_reweighted_5 = (TH1D*)unfold_reco_reweighted_5.Hunfold(RooUnfold::kErrors); h_unfold_reco_reweighted_5->SetName("h_unfold_reco_reweighted_5");
  RooUnfoldBayes unfold_reco_reweighted_6(resp_reco_reweighted, h_recojet_pt, 6); TH1D* h_unfold_reco_reweighted_6 = (TH1D*)unfold_reco_reweighted_6.Hunfold(RooUnfold::kErrors); h_unfold_reco_reweighted_6->SetName("h_unfold_reco_reweighted_6");

  RooUnfoldBayes unfold_calib_reweighted_1(resp_calib_reweighted, h_calibjet_pt, 1); TH1D* h_unfold_calib_reweighted_1 = (TH1D*)unfold_calib_reweighted_1.Hunfold(RooUnfold::kErrors); h_unfold_calib_reweighted_1->SetName("h_unfold_calib_reweighted_1");
  RooUnfoldBayes unfold_calib_reweighted_2(resp_calib_reweighted, h_calibjet_pt, 2); TH1D* h_unfold_calib_reweighted_2 = (TH1D*)unfold_calib_reweighted_2.Hunfold(RooUnfold::kErrors); h_unfold_calib_reweighted_2->SetName("h_unfold_calib_reweighted_2");
  RooUnfoldBayes unfold_calib_reweighted_3(resp_calib_reweighted, h_calibjet_pt, 3); TH1D* h_unfold_calib_reweighted_3 = (TH1D*)unfold_calib_reweighted_3.Hunfold(RooUnfold::kErrors); h_unfold_calib_reweighted_3->SetName("h_unfold_calib_reweighted_3");
  RooUnfoldBayes unfold_calib_reweighted_4(resp_calib_reweighted, h_calibjet_pt, 4); TH1D* h_unfold_calib_reweighted_4 = (TH1D*)unfold_calib_reweighted_4.Hunfold(RooUnfold::kErrors); h_unfold_calib_reweighted_4->SetName("h_unfold_calib_reweighted_4");
  RooUnfoldBayes unfold_calib_reweighted_5(resp_calib_reweighted, h_calibjet_pt, 5); TH1D* h_unfold_calib_reweighted_5 = (TH1D*)unfold_calib_reweighted_5.Hunfold(RooUnfold::kErrors); h_unfold_calib_reweighted_5->SetName("h_unfold_calib_reweighted_5");
  RooUnfoldBayes unfold_calib_reweighted_6(resp_calib_reweighted, h_calibjet_pt, 6); TH1D* h_unfold_calib_reweighted_6 = (TH1D*)unfold_calib_reweighted_6.Hunfold(RooUnfold::kErrors); h_unfold_calib_reweighted_6->SetName("h_unfold_calib_reweighted_6");

  //********** Efficiency Correction **********//
  for (int ib = 1; ib <= h_efficiency_reco->GetNbinsX(); ++ib) {
    double reco_efficiency = h_efficiency_reco->GetBinContent(ib);

    h_unfold_reco_1->SetBinContent(ib, h_unfold_reco_1->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_1->SetBinError(ib, h_unfold_reco_1->GetBinError(ib) / (double)reco_efficiency);
    h_unfold_reco_2->SetBinContent(ib, h_unfold_reco_2->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_2->SetBinError(ib, h_unfold_reco_2->GetBinError(ib) / (double)reco_efficiency);
    h_unfold_reco_3->SetBinContent(ib, h_unfold_reco_3->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_3->SetBinError(ib, h_unfold_reco_3->GetBinError(ib) / (double)reco_efficiency);
    h_unfold_reco_4->SetBinContent(ib, h_unfold_reco_4->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_4->SetBinError(ib, h_unfold_reco_4->GetBinError(ib) / (double)reco_efficiency);
    h_unfold_reco_5->SetBinContent(ib, h_unfold_reco_5->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_5->SetBinError(ib, h_unfold_reco_5->GetBinError(ib) / (double)reco_efficiency);
    h_unfold_reco_6->SetBinContent(ib, h_unfold_reco_6->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_6->SetBinError(ib, h_unfold_reco_6->GetBinError(ib) / (double)reco_efficiency);

    h_unfold_reco_reweighted_1->SetBinContent(ib, h_unfold_reco_reweighted_1->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_reweighted_1->SetBinError(ib, h_unfold_reco_reweighted_1->GetBinError(ib) / (double)reco_efficiency);
    h_unfold_reco_reweighted_2->SetBinContent(ib, h_unfold_reco_reweighted_2->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_reweighted_2->SetBinError(ib, h_unfold_reco_reweighted_2->GetBinError(ib) / (double)reco_efficiency);
    h_unfold_reco_reweighted_3->SetBinContent(ib, h_unfold_reco_reweighted_3->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_reweighted_3->SetBinError(ib, h_unfold_reco_reweighted_3->GetBinError(ib) / (double)reco_efficiency);
    h_unfold_reco_reweighted_4->SetBinContent(ib, h_unfold_reco_reweighted_4->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_reweighted_4->SetBinError(ib, h_unfold_reco_reweighted_4->GetBinError(ib) / (double)reco_efficiency);
    h_unfold_reco_reweighted_5->SetBinContent(ib, h_unfold_reco_reweighted_5->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_reweighted_5->SetBinError(ib, h_unfold_reco_reweighted_5->GetBinError(ib) / (double)reco_efficiency);
    h_unfold_reco_reweighted_6->SetBinContent(ib, h_unfold_reco_reweighted_6->GetBinContent(ib) / (double)reco_efficiency); h_unfold_reco_reweighted_6->SetBinError(ib, h_unfold_reco_reweighted_6->GetBinError(ib) / (double)reco_efficiency);

    double calib_efficiency = h_efficiency_calib->GetBinContent(ib);

    h_unfold_calib_1->SetBinContent(ib, h_unfold_calib_1->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_1->SetBinError(ib, h_unfold_calib_1->GetBinError(ib) / (double)calib_efficiency);
    h_unfold_calib_2->SetBinContent(ib, h_unfold_calib_2->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_2->SetBinError(ib, h_unfold_calib_2->GetBinError(ib) / (double)calib_efficiency);
    h_unfold_calib_3->SetBinContent(ib, h_unfold_calib_3->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_3->SetBinError(ib, h_unfold_calib_3->GetBinError(ib) / (double)calib_efficiency);
    h_unfold_calib_4->SetBinContent(ib, h_unfold_calib_4->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_4->SetBinError(ib, h_unfold_calib_4->GetBinError(ib) / (double)calib_efficiency);
    h_unfold_calib_5->SetBinContent(ib, h_unfold_calib_5->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_5->SetBinError(ib, h_unfold_calib_5->GetBinError(ib) / (double)calib_efficiency);
    h_unfold_calib_6->SetBinContent(ib, h_unfold_calib_6->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_6->SetBinError(ib, h_unfold_calib_6->GetBinError(ib) / (double)calib_efficiency);

    h_unfold_calib_reweighted_1->SetBinContent(ib, h_unfold_calib_reweighted_1->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_reweighted_1->SetBinError(ib, h_unfold_calib_reweighted_1->GetBinError(ib) / (double)calib_efficiency);
    h_unfold_calib_reweighted_2->SetBinContent(ib, h_unfold_calib_reweighted_2->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_reweighted_2->SetBinError(ib, h_unfold_calib_reweighted_2->GetBinError(ib) / (double)calib_efficiency);
    h_unfold_calib_reweighted_3->SetBinContent(ib, h_unfold_calib_reweighted_3->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_reweighted_3->SetBinError(ib, h_unfold_calib_reweighted_3->GetBinError(ib) / (double)calib_efficiency);
    h_unfold_calib_reweighted_4->SetBinContent(ib, h_unfold_calib_reweighted_4->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_reweighted_4->SetBinError(ib, h_unfold_calib_reweighted_4->GetBinError(ib) / (double)calib_efficiency);
    h_unfold_calib_reweighted_5->SetBinContent(ib, h_unfold_calib_reweighted_5->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_reweighted_5->SetBinError(ib, h_unfold_calib_reweighted_5->GetBinError(ib) / (double)calib_efficiency);
    h_unfold_calib_reweighted_6->SetBinContent(ib, h_unfold_calib_reweighted_6->GetBinContent(ib) / (double)calib_efficiency); h_unfold_calib_reweighted_6->SetBinError(ib, h_unfold_calib_reweighted_6->GetBinError(ib) / (double)calib_efficiency);
  }
 
  //********** Writing **********//
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();
  h_unfold_reco_1->Write();
  h_unfold_reco_2->Write();
  h_unfold_reco_3->Write();
  h_unfold_reco_4->Write();
  h_unfold_reco_5->Write();
  h_unfold_reco_6->Write();

  h_unfold_calib_1->Write();
  h_unfold_calib_2->Write();
  h_unfold_calib_3->Write();
  h_unfold_calib_4->Write();
  h_unfold_calib_5->Write();
  h_unfold_calib_6->Write();

  h_unfold_reco_reweighted_1->Write();
  h_unfold_reco_reweighted_2->Write();
  h_unfold_reco_reweighted_3->Write();
  h_unfold_reco_reweighted_4->Write();
  h_unfold_reco_reweighted_5->Write();
  h_unfold_reco_reweighted_6->Write();

  h_unfold_calib_reweighted_1->Write();
  h_unfold_calib_reweighted_2->Write();
  h_unfold_calib_reweighted_3->Write();
  h_unfold_calib_reweighted_4->Write();
  h_unfold_calib_reweighted_5->Write();
  h_unfold_calib_reweighted_6->Write();
  f_out->Close();
  std::cout << "All done!" << std::endl;
}