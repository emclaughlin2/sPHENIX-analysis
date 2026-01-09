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
#include "/sphenix/user/hanpuj/CaloDataAna24_skimmed/src/draw_template.C"  

void toy_errors(RooUnfoldBayes& unfold, TH2D* h_unfolded, int ntoys = 1000) {
    // Vectors to hold toy results
    std::vector<TVectorD> vx;   // unfolded values for each toy
    std::vector<TVectorD> vxe;  // errors for each toy (optional)
    std::vector<double> chi2;   // chi2 per toy (optional)

    // Run toys
    unfold.SetNToys(ntoys);
    unfold.RunToys(ntoys, vx, vxe, chi2);

    // Check toy results
    if (vx.empty()) {
        std::cerr << "No toy results returned!" << std::endl;
        return;
    }

    int nbins = vx[0].GetNoElements();

    // Calculate mean per bin over toys
    TVectorD mean(nbins);
    for (int i = 0; i < nbins; ++i) {
        double sum = 0;
        for (int t = 0; t < ntoys; ++t) {
            sum += vx[t][i];
        }
        mean[i] = sum / ntoys;
    }

    // Calculate variance per bin
    TVectorD variance(nbins);
    for (int i = 0; i < nbins; ++i) {
        double var = 0;
        for (int t = 0; t < ntoys; ++t) {
            double diff = vx[t][i] - mean[i];
            var += diff * diff;
        }
        variance[i] = var / (ntoys - 1);
    }

    // Apply errors to histogram
    int nx = h_unfolded->GetNbinsX();
    int ny = h_unfolded->GetNbinsY();
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            int flatIndex = (j - 1) * nx + (i - 1);  
            // statistical uncertainty is quad sum of unfolded data errors from kError 
            // and response matrix errors determined from variane of toys sampling
            h_unfolded->SetBinError(i, j, std::sqrt(h_unfolded->GetBinError(i,j)*h_unfolded->GetBinError(i,j) + variance[flatIndex]));
        }
    }
}


void do_unfolding() {
  SetAtlasStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  std::cout << gSystem->GetLibraries() << std::endl;
  std::cout << RooUnfoldResponse::Class()->GetImplFileName() << std::endl;

  //********** General Set up **********//
  const float PI = TMath::Pi();
  const float jet_radius = 0.4;
 
  //********** Files **********//
  TFile *f_sim = new TFile("analysis_sim_output/output_dijet_jet10_50_4bin_newetbin_trim.root", "READ");
  TFile *f_data = new TFile("analysis_data_output/output_0mrad_1.5mrad_4bin_newetbin_dijet.root", "READ");
  TFile *f_out = new TFile("output_unfolded_data_4bin_newetbin_dijet_1000toys.root", "RECREATE");

  // unfolding for full response matrix, trim 5 and trim 10 entries response matrices
  RooUnfoldResponse* h_respmatrix_calib = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet");
  RooUnfoldResponse* h_respmatrix_calib_trim_5 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_trim_5");
  RooUnfoldResponse* h_respmatrix_calib_trim_10 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_trim_10");
  TH2D* h_truth_calib_all_record = (TH2D*)f_sim->Get("h_truth_calib_dijet");
  TH2D* h_measure_calib_all_record = (TH2D*)f_data->Get("h_calibjet_pt_dijet_eff");

  int ntoys = 1000;
  RooUnfoldBayes unfold_calib_all_1(h_respmatrix_calib, h_measure_calib_all_record, 1, false, true);
  TH2D* h_unfold_calib_all_1 = (TH2D*)unfold_calib_all_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_1->SetName("h_unfold_calib_all_1");
  toy_errors(unfold_calib_all_1, h_unfold_calib_all_1, ntoys);
  RooUnfoldBayes unfold_calib_all_2(h_respmatrix_calib, h_measure_calib_all_record, 2, false, true);
  TH2D* h_unfold_calib_all_2 = (TH2D*)unfold_calib_all_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_2->SetName("h_unfold_calib_all_2");
  toy_errors(unfold_calib_all_2, h_unfold_calib_all_2, ntoys);
  RooUnfoldBayes unfold_calib_all_3(h_respmatrix_calib, h_measure_calib_all_record, 3, false, true);
  TH2D* h_unfold_calib_all_3 = (TH2D*)unfold_calib_all_3.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_3->SetName("h_unfold_calib_all_3");
  toy_errors(unfold_calib_all_3, h_unfold_calib_all_3, ntoys);
  RooUnfoldBayes unfold_calib_all_4(h_respmatrix_calib, h_measure_calib_all_record, 4, false, true);
  TH2D* h_unfold_calib_all_4 = (TH2D*)unfold_calib_all_4.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_4->SetName("h_unfold_calib_all_4");
  toy_errors(unfold_calib_all_4, h_unfold_calib_all_4, ntoys);
  RooUnfoldBayes unfold_calib_all_5(h_respmatrix_calib, h_measure_calib_all_record, 5, false, true);
  TH2D* h_unfold_calib_all_5 = (TH2D*)unfold_calib_all_5.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_5->SetName("h_unfold_calib_all_5");
  toy_errors(unfold_calib_all_5, h_unfold_calib_all_5, ntoys);
  RooUnfoldBayes unfold_calib_all_6(h_respmatrix_calib, h_measure_calib_all_record, 6, false, true);
  TH2D* h_unfold_calib_all_6 = (TH2D*)unfold_calib_all_6.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_6->SetName("h_unfold_calib_all_6");
  toy_errors(unfold_calib_all_6, h_unfold_calib_all_6, ntoys);
  RooUnfoldBayes unfold_calib_all_7(h_respmatrix_calib, h_measure_calib_all_record, 7, false, true);
  TH2D* h_unfold_calib_all_7 = (TH2D*)unfold_calib_all_7.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_7->SetName("h_unfold_calib_all_7");
  toy_errors(unfold_calib_all_7, h_unfold_calib_all_7, ntoys);
  RooUnfoldBayes unfold_calib_all_8(h_respmatrix_calib, h_measure_calib_all_record, 8, false, true);
  TH2D* h_unfold_calib_all_8 = (TH2D*)unfold_calib_all_8.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_8->SetName("h_unfold_calib_all_8");
  toy_errors(unfold_calib_all_8, h_unfold_calib_all_8, ntoys);
  RooUnfoldBayes unfold_calib_all_9(h_respmatrix_calib, h_measure_calib_all_record, 9, false, true);
  TH2D* h_unfold_calib_all_9 = (TH2D*)unfold_calib_all_9.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_9->SetName("h_unfold_calib_all_9");
  toy_errors(unfold_calib_all_9, h_unfold_calib_all_9, ntoys);
  RooUnfoldBayes unfold_calib_all_10(h_respmatrix_calib, h_measure_calib_all_record, 10, false, true);
  TH2D* h_unfold_calib_all_10 = (TH2D*)unfold_calib_all_10.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_10->SetName("h_unfold_calib_all_10");
  toy_errors(unfold_calib_all_10, h_unfold_calib_all_10, ntoys);
  RooUnfoldBayes unfold_calib_all_11(h_respmatrix_calib, h_measure_calib_all_record, 11, false, true);
  TH2D* h_unfold_calib_all_11 = (TH2D*)unfold_calib_all_11.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_11->SetName("h_unfold_calib_all_11");
  toy_errors(unfold_calib_all_11, h_unfold_calib_all_11, ntoys);
  RooUnfoldBayes unfold_calib_all_12(h_respmatrix_calib, h_measure_calib_all_record, 12, false, true);
  TH2D* h_unfold_calib_all_12 = (TH2D*)unfold_calib_all_12.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_12->SetName("h_unfold_calib_all_12");
  toy_errors(unfold_calib_all_12, h_unfold_calib_all_12, ntoys);
  RooUnfoldBayes unfold_calib_all_13(h_respmatrix_calib, h_measure_calib_all_record, 13, false, true);
  TH2D* h_unfold_calib_all_13 = (TH2D*)unfold_calib_all_13.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_13->SetName("h_unfold_calib_all_13");
  toy_errors(unfold_calib_all_13, h_unfold_calib_all_13, ntoys);
  RooUnfoldBayes unfold_calib_all_14(h_respmatrix_calib, h_measure_calib_all_record, 14, false, true);
  TH2D* h_unfold_calib_all_14 = (TH2D*)unfold_calib_all_14.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_14->SetName("h_unfold_calib_all_14");
  toy_errors(unfold_calib_all_14, h_unfold_calib_all_14, ntoys);
  RooUnfoldBayes unfold_calib_all_15(h_respmatrix_calib, h_measure_calib_all_record, 15, false, true);
  TH2D* h_unfold_calib_all_15 = (TH2D*)unfold_calib_all_15.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_15->SetName("h_unfold_calib_all_15");
  toy_errors(unfold_calib_all_15, h_unfold_calib_all_15, ntoys);
  RooUnfoldBayes unfold_calib_all_16(h_respmatrix_calib, h_measure_calib_all_record, 16, false, true);
  TH2D* h_unfold_calib_all_16 = (TH2D*)unfold_calib_all_16.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_16->SetName("h_unfold_calib_all_16");
  toy_errors(unfold_calib_all_16, h_unfold_calib_all_16, ntoys);
  RooUnfoldBayes unfold_calib_all_17(h_respmatrix_calib, h_measure_calib_all_record, 17, false, true);
  TH2D* h_unfold_calib_all_17 = (TH2D*)unfold_calib_all_17.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_17->SetName("h_unfold_calib_all_17");
  toy_errors(unfold_calib_all_17, h_unfold_calib_all_17, ntoys);
  RooUnfoldBayes unfold_calib_all_18(h_respmatrix_calib, h_measure_calib_all_record, 18, false, true);
  TH2D* h_unfold_calib_all_18 = (TH2D*)unfold_calib_all_18.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_18->SetName("h_unfold_calib_all_18");
  toy_errors(unfold_calib_all_18, h_unfold_calib_all_18, ntoys);
  RooUnfoldBayes unfold_calib_all_19(h_respmatrix_calib, h_measure_calib_all_record, 19, false, true);
  TH2D* h_unfold_calib_all_19 = (TH2D*)unfold_calib_all_19.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_19->SetName("h_unfold_calib_all_19");
  toy_errors(unfold_calib_all_19, h_unfold_calib_all_19, ntoys);
  RooUnfoldBayes unfold_calib_all_20(h_respmatrix_calib, h_measure_calib_all_record, 20, false, true);
  TH2D* h_unfold_calib_all_20 = (TH2D*)unfold_calib_all_20.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_all_20->SetName("h_unfold_calib_all_20");
  toy_errors(unfold_calib_all_20, h_unfold_calib_all_20, ntoys);

  RooUnfoldBayes unfold_calib_trim_5_1(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 1, false, true);
  TH2D* h_unfold_calib_trim_5_1 = (TH2D*)unfold_calib_trim_5_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_1->SetName("h_unfold_calib_trim_5_1");
  toy_errors(unfold_calib_trim_5_1, h_unfold_calib_trim_5_1, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_2(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 2, false, true);
  TH2D* h_unfold_calib_trim_5_2 = (TH2D*)unfold_calib_trim_5_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_2->SetName("h_unfold_calib_trim_5_2");
  toy_errors(unfold_calib_trim_5_2, h_unfold_calib_trim_5_2, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_3(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 3, false, true);
  TH2D* h_unfold_calib_trim_5_3 = (TH2D*)unfold_calib_trim_5_3.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_3->SetName("h_unfold_calib_trim_5_3");
  toy_errors(unfold_calib_trim_5_3, h_unfold_calib_trim_5_3, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_4(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 4, false, true);
  TH2D* h_unfold_calib_trim_5_4 = (TH2D*)unfold_calib_trim_5_4.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_4->SetName("h_unfold_calib_trim_5_4");
  toy_errors(unfold_calib_trim_5_4, h_unfold_calib_trim_5_4, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_5(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 5, false, true);
  TH2D* h_unfold_calib_trim_5_5 = (TH2D*)unfold_calib_trim_5_5.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_5->SetName("h_unfold_calib_trim_5_5");
  toy_errors(unfold_calib_trim_5_5, h_unfold_calib_trim_5_5, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_6(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 6, false, true);
  TH2D* h_unfold_calib_trim_5_6 = (TH2D*)unfold_calib_trim_5_6.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_6->SetName("h_unfold_calib_trim_5_6");
  toy_errors(unfold_calib_trim_5_6, h_unfold_calib_trim_5_6, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_7(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 7, false, true);
  TH2D* h_unfold_calib_trim_5_7 = (TH2D*)unfold_calib_trim_5_7.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_7->SetName("h_unfold_calib_trim_5_7");
  toy_errors(unfold_calib_trim_5_7, h_unfold_calib_trim_5_7, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_8(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 8, false, true);
  TH2D* h_unfold_calib_trim_5_8 = (TH2D*)unfold_calib_trim_5_8.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_8->SetName("h_unfold_calib_trim_5_8");
  toy_errors(unfold_calib_trim_5_8, h_unfold_calib_trim_5_8, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_9(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 9, false, true);
  TH2D* h_unfold_calib_trim_5_9 = (TH2D*)unfold_calib_trim_5_9.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_9->SetName("h_unfold_calib_trim_5_9");
  toy_errors(unfold_calib_trim_5_9, h_unfold_calib_trim_5_9, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_10(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 10, false, true);
  TH2D* h_unfold_calib_trim_5_10 = (TH2D*)unfold_calib_trim_5_10.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_10->SetName("h_unfold_calib_trim_5_10");
  toy_errors(unfold_calib_trim_5_10, h_unfold_calib_trim_5_10, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_11(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 11, false, true);
  TH2D* h_unfold_calib_trim_5_11 = (TH2D*)unfold_calib_trim_5_11.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_11->SetName("h_unfold_calib_trim_5_11");
  toy_errors(unfold_calib_trim_5_11, h_unfold_calib_trim_5_11, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_12(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 12, false, true);
  TH2D* h_unfold_calib_trim_5_12 = (TH2D*)unfold_calib_trim_5_12.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_12->SetName("h_unfold_calib_trim_5_12");
  toy_errors(unfold_calib_trim_5_12, h_unfold_calib_trim_5_12, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_13(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 13, false, true);
  TH2D* h_unfold_calib_trim_5_13 = (TH2D*)unfold_calib_trim_5_13.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_13->SetName("h_unfold_calib_trim_5_13");
  toy_errors(unfold_calib_trim_5_13, h_unfold_calib_trim_5_13, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_14(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 14, false, true);
  TH2D* h_unfold_calib_trim_5_14 = (TH2D*)unfold_calib_trim_5_14.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_14->SetName("h_unfold_calib_trim_5_14");
  toy_errors(unfold_calib_trim_5_14, h_unfold_calib_trim_5_14, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_15(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 15, false, true);
  TH2D* h_unfold_calib_trim_5_15 = (TH2D*)unfold_calib_trim_5_15.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_15->SetName("h_unfold_calib_trim_5_15");
  toy_errors(unfold_calib_trim_5_15, h_unfold_calib_trim_5_15, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_16(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 16, false, true);
  TH2D* h_unfold_calib_trim_5_16 = (TH2D*)unfold_calib_trim_5_16.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_16->SetName("h_unfold_calib_trim_5_16");
  toy_errors(unfold_calib_trim_5_16, h_unfold_calib_trim_5_16, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_17(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 17, false, true);
  TH2D* h_unfold_calib_trim_5_17 = (TH2D*)unfold_calib_trim_5_17.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_17->SetName("h_unfold_calib_trim_5_17");
  toy_errors(unfold_calib_trim_5_17, h_unfold_calib_trim_5_17, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_18(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 18, false, true);
  TH2D* h_unfold_calib_trim_5_18 = (TH2D*)unfold_calib_trim_5_18.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_18->SetName("h_unfold_calib_trim_5_18");
  toy_errors(unfold_calib_trim_5_18, h_unfold_calib_trim_5_18, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_19(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 19, false, true);
  TH2D* h_unfold_calib_trim_5_19 = (TH2D*)unfold_calib_trim_5_19.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_19->SetName("h_unfold_calib_trim_5_19");
  toy_errors(unfold_calib_trim_5_19, h_unfold_calib_trim_5_19, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_20(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 20, false, true);
  TH2D* h_unfold_calib_trim_5_20 = (TH2D*)unfold_calib_trim_5_20.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_20->SetName("h_unfold_calib_trim_5_20");
  toy_errors(unfold_calib_trim_5_20, h_unfold_calib_trim_5_20, ntoys);

  RooUnfoldBayes unfold_calib_trim_10_1(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 1, false, true);
  TH2D* h_unfold_calib_trim_10_1 = (TH2D*)unfold_calib_trim_10_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_1->SetName("h_unfold_calib_trim_10_1");
  toy_errors(unfold_calib_trim_10_1, h_unfold_calib_trim_10_1, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_2(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 2, false, true);
  TH2D* h_unfold_calib_trim_10_2 = (TH2D*)unfold_calib_trim_10_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_2->SetName("h_unfold_calib_trim_10_2");
  toy_errors(unfold_calib_trim_10_2, h_unfold_calib_trim_10_2, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_3(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 3, false, true);
  TH2D* h_unfold_calib_trim_10_3 = (TH2D*)unfold_calib_trim_10_3.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_3->SetName("h_unfold_calib_trim_10_3");
  toy_errors(unfold_calib_trim_10_3, h_unfold_calib_trim_10_3, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_4(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 4, false, true);
  TH2D* h_unfold_calib_trim_10_4 = (TH2D*)unfold_calib_trim_10_4.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_4->SetName("h_unfold_calib_trim_10_4");
  toy_errors(unfold_calib_trim_10_4, h_unfold_calib_trim_10_4, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_5(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 5, false, true);
  TH2D* h_unfold_calib_trim_10_5 = (TH2D*)unfold_calib_trim_10_5.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_5->SetName("h_unfold_calib_trim_10_5");
  toy_errors(unfold_calib_trim_10_5, h_unfold_calib_trim_10_5, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_6(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 6, false, true);
  TH2D* h_unfold_calib_trim_10_6 = (TH2D*)unfold_calib_trim_10_6.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_6->SetName("h_unfold_calib_trim_10_6");
  toy_errors(unfold_calib_trim_10_6, h_unfold_calib_trim_10_6, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_7(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 7, false, true);
  TH2D* h_unfold_calib_trim_10_7 = (TH2D*)unfold_calib_trim_10_7.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_7->SetName("h_unfold_calib_trim_10_7");
  toy_errors(unfold_calib_trim_10_7, h_unfold_calib_trim_10_7, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_8(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 8, false, true);
  TH2D* h_unfold_calib_trim_10_8 = (TH2D*)unfold_calib_trim_10_8.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_8->SetName("h_unfold_calib_trim_10_8");
  toy_errors(unfold_calib_trim_10_8, h_unfold_calib_trim_10_8, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_9(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 9, false, true);
  TH2D* h_unfold_calib_trim_10_9 = (TH2D*)unfold_calib_trim_10_9.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_9->SetName("h_unfold_calib_trim_10_9");
  toy_errors(unfold_calib_trim_10_9, h_unfold_calib_trim_10_9, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_10(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 10, false, true);
  TH2D* h_unfold_calib_trim_10_10 = (TH2D*)unfold_calib_trim_10_10.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_10->SetName("h_unfold_calib_trim_10_10");
  toy_errors(unfold_calib_trim_10_10, h_unfold_calib_trim_10_10, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_11(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 11, false, true);
  TH2D* h_unfold_calib_trim_10_11 = (TH2D*)unfold_calib_trim_10_11.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_11->SetName("h_unfold_calib_trim_10_11");
  toy_errors(unfold_calib_trim_10_11, h_unfold_calib_trim_10_11, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_12(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 12, false, true);
  TH2D* h_unfold_calib_trim_10_12 = (TH2D*)unfold_calib_trim_10_12.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_12->SetName("h_unfold_calib_trim_10_12");
  toy_errors(unfold_calib_trim_10_12, h_unfold_calib_trim_10_12, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_13(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 13, false, true);
  TH2D* h_unfold_calib_trim_10_13 = (TH2D*)unfold_calib_trim_10_13.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_13->SetName("h_unfold_calib_trim_10_13");
  toy_errors(unfold_calib_trim_10_13, h_unfold_calib_trim_10_13, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_14(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 14, false, true);
  TH2D* h_unfold_calib_trim_10_14 = (TH2D*)unfold_calib_trim_10_14.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_14->SetName("h_unfold_calib_trim_10_14");
  toy_errors(unfold_calib_trim_10_14, h_unfold_calib_trim_10_14, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_15(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 15, false, true);
  TH2D* h_unfold_calib_trim_10_15 = (TH2D*)unfold_calib_trim_10_15.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_15->SetName("h_unfold_calib_trim_10_15");
  toy_errors(unfold_calib_trim_10_15, h_unfold_calib_trim_10_15, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_16(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 16, false, true);
  TH2D* h_unfold_calib_trim_10_16 = (TH2D*)unfold_calib_trim_10_16.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_16->SetName("h_unfold_calib_trim_10_16");
  toy_errors(unfold_calib_trim_10_16, h_unfold_calib_trim_10_16, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_17(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 17, false, true);
  TH2D* h_unfold_calib_trim_10_17 = (TH2D*)unfold_calib_trim_10_17.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_17->SetName("h_unfold_calib_trim_10_17");
  toy_errors(unfold_calib_trim_10_17, h_unfold_calib_trim_10_17, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_18(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 18, false, true);
  TH2D* h_unfold_calib_trim_10_18 = (TH2D*)unfold_calib_trim_10_18.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_18->SetName("h_unfold_calib_trim_10_18");
  toy_errors(unfold_calib_trim_10_18, h_unfold_calib_trim_10_18, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_19(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 19, false, true);
  TH2D* h_unfold_calib_trim_10_19 = (TH2D*)unfold_calib_trim_10_19.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_19->SetName("h_unfold_calib_trim_10_19");
  toy_errors(unfold_calib_trim_10_19, h_unfold_calib_trim_10_19, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_20(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 20, false, true);
  TH2D* h_unfold_calib_trim_10_20 = (TH2D*)unfold_calib_trim_10_20.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_20->SetName("h_unfold_calib_trim_10_20");
  toy_errors(unfold_calib_trim_10_20, h_unfold_calib_trim_10_20, ntoys);

  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();
  h_truth_calib_all_record->Write();
  h_measure_calib_all_record->Write();
  h_unfold_calib_all_1->Write(); h_unfold_calib_all_2->Write(); h_unfold_calib_all_3->Write(); h_unfold_calib_all_4->Write(); h_unfold_calib_all_5->Write();
  h_unfold_calib_all_6->Write(); h_unfold_calib_all_7->Write(); h_unfold_calib_all_8->Write(); h_unfold_calib_all_9->Write(); h_unfold_calib_all_10->Write();
  h_unfold_calib_all_11->Write(); h_unfold_calib_all_12->Write(); h_unfold_calib_all_13->Write(); h_unfold_calib_all_14->Write(); h_unfold_calib_all_15->Write();
  h_unfold_calib_all_16->Write(); h_unfold_calib_all_17->Write(); h_unfold_calib_all_18->Write(); h_unfold_calib_all_19->Write(); h_unfold_calib_all_20->Write();
  h_unfold_calib_trim_5_1->Write(); h_unfold_calib_trim_5_2->Write(); h_unfold_calib_trim_5_3->Write(); h_unfold_calib_trim_5_4->Write(); h_unfold_calib_trim_5_5->Write();
  h_unfold_calib_trim_5_6->Write(); h_unfold_calib_trim_5_7->Write(); h_unfold_calib_trim_5_8->Write(); h_unfold_calib_trim_5_9->Write(); h_unfold_calib_trim_5_10->Write();
  h_unfold_calib_trim_5_11->Write(); h_unfold_calib_trim_5_12->Write(); h_unfold_calib_trim_5_13->Write(); h_unfold_calib_trim_5_14->Write(); h_unfold_calib_trim_5_15->Write();
  h_unfold_calib_trim_5_16->Write(); h_unfold_calib_trim_5_17->Write(); h_unfold_calib_trim_5_18->Write(); h_unfold_calib_trim_5_19->Write(); h_unfold_calib_trim_5_20->Write();
  h_unfold_calib_trim_10_1->Write(); h_unfold_calib_trim_10_2->Write(); h_unfold_calib_trim_10_3->Write(); h_unfold_calib_trim_10_4->Write(); h_unfold_calib_trim_10_5->Write();
  h_unfold_calib_trim_10_6->Write(); h_unfold_calib_trim_10_7->Write(); h_unfold_calib_trim_10_8->Write(); h_unfold_calib_trim_10_9->Write(); h_unfold_calib_trim_10_10->Write();
  h_unfold_calib_trim_10_11->Write(); h_unfold_calib_trim_10_12->Write(); h_unfold_calib_trim_10_13->Write(); h_unfold_calib_trim_10_14->Write(); h_unfold_calib_trim_10_15->Write();
  h_unfold_calib_trim_10_16->Write(); h_unfold_calib_trim_10_17->Write(); h_unfold_calib_trim_10_18->Write(); h_unfold_calib_trim_10_19->Write(); h_unfold_calib_trim_10_20->Write();

  f_out->Close();
  std::cout << "All done!" << std::endl;
}