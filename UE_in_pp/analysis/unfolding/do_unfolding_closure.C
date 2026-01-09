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

void do_unfolding_closure() {
  SetAtlasStyle();
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  std::cout << gSystem->GetLibraries() << std::endl;
  std::cout << RooUnfoldResponse::Class()->GetImplFileName() << std::endl;

  //********** General Set up **********//
  const float PI = TMath::Pi();
  const float jet_radius = 0.4;

  std::vector<TH1F*> h_input;
  std::vector<int> color;
  std::vector<int> markerstyle;
  std::vector<std::string> text;
  std::vector<std::string> legend;
 
  //********** Files **********//
  TFile *f_sim = new TFile("analysis_sim_output/output_dijet_jet10_50_4bin_newetbin_trim.root", "READ");
  TFile *f_out = new TFile("output_closure_sim_4bin_newetbin_dijet_1000toys.root", "RECREATE");

  RooUnfoldResponse* h_respmatrix_calib = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet");
  RooUnfoldResponse* h_respmatrix_calib_trim_5 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_trim_5");
  RooUnfoldResponse* h_respmatrix_calib_trim_10 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_trim_10");
  TH2D* h_meas_calib_all = (TH2D*)h_respmatrix_calib->Hmeasured(); h_meas_calib_all->SetName("h_meas_calib_all");
  TH2D* h_truth_calib_all = (TH2D*)h_respmatrix_calib->Htruth(); h_truth_calib_all->SetName("h_truth_calib_all");
  TH2D* h_truth_calib_all_record = (TH2D*)f_sim->Get("h_truth_calib_dijet");
  TH2D* h_measure_calib_all_record = (TH2D*)f_sim->Get("h_measure_calib_dijet");

  int ntoys = 1000;
  RooUnfoldBayes unfold_calib_full_1(h_respmatrix_calib, h_meas_calib_all, 1, false, true);
  TH2D* h_unfold_calib_full_1 = (TH2D*)unfold_calib_full_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_full_1->SetName("h_unfold_calib_full_1");
  toy_errors(unfold_calib_full_1, h_unfold_calib_full_1, ntoys);
  RooUnfoldBayes unfold_calib_full_2(h_respmatrix_calib, h_meas_calib_all, 2, false, true);
  TH2D* h_unfold_calib_full_2 = (TH2D*)unfold_calib_full_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_full_2->SetName("h_unfold_calib_full_2");
  toy_errors(unfold_calib_full_2, h_unfold_calib_full_2, ntoys);
  RooUnfoldBayes unfold_calib_full_extra_1(h_respmatrix_calib, h_measure_calib_all_record, 1, false, true);
  TH2D* h_unfold_calib_full_extra_1 = (TH2D*)unfold_calib_full_extra_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_full_extra_1->SetName("h_unfold_calib_full_extra_1");
  toy_errors(unfold_calib_full_extra_1, h_unfold_calib_full_extra_1, ntoys);
  RooUnfoldBayes unfold_calib_full_extra_2(h_respmatrix_calib, h_measure_calib_all_record, 2, false, true);
  TH2D* h_unfold_calib_full_extra_2 = (TH2D*)unfold_calib_full_extra_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_full_extra_2->SetName("h_unfold_calib_full_extra_2");
  toy_errors(unfold_calib_full_extra_2, h_unfold_calib_full_extra_2, ntoys);

  RooUnfoldBayes unfold_calib_trim_5_full_1(h_respmatrix_calib_trim_5, h_meas_calib_all, 1, false, true);
  TH2D* h_unfold_calib_trim_5_full_1 = (TH2D*)unfold_calib_trim_5_full_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_full_1->SetName("h_unfold_calib_trim_5_full_1");
  toy_errors(unfold_calib_trim_5_full_1, h_unfold_calib_trim_5_full_1, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_full_2(h_respmatrix_calib_trim_5, h_meas_calib_all, 2, false, true);
  TH2D* h_unfold_calib_trim_5_full_2 = (TH2D*)unfold_calib_trim_5_full_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_full_2->SetName("h_unfold_calib_trim_5_full_2");
  toy_errors(unfold_calib_trim_5_full_2, h_unfold_calib_trim_5_full_2, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_full_extra_1(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 1, false, true);
  TH2D* h_unfold_calib_trim_5_full_extra_1 = (TH2D*)unfold_calib_trim_5_full_extra_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_full_extra_1->SetName("h_unfold_calib_trim_5_full_extra_1");
  toy_errors(unfold_calib_trim_5_full_extra_1, h_unfold_calib_trim_5_full_extra_1, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_full_extra_2(h_respmatrix_calib_trim_5, h_measure_calib_all_record, 2, false, true);
  TH2D* h_unfold_calib_trim_5_full_extra_2 = (TH2D*)unfold_calib_trim_5_full_extra_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_full_extra_2->SetName("h_unfold_calib_trim_5_full_extra_2");
  toy_errors(unfold_calib_trim_5_full_extra_2, h_unfold_calib_trim_5_full_extra_2, ntoys);

  RooUnfoldBayes unfold_calib_trim_10_full_1(h_respmatrix_calib_trim_10, h_meas_calib_all, 1, false, true);
  TH2D* h_unfold_calib_trim_10_full_1 = (TH2D*)unfold_calib_trim_10_full_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_full_1->SetName("h_unfold_calib_trim_10_full_1");
  toy_errors(unfold_calib_trim_10_full_1, h_unfold_calib_trim_10_full_1, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_full_2(h_respmatrix_calib_trim_10, h_meas_calib_all, 2, false, true);
  TH2D* h_unfold_calib_trim_10_full_2 = (TH2D*)unfold_calib_trim_10_full_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_full_2->SetName("h_unfold_calib_trim_10_full_2");
  toy_errors(unfold_calib_trim_10_full_2, h_unfold_calib_trim_10_full_2, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_full_extra_1(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 1, false, true);
  TH2D* h_unfold_calib_trim_10_full_extra_1 = (TH2D*)unfold_calib_trim_10_full_extra_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_full_extra_1->SetName("h_unfold_calib_trim_10_full_extra_1");
  toy_errors(unfold_calib_trim_10_full_extra_1, h_unfold_calib_trim_10_full_extra_1, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_full_extra_2(h_respmatrix_calib_trim_10, h_measure_calib_all_record, 2, false, true);
  TH2D* h_unfold_calib_trim_10_full_extra_2 = (TH2D*)unfold_calib_trim_10_full_extra_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_full_extra_2->SetName("h_unfold_calib_trim_10_full_extra_2");
  toy_errors(unfold_calib_trim_10_full_extra_2, h_unfold_calib_trim_10_full_extra_2, ntoys);
  
  // Half closure test on calib //
  RooUnfoldResponse* h_respmatrix_calib_half1 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_half1");
  RooUnfoldResponse* h_respmatrix_calib_half1_trim_5 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_half1_trim_5");
  RooUnfoldResponse* h_respmatrix_calib_half1_trim_10 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_half1_trim_10");
  TH2D* h_truth_calib_dijet_half1_record = (TH2D*)f_sim->Get("h_truth_calib_dijet_half1");
  TH2D* h_measure_calib_dijet_half1_record = (TH2D*)f_sim->Get("h_measure_calib_dijet_half1");
  RooUnfoldResponse* h_respmatrix_calib_half2 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_half2");
  RooUnfoldResponse* h_respmatrix_calib_half2_trim_5 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_half2_trim_5");
  RooUnfoldResponse* h_respmatrix_calib_half2_trim_10 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_half2_trim_10");
  TH2D* h_truth_calib_dijet_half2_record = (TH2D*)f_sim->Get("h_truth_calib_dijet_half2");
  TH2D* h_measure_calib_dijet_half2_record = (TH2D*)f_sim->Get("h_measure_calib_dijet_half2");
  
  RooUnfoldBayes unfold_calib_half_1(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 1, false, true);
  TH2D* h_unfold_calib_half_1 = (TH2D*)unfold_calib_half_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_1->SetName("h_unfold_calib_half_1");
  toy_errors(unfold_calib_half_1, h_unfold_calib_half_1, ntoys);
  RooUnfoldBayes unfold_calib_half_2(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 2, false, true);
  TH2D* h_unfold_calib_half_2 = (TH2D*)unfold_calib_half_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_2->SetName("h_unfold_calib_half_2");
  toy_errors(unfold_calib_half_2, h_unfold_calib_half_2, ntoys);
  RooUnfoldBayes unfold_calib_half_3(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 3, false, true);
  TH2D* h_unfold_calib_half_3 = (TH2D*)unfold_calib_half_3.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_3->SetName("h_unfold_calib_half_3");
  toy_errors(unfold_calib_half_3, h_unfold_calib_half_3, ntoys);
  RooUnfoldBayes unfold_calib_half_4(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 4, false, true);
  TH2D* h_unfold_calib_half_4 = (TH2D*)unfold_calib_half_4.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_4->SetName("h_unfold_calib_half_4");
  toy_errors(unfold_calib_half_4, h_unfold_calib_half_4, ntoys);
  RooUnfoldBayes unfold_calib_half_5(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 5, false, true);
  TH2D* h_unfold_calib_half_5 = (TH2D*)unfold_calib_half_5.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_5->SetName("h_unfold_calib_half_5");
  toy_errors(unfold_calib_half_5, h_unfold_calib_half_5, ntoys);
  RooUnfoldBayes unfold_calib_half_6(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 6, false, true);
  TH2D* h_unfold_calib_half_6 = (TH2D*)unfold_calib_half_6.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_6->SetName("h_unfold_calib_half_6");
  toy_errors(unfold_calib_half_6, h_unfold_calib_half_6, ntoys);
  RooUnfoldBayes unfold_calib_half_7(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 7, false, true);
  TH2D* h_unfold_calib_half_7 = (TH2D*)unfold_calib_half_7.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_7->SetName("h_unfold_calib_half_7");
  toy_errors(unfold_calib_half_7, h_unfold_calib_half_7, ntoys);
  RooUnfoldBayes unfold_calib_half_8(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 8, false, true);
  TH2D* h_unfold_calib_half_8 = (TH2D*)unfold_calib_half_8.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_8->SetName("h_unfold_calib_half_8");
  toy_errors(unfold_calib_half_8, h_unfold_calib_half_8, ntoys);
  RooUnfoldBayes unfold_calib_half_9(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 9, false, true);
  TH2D* h_unfold_calib_half_9 = (TH2D*)unfold_calib_half_9.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_9->SetName("h_unfold_calib_half_9");
  toy_errors(unfold_calib_half_9, h_unfold_calib_half_9, ntoys);
  RooUnfoldBayes unfold_calib_half_10(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 10, false, true);
  TH2D* h_unfold_calib_half_10 = (TH2D*)unfold_calib_half_10.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_10->SetName("h_unfold_calib_half_10");
  toy_errors(unfold_calib_half_10, h_unfold_calib_half_10, ntoys);
  RooUnfoldBayes unfold_calib_half_11(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 11, false, true);
  TH2D* h_unfold_calib_half_11 = (TH2D*)unfold_calib_half_11.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_11->SetName("h_unfold_calib_half_11");
  toy_errors(unfold_calib_half_11, h_unfold_calib_half_11, ntoys);
  RooUnfoldBayes unfold_calib_half_12(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 12, false, true);
  TH2D* h_unfold_calib_half_12 = (TH2D*)unfold_calib_half_12.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_12->SetName("h_unfold_calib_half_12");
  toy_errors(unfold_calib_half_12, h_unfold_calib_half_12, ntoys);
  RooUnfoldBayes unfold_calib_half_13(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 13, false, true);
  TH2D* h_unfold_calib_half_13 = (TH2D*)unfold_calib_half_13.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_13->SetName("h_unfold_calib_half_13");
  toy_errors(unfold_calib_half_13, h_unfold_calib_half_13, ntoys);
  RooUnfoldBayes unfold_calib_half_14(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 14, false, true);
  TH2D* h_unfold_calib_half_14 = (TH2D*)unfold_calib_half_14.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_14->SetName("h_unfold_calib_half_14");
  toy_errors(unfold_calib_half_14, h_unfold_calib_half_14, ntoys);
  RooUnfoldBayes unfold_calib_half_15(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 15, false, true);
  TH2D* h_unfold_calib_half_15 = (TH2D*)unfold_calib_half_15.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_15->SetName("h_unfold_calib_half_15");
  toy_errors(unfold_calib_half_15, h_unfold_calib_half_15, ntoys);
  RooUnfoldBayes unfold_calib_half_16(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 16, false, true);
  TH2D* h_unfold_calib_half_16 = (TH2D*)unfold_calib_half_16.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_16->SetName("h_unfold_calib_half_16");
  toy_errors(unfold_calib_half_16, h_unfold_calib_half_16, ntoys);
  RooUnfoldBayes unfold_calib_half_17(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 17, false, true);
  TH2D* h_unfold_calib_half_17 = (TH2D*)unfold_calib_half_17.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_17->SetName("h_unfold_calib_half_17");
  toy_errors(unfold_calib_half_17, h_unfold_calib_half_17, ntoys);
  RooUnfoldBayes unfold_calib_half_18(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 18, false, true);
  TH2D* h_unfold_calib_half_18 = (TH2D*)unfold_calib_half_18.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_18->SetName("h_unfold_calib_half_18");
  toy_errors(unfold_calib_half_18, h_unfold_calib_half_18, ntoys);
  RooUnfoldBayes unfold_calib_half_19(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 19, false, true);
  TH2D* h_unfold_calib_half_19 = (TH2D*)unfold_calib_half_19.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_19->SetName("h_unfold_calib_half_19");
  toy_errors(unfold_calib_half_19, h_unfold_calib_half_19, ntoys);
  RooUnfoldBayes unfold_calib_half_20(h_respmatrix_calib_half1, h_measure_calib_dijet_half2_record, 20, false, true);
  TH2D* h_unfold_calib_half_20 = (TH2D*)unfold_calib_half_20.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_half_20->SetName("h_unfold_calib_half_20");
  toy_errors(unfold_calib_half_20, h_unfold_calib_half_20, ntoys);

  RooUnfoldBayes unfold_calib_trim_5_half_1(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 1, false, true);
  TH2D* h_unfold_calib_trim_5_half_1 = (TH2D*)unfold_calib_trim_5_half_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_1->SetName("h_unfold_calib_trim_5_half_1");
  toy_errors(unfold_calib_trim_5_half_1, h_unfold_calib_trim_5_half_1, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_2(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 2, false, true);
  TH2D* h_unfold_calib_trim_5_half_2 = (TH2D*)unfold_calib_trim_5_half_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_2->SetName("h_unfold_calib_trim_5_half_2");
  toy_errors(unfold_calib_trim_5_half_2, h_unfold_calib_trim_5_half_2, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_3(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 3, false, true);
  TH2D* h_unfold_calib_trim_5_half_3 = (TH2D*)unfold_calib_trim_5_half_3.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_3->SetName("h_unfold_calib_trim_5_half_3");
  toy_errors(unfold_calib_trim_5_half_3, h_unfold_calib_trim_5_half_3, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_4(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 4, false, true);
  TH2D* h_unfold_calib_trim_5_half_4 = (TH2D*)unfold_calib_trim_5_half_4.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_4->SetName("h_unfold_calib_trim_5_half_4");
  toy_errors(unfold_calib_trim_5_half_4, h_unfold_calib_trim_5_half_4, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_5(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 5, false, true);
  TH2D* h_unfold_calib_trim_5_half_5 = (TH2D*)unfold_calib_trim_5_half_5.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_5->SetName("h_unfold_calib_trim_5_half_5");
  toy_errors(unfold_calib_trim_5_half_5, h_unfold_calib_trim_5_half_5, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_6(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 6, false, true);
  TH2D* h_unfold_calib_trim_5_half_6 = (TH2D*)unfold_calib_trim_5_half_6.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_6->SetName("h_unfold_calib_trim_5_half_6");
  toy_errors(unfold_calib_trim_5_half_6, h_unfold_calib_trim_5_half_6, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_7(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 7, false, true);
  TH2D* h_unfold_calib_trim_5_half_7 = (TH2D*)unfold_calib_trim_5_half_7.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_7->SetName("h_unfold_calib_trim_5_half_7");
  toy_errors(unfold_calib_trim_5_half_7, h_unfold_calib_trim_5_half_7, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_8(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 8, false, true);
  TH2D* h_unfold_calib_trim_5_half_8 = (TH2D*)unfold_calib_trim_5_half_8.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_8->SetName("h_unfold_calib_trim_5_half_8");
  toy_errors(unfold_calib_trim_5_half_8, h_unfold_calib_trim_5_half_8, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_9(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 9, false, true);
  TH2D* h_unfold_calib_trim_5_half_9 = (TH2D*)unfold_calib_trim_5_half_9.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_9->SetName("h_unfold_calib_trim_5_half_9");
  toy_errors(unfold_calib_trim_5_half_9, h_unfold_calib_trim_5_half_9, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_10(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 10, false, true);
  TH2D* h_unfold_calib_trim_5_half_10 = (TH2D*)unfold_calib_trim_5_half_10.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_10->SetName("h_unfold_calib_trim_5_half_10");
  toy_errors(unfold_calib_trim_5_half_10, h_unfold_calib_trim_5_half_10, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_11(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 11, false, true);
  TH2D* h_unfold_calib_trim_5_half_11 = (TH2D*)unfold_calib_trim_5_half_11.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_11->SetName("h_unfold_calib_trim_5_half_11");
  toy_errors(unfold_calib_trim_5_half_11, h_unfold_calib_trim_5_half_11, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_12(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 12, false, true);
  TH2D* h_unfold_calib_trim_5_half_12 = (TH2D*)unfold_calib_trim_5_half_12.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_12->SetName("h_unfold_calib_trim_5_half_12");
  toy_errors(unfold_calib_trim_5_half_12, h_unfold_calib_trim_5_half_12, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_13(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 13, false, true);
  TH2D* h_unfold_calib_trim_5_half_13 = (TH2D*)unfold_calib_trim_5_half_13.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_13->SetName("h_unfold_calib_trim_5_half_13");
  toy_errors(unfold_calib_trim_5_half_13, h_unfold_calib_trim_5_half_13, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_14(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 14, false, true);
  TH2D* h_unfold_calib_trim_5_half_14 = (TH2D*)unfold_calib_trim_5_half_14.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_14->SetName("h_unfold_calib_trim_5_half_14");
  toy_errors(unfold_calib_trim_5_half_14, h_unfold_calib_trim_5_half_14, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_15(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 15, false, true);
  TH2D* h_unfold_calib_trim_5_half_15 = (TH2D*)unfold_calib_trim_5_half_15.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_15->SetName("h_unfold_calib_trim_5_half_15");
  toy_errors(unfold_calib_trim_5_half_15, h_unfold_calib_trim_5_half_15, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_16(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 16, false, true);
  TH2D* h_unfold_calib_trim_5_half_16 = (TH2D*)unfold_calib_trim_5_half_16.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_16->SetName("h_unfold_calib_trim_5_half_16");
  toy_errors(unfold_calib_trim_5_half_16, h_unfold_calib_trim_5_half_16, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_17(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 17, false, true);
  TH2D* h_unfold_calib_trim_5_half_17 = (TH2D*)unfold_calib_trim_5_half_17.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_17->SetName("h_unfold_calib_trim_5_half_17");
  toy_errors(unfold_calib_trim_5_half_17, h_unfold_calib_trim_5_half_17, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_18(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 18, false, true);
  TH2D* h_unfold_calib_trim_5_half_18 = (TH2D*)unfold_calib_trim_5_half_18.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_18->SetName("h_unfold_calib_trim_5_half_18");
  toy_errors(unfold_calib_trim_5_half_18, h_unfold_calib_trim_5_half_18, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_19(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 19, false, true);
  TH2D* h_unfold_calib_trim_5_half_19 = (TH2D*)unfold_calib_trim_5_half_19.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_19->SetName("h_unfold_calib_trim_5_half_19");
  toy_errors(unfold_calib_trim_5_half_19, h_unfold_calib_trim_5_half_19, ntoys);
  RooUnfoldBayes unfold_calib_trim_5_half_20(h_respmatrix_calib_half1_trim_5, h_measure_calib_dijet_half2_record, 20, false, true);
  TH2D* h_unfold_calib_trim_5_half_20 = (TH2D*)unfold_calib_trim_5_half_20.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_5_half_20->SetName("h_unfold_calib_trim_5_half_20");
  toy_errors(unfold_calib_trim_5_half_20, h_unfold_calib_trim_5_half_20, ntoys);

  RooUnfoldBayes unfold_calib_trim_10_half_1(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 1, false, true);
  TH2D* h_unfold_calib_trim_10_half_1 = (TH2D*)unfold_calib_trim_10_half_1.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_1->SetName("h_unfold_calib_trim_10_half_1");
  toy_errors(unfold_calib_trim_10_half_1, h_unfold_calib_trim_10_half_1, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_2(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 2, false, true);
  TH2D* h_unfold_calib_trim_10_half_2 = (TH2D*)unfold_calib_trim_10_half_2.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_2->SetName("h_unfold_calib_trim_10_half_2");
  toy_errors(unfold_calib_trim_10_half_2, h_unfold_calib_trim_10_half_2, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_3(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 3, false, true);
  TH2D* h_unfold_calib_trim_10_half_3 = (TH2D*)unfold_calib_trim_10_half_3.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_3->SetName("h_unfold_calib_trim_10_half_3");
  toy_errors(unfold_calib_trim_10_half_3, h_unfold_calib_trim_10_half_3, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_4(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 4, false, true);
  TH2D* h_unfold_calib_trim_10_half_4 = (TH2D*)unfold_calib_trim_10_half_4.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_4->SetName("h_unfold_calib_trim_10_half_4");
  toy_errors(unfold_calib_trim_10_half_4, h_unfold_calib_trim_10_half_4, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_5(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 5, false, true);
  TH2D* h_unfold_calib_trim_10_half_5 = (TH2D*)unfold_calib_trim_10_half_5.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_5->SetName("h_unfold_calib_trim_10_half_5");
  toy_errors(unfold_calib_trim_10_half_5, h_unfold_calib_trim_10_half_5, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_6(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 6, false, true);
  TH2D* h_unfold_calib_trim_10_half_6 = (TH2D*)unfold_calib_trim_10_half_6.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_6->SetName("h_unfold_calib_trim_10_half_6");
  toy_errors(unfold_calib_trim_10_half_6, h_unfold_calib_trim_10_half_6, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_7(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 7, false, true);
  TH2D* h_unfold_calib_trim_10_half_7 = (TH2D*)unfold_calib_trim_10_half_7.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_7->SetName("h_unfold_calib_trim_10_half_7");
  toy_errors(unfold_calib_trim_10_half_7, h_unfold_calib_trim_10_half_7, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_8(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 8, false, true);
  TH2D* h_unfold_calib_trim_10_half_8 = (TH2D*)unfold_calib_trim_10_half_8.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_8->SetName("h_unfold_calib_trim_10_half_8");
  toy_errors(unfold_calib_trim_10_half_8, h_unfold_calib_trim_10_half_8, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_9(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 9, false, true);
  TH2D* h_unfold_calib_trim_10_half_9 = (TH2D*)unfold_calib_trim_10_half_9.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_9->SetName("h_unfold_calib_trim_10_half_9");
  toy_errors(unfold_calib_trim_10_half_9, h_unfold_calib_trim_10_half_9, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_10(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 10, false, true);
  TH2D* h_unfold_calib_trim_10_half_10 = (TH2D*)unfold_calib_trim_10_half_10.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_10->SetName("h_unfold_calib_trim_10_half_10");
  toy_errors(unfold_calib_trim_10_half_10, h_unfold_calib_trim_10_half_10, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_11(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 11, false, true);
  TH2D* h_unfold_calib_trim_10_half_11 = (TH2D*)unfold_calib_trim_10_half_11.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_11->SetName("h_unfold_calib_trim_10_half_11");
  toy_errors(unfold_calib_trim_10_half_11, h_unfold_calib_trim_10_half_11, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_12(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 12, false, true);
  TH2D* h_unfold_calib_trim_10_half_12 = (TH2D*)unfold_calib_trim_10_half_12.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_12->SetName("h_unfold_calib_trim_10_half_12");
  toy_errors(unfold_calib_trim_10_half_12, h_unfold_calib_trim_10_half_12, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_13(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 13, false, true);
  TH2D* h_unfold_calib_trim_10_half_13 = (TH2D*)unfold_calib_trim_10_half_13.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_13->SetName("h_unfold_calib_trim_10_half_13");
  toy_errors(unfold_calib_trim_10_half_13, h_unfold_calib_trim_10_half_13, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_14(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 14, false, true);
  TH2D* h_unfold_calib_trim_10_half_14 = (TH2D*)unfold_calib_trim_10_half_14.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_14->SetName("h_unfold_calib_trim_10_half_14");
  toy_errors(unfold_calib_trim_10_half_14, h_unfold_calib_trim_10_half_14, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_15(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 15, false, true);
  TH2D* h_unfold_calib_trim_10_half_15 = (TH2D*)unfold_calib_trim_10_half_15.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_15->SetName("h_unfold_calib_trim_10_half_15");
  toy_errors(unfold_calib_trim_10_half_15, h_unfold_calib_trim_10_half_15, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_16(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 16, false, true);
  TH2D* h_unfold_calib_trim_10_half_16 = (TH2D*)unfold_calib_trim_10_half_16.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_16->SetName("h_unfold_calib_trim_10_half_16");
  toy_errors(unfold_calib_trim_10_half_16, h_unfold_calib_trim_10_half_16, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_17(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 17, false, true);
  TH2D* h_unfold_calib_trim_10_half_17 = (TH2D*)unfold_calib_trim_10_half_17.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_17->SetName("h_unfold_calib_trim_10_half_17");
  toy_errors(unfold_calib_trim_10_half_17, h_unfold_calib_trim_10_half_17, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_18(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 18, false, true);
  TH2D* h_unfold_calib_trim_10_half_18 = (TH2D*)unfold_calib_trim_10_half_18.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_18->SetName("h_unfold_calib_trim_10_half_18");
  toy_errors(unfold_calib_trim_10_half_18, h_unfold_calib_trim_10_half_18, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_19(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 19, false, true);
  TH2D* h_unfold_calib_trim_10_half_19 = (TH2D*)unfold_calib_trim_10_half_19.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_19->SetName("h_unfold_calib_trim_10_half_19");
  toy_errors(unfold_calib_trim_10_half_19, h_unfold_calib_trim_10_half_19, ntoys);
  RooUnfoldBayes unfold_calib_trim_10_half_20(h_respmatrix_calib_half1_trim_10, h_measure_calib_dijet_half2_record, 20, false, true);
  TH2D* h_unfold_calib_trim_10_half_20 = (TH2D*)unfold_calib_trim_10_half_20.Hunfold(RooUnfolding::kErrors); 
  h_unfold_calib_trim_10_half_20->SetName("h_unfold_calib_trim_10_half_20");
  toy_errors(unfold_calib_trim_10_half_20, h_unfold_calib_trim_10_half_20, ntoys);

  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();
  h_truth_calib_all_record->Write(); h_measure_calib_all_record->Write(); h_truth_calib_all->Write(); h_meas_calib_all->Write();
  h_unfold_calib_full_1->Write(); h_unfold_calib_full_2->Write(); h_unfold_calib_full_extra_1->Write(); h_unfold_calib_full_extra_2->Write();
  h_unfold_calib_trim_5_full_1->Write(); h_unfold_calib_trim_5_full_2->Write(); h_unfold_calib_trim_5_full_extra_1->Write(); h_unfold_calib_trim_5_full_extra_2->Write();
  h_unfold_calib_trim_10_full_1->Write(); h_unfold_calib_trim_10_full_2->Write(); h_unfold_calib_trim_10_full_extra_1->Write(); h_unfold_calib_trim_10_full_extra_2->Write();
  
  h_truth_calib_dijet_half1_record->Write(); h_measure_calib_dijet_half1_record->Write(); h_truth_calib_dijet_half2_record->Write(); h_measure_calib_dijet_half2_record->Write();
  h_unfold_calib_half_1->Write(); h_unfold_calib_half_2->Write(); h_unfold_calib_half_3->Write(); h_unfold_calib_half_4->Write(); h_unfold_calib_half_5->Write();
  h_unfold_calib_half_6->Write(); h_unfold_calib_half_7->Write(); h_unfold_calib_half_8->Write(); h_unfold_calib_half_9->Write(); h_unfold_calib_half_10->Write();
  h_unfold_calib_half_11->Write(); h_unfold_calib_half_12->Write(); h_unfold_calib_half_13->Write(); h_unfold_calib_half_14->Write(); h_unfold_calib_half_15->Write();
  h_unfold_calib_half_16->Write(); h_unfold_calib_half_17->Write(); h_unfold_calib_half_18->Write(); h_unfold_calib_half_19->Write(); h_unfold_calib_half_20->Write();
  h_unfold_calib_trim_5_half_1->Write(); h_unfold_calib_trim_5_half_2->Write(); h_unfold_calib_trim_5_half_3->Write(); h_unfold_calib_trim_5_half_4->Write(); h_unfold_calib_trim_5_half_5->Write();
  h_unfold_calib_trim_5_half_6->Write(); h_unfold_calib_trim_5_half_7->Write(); h_unfold_calib_trim_5_half_8->Write(); h_unfold_calib_trim_5_half_9->Write(); h_unfold_calib_trim_5_half_10->Write();
  h_unfold_calib_trim_5_half_11->Write(); h_unfold_calib_trim_5_half_12->Write(); h_unfold_calib_trim_5_half_13->Write(); h_unfold_calib_trim_5_half_14->Write(); h_unfold_calib_trim_5_half_15->Write();
  h_unfold_calib_trim_5_half_16->Write(); h_unfold_calib_trim_5_half_17->Write(); h_unfold_calib_trim_5_half_18->Write(); h_unfold_calib_trim_5_half_19->Write(); h_unfold_calib_trim_5_half_20->Write();
  h_unfold_calib_trim_10_half_1->Write(); h_unfold_calib_trim_10_half_2->Write(); h_unfold_calib_trim_10_half_3->Write(); h_unfold_calib_trim_10_half_4->Write(); h_unfold_calib_trim_10_half_5->Write();
  h_unfold_calib_trim_10_half_6->Write(); h_unfold_calib_trim_10_half_7->Write(); h_unfold_calib_trim_10_half_8->Write(); h_unfold_calib_trim_10_half_9->Write(); h_unfold_calib_trim_10_half_10->Write();
  h_unfold_calib_trim_10_half_11->Write(); h_unfold_calib_trim_10_half_12->Write(); h_unfold_calib_trim_10_half_13->Write(); h_unfold_calib_trim_10_half_14->Write(); h_unfold_calib_trim_10_half_15->Write();
  h_unfold_calib_trim_10_half_16->Write(); h_unfold_calib_trim_10_half_17->Write(); h_unfold_calib_trim_10_half_18->Write(); h_unfold_calib_trim_10_half_19->Write(); h_unfold_calib_trim_10_half_20->Write();
  
  f_out->Close();
  std::cout << "All done!" << std::endl;
}