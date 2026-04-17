#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>

void get_reweightmatrix(const char* simfile = "analysis_sim_run21_output/output_dijet_sim_iter_2.root", 
                        const char* datafile = "analysis_data_run21_output/output_pu_correct_data_dijet.root", 
                        const char* reweightfile = "output_reweighted_respmatrix_run21_iter_2.root") {
  // Read Files
  //********** Files **********//
  TFile *f_sim = new TFile(simfile, "READ");
  TFile *f_data = new TFile(datafile, "READ");

  std::vector<std::string> syst = {"calib_dijet","calib_dijet_jesdown","calib_dijet_jesup","calib_dijet_jerdown","calib_dijet_jerup",
                                  "calib_dijet_clus_smear","calib_dijet_ohcal_mc_data_var","calib_dijet_2sigma_noise","calib_dijet_4sigma_noise",
                                  "calib_dijet_half1","calib_dijet_half2"};
  std::vector<std::string> trim = {"","_trim_5","_trim_10"};

  RooUnfoldResponse* h_respmatrix[11][3];
  TH2D* h_truth_sim[11];
  for (int i = 0; i < 11; i++) {
    h_truth_sim[i] = (TH2D*)f_sim->Get(("h_truth_"+syst[i]).c_str());
    for (int j = 0; j < 3; j++) {
      h_respmatrix[i][j] = (RooUnfoldResponse*)f_sim->Get(("h_respmatrix_"+syst[i]+trim[j]).c_str());
    }
  }
  TH2D* h_measure_data = (TH2D*)f_data->Get("h_calibjet_pt_dijet_eff");
  //TH2D* h_measure_data = (TH2D*)f_data->Get("h_calibjet_pt_dijet_pu_correct_et");
  TH2D* h_measure_data_2sigma = (TH2D*)f_data->Get("h_calibjet_pt_2sigma_noise");
  TH2D* h_measure_data_4sigma = (TH2D*)f_data->Get("h_calibjet_pt_4sigma_noise");

  TH2D* h_norm_truth_sim[11];
  RooUnfoldBayes pre_unfold[11][3];
  TH2D* pseudo_truth[11][3];
  TH2D* weights[11][3];
  TH2D* h_rw_truth[11][3];
  for (int i = 0; i < 11; i++) {
    h_norm_truth_sim[i] = dynamic_cast<TH2D*>(h_truth_sim[i]->Clone(("h_norm_truth_"+syst[i]).c_str()));
    h_norm_truth_sim[i]->Scale(1.0/h_norm_truth_sim[i]->Integral());
    for (int j = 0; j < 3; j++) {
      if (i == 7) {
        pre_unfold[i][j] = RooUnfoldBayes(h_respmatrix[i][j], h_measure_data_2sigma, 1, false, true);
      } else if (i == 8) {
        pre_unfold[i][j] = RooUnfoldBayes(h_respmatrix[i][j], h_measure_data_4sigma, 1, false, true);
      } else {
        pre_unfold[i][j] = RooUnfoldBayes(h_respmatrix[i][j], h_measure_data, 1, false, true);
      }
      pseudo_truth[i][j] = dynamic_cast<TH2D*>(pre_unfold[i][j].Hunfold(RooUnfolding::kErrors));
      pseudo_truth[i][j]->SetName(("pseudo_truth_"+syst[i]+trim[j]).c_str());

      weights[i][j] = dynamic_cast<TH2D*>(pseudo_truth[i][j]->Clone(("weights_"+syst[i]+trim[j]).c_str()));
      weights[i][j]->Scale(1.0/weights[i][j]->Integral());
      weights[i][j]->Divide(h_norm_truth_sim[i]);

      h_rw_truth[i][j] = dynamic_cast<TH2D*>(h_truth_sim[i]->Clone(("h_rw_truth_"+syst[i]+trim[j]).c_str()));
      h_rw_truth[i][j]->Multiply(weights[i][j]);
    }
  }

  // 1D unfoldings 
  std::vector<std::string> syst_1D = {"jetpt","caloet"};
  RooUnfoldResponse* h_1D_respmatrix[2][3];
  TH1D* h_1D_truth_sim[2]; TH1D* h_1D_measure_data[2];
  for (int j = 0; j < 3; j++) {
    h_1D_respmatrix[0][j] = (RooUnfoldResponse*)f_sim->Get(("h_jetpt_respmatrix"+trim[j]).c_str());
    h_1D_respmatrix[1][j] = (RooUnfoldResponse*)f_sim->Get(("h_caloet_respmatrix"+trim[j]).c_str());
  } 
  h_1D_truth_sim[0] = (TH1D*)h_1D_respmatrix[0][0]->Htruth(); h_1D_truth_sim[0]->SetName("h_1D_jetpt_truth");
  h_1D_truth_sim[1] = (TH1D*)h_1D_respmatrix[1][0]->Htruth(); h_1D_truth_sim[1]->SetName("h_1D_caloet_truth");
  h_1D_measure_data[0] = (TH1D*)f_data->Get("h_jetpt");
  h_1D_measure_data[1] = (TH1D*)f_data->Get("h_caloet");

  TH2D* test_respmatrix = dynamic_cast<TH2D*>(h_1D_respmatrix[0][0]->Hresponse()); test_respmatrix->SetName("jetpt_test_resp");
  RooUnfoldBayes pre_unfold_test = RooUnfoldBayes(h_1D_respmatrix[0][0], h_1D_measure_data[0], 1, false, true);
  TH1D* pseudo_truth_test = dynamic_cast<TH1D*>(pre_unfold_test.Hunfold(RooUnfolding::kErrors));
  pseudo_truth_test->SetName("jetpt_test");

  TH1D* h_1D_norm_truth_sim[2];
  RooUnfoldBayes pre_unfold_1D[2][3];
  TH1D* pseudo_truth_1D[2][3];
  TH1D* weights_1D[2][3];
  TH1D* h_1D_rw_truth[2][3];
  for (int i = 0; i < 2; i++) {
    h_1D_norm_truth_sim[i] = dynamic_cast<TH1D*>(h_1D_truth_sim[i]->Clone(("h_1D_norm_truth_"+syst_1D[i]).c_str()));
    h_1D_norm_truth_sim[i]->Scale(1.0/h_1D_norm_truth_sim[i]->Integral());
    for (int j = 0; j < 3; j++) {
      pre_unfold_1D[i][j] = RooUnfoldBayes(h_1D_respmatrix[i][j], h_1D_measure_data[i], 1, false, true);
      pseudo_truth_1D[i][j] = dynamic_cast<TH1D*>(pre_unfold_1D[i][j].Hunfold(RooUnfolding::kErrors));
      pseudo_truth_1D[i][j]->SetName(("pseudo_truth_"+syst_1D[i]+trim[j]).c_str());

      weights_1D[i][j] = dynamic_cast<TH1D*>(pseudo_truth_1D[i][j]->Clone(("weights_"+syst_1D[i]+trim[j]).c_str()));
      weights_1D[i][j]->Scale(1.0/weights_1D[i][j]->Integral());
      weights_1D[i][j]->Divide(h_1D_norm_truth_sim[i]);

      h_1D_rw_truth[i][j] = dynamic_cast<TH1D*>(h_1D_truth_sim[i]->Clone(("h_1D_rw_truth_"+syst[i]+trim[j]).c_str()));
      h_1D_rw_truth[i][j]->Multiply(weights_1D[i][j]);
    }
  }
  
  /*
  // unfolding for full response matrix, trim 5 and trim 10 entries response matrices
  RooUnfoldResponse* h_respmatrix_calib = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet");
  RooUnfoldResponse* h_respmatrix_calib_trim_5 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_trim_5");
  RooUnfoldResponse* h_respmatrix_calib_trim_10 = (RooUnfoldResponse*)f_sim->Get("h_respmatrix_calib_dijet_trim_10");
  TH2D* h_truth_sim = (TH2D*)f_sim->Get("h_truth_calib_dijet");
  TH2D* h_measure_sim = (TH2D*)f_sim->Get("h_measure_calib_dijet");
  TH2D* h_measure_data = (TH2D*)f_data->Get("h_calibjet_pt_dijet_eff");

  RooUnfoldBayes pre_unfold(h_respmatrix_calib, h_measure_data, 1, false, true);
  TH2D* pseudo_truth = dynamic_cast<TH2D*>(pre_unfold.Hunfold(RooUnfolding::kErrors));
  pseudo_truth->SetName("pseudo_truth");
  RooUnfoldBayes pre_unfold_5(h_respmatrix_calib_trim_5, h_measure_data, 1, false, true);
  TH2D* pseudo_truth_trim_5 = dynamic_cast<TH2D*>(pre_unfold_5.Hunfold(RooUnfolding::kErrors));
  pseudo_truth_trim_5->SetName("pseudo_truth_trim_5");
  RooUnfoldBayes pre_unfold_10(h_respmatrix_calib_trim_10, h_measure_data, 1, false, true);
  TH2D* pseudo_truth_trim_10 = dynamic_cast<TH2D*>(pre_unfold_10.Hunfold(RooUnfolding::kErrors));
  pseudo_truth_trim_10->SetName("pseudo_truth_trim_10");

  TH2D* weights = dynamic_cast<TH2D*>(pseudo_truth->Clone("weights"));
  weights->Scale(1.0/weights->Integral());
  weights->Divide(h_norm_truth_sim);
  
  TH2D* weights_trim_5 = dynamic_cast<TH2D*>(pseudo_truth_trim_5->Clone("weights_trim_5"));
  weights_trim_5->Scale(1.0/weights_trim_5->Integral());
  weights_trim_5->Divide(h_norm_truth_sim);
  
  TH2D* weights_trim_10 = dynamic_cast<TH2D*>(pseudo_truth_trim_10->Clone("weights_trim_10"));
  weights_trim_10->Scale(1.0/weights_trim_10->Integral());
  weights_trim_10->Divide(h_norm_truth_sim);

  TH2D* h_rw_truth = dynamic_cast<TH2D*>(h_truth_sim->Clone("h_rw_truth"));
  TH2D* h_rw_truth_trim_5 = dynamic_cast<TH2D*>(h_truth_sim->Clone("h_rw_truth_trim_5"));
  TH2D* h_rw_truth_trim_10 = dynamic_cast<TH2D*>(h_truth_sim->Clone("h_rw_truth_trim_10"));

  h_rw_truth->Multiply(weights);
  h_rw_truth_trim_5->Multiply(weights_trim_5);
  h_rw_truth_trim_10->Multiply(weights_trim_10);

  TFile *f_out = new TFile("output_reweighted_respmatrix_4bin_newetbin_dijet.root", "RECREATE");
  f_out->cd();
  h_truth_sim->Write();
  pseudo_truth->Write();
  pseudo_truth_trim_5->Write();
  pseudo_truth_trim_10->Write();
  h_norm_truth_sim->Write();
  weights->Write();
  weights_trim_5->Write();
  weights_trim_10->Write();
  h_rw_truth->Write();
  h_rw_truth_trim_5->Write();
  h_rw_truth_trim_10->Write();
  f_out->Close();
  */


  TFile *f_out = new TFile(reweightfile, "RECREATE");
  f_out->cd();
  h_measure_data->Write();
  h_measure_data_2sigma->Write();
  h_measure_data_4sigma->Write();
  for (int i = 0; i < 11; i++) {
    h_truth_sim[i]->Write();
    h_norm_truth_sim[i]->Write();
    for (int j = 0; j < 3; j++) {
      pseudo_truth[i][j]->Write();
      weights[i][j]->Write();
      h_rw_truth[i][j]->Write();
    }
  }
  for (int i = 0; i < 2; i++) {
    h_1D_measure_data[i]->Write();
    h_1D_truth_sim[i]->Write();
    h_1D_norm_truth_sim[i]->Write();
    for (int j = 0; j < 3; j++) {
      pseudo_truth_1D[i][j]->Write();
      weights_1D[i][j]->Write();
      h_1D_rw_truth[i][j]->Write();
    }
  }

  pseudo_truth_test->Write();
  test_respmatrix->Write();
  f_out->Close();

}