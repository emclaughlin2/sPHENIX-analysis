#include <iostream>
#include <TFile.h>
#include <TH1D.h>
#include "/sphenix/user/hanpuj/CaloDataAna24_skimmed/src/draw_template.C"
#include "unfold_Def.h"

void get_purityefficiency(string infilename = "analysis_sim_run21_output/output_efrac_bkg_cut_sim_iter_3.root", string outfilename = "output_purityefficiency_efrac_bkg_cut.root") {
  // Read Files
  TFile *f_in = new TFile(infilename.c_str(), "READ");
  if (!f_in) {
    std::cout << "Error: cannot open output_sim.root" << std::endl;
    return;
  }

  // Read histograms
  TH2D* h_truth_calib_dijet = (TH2D*)f_in->Get("h_truth_calib_dijet");
  TH2D* h_measure_calib_dijet = (TH2D*)f_in->Get("h_measure_calib_dijet");
  RooUnfoldResponse* h_respmatrix_calib_dijet = (RooUnfoldResponse*)f_in->Get("h_respmatrix_calib_dijet");
  TH2D* h_truth_calib_dijet_jesdown = (TH2D*)f_in->Get("h_truth_calib_dijet_jesdown");
  TH2D* h_measure_calib_dijet_jesdown = (TH2D*)f_in->Get("h_measure_calib_dijet_jesdown");
  RooUnfoldResponse* h_respmatrix_calib_dijet_jesdown = (RooUnfoldResponse*)f_in->Get("h_respmatrix_calib_dijet_jesdown");
  TH2D* h_truth_calib_dijet_jesup = (TH2D*)f_in->Get("h_truth_calib_dijet_jesup");
  TH2D* h_measure_calib_dijet_jesup = (TH2D*)f_in->Get("h_measure_calib_dijet_jesup");
  RooUnfoldResponse* h_respmatrix_calib_dijet_jesup = (RooUnfoldResponse*)f_in->Get("h_respmatrix_calib_dijet_jesup");
  TH2D* h_truth_calib_dijet_jerdown = (TH2D*)f_in->Get("h_truth_calib_dijet_jerdown");
  TH2D* h_measure_calib_dijet_jerdown = (TH2D*)f_in->Get("h_measure_calib_dijet_jerdown");
  RooUnfoldResponse* h_respmatrix_calib_dijet_jerdown = (RooUnfoldResponse*)f_in->Get("h_respmatrix_calib_dijet_jerdown");
  TH2D* h_truth_calib_dijet_jerup = (TH2D*)f_in->Get("h_truth_calib_dijet_jerup");
  TH2D* h_measure_calib_dijet_jerup = (TH2D*)f_in->Get("h_measure_calib_dijet_jerup");
  RooUnfoldResponse* h_respmatrix_calib_dijet_jerup = (RooUnfoldResponse*)f_in->Get("h_respmatrix_calib_dijet_jerup");
  TH2D* h_miss_calib_dijet = (TH2D*)f_in->Get("h_miss_calib_dijet");
  TH2D* h_fake_calib_dijet = (TH2D*)f_in->Get("h_fake_calib_dijet");
  TH2D* h_miss_calib_dijet_jesdown = (TH2D*)f_in->Get("h_miss_calib_dijet_jesdown");
  TH2D* h_fake_calib_dijet_jesdown = (TH2D*)f_in->Get("h_fake_calib_dijet_jesdown");
  TH2D* h_miss_calib_dijet_jesup = (TH2D*)f_in->Get("h_miss_calib_dijet_jesup");
  TH2D* h_fake_calib_dijet_jesup = (TH2D*)f_in->Get("h_fake_calib_dijet_jesup");
  TH2D* h_miss_calib_dijet_jerdown = (TH2D*)f_in->Get("h_miss_calib_dijet_jerdown");
  TH2D* h_fake_calib_dijet_jerdown = (TH2D*)f_in->Get("h_fake_calib_dijet_jerdown");
  TH2D* h_miss_calib_dijet_jerup = (TH2D*)f_in->Get("h_miss_calib_dijet_jerup");
  TH2D* h_fake_calib_dijet_jerup = (TH2D*)f_in->Get("h_fake_calib_dijet_jerup");
  

  // Get ET and PT projections 
  TH1D* h_measure_jet_pt_calib_dijet = (TH1D*)h_measure_calib_dijet->ProjectionX("h_measure_jet_pt_calib_dijet");
  TH1D* h_truth_jet_pt_calib_dijet = (TH1D*)h_truth_calib_dijet->ProjectionX("h_truth_jet_pt_calib_dijet");
  TH1D* h_measure_jet_pt_calib_dijet_jesdown = (TH1D*)h_measure_calib_dijet_jesdown->ProjectionX("h_measure_jet_pt_calib_dijet_jesdown");
  TH1D* h_truth_jet_pt_calib_dijet_jesdown = (TH1D*)h_truth_calib_dijet_jesdown->ProjectionX("h_truth_jet_pt_calib_dijet_jesdown");
  TH1D* h_measure_jet_pt_calib_dijet_jesup = (TH1D*)h_measure_calib_dijet_jesup->ProjectionX("h_measure_jet_pt_calib_dijet_jesup");
  TH1D* h_truth_jet_pt_calib_dijet_jesup = (TH1D*)h_truth_calib_dijet_jesup->ProjectionX("h_truth_jet_pt_calib_dijet_jesup");
  TH1D* h_measure_jet_pt_calib_dijet_jerdown = (TH1D*)h_measure_calib_dijet_jerdown->ProjectionX("h_measure_jet_pt_calib_dijet_jerdown");
  TH1D* h_truth_jet_pt_calib_dijet_jerdown = (TH1D*)h_truth_calib_dijet_jerdown->ProjectionX("h_truth_jet_pt_calib_dijet_jerdown");
  TH1D* h_measure_jet_pt_calib_dijet_jerup = (TH1D*)h_measure_calib_dijet_jerup->ProjectionX("h_measure_jet_pt_calib_dijet_jerup");
  TH1D* h_truth_jet_pt_calib_dijet_jerup = (TH1D*)h_truth_calib_dijet_jerup->ProjectionX("h_truth_jet_pt_calib_dijet_jerup");

  TH1D* h_fake_jet_pt_calib_dijet = (TH1D*)h_fake_calib_dijet->ProjectionX("h_fake_jet_pt_calib_dijet");
  TH1D* h_miss_jet_pt_calib_dijet = (TH1D*)h_miss_calib_dijet->ProjectionX("h_miss_jet_pt_calib_dijet");
  TH1D* h_fake_jet_pt_calib_dijet_jesdown = (TH1D*)h_fake_calib_dijet_jesdown->ProjectionX("h_fake_jet_pt_calib_dijet_jesdown");
  TH1D* h_miss_jet_pt_calib_dijet_jesdown = (TH1D*)h_miss_calib_dijet_jesdown->ProjectionX("h_miss_jet_pt_calib_dijet_jesdown");
  TH1D* h_fake_jet_pt_calib_dijet_jesup = (TH1D*)h_fake_calib_dijet_jesup->ProjectionX("h_fake_jet_pt_calib_dijet_jesup");
  TH1D* h_miss_jet_pt_calib_dijet_jesup = (TH1D*)h_miss_calib_dijet_jesup->ProjectionX("h_miss_jet_pt_calib_dijet_jesup");
  TH1D* h_fake_jet_pt_calib_dijet_jerdown = (TH1D*)h_fake_calib_dijet_jerdown->ProjectionX("h_fake_jet_pt_calib_dijet_jerdown");
  TH1D* h_miss_jet_pt_calib_dijet_jerdown = (TH1D*)h_miss_calib_dijet_jerdown->ProjectionX("h_miss_jet_pt_calib_dijet_jerdown");
  TH1D* h_fake_jet_pt_calib_dijet_jerup = (TH1D*)h_fake_calib_dijet_jerup->ProjectionX("h_fake_jet_pt_calib_dijet_jerup");
  TH1D* h_miss_jet_pt_calib_dijet_jerup = (TH1D*)h_miss_calib_dijet_jerup->ProjectionX("h_miss_jet_pt_calib_dijet_jerup");
  
  TH1D* h_measure_trans_et_calib_dijet = (TH1D*)h_measure_calib_dijet->ProjectionY("h_measure_trans_et_calib_dijet");
  TH1D* h_truth_trans_et_calib_dijet = (TH1D*)h_truth_calib_dijet->ProjectionY("h_truth_trans_et_calib_dijet");
  TH1D* h_measure_trans_et_calib_dijet_jesdown = (TH1D*)h_measure_calib_dijet_jesdown->ProjectionY("h_measure_trans_et_calib_dijet_jesdown");
  TH1D* h_truth_trans_et_calib_dijet_jesdown = (TH1D*)h_truth_calib_dijet_jesdown->ProjectionY("h_truth_trans_et_calib_dijet_jesdown");
  TH1D* h_measure_trans_et_calib_dijet_jesup = (TH1D*)h_measure_calib_dijet_jesup->ProjectionY("h_measure_trans_et_calib_dijet_jesup");
  TH1D* h_truth_trans_et_calib_dijet_jesup = (TH1D*)h_truth_calib_dijet_jesup->ProjectionY("h_truth_trans_et_calib_dijet_jesup");
  TH1D* h_measure_trans_et_calib_dijet_jerdown = (TH1D*)h_measure_calib_dijet_jerdown->ProjectionY("h_measure_trans_et_calib_dijet_jerdown");
  TH1D* h_truth_trans_et_calib_dijet_jerdown = (TH1D*)h_truth_calib_dijet_jerdown->ProjectionY("h_truth_trans_et_calib_dijet_jerdown");
  TH1D* h_measure_trans_et_calib_dijet_jerup = (TH1D*)h_measure_calib_dijet_jerup->ProjectionY("h_measure_trans_et_calib_dijet_jerup");
  TH1D* h_truth_trans_et_calib_dijet_jerup = (TH1D*)h_truth_calib_dijet_jerup->ProjectionY("h_truth_trans_et_calib_dijet_jerup");

  TH1D* h_fake_trans_et_calib_dijet = (TH1D*)h_fake_calib_dijet->ProjectionY("h_fake_trans_et_calib_dijet");
  TH1D* h_miss_trans_et_calib_dijet = (TH1D*)h_miss_calib_dijet->ProjectionY("h_miss_trans_et_calib_dijet");
  TH1D* h_fake_trans_et_calib_dijet_jesdown = (TH1D*)h_fake_calib_dijet_jesdown->ProjectionY("h_fake_trans_et_calib_dijet_jesdown");
  TH1D* h_miss_trans_et_calib_dijet_jesdown = (TH1D*)h_miss_calib_dijet_jesdown->ProjectionY("h_miss_trans_et_calib_dijet_jesdown");
  TH1D* h_fake_trans_et_calib_dijet_jesup = (TH1D*)h_fake_calib_dijet_jesup->ProjectionY("h_fake_trans_et_calib_dijet_jesup");
  TH1D* h_miss_trans_et_calib_dijet_jesup = (TH1D*)h_miss_calib_dijet_jesup->ProjectionY("h_miss_trans_et_calib_dijet_jesup");
  TH1D* h_fake_trans_et_calib_dijet_jerdown = (TH1D*)h_fake_calib_dijet_jerdown->ProjectionY("h_fake_trans_et_calib_dijet_jerdown");
  TH1D* h_miss_trans_et_calib_dijet_jerdown = (TH1D*)h_miss_calib_dijet_jerdown->ProjectionY("h_miss_trans_et_calib_dijet_jerdown");
  TH1D* h_fake_trans_et_calib_dijet_jerup = (TH1D*)h_fake_calib_dijet_jerup->ProjectionY("h_fake_trans_et_calib_dijet_jerup");
  TH1D* h_miss_trans_et_calib_dijet_jerup = (TH1D*)h_miss_calib_dijet_jerup->ProjectionY("h_miss_trans_et_calib_dijet_jerup");

  // purity is (measured - fake)/measured, efficiency is (truth - miss)/truth
  
    // Form purity and efficiency histograms
  TH2D* h_purity_calib_dijet = (TH2D*)h_respmatrix_calib_dijet->Hmeasured();
  TH2D* h_efficiency_calib_dijet = (TH2D*)h_respmatrix_calib_dijet->Htruth();
  TH2D* h_purity_calib_dijet_jesdown = (TH2D*)h_respmatrix_calib_dijet_jesdown->Hmeasured();
  TH2D* h_efficiency_calib_dijet_jesdown = (TH2D*)h_respmatrix_calib_dijet_jesdown->Htruth();
  TH2D* h_purity_calib_dijet_jesup = (TH2D*)h_respmatrix_calib_dijet_jesup->Hmeasured();
  TH2D* h_efficiency_calib_dijet_jesup = (TH2D*)h_respmatrix_calib_dijet_jesup->Htruth();
  TH2D* h_purity_calib_dijet_jerdown = (TH2D*)h_respmatrix_calib_dijet_jerdown->Hmeasured();
  TH2D* h_efficiency_calib_dijet_jerdown = (TH2D*)h_respmatrix_calib_dijet_jerdown->Htruth();
  TH2D* h_purity_calib_dijet_jerup = (TH2D*)h_respmatrix_calib_dijet_jerup->Hmeasured();
  TH2D* h_efficiency_calib_dijet_jerup = (TH2D*)h_respmatrix_calib_dijet_jerup->Htruth();
  
  h_purity_calib_dijet->SetName("h_purity_calib_dijet");
  h_efficiency_calib_dijet->SetName("h_efficiency_calib_dijet");
  h_purity_calib_dijet_jesdown->SetName("h_purity_calib_dijet_jesdown");
  h_efficiency_calib_dijet_jesdown->SetName("h_efficiency_calib_dijet_jesdown");
  h_purity_calib_dijet_jesup->SetName("h_purity_calib_dijet_jesup");
  h_efficiency_calib_dijet_jesup->SetName("h_efficiency_calib_dijet_jesup");
  h_purity_calib_dijet_jerdown->SetName("h_purity_calib_dijet_jerdown");
  h_efficiency_calib_dijet_jerdown->SetName("h_efficiency_calib_dijet_jerdown");
  h_purity_calib_dijet_jerup->SetName("h_purity_calib_dijet_jerup");
  h_efficiency_calib_dijet_jerup->SetName("h_efficiency_calib_dijet_jerup");

  for (int ix = 1; ix < h_purity_calib_dijet->GetNbinsX() + 1; ix++) {
    for (int iy = 1; iy < h_purity_calib_dijet->GetNbinsY() + 1; iy++) {
      h_purity_calib_dijet->SetBinContent(ix, iy, h_purity_calib_dijet->GetBinContent(ix, iy) - h_fake_calib_dijet->GetBinContent(ix, iy));
      h_purity_calib_dijet->SetBinError(ix , iy, sqrt(h_purity_calib_dijet->GetBinError(ix, iy)*h_purity_calib_dijet->GetBinError(ix, iy) + h_fake_calib_dijet->GetBinError(ix, iy)*h_fake_calib_dijet->GetBinError(ix, iy)));
      h_purity_calib_dijet_jesdown->SetBinContent(ix, iy, h_purity_calib_dijet_jesdown->GetBinContent(ix, iy) - h_fake_calib_dijet_jesdown->GetBinContent(ix, iy));
      h_purity_calib_dijet_jesdown->SetBinError(ix , iy, sqrt(h_purity_calib_dijet_jesdown->GetBinError(ix, iy)*h_purity_calib_dijet_jesdown->GetBinError(ix, iy) + h_fake_calib_dijet_jesdown->GetBinError(ix, iy)*h_fake_calib_dijet_jesdown->GetBinError(ix, iy)));
      h_purity_calib_dijet_jesup->SetBinContent(ix, iy, h_purity_calib_dijet_jesup->GetBinContent(ix, iy) - h_fake_calib_dijet_jesup->GetBinContent(ix, iy));
      h_purity_calib_dijet_jesup->SetBinError(ix , iy, sqrt(h_purity_calib_dijet_jesup->GetBinError(ix, iy)*h_purity_calib_dijet_jesup->GetBinError(ix, iy) + h_fake_calib_dijet_jesup->GetBinError(ix, iy)*h_fake_calib_dijet_jesup->GetBinError(ix, iy)));
      h_purity_calib_dijet_jerdown->SetBinContent(ix, iy, h_purity_calib_dijet_jerdown->GetBinContent(ix, iy) - h_fake_calib_dijet_jerdown->GetBinContent(ix, iy));
      h_purity_calib_dijet_jerdown->SetBinError(ix , iy, sqrt(h_purity_calib_dijet_jerdown->GetBinError(ix, iy)*h_purity_calib_dijet_jerdown->GetBinError(ix, iy) + h_fake_calib_dijet_jerdown->GetBinError(ix, iy)*h_fake_calib_dijet_jerdown->GetBinError(ix, iy)));
      h_purity_calib_dijet_jerup->SetBinContent(ix, iy, h_purity_calib_dijet_jerup->GetBinContent(ix, iy) - h_fake_calib_dijet_jerup->GetBinContent(ix, iy));
      h_purity_calib_dijet_jerup->SetBinError(ix , iy, sqrt(h_purity_calib_dijet_jerup->GetBinError(ix, iy)*h_purity_calib_dijet_jerup->GetBinError(ix, iy) + h_fake_calib_dijet_jerup->GetBinError(ix, iy)*h_fake_calib_dijet_jerup->GetBinError(ix, iy)));
    }
  }

  for (int ix = 1; ix < h_efficiency_calib_dijet->GetNbinsX() + 1; ix++) {
    for (int iy = 1; iy < h_efficiency_calib_dijet->GetNbinsY() + 1; iy++) {
      h_efficiency_calib_dijet->SetBinContent(ix, iy, h_efficiency_calib_dijet->GetBinContent(ix, iy) - h_miss_calib_dijet->GetBinContent(ix, iy));
      h_efficiency_calib_dijet->SetBinError(ix, iy, sqrt(h_efficiency_calib_dijet->GetBinError(ix,iy)*h_efficiency_calib_dijet->GetBinError(ix,iy) + h_miss_calib_dijet->GetBinError(ix, iy)*h_miss_calib_dijet->GetBinError(ix,iy)));
      h_efficiency_calib_dijet_jesdown->SetBinContent(ix, iy, h_efficiency_calib_dijet_jesdown->GetBinContent(ix, iy) - h_miss_calib_dijet_jesdown->GetBinContent(ix, iy));
      h_efficiency_calib_dijet_jesdown->SetBinError(ix, iy, sqrt(h_efficiency_calib_dijet_jesdown->GetBinError(ix,iy)*h_efficiency_calib_dijet_jesdown->GetBinError(ix,iy) + h_miss_calib_dijet_jesdown->GetBinError(ix, iy)*h_miss_calib_dijet_jesdown->GetBinError(ix,iy)));
      h_efficiency_calib_dijet_jesup->SetBinContent(ix, iy, h_efficiency_calib_dijet_jesup->GetBinContent(ix, iy) - h_miss_calib_dijet_jesup->GetBinContent(ix, iy));
      h_efficiency_calib_dijet_jesup->SetBinError(ix, iy, sqrt(h_efficiency_calib_dijet_jesup->GetBinError(ix,iy)*h_efficiency_calib_dijet_jesup->GetBinError(ix,iy) + h_miss_calib_dijet_jesup->GetBinError(ix, iy)*h_miss_calib_dijet_jesup->GetBinError(ix,iy)));
      h_efficiency_calib_dijet_jerdown->SetBinContent(ix, iy, h_efficiency_calib_dijet_jerdown->GetBinContent(ix, iy) - h_miss_calib_dijet_jerdown->GetBinContent(ix, iy));
      h_efficiency_calib_dijet_jerdown->SetBinError(ix, iy, sqrt(h_efficiency_calib_dijet_jerdown->GetBinError(ix,iy)*h_efficiency_calib_dijet_jerdown->GetBinError(ix,iy) + h_miss_calib_dijet_jerdown->GetBinError(ix, iy)*h_miss_calib_dijet_jerdown->GetBinError(ix,iy)));
      h_efficiency_calib_dijet_jerup->SetBinContent(ix, iy, h_efficiency_calib_dijet_jerup->GetBinContent(ix, iy) - h_miss_calib_dijet_jerup->GetBinContent(ix, iy));
      h_efficiency_calib_dijet_jerup->SetBinError(ix, iy, sqrt(h_efficiency_calib_dijet_jerup->GetBinError(ix,iy)*h_efficiency_calib_dijet_jerup->GetBinError(ix,iy) + h_miss_calib_dijet_jerup->GetBinError(ix, iy)*h_miss_calib_dijet_jerup->GetBinError(ix,iy)));
    }
  }

  // Get ET and PT projections 
  TH1D* h_purity_jet_pt_calib_dijet = (TH1D*)h_purity_calib_dijet->ProjectionX("h_purity_jet_pt_calib_dijet");
  TH1D* h_efficiency_jet_pt_calib_dijet = (TH1D*)h_efficiency_calib_dijet->ProjectionX("h_efficiency_jet_pt_calib_dijet");
  TH1D* h_purity_jet_pt_calib_dijet_jesdown = (TH1D*)h_purity_calib_dijet_jesdown->ProjectionX("h_purity_jet_pt_calib_dijet_jesdown");
  TH1D* h_efficiency_jet_pt_calib_dijet_jesdown = (TH1D*)h_efficiency_calib_dijet_jesdown->ProjectionX("h_efficiency_jet_pt_calib_dijet_jesdown");
  TH1D* h_purity_jet_pt_calib_dijet_jesup = (TH1D*)h_purity_calib_dijet_jesup->ProjectionX("h_purity_jet_pt_calib_dijet_jesup");
  TH1D* h_efficiency_jet_pt_calib_dijet_jesup = (TH1D*)h_efficiency_calib_dijet_jesup->ProjectionX("h_efficiency_jet_pt_calib_dijet_jesup");
  TH1D* h_purity_jet_pt_calib_dijet_jerdown = (TH1D*)h_purity_calib_dijet_jerdown->ProjectionX("h_purity_jet_pt_calib_dijet_jerdown");
  TH1D* h_efficiency_jet_pt_calib_dijet_jerdown = (TH1D*)h_efficiency_calib_dijet_jerdown->ProjectionX("h_efficiency_jet_pt_calib_dijet_jerdown");
  TH1D* h_purity_jet_pt_calib_dijet_jerup = (TH1D*)h_purity_calib_dijet_jerup->ProjectionX("h_purity_jet_pt_calib_dijet_jerup");
  TH1D* h_efficiency_jet_pt_calib_dijet_jerup = (TH1D*)h_efficiency_calib_dijet_jerup->ProjectionX("h_efficiency_jet_pt_calib_dijet_jerup");
  
  TH1D* h_purity_trans_et_calib_dijet = (TH1D*)h_purity_calib_dijet->ProjectionY("h_purity_trans_et_calib_dijet");
  TH1D* h_efficiency_trans_et_calib_dijet = (TH1D*)h_efficiency_calib_dijet->ProjectionY("h_efficiency_trans_et_calib_dijet");
  TH1D* h_purity_trans_et_calib_dijet_jesdown = (TH1D*)h_purity_calib_dijet_jesdown->ProjectionY("h_purity_trans_et_calib_dijet_jesdown");
  TH1D* h_efficiency_trans_et_calib_dijet_jesdown = (TH1D*)h_efficiency_calib_dijet_jesdown->ProjectionY("h_efficiency_trans_et_calib_dijet_jesdown");
  TH1D* h_purity_trans_et_calib_dijet_jesup = (TH1D*)h_purity_calib_dijet_jesup->ProjectionY("h_purity_trans_et_calib_dijet_jesup");
  TH1D* h_efficiency_trans_et_calib_dijet_jesup = (TH1D*)h_efficiency_calib_dijet_jesup->ProjectionY("h_efficiency_trans_et_calib_dijet_jesup");
  TH1D* h_purity_trans_et_calib_dijet_jerdown = (TH1D*)h_purity_calib_dijet_jerdown->ProjectionY("h_purity_trans_et_calib_dijet_jerdown");
  TH1D* h_efficiency_trans_et_calib_dijet_jerdown = (TH1D*)h_efficiency_calib_dijet_jerdown->ProjectionY("h_efficiency_trans_et_calib_dijet_jerdown");
  TH1D* h_purity_trans_et_calib_dijet_jerup = (TH1D*)h_purity_calib_dijet_jerup->ProjectionY("h_purity_trans_et_calib_dijet_jerup");
  TH1D* h_efficiency_trans_et_calib_dijet_jerup = (TH1D*)h_efficiency_calib_dijet_jerup->ProjectionY("h_efficiency_trans_et_calib_dijet_jerup");
  
  // Calculate purity and efficiency
  h_purity_calib_dijet->Divide(h_purity_calib_dijet, h_measure_calib_dijet, 1, 1, "B");
  h_efficiency_calib_dijet->Divide(h_efficiency_calib_dijet, h_truth_calib_dijet, 1, 1, "B");
  h_purity_calib_dijet_jesdown->Divide(h_purity_calib_dijet_jesdown, h_measure_calib_dijet_jesdown, 1, 1, "B");
  h_efficiency_calib_dijet_jesdown->Divide(h_efficiency_calib_dijet_jesdown, h_truth_calib_dijet_jesdown, 1, 1, "B");
  h_purity_calib_dijet_jesup->Divide(h_purity_calib_dijet_jesup, h_measure_calib_dijet_jesup, 1, 1, "B");
  h_efficiency_calib_dijet_jesup->Divide(h_efficiency_calib_dijet_jesup, h_truth_calib_dijet_jesup, 1, 1, "B");
  h_purity_calib_dijet_jerdown->Divide(h_purity_calib_dijet_jerdown, h_measure_calib_dijet_jerdown, 1, 1, "B");
  h_efficiency_calib_dijet_jerdown->Divide(h_efficiency_calib_dijet_jerdown, h_truth_calib_dijet_jerdown, 1, 1, "B");
  h_purity_calib_dijet_jerup->Divide(h_purity_calib_dijet_jerup, h_measure_calib_dijet_jerup, 1, 1, "B");
  h_efficiency_calib_dijet_jerup->Divide(h_efficiency_calib_dijet_jerup, h_truth_calib_dijet_jerup, 1, 1, "B");

  // Write histograms
  TFile *f_out = new TFile("output_purityefficiency.root", "RECREATE");
  f_out->cd();
  h_purity_calib_dijet->Write();
  h_efficiency_calib_dijet->Write();
  h_purity_calib_dijet_jesdown->Write();
  h_efficiency_calib_dijet_jesdown->Write();
  h_purity_calib_dijet_jesup->Write();
  h_efficiency_calib_dijet_jesup->Write();
  h_purity_calib_dijet_jerdown->Write();
  h_efficiency_calib_dijet_jerdown->Write();
  h_purity_calib_dijet_jerup->Write();
  h_efficiency_calib_dijet_jerup->Write();
  f_out->Close();
  /*
  // Apply purity correction.
  TFile* f_data = new TFile("output_data.root", "READ");
  TFile* f_dataout = new TFile("output_data_puritycorr.root", "RECREATE");

  // Read histograms
  TH2D* h_calibjet_pt_dijet_eff = (TH2D*)f_data->Get("h_calibjet_pt_dijet_eff");
  TH2D* h_calibjet_pt_dijet_effdown = (TH2D*)f_data->Get("h_calibjet_pt_dijet_effdown");
  TH2D* h_calibjet_pt_dijet_effup = (TH2D*)f_data->Get("h_calibjet_pt_dijet_effup");
  TH2D* h_calibjet_pt_frac_eff = (TH2D*)f_data->Get("h_calibjet_pt_frac_eff");
  TH2D* h_calibjet_pt_frac_effdown = (TH2D*)f_data->Get("h_calibjet_pt_frac_effdown");
  TH2D* h_calibjet_pt_frac_effup = (TH2D*)f_data->Get("h_calibjet_pt_frac_effup");

  // Histograms for purity-corrected jet pT.
  TH2D* h_calibjet_pt_puritycorr_dijet = (TH2D*)h_calibjet_pt_dijet_eff->Clone("h_calibjet_pt_puritycorr_dijet");
  TH2D* h_calibjet_pt_puritycorr_dijet_effdown = (TH2D*)h_calibjet_pt_dijet_effdown->Clone("h_calibjet_pt_puritycorr_dijet_effdown");
  TH2D* h_calibjet_pt_puritycorr_dijet_effup = (TH2D*)h_calibjet_pt_dijet_effup->Clone("h_calibjet_pt_puritycorr_dijet_effup");
  TH2D* h_calibjet_pt_puritycorr_dijet_jesdown = (TH2D*)h_calibjet_pt_dijet_eff->Clone("h_calibjet_pt_puritycorr_dijet_jesdown");
  TH2D* h_calibjet_pt_puritycorr_dijet_jesup = (TH2D*)h_calibjet_pt_dijet_eff->Clone("h_calibjet_pt_puritycorr_dijet_jesup");
  TH2D* h_calibjet_pt_puritycorr_dijet_jerdown = (TH2D*)h_calibjet_pt_dijet_eff->Clone("h_calibjet_pt_puritycorr_dijet_jerdown");
  TH2D* h_calibjet_pt_puritycorr_dijet_jerup = (TH2D*)h_calibjet_pt_dijet_eff->Clone("h_calibjet_pt_puritycorr_dijet_jerup");

  TH2D* h_calibjet_pt_puritycorr_frac = (TH2D*)h_calibjet_pt_frac_eff->Clone("h_calibjet_pt_puritycorr_frac");
  TH2D* h_calibjet_pt_puritycorr_frac_effdown = (TH2D*)h_calibjet_pt_frac_effdown->Clone("h_calibjet_pt_puritycorr_frac_effdown");
  TH2D* h_calibjet_pt_puritycorr_frac_effup = (TH2D*)h_calibjet_pt_frac_effup->Clone("h_calibjet_pt_puritycorr_frac_effup");
  TH2D* h_calibjet_pt_puritycorr_frac_jesdown = (TH2D*)h_calibjet_pt_frac_eff->Clone("h_calibjet_pt_puritycorr_frac_jesdown");
  TH2D* h_calibjet_pt_puritycorr_frac_jesup = (TH2D*)h_calibjet_pt_frac_eff->Clone("h_calibjet_pt_puritycorr_frac_jesup");
  TH2D* h_calibjet_pt_puritycorr_frac_jerdown = (TH2D*)h_calibjet_pt_frac_eff->Clone("h_calibjet_pt_puritycorr_frac_jerdown");
  TH2D* h_calibjet_pt_puritycorr_frac_jerup = (TH2D*)h_calibjet_pt_frac_eff->Clone("h_calibjet_pt_puritycorr_frac_jerup");

  h_calibjet_pt_puritycorr_dijet->Divide(h_purity_calib_dijet);
  h_calibjet_pt_puritycorr_dijet_effdown->Divide(h_purity_calib_dijet);
  h_calibjet_pt_puritycorr_dijet_effup->Divide(h_purity_calib_dijet);
  h_calibjet_pt_puritycorr_dijet_jesdown->Divide(h_purity_calib_dijet_jesdown);
  h_calibjet_pt_puritycorr_dijet_jesup->Divide(h_purity_calib_dijet_jesup);
  h_calibjet_pt_puritycorr_dijet_jerdown->Divide(h_purity_calib_dijet_jerdown);
  h_calibjet_pt_puritycorr_dijet_jerup->Divide(h_purity_calib_dijet_jerup);

  h_calibjet_pt_puritycorr_frac->Divide(h_purity_calib_frac);
  h_calibjet_pt_puritycorr_frac_effdown->Divide(h_purity_calib_frac);
  h_calibjet_pt_puritycorr_frac_effup->Divide(h_purity_calib_frac);
  h_calibjet_pt_puritycorr_frac_jesdown->Divide(h_purity_calib_frac_jesdown);
  h_calibjet_pt_puritycorr_frac_jesup->Divide(h_purity_calib_frac_jesup);
  h_calibjet_pt_puritycorr_frac_jerdown->Divide(h_purity_calib_frac_jerdown);
  h_calibjet_pt_puritycorr_frac_jerup->Divide(h_purity_calib_frac_jerup);

  // Write histograms
  f_dataout->cd();
  h_calibjet_pt_puritycorr_dijet->Write();
  h_calibjet_pt_puritycorr_dijet_effdown->Write();
  h_calibjet_pt_puritycorr_dijet_effup->Write();
  h_calibjet_pt_puritycorr_dijet_jesdown->Write();
  h_calibjet_pt_puritycorr_dijet_jesup->Write();
  h_calibjet_pt_puritycorr_dijet_jerdown->Write();
  h_calibjet_pt_puritycorr_dijet_jerup->Write();

  h_calibjet_pt_puritycorr_frac->Write();
  h_calibjet_pt_puritycorr_frac_effdown->Write();
  h_calibjet_pt_puritycorr_frac_effup->Write();
  h_calibjet_pt_puritycorr_frac_jesdown->Write();
  h_calibjet_pt_puritycorr_frac_jesup->Write();
  h_calibjet_pt_puritycorr_frac_jerdown->Write();
  h_calibjet_pt_puritycorr_frac_jerup->Write();
  f_dataout->Close();
  */

  /*
  // Draw purity and efficiency.
  TH1D* h_purity_calib_dijet_proj = (TH1D*)h_respmatrix_calib_dijet->Hmeasured();
  TH1D* h_efficiency_calib_dijet_proj = (TH1D*)h_respmatrix_calib_dijet->Htruth();

  TH1D* h_purity_calib_dijet_jesdown_proj = (TH1D*)h_respmatrix_calib_dijet_jesdown->Hmeasured();
  TH1D* h_efficiency_calib_dijet_jesdown_proj = (TH1D*)h_respmatrix_calib_dijet_jesdown->Htruth();

  TH1D* h_purity_calib_dijet_jesup_proj = (TH1D*)h_respmatrix_calib_dijet_jesup->Hmeasured();
  TH1D* h_efficiency_calib_dijet_jesup_proj = (TH1D*)h_respmatrix_calib_dijet_jesup->Htruth();

  TH1D* h_purity_calib_dijet_jerdown_proj = (TH1D*)h_respmatrix_calib_dijet_jerdown->Hmeasured();
  TH1D* h_efficiency_calib_dijet_jerdown_proj = (TH1D*)h_respmatrix_calib_dijet_jerdown->Htruth();

  TH1D* h_purity_calib_dijet_jerup_proj = (TH1D*)h_respmatrix_calib_dijet_jerup->Hmeasured();
  TH1D* h_efficiency_calib_dijet_jerup_proj = (TH1D*)h_respmatrix_calib_dijet_jerup->Htruth();

  TH1D* h_purity_calib_frac_proj = (TH1D*)h_respmatrix_calib_frac->Hmeasured();
  TH1D* h_efficiency_calib_frac_proj = (TH1D*)h_respmatrix_calib_frac->Htruth();

  TH1D* h_purity_calib_frac_jesdown_proj = (TH1D*)h_respmatrix_calib_frac_jesdown->Hmeasured();
  TH1D* h_efficiency_calib_frac_jesdown_proj = (TH1D*)h_respmatrix_calib_frac_jesdown->Htruth();

  TH1D* h_purity_calib_frac_jesup_proj = (TH1D*)h_respmatrix_calib_frac_jesup->Hmeasured();
  TH1D* h_efficiency_calib_frac_jesup_proj = (TH1D*)h_respmatrix_calib_frac_jesup->Htruth();

  TH1D* h_purity_calib_frac_jerdown_proj = (TH1D*)h_respmatrix_calib_frac_jerdown->Hmeasured();
  TH1D* h_efficiency_calib_frac_jerdown_proj = (TH1D*)h_respmatrix_calib_frac_jerdown->Htruth();

  TH1D* h_purity_calib_frac_jerup_proj = (TH1D*)h_respmatrix_calib_frac_jerup->Hmeasured();
  TH1D* h_efficiency_calib_frac_jerup_proj = (TH1D*)h_respmatrix_calib_frac_jerup->Htruth();
  */

  std::vector<TH1F*> h_input;
  std::vector<int> color;
  std::vector<int> markerstyle;
  std::vector<std::string> text;
  std::vector<std::string> legend;

  h_input.push_back((TH1F*)h_measure_jet_pt_calib_dijet);
  h_input.push_back((TH1F*)h_purity_jet_pt_calib_dijet);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Simulation");
  text.push_back("PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  //text.push_back("#DeltaR_{matching} < 0.3");
  legend.push_back("Reco jet spectrum");
  legend.push_back("Matched reco jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              false, 0, 0.5, true,
                              true, 0.0, 1.1,
                              true, "p_{T}^{reco jet} [GeV]", "Arbitrary Unit", "Purity", 1,
                              true, text, 0.45, 0.9, 0.05,
                              true, legend, 0.25, 0.2, 0.05,
                              "test_figure/purity_jet_pt_calib_dijet.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_measure_jet_pt_calib_dijet_jesdown);
  h_input.push_back((TH1F*)h_purity_jet_pt_calib_dijet_jesdown);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Simulation");
  text.push_back("PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  //text.push_back("#DeltaR_{matching} < 0.3");
  //text.push_back("With JES jet_pt_calibration");
  legend.push_back("Reco jet spectrum");
  legend.push_back("Matched reco jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              false, 0, 0.5, true,
                              true, 0.0, 1.1,
                              true, "p_{T}^{reco jet} [GeV]", "Arbitrary Unit", "Purity", 1,
                              true, text, 0.45, 0.9, 0.05,
                              true, legend, 0.25, 0.2, 0.05,
                              "test_figure/purity_jet_pt_calib_dijet_jesdown.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_measure_jet_pt_calib_dijet_jesup);
  h_input.push_back((TH1F*)h_purity_jet_pt_calib_dijet_jesup);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Simulation");
  text.push_back("PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  //text.push_back("#DeltaR_{matching} < 0.3");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Matched truth jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              false, 0, 0.5, true,
                              true, 0.0, 1.1,
                              true, "p_{T}^{truth jet} [GeV]", "Arbitrary Unit", "Purity", 1,
                              true, text, 0.45, 0.9, 0.05,
                              true, legend, 0.25, 0.2, 0.05,
                              "test_figure/purity_jet_pt_calib_dijet_jesup.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_measure_jet_pt_calib_dijet_jerdown);
  h_input.push_back((TH1F*)h_purity_jet_pt_calib_dijet_jerdown);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Simulation");
  text.push_back("PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  //text.push_back("#DeltaR_{matching} < 0.3");
  //text.push_back("With JES jet_pt_calibration");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Matched truth jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              false, 0, 0.5, true,
                              true, 0.0, 1.1,
                              true, "p_{T}^{truth jet} [GeV]", "Arbitrary Unit", "Purity", 1,
                              true, text, 0.45, 0.9, 0.05,
                              true, legend, 0.25, 0.2, 0.05,
                              "test_figure/purity_jet_pt_calib_dijet_jerdown.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_measure_jet_pt_calib_dijet_jerup);
  h_input.push_back((TH1F*)h_purity_jet_pt_calib_dijet_jerup);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Simulation");
  text.push_back("PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  //text.push_back("#DeltaR_{matching} < 0.3");
  //text.push_back("With JES jet_pt_calibration");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Matched truth jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              false, 0, 0.5, true,
                              true, 0.0, 1.1,
                              true, "p_{T}^{truth jet} [GeV]", "Arbitrary Unit", "Purity", 1,
                              true, text, 0.45, 0.9, 0.05,
                              true, legend, 0.25, 0.2, 0.05,
                              "test_figure/purity_jet_pt_calib_dijet_jerup.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_truth_jet_pt_calib_dijet);
  h_input.push_back((TH1F*)h_efficiency_jet_pt_calib_dijet);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Simulation");
  text.push_back("PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  //text.push_back("#DeltaR_{matching} < 0.3");
  legend.push_back("Reco jet spectrum");
  legend.push_back("Matched reco jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              true, 1e-9, 1e-1, true,
                              true, 0., 1.2,
                              true, "p_{T}^{reco jet} [GeV]", "Arbitrary Unit", "Efficiency", 1,
                              true, text, 0.45, 0.9, 0.05,
                              true, legend, 0.25, 0.2, 0.05,
                              "test_figure/efficiency_jet_pt_calib_dijet.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_truth_jet_pt_calib_dijet_jesdown);
  h_input.push_back((TH1F*)h_efficiency_jet_pt_calib_dijet_jesdown);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Simulation");
  text.push_back("PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  //text.push_back("#DeltaR_{matching} < 0.3");
  //text.push_back("With JES jet_pt_calibration");
  legend.push_back("Reco jet spectrum");
  legend.push_back("Matched reco jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              true, 1e-9, 1e-1, true,
                              true, 0., 1.2,
                              true, "p_{T}^{reco jet} [GeV]", "Arbitrary Unit", "Efficiency", 1,
                              true, text, 0.45, 0.9, 0.05,
                              true, legend, 0.25, 0.2, 0.05,
                              "test_figure/efficiency_jet_pt_calib_dijet_jesdown.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_truth_jet_pt_calib_dijet_jesup);
  h_input.push_back((TH1F*)h_efficiency_jet_pt_calib_dijet_jesup);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Simulation");
  text.push_back("PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  //text.push_back("#DeltaR_{matching} < 0.3");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Matched truth jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              true, 1e-9, 1e-1, true,
                              true, 0, 1.2,
                              true, "p_{T}^{truth jet} [GeV]", "Arbitrary Unit", "Efficiency", 1,
                              true, text, 0.45, 0.9, 0.05,
                              true, legend, 0.25, 0.2, 0.05,
                              "test_figure/efficiency_jet_pt_calib_dijet_jesup.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_truth_jet_pt_calib_dijet_jerdown);
  h_input.push_back((TH1F*)h_efficiency_jet_pt_calib_dijet_jerdown);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Simulation");
  text.push_back("PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  //text.push_back("#DeltaR_{matching} < 0.3");
  //text.push_back("With JES jet_pt_calibration");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Matched truth jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              true, 1e-9, 1e-1, true,
                              true, 0, 1.2,
                              true, "p_{T}^{truth jet} [GeV]", "Arbitrary Unit", "Efficiency", 1,
                              true, text, 0.45, 0.9, 0.05,
                              true, legend, 0.25, 0.2, 0.05,
                              "test_figure/efficiency_jet_pt_calib_dijet_jerdown.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();

  h_input.push_back((TH1F*)h_truth_jet_pt_calib_dijet_jerup);
  h_input.push_back((TH1F*)h_efficiency_jet_pt_calib_dijet_jerup);
  color.push_back(kBlack);
  color.push_back(kRed);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Simulation");
  text.push_back("PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  //text.push_back("#DeltaR_{matching} < 0.3");
  //text.push_back("With JES jet_pt_calibration");
  legend.push_back("Truth jet spectrum");
  legend.push_back("Matched truth jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              true, 1e-9, 1e-1, true,
                              true, 0, 1.2,
                              true, "p_{T}^{truth jet} [GeV]", "Arbitrary Unit", "Efficiency", 1,
                              true, text, 0.45, 0.9, 0.05,
                              true, legend, 0.25, 0.2, 0.05,
                              "test_figure/efficiency_jet_pt_calib_dijet_jerup.pdf");
  h_input.clear();
  color.clear();
  markerstyle.clear();
  text.clear();
  legend.clear();
}