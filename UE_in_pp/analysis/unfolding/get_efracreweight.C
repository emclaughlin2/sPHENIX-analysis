#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>

void get_efracreweight(const char* simfile = "sphenix_primary_analysis_sim_run28_output/output_dijet_bkg_cut_sim_iter_3.root", 
                        const char* datafile = "analysis_data_run28_output/output_pu_correct_data_dijet_bkg_cut.root", 
                        const char* reweightfile = "output_efrac_reweight_run28_dijet_bkg_cut_iter_1.root") {
  // Read Files
  //********** Files **********//
  TFile *f_sim  = new TFile(simfile,  "READ");
  TFile *f_data = new TFile(datafile, "READ");

  TH2D* h_2D_efrac_sim  = (TH2D*)f_sim ->Get("h_efrac_lead_pt_dijet");
  h_2D_efrac_sim ->SetName("h_2D_efrac_sim");
  TH2D* h_2D_efrac_data = (TH2D*)f_data->Get("h_efrac_lead_pt_dijet");
  h_2D_efrac_data->SetName("h_2D_efrac_data");

  // Clone structure for pT-slice-normalized versions
  TH2D* h_sim_norm  = (TH2D*)h_2D_efrac_sim ->Clone("h_sim_norm");
  TH2D* h_data_norm = (TH2D*)h_2D_efrac_data->Clone("h_data_norm");
  h_sim_norm ->Reset();
  h_data_norm->Reset();

  int nBinsX = h_2D_efrac_sim->GetNbinsX();
  int nBinsY = h_2D_efrac_sim->GetNbinsY();

  // For each pT bin, normalize the efrac (Y) distribution over iy 1-10 only
  for (int ix = 1; ix <= nBinsX; ix++) {

      double sum_sim  = 0.0;
      double sum_data = 0.0;
      for (int iy = 1; iy <= 10; iy++) {
          sum_sim  += h_2D_efrac_sim ->GetBinContent(ix, iy);
          sum_data += h_2D_efrac_data->GetBinContent(ix, iy);
      }

      if (sum_sim == 0 || sum_data == 0) {
          std::cout << "WARNING: empty pT bin ix=" << ix
                    << "  sum_sim=" << sum_sim
                    << "  sum_data=" << sum_data
                    << " -- skipping (weights left at 0)" << std::endl;
          continue;
      }

      for (int iy = 1; iy <= nBinsY; iy++) {

          // iy 11-15: no reweighting, set both to 1 so ratio = 1
          if (iy > 10) {
              h_sim_norm ->SetBinContent(ix, iy, 1.0);
              h_sim_norm ->SetBinError  (ix, iy, 0.0);
              h_data_norm->SetBinContent(ix, iy, 1.0);
              h_data_norm->SetBinError  (ix, iy, 0.0);
              continue;
          }

          h_sim_norm ->SetBinContent(ix, iy, h_2D_efrac_sim ->GetBinContent(ix, iy) / sum_sim);
          h_sim_norm ->SetBinError  (ix, iy, h_2D_efrac_sim ->GetBinError  (ix, iy) / sum_sim);
          h_data_norm->SetBinContent(ix, iy, h_2D_efrac_data->GetBinContent(ix, iy) / sum_data);
          h_data_norm->SetBinError  (ix, iy, h_2D_efrac_data->GetBinError  (ix, iy) / sum_data);
      }
  }

  // data / MC → 2D weight histogram
  TH2D* efrac_weight = (TH2D*)h_data_norm->Clone("efrac_weight");
  efrac_weight->Divide(h_sim_norm);

  // Safety renormalization over iy 1-10 only
  // norm should be ~1; any deviation indicates empty bins causing yield drift
  for (int ix = 1; ix <= nBinsX; ix++) {

      double sum_sim      = 0.0;
      double weighted_sum = 0.0;
      for (int iy = 1; iy <= 10; iy++) {
          double s = h_2D_efrac_sim->GetBinContent(ix, iy);
          sum_sim      += s;
          weighted_sum += efrac_weight->GetBinContent(ix, iy) * s;
      }

      if (sum_sim == 0 || weighted_sum == 0) continue;

      double norm = weighted_sum / sum_sim;
      if (std::abs(norm - 1.0) > 0.01) {
          std::cout << "WARNING: pT bin ix=" << ix
                    << " has weighted yield norm=" << norm
                    << " (deviates >1% from 1) -- correcting" << std::endl;
      }

      for (int iy = 1; iy <= 10; iy++) {
          efrac_weight->SetBinContent(ix, iy, efrac_weight->GetBinContent(ix, iy) / norm);
          efrac_weight->SetBinError  (ix, iy, efrac_weight->GetBinError  (ix, iy) / norm);
      }
      // iy 11-15 are exactly 1.0 and unaffected
  }

  // Verify final weighted yields match original MC yields per pT bin
  std::cout << "\n--- Normalization check: weighted MC yield vs original per pT bin ---" << std::endl;
  for (int ix = 1; ix <= nBinsX; ix++) {
      double orig_sum     = 0.0;
      double weighted_sum = 0.0;
      for (int iy = 1; iy <= 10; iy++) {
          double s = h_2D_efrac_sim->GetBinContent(ix, iy);
          orig_sum     += s;
          weighted_sum += efrac_weight->GetBinContent(ix, iy) * s;
      }
      double pt_lo = h_2D_efrac_sim->GetXaxis()->GetBinLowEdge(ix);
      double pt_hi = h_2D_efrac_sim->GetXaxis()->GetBinUpEdge(ix);
      std::cout << "  pT bin ix=" << ix
                << " [" << pt_lo << ", " << pt_hi << "]"
                << "  orig=" << orig_sum
                << "  weighted=" << weighted_sum;
      if (orig_sum > 0)
          std::cout << "  ratio=" << weighted_sum / orig_sum;
      std::cout << std::endl;
  }

  TFile *f_out = new TFile(reweightfile, "RECREATE");
  f_out->cd();
  h_2D_efrac_sim ->Write();
  h_2D_efrac_data->Write();
  h_sim_norm     ->Write();
  h_data_norm    ->Write();
  efrac_weight   ->Write();
  f_out->Close();

}