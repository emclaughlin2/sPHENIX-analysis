#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>

void get_vertexreweight(const char* simfile = "analysis_sim_run28_output/output_none_bkg_cut_sim_iter_3.root", 
                        const char* datafile = "analysis_data_run28_output/output_pu_correct_data_efrac_bkg_cut.root", 
                        const char* reweightfile = "output_vertex_reweight_run28_none_bkg_cut_iter_1.root") {
  // Read Files
  //********** Files **********//
  TFile *f_sim = new TFile(simfile, "READ");
  TFile *f_data = new TFile(datafile, "READ");

  TH1D* h_zvertex_sim = (TH1D*)f_sim->Get("h_zvertex");
  h_zvertex_sim->SetName("h_zvertex_sim");
  TH1D* h_zvertex_data = (TH1D*)f_data->Get("h_zvertex");
  h_zvertex_data->SetName("h_zvertex_data");

  TH1D* h_zvertex_sim_resize = new TH1D("h_zvertex_sim_resize","",120,-60,60);
  for (int i = 1; i < h_zvertex_sim_resize->GetNbinsX() + 1; i++) {
    h_zvertex_sim_resize->SetBinContent(i, h_zvertex_sim->GetBinContent(140+i));
    h_zvertex_sim_resize->SetBinError(i, h_zvertex_sim->GetBinError(140+i));
  }

  h_zvertex_sim_resize->Scale(1.0/ h_zvertex_sim_resize->Integral());
  h_zvertex_sim->Scale(1.0/ h_zvertex_sim->Integral());
  h_zvertex_data->Scale(1.0/ h_zvertex_data->Integral());

  TH1D* vertex_weight = (TH1D*)h_zvertex_data->Clone("vertex_weight");
  vertex_weight->Divide(h_zvertex_sim_resize);

  TFile *f_out = new TFile(reweightfile, "RECREATE");
  f_out->cd();
  h_zvertex_sim->Write();
  h_zvertex_sim_resize->Write();
  h_zvertex_data->Write();
  vertex_weight->Write();
  f_out->Close();

}