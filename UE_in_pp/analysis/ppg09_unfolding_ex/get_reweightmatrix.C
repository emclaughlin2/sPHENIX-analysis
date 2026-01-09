#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <iostream>
#include "/sphenix/user/hanpuj/CaloDataAna24_skimmed/src/draw_template.C"

void draw_2D_plot2(TH2F *h_input,
                   bool do_rebin, int rebin_xfactor, int rebin_yfactor, bool do_normalize,
                   bool set_xrange, float xlow, float xhigh, bool set_logx,
                   bool set_yrange, float ylow, float yhigh, bool set_logy,
                   bool set_zrange, float zlow, float zhigh, bool set_logz,
                   bool set_title, std::string xtitle, std::string ytitle, std::string ztitle,
                   bool set_text, std::vector<std::string> text, std::string text2, float xstart, float ystart, float size,
                   std::string output_name) {
  TH2F *h = (TH2F*)h_input->Clone("h");
  TCanvas *can = new TCanvas("can", "", 800, 700);
  gStyle->SetPalette(57);
  can->SetTopMargin(0.12);
  can->SetLeftMargin(0.15);
  can->SetBottomMargin(0.15);
  can->SetRightMargin(0.2);
  if (do_rebin) h->Rebin2D(rebin_xfactor, rebin_yfactor);
  if (do_normalize) {
    for (int ixb = 1; ixb <= h->GetNbinsX(); ++ixb) {
      for (int iyb = 1; iyb <= h->GetNbinsY(); ++iyb) {
        h->SetBinContent(ixb, iyb, h->GetBinContent(ixb, iyb)/(float)(h->GetBinWidth(ixb)*h->GetBinWidth(iyb)));
        h->SetBinError(ixb, iyb, h->GetBinError(ixb, iyb)/(float)(h->GetBinWidth(ixb)*h->GetBinWidth(iyb)));
      }
    }
  }
  if (set_xrange) h->GetXaxis()->SetRangeUser(xlow, xhigh);
  if (set_logx) can->SetLogx();
  if (set_yrange) h->GetYaxis()->SetRangeUser(ylow, yhigh);
  if (set_logy) can->SetLogy();
  if (set_zrange) h->GetZaxis()->SetRangeUser(zlow, zhigh);
  if (set_logz) can->SetLogz();
  if (set_title) {
    h->GetXaxis()->SetTitle(xtitle.c_str());
    h->GetYaxis()->SetTitle(ytitle.c_str());
    h->GetZaxis()->SetTitle(ztitle.c_str());
  }
  h->GetXaxis()->SetTitleSize(0.065);
  h->GetYaxis()->SetTitleSize(0.065);
  h->GetZaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetTitleOffset(1.04);
  h->GetYaxis()->SetTitleOffset(1.06);
  h->GetZaxis()->SetTitleOffset(1.22);
  h->GetXaxis()->SetLabelSize(0.055);
  h->GetYaxis()->SetLabelSize(0.055);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetZaxis()->CenterTitle();
  h->Draw("colz");
  myText(0.03, 0.97, 1, "#bf{#it{sPHENIX}} Simulation    PYTHIA8 p+p#sqrt{s} = 200 GeV", 0.05);
  myText(0.03, 0.91, 1, text2.c_str(), 0.05);
  //myText(0.5, 0.97, 1, "c", 0.05);
  if (set_text) {
    for (int i = 0; i < text.size(); i++) {
      myText(xstart, ystart-i*(size+0.01), 1, text[i].c_str(), size);
    }
  }
  can->SaveAs(output_name.c_str());
  delete can;
}


void get_reweight_hist(TFile* f_data, TH2D* h_respmatrix, TH1D*& h_reweight, std::string input_histname, std::string output_histname) { 
  TH1D* h_data_reco = (TH1D*)f_data->Get(input_histname.c_str());
  TH1D* h_sim_reco = (TH1D*)h_respmatrix->ProjectionX("h_sim_reco");
  h_reweight = (TH1D*)h_data_reco->Clone(output_histname.c_str());
  double normalization_reco = h_data_reco->Integral() / h_sim_reco->Integral();
  h_sim_reco->Scale(normalization_reco);
  h_reweight->Divide(h_sim_reco);

  std::string prefix_to_remove = "h_reweight_";
  std::string new_prefix = "reweightfactor_";
  std::string plot_name = new_prefix + output_histname.substr(prefix_to_remove.length());

  std::vector<TH1F*> h_input;
  std::vector<int> color;
  std::vector<int> markerstyle;
  std::vector<std::string> text;
  std::vector<std::string> legend;

  h_input.push_back((TH1F*)h_sim_reco);
  h_input.push_back((TH1F*)h_data_reco);
  color.push_back(kRed);
  color.push_back(kBlack);
  markerstyle.push_back(20);
  markerstyle.push_back(20);
  text.push_back("#bf{#it{sPHENIX}} Internal");
  text.push_back("Data & PYTHIA8 p+p#sqrt{s} = 200 GeV");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4");
  text.push_back("|#eta^{jet}| < 0.7");
  //text.push_back("#DeltaR_{matching} < 0.3");
  legend.push_back("Simulation jet spectrum");
  legend.push_back("Data jet spectrum");
  draw_1D_multiple_plot_ratio(h_input, color, markerstyle,
                              false, 10, true,
                              true, 15, 72, false,
                              false, 0, 0.5, true,
                              true, 0., 2.,
                              true, "p_{T}^{reco jet} [GeV]", "Arbitrary Unit", "Reweight factor", 0,
                              true, text, 0.4, 0.9, 0.05,
                              true, legend, 0.25, 0.2, 0.05,
                              Form("figure/%s.png", plot_name.c_str()));
  h_input.clear();
  color.clear();
  markerstyle.clear();
   
  legend.clear();
}

void get_reweightmatrix() {
  // Read Files
  TFile* f_data = new TFile("output_data_puritycorr.root", "READ");
  TFile* f_sim = new TFile("output_sim.root", "READ");
  TFile* f_out = new TFile("output_reweightmatrix.root", "RECREATE");

  TH2D* h_respmatrix_calib_dijet = (TH2D*)f_sim->Get("h_respmatrix_calib_dijet");
  TH2D* h_respmatrix_calib_dijet_jesdown = (TH2D*)f_sim->Get("h_respmatrix_calib_dijet_jesdown");
  TH2D* h_respmatrix_calib_dijet_jesup = (TH2D*)f_sim->Get("h_respmatrix_calib_dijet_jesup");
  TH2D* h_respmatrix_calib_dijet_jerdown = (TH2D*)f_sim->Get("h_respmatrix_calib_dijet_jerdown");
  TH2D* h_respmatrix_calib_dijet_jerup = (TH2D*)f_sim->Get("h_respmatrix_calib_dijet_jerup");
  
  TH2D* h_respmatrix_calib_frac = (TH2D*)f_sim->Get("h_respmatrix_calib_frac");
  TH2D* h_respmatrix_calib_frac_jesdown = (TH2D*)f_sim->Get("h_respmatrix_calib_frac_jesdown");
  TH2D* h_respmatrix_calib_frac_jesup = (TH2D*)f_sim->Get("h_respmatrix_calib_frac_jesup");
  TH2D* h_respmatrix_calib_frac_jerdown = (TH2D*)f_sim->Get("h_respmatrix_calib_frac_jerdown");
  TH2D* h_respmatrix_calib_frac_jerup = (TH2D*)f_sim->Get("h_respmatrix_calib_frac_jerup");

  // Form reweight histograms.
  TH1D* h_reweight_calib_dijet; get_reweight_hist(f_data, h_respmatrix_calib_dijet, h_reweight_calib_dijet, "h_calibjet_pt_puritycorr_dijet", "h_reweight_calib_dijet");
  TH1D* h_reweight_calib_dijet_effdown; get_reweight_hist(f_data, h_respmatrix_calib_dijet, h_reweight_calib_dijet_effdown, "h_calibjet_pt_puritycorr_dijet_effdown", "h_reweight_calib_dijet_effdown");
  TH1D* h_reweight_calib_dijet_effup; get_reweight_hist(f_data, h_respmatrix_calib_dijet, h_reweight_calib_dijet_effup, "h_calibjet_pt_puritycorr_dijet_effup", "h_reweight_calib_dijet_effup");
  TH1D* h_reweight_calib_dijet_jesdown; get_reweight_hist(f_data, h_respmatrix_calib_dijet_jesdown, h_reweight_calib_dijet_jesdown, "h_calibjet_pt_puritycorr_dijet_jesdown", "h_reweight_calib_dijet_jesdown");
  TH1D* h_reweight_calib_dijet_jesup; get_reweight_hist(f_data, h_respmatrix_calib_dijet_jesup, h_reweight_calib_dijet_jesup, "h_calibjet_pt_puritycorr_dijet_jesup", "h_reweight_calib_dijet_jesup");
  TH1D* h_reweight_calib_dijet_jerdown; get_reweight_hist(f_data, h_respmatrix_calib_dijet_jerdown, h_reweight_calib_dijet_jerdown, "h_calibjet_pt_puritycorr_dijet_jerdown", "h_reweight_calib_dijet_jerdown");
  TH1D* h_reweight_calib_dijet_jerup; get_reweight_hist(f_data, h_respmatrix_calib_dijet_jerup, h_reweight_calib_dijet_jerup, "h_calibjet_pt_puritycorr_dijet_jerup", "h_reweight_calib_dijet_jerup");

  TH1D* h_reweight_calib_frac; get_reweight_hist(f_data, h_respmatrix_calib_frac, h_reweight_calib_frac, "h_calibjet_pt_puritycorr_frac", "h_reweight_calib_frac");
  TH1D* h_reweight_calib_frac_effdown; get_reweight_hist(f_data, h_respmatrix_calib_frac, h_reweight_calib_frac_effdown, "h_calibjet_pt_puritycorr_frac_effdown", "h_reweight_calib_frac_effdown");
  TH1D* h_reweight_calib_frac_effup; get_reweight_hist(f_data, h_respmatrix_calib_frac, h_reweight_calib_frac_effup, "h_calibjet_pt_puritycorr_frac_effup", "h_reweight_calib_frac_effup");
  TH1D* h_reweight_calib_frac_jesdown; get_reweight_hist(f_data, h_respmatrix_calib_frac_jesdown, h_reweight_calib_frac_jesdown, "h_calibjet_pt_puritycorr_frac_jesdown", "h_reweight_calib_frac_jesdown");
  TH1D* h_reweight_calib_frac_jesup; get_reweight_hist(f_data, h_respmatrix_calib_frac_jesup, h_reweight_calib_frac_jesup, "h_calibjet_pt_puritycorr_frac_jesup", "h_reweight_calib_frac_jesup");
  TH1D* h_reweight_calib_frac_jerdown; get_reweight_hist(f_data, h_respmatrix_calib_frac_jerdown, h_reweight_calib_frac_jerdown, "h_calibjet_pt_puritycorr_frac_jerdown", "h_reweight_calib_frac_jerdown");
  TH1D* h_reweight_calib_frac_jerup; get_reweight_hist(f_data, h_respmatrix_calib_frac_jerup, h_reweight_calib_frac_jerup, "h_calibjet_pt_puritycorr_frac_jerup", "h_reweight_calib_frac_jerup");

  // Form reweighted response matrix.
  TH2D* h_respmatrix_calib_dijet_reweighted = (TH2D*)h_respmatrix_calib_dijet->Clone("h_respmatrix_calib_dijet_reweighted");
  TH2D* h_respmatrix_calib_dijet_effdown_reweighted = (TH2D*)h_respmatrix_calib_dijet->Clone("h_respmatrix_calib_dijet_effdown_reweighted");
  TH2D* h_respmatrix_calib_dijet_effup_reweighted = (TH2D*)h_respmatrix_calib_dijet->Clone("h_respmatrix_calib_dijet_effup_reweighted");
  TH2D* h_respmatrix_calib_dijet_jesdown_reweighted = (TH2D*)h_respmatrix_calib_dijet_jesdown->Clone("h_respmatrix_calib_dijet_jesdown_reweighted");
  TH2D* h_respmatrix_calib_dijet_jesup_reweighted = (TH2D*)h_respmatrix_calib_dijet_jesup->Clone("h_respmatrix_calib_dijet_jesup_reweighted");
  TH2D* h_respmatrix_calib_dijet_jerdown_reweighted = (TH2D*)h_respmatrix_calib_dijet_jerdown->Clone("h_respmatrix_calib_dijet_jerdown_reweighted");
  TH2D* h_respmatrix_calib_dijet_jerup_reweighted = (TH2D*)h_respmatrix_calib_dijet_jerup->Clone("h_respmatrix_calib_dijet_jerup_reweighted");

  TH2D* h_respmatrix_calib_frac_reweighted = (TH2D*)h_respmatrix_calib_frac->Clone("h_respmatrix_calib_frac_reweighted");
  TH2D* h_respmatrix_calib_frac_effdown_reweighted = (TH2D*)h_respmatrix_calib_frac->Clone("h_respmatrix_calib_frac_effdown_reweighted");
  TH2D* h_respmatrix_calib_frac_effup_reweighted = (TH2D*)h_respmatrix_calib_frac->Clone("h_respmatrix_calib_frac_effup_reweighted");
  TH2D* h_respmatrix_calib_frac_jesdown_reweighted = (TH2D*)h_respmatrix_calib_frac_jesdown->Clone("h_respmatrix_calib_frac_jesdown_reweighted");
  TH2D* h_respmatrix_calib_frac_jesup_reweighted = (TH2D*)h_respmatrix_calib_frac_jesup->Clone("h_respmatrix_calib_frac_jesup_reweighted");
  TH2D* h_respmatrix_calib_frac_jerdown_reweighted = (TH2D*)h_respmatrix_calib_frac_jerdown->Clone("h_respmatrix_calib_frac_jerdown_reweighted");
  TH2D* h_respmatrix_calib_frac_jerup_reweighted = (TH2D*)h_respmatrix_calib_frac_jerup->Clone("h_respmatrix_calib_frac_jerup_reweighted");

  int nbin_meas = h_respmatrix_calib_dijet_reweighted->GetNbinsX();
  int nbin_truth = h_respmatrix_calib_dijet_reweighted->GetNbinsY();
  for (int ibm = 1; ibm <= nbin_meas; ++ibm) {
    for (int ibt = 1; ibt <= nbin_truth; ++ibt) {
      h_respmatrix_calib_dijet_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_dijet_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_dijet->GetBinContent(ibm));
      h_respmatrix_calib_dijet_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_dijet_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_dijet->GetBinContent(ibm));

      h_respmatrix_calib_dijet_effdown_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_dijet_effdown_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_dijet_effdown->GetBinContent(ibm));
      h_respmatrix_calib_dijet_effdown_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_dijet_effdown_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_dijet_effdown->GetBinContent(ibm));

      h_respmatrix_calib_dijet_effup_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_dijet_effup_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_dijet_effup->GetBinContent(ibm));
      h_respmatrix_calib_dijet_effup_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_dijet_effup_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_dijet_effup->GetBinContent(ibm));

      h_respmatrix_calib_dijet_jesdown_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_dijet_jesdown_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_dijet_jesdown->GetBinContent(ibm));
      h_respmatrix_calib_dijet_jesdown_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_dijet_jesdown_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_dijet_jesdown->GetBinContent(ibm));

      h_respmatrix_calib_dijet_jesup_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_dijet_jesup_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_dijet_jesup->GetBinContent(ibm));
      h_respmatrix_calib_dijet_jesup_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_dijet_jesup_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_dijet_jesup->GetBinContent(ibm));

      h_respmatrix_calib_dijet_jerdown_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_dijet_jerdown_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_dijet_jerdown->GetBinContent(ibm));
      h_respmatrix_calib_dijet_jerdown_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_dijet_jerdown_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_dijet_jerdown->GetBinContent(ibm));

      h_respmatrix_calib_dijet_jerup_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_dijet_jerup_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_dijet_jerup->GetBinContent(ibm));
      h_respmatrix_calib_dijet_jerup_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_dijet_jerup_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_dijet_jerup->GetBinContent(ibm));

      h_respmatrix_calib_frac_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_frac_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_frac->GetBinContent(ibm));
      h_respmatrix_calib_frac_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_frac_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_frac->GetBinContent(ibm));

      h_respmatrix_calib_frac_effdown_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_frac_effdown_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_frac_effdown->GetBinContent(ibm));
      h_respmatrix_calib_frac_effdown_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_frac_effdown_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_frac_effdown->GetBinContent(ibm));

      h_respmatrix_calib_frac_effup_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_frac_effup_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_frac_effup->GetBinContent(ibm));
      h_respmatrix_calib_frac_effup_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_frac_effup_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_frac_effup->GetBinContent(ibm));

      h_respmatrix_calib_frac_jesdown_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_frac_jesdown_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_frac_jesdown->GetBinContent(ibm));
      h_respmatrix_calib_frac_jesdown_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_frac_jesdown_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_frac_jesdown->GetBinContent(ibm));

      h_respmatrix_calib_frac_jesup_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_frac_jesup_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_frac_jesup->GetBinContent(ibm));
      h_respmatrix_calib_frac_jesup_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_frac_jesup_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_frac_jesup->GetBinContent(ibm));

      h_respmatrix_calib_frac_jerdown_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_frac_jerdown_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_frac_jerdown->GetBinContent(ibm));
      h_respmatrix_calib_frac_jerdown_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_frac_jerdown_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_frac_jerdown->GetBinContent(ibm));

      h_respmatrix_calib_frac_jerup_reweighted->SetBinContent(ibm, ibt, h_respmatrix_calib_frac_jerup_reweighted->GetBinContent(ibm, ibt) * h_reweight_calib_frac_jerup->GetBinContent(ibm));
      h_respmatrix_calib_frac_jerup_reweighted->SetBinError(ibm, ibt, h_respmatrix_calib_frac_jerup_reweighted->GetBinError(ibm, ibt) * h_reweight_calib_frac_jerup->GetBinContent(ibm));
    }
  }

  // Draw response matrix.
  std::vector<std::string> text;
  std::string text2;
  //text.push_back("");
  //text.push_back("#DeltaR_{matching} < 0.3");
  text.push_back("anti-k_{t} #kern[-0.5]{#it{R}} = 0.4, |#eta^{jet}| < 0.7");

  text2 = "Dijet requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_calib_dijet.png");
   

  text2 = "Dijet requirement";            
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet_jesdown,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_calib_dijet_jesdown.png");
   

  text2 = "Dijet requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet_jesup,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_calib_dijet_jesup.png");
   

  text2 = "Dijet requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet_jerdown,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_calib_dijet_jerdown.png");
   

  text2 = "Dijet requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet_jerup,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_calib_dijet_jerup.png");
   

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_calib_frac.png");
   

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac_jesdown,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_calib_frac_jesdown.png");
   

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac_jesup,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_calib_frac_jesup.png");
   

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac_jerdown,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_calib_frac_jerdown.png");
   

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac_jerup,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_calib_frac_jerup.png");
   

  text2 = "Dijet requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_dijet.png");
   

  text2 = "Dijet requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet_effdown_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_dijet_effdown.png");
   

  text2 = "Dijet requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet_effup_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_dijet_effup.png");
   

  text2 = "Dijet requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet_jesdown_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_dijet_jesdown.png");
   

  text2 = "Dijet requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet_jesup_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_dijet_jesup.png");
   

  text2 = "Dijet requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet_jerdown_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_dijet_jerdown.png");
   

  text2 = "Dijet requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_dijet_jerup_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_dijet_jerup.png");
   

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_frac.png");
   

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac_effdown_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_frac_effdown.png");
   

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac_effup_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_frac_effup.png");
   

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac_jesdown_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_frac_jesdown.png");
                

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac_jesup_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_frac_jesup.png");
   

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac_jerdown_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_frac_jerdown.png");
   

  text2 = "Energy fraction requirement";
  draw_2D_plot2((TH2F*)h_respmatrix_calib_frac_jerup_reweighted,
               false, 1, 1, false,
               false, 0, 40, false,
               false, 0, 40, false,
               true, 1e-9, 2, true,
               true, "p_{T}^{reco} [GeV]", "p_{T}^{truth} [GeV]", "",
               true, text, text2, 0.18, 0.82, 0.05,
               "figure/respmatrix_reweighted_calib_frac_jerup.png");
   

  f_out->cd();
  h_respmatrix_calib_dijet->Write();
  h_respmatrix_calib_dijet_jesdown->Write();
  h_respmatrix_calib_dijet_jesup->Write();
  h_respmatrix_calib_dijet_jerdown->Write();
  h_respmatrix_calib_dijet_jerup->Write();
  h_respmatrix_calib_frac->Write();
  h_respmatrix_calib_frac_jesdown->Write();
  h_respmatrix_calib_frac_jesup->Write();
  h_respmatrix_calib_frac_jerdown->Write();
  h_respmatrix_calib_frac_jerup->Write();

  h_reweight_calib_dijet->Write();
  h_reweight_calib_dijet_effdown->Write();
  h_reweight_calib_dijet_effup->Write();
  h_reweight_calib_dijet_jesdown->Write();
  h_reweight_calib_dijet_jesup->Write();
  h_reweight_calib_dijet_jerdown->Write();
  h_reweight_calib_dijet_jerup->Write();
  h_reweight_calib_frac->Write();
  h_reweight_calib_frac_effdown->Write();
  h_reweight_calib_frac_effup->Write();
  h_reweight_calib_frac_jesdown->Write();
  h_reweight_calib_frac_jesup->Write();
  h_reweight_calib_frac_jerdown->Write();
  h_reweight_calib_frac_jerup->Write();

  h_respmatrix_calib_dijet_reweighted->Write();
  h_respmatrix_calib_dijet_effdown_reweighted->Write();
  h_respmatrix_calib_dijet_effup_reweighted->Write();
  h_respmatrix_calib_dijet_jesdown_reweighted->Write();
  h_respmatrix_calib_dijet_jesup_reweighted->Write();
  h_respmatrix_calib_dijet_jerdown_reweighted->Write();
  h_respmatrix_calib_dijet_jerup_reweighted->Write();
  h_respmatrix_calib_frac_reweighted->Write();
  h_respmatrix_calib_frac_effdown_reweighted->Write();
  h_respmatrix_calib_frac_effup_reweighted->Write();
  h_respmatrix_calib_frac_jesdown_reweighted->Write();
  h_respmatrix_calib_frac_jesup_reweighted->Write();
  h_respmatrix_calib_frac_jerdown_reweighted->Write();
  h_respmatrix_calib_frac_jerup_reweighted->Write();
  f_out->Close();
}