#include <iostream>
#include <string>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

void draw_one(
    const std::string &infile,
    const std::string &hist_name,
    const std::string &fit_name,
    const std::string &canvas_name,
    const std::string &particle_label,
    const std::string &out_prefix)
{
  TFile *f = TFile::Open(infile.c_str(), "READ");
  if (!f || f->IsZombie())
  {
    std::cerr << "ERROR: cannot open " << infile << std::endl;
    return;
  }

  // Get latest cycle by requesting object name directly.
  TH1 *h = nullptr;
  TF1 *fit = nullptr;
  f->GetObject(hist_name.c_str(), h);
  f->GetObject(fit_name.c_str(), fit);

  if (!h)
  {
    std::cerr << "ERROR: missing histogram " << hist_name
              << " in " << infile << std::endl;
    f->Close();
    return;
  }
  if (!fit)
  {
    std::cerr << "WARNING: missing fit " << fit_name
              << " in " << infile << std::endl;
  }

  // Detach from file so objects survive file close.
  TH1 *hclone = (TH1 *) h->Clone((hist_name + "_clone").c_str());
  hclone->SetDirectory(nullptr);

  TF1 *fclone = nullptr;
  if (fit)
  {
    fclone = (TF1 *) fit->Clone((fit_name + "_" + particle_label + "_clone").c_str());
  }
  f->Close();

  TCanvas *c = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 900, 700);
  c->cd();

  // Histogram style tuned for sPHENIX-like look.
  hclone->SetLineColor(kBlack);
  hclone->SetLineWidth(2);
  hclone->SetMarkerStyle(20);
  hclone->SetMarkerSize(1.0);
  hclone->SetMarkerColor(kBlack);
  hclone->SetStats(0);

  hclone->GetXaxis()->SetTitle("E/p");
  hclone->GetYaxis()->SetTitle("Counts");
  hclone->GetXaxis()->SetTitleOffset(1.1);
  if (particle_label == "Isolated EM") hclone->GetYaxis()->SetTitleOffset(1.3);
  else hclone->GetYaxis()->SetTitleOffset(1.2);
  if (particle_label == "Isolated hadron") hclone->GetXaxis()->SetRangeUser(-0.5, 1.5);
  else hclone->GetXaxis()->SetRangeUser(0.0, 2.0);
  if (particle_label == "Isolated hadron") hclone->GetYaxis()->SetRangeUser(0.0, 450);
  else hclone->GetYaxis()->SetRangeUser(0.0, 5000);

  hclone->Draw("E1");

  if (fclone)
  {
    //fclone->SetLineColor(kRed + 1);
    fclone->SetLineColor(kBlue + 1);
    fclone->SetLineWidth(2);
    fclone->SetLineStyle(1);
    hclone->GetListOfFunctions()->Add(fclone);  // keep fit with histogram
    fclone->Draw("SAME");
  }

  // Legend style (minimal box, sPHENIX-friendly).
  TLegend *leg = new TLegend(0.56, 0.76, 0.86, 0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  if (particle_label == "Isolated hadron") leg->AddEntry(hclone, "3GeV Isolated #pi^{+/-}", "lep");
  else leg->AddEntry(hclone, "1GeV Isolated e^{+/-} + #gamma", "lep");
  leg->AddEntry("", "Matched Topoclusters", "");
  if (fclone) leg->AddEntry(fclone, "Gaussian Fit", "l");
  leg->Draw();

  // sPHENIX label style.
  TLatex lat;
  lat.SetNDC();
  lat.SetTextFont(42);
  lat.SetTextSize(0.04);
  lat.DrawLatex(0.2, 0.87, "#it{#bf{sPHENIX}} Simulation Internal");
  lat.DrawLatex(0.2, 0.82, "200 GeV Pythia8 p+p");

  c->Modified();
  c->Update();

  c->SaveAs((out_prefix + ".pdf").c_str());
  c->SaveAs((out_prefix + ".png").c_str());
}

void plot_isolated_resolution_sphenix()
{
  // Load sPHENIX style exactly as requested.
  gROOT->LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
  gROOT->ProcessLine("SetsPhenixStyle()");

  // Global style touch-ups (optional but common).
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // Hadron from Hadron_Isolated_Resolution_3Plots_fit.root
  draw_one(
      "Hadron_Isolated_Resolution_3Plots_fit.root",
      "hEP_hadron_3s",
      "gausFit_3s",
      "c_hadron_3s_sphenix",
      "Isolated hadron",
      "Hadron_Isolated_Resolution_3Plots_sphenix");

  // EM from EM_Isolated_Resolution_3Plots_fit.root
  draw_one(
      "EM_Isolated_Resolution_3Plots_fit.root",
      "hEP_em_3s",
      "gausFit_3s",
      "c_em_3s_sphenix",
      "Isolated EM",
      "EM_Isolated_Resolution_3Plots_sphenix");
}