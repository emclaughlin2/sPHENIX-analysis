#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TKey.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TError.h"

namespace {

TH1D* BuildEfficiency(const TH1* pass, const TH1* total, const char* name)
{
  if (!pass || !total) return nullptr;
  TH1D* eff = dynamic_cast<TH1D*>(pass->Clone(name));
  if (!eff) return nullptr;
  eff->SetDirectory(nullptr);
  eff->Reset("ICESM");
  eff->Divide(pass, total, 1.0, 1.0, "B");
  return eff;
}

TH1D* BuildRatio(const TH1* numerator, const TH1* denominator, const char* name)
{
  if (!numerator || !denominator) return nullptr;
  TH1D* ratio = dynamic_cast<TH1D*>(numerator->Clone(name));
  if (!ratio) return nullptr;
  ratio->SetDirectory(nullptr);
  ratio->Divide(denominator);
  return ratio;
}

void DrawEfficiencyWithRatio(TH1D* effData, TH1D* effMC, const std::string& outName, const std::string& xTitle)
{
  if (!effData || !effMC) return;

  TH1D* ratio = BuildRatio(effData, effMC, (outName + "_ratio").c_str());
  if (!ratio) return;

  TCanvas* c = new TCanvas((outName + "_c").c_str(), "", 900, 900);
  c->cd();

  TPad* pTop = new TPad((outName + "_top").c_str(), "", 0.0, 0.30, 1.0, 1.0);
  TPad* pBot = new TPad((outName + "_bot").c_str(), "", 0.0, 0.0, 1.0, 0.30);
  pTop->SetBottomMargin(0.02);
  pBot->SetTopMargin(0.05);
  pBot->SetBottomMargin(0.28);
  pTop->Draw();
  pBot->Draw();

  pTop->cd();
  effData->SetMarkerStyle(20);
  effData->SetMarkerColor(kBlack);
  effData->SetLineColor(kBlack);
  effMC->SetMarkerStyle(24);
  effMC->SetMarkerColor(kRed + 1);
  effMC->SetLineColor(kRed + 1);

  effData->GetYaxis()->SetTitle("Efficiency");
  effData->GetYaxis()->SetRangeUser(0.0, 1.15);
  effData->GetXaxis()->SetLabelSize(0.0);
  effData->Draw("E1");
  effMC->Draw("E1 SAME");

  TLegend* leg = new TLegend(0.62, 0.72, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(effData, "Data", "lep");
  leg->AddEntry(effMC, "MC", "lep");
  leg->Draw();

  pBot->cd();
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(kBlue + 1);
  ratio->SetLineColor(kBlue + 1);
  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetTitleSize(0.09);
  ratio->GetYaxis()->SetTitleOffset(0.45);
  ratio->GetYaxis()->SetLabelSize(0.08);
  ratio->GetYaxis()->SetRangeUser(0.8, 1.2);
  ratio->GetXaxis()->SetTitle(xTitle.c_str());
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.05);
  ratio->GetXaxis()->SetLabelSize(0.085);
  ratio->Draw("E1");

  TLine* line1 = new TLine(ratio->GetXaxis()->GetXmin(), 1.0, ratio->GetXaxis()->GetXmax(), 1.0);
  line1->SetLineStyle(2);
  line1->Draw("SAME");

  c->SaveAs((outName + ".pdf").c_str());
  c->SaveAs((outName + ".png").c_str());
}

void DrawSimpleComparison(TH1* hData, TH1* hMC, const std::string& outName)
{
  if (!hData || !hMC) return;

  if (hData->InheritsFrom(TH2::Class())) {
    TCanvas* c2 = new TCanvas((outName + "_c2").c_str(), "", 1200, 500);
    c2->Divide(2, 1);
    c2->cd(1);
    hData->SetTitle((std::string(hData->GetName()) + " (Data)").c_str());
    hData->Draw("COLZ");
    c2->cd(2);
    hMC->SetTitle((std::string(hMC->GetName()) + " (MC)").c_str());
    hMC->Draw("COLZ");
    c2->SaveAs((outName + ".pdf").c_str());
    c2->SaveAs((outName + ".png").c_str());
    return;
  }

  TH1* hDataDraw = dynamic_cast<TH1*>(hData->Clone((std::string(hData->GetName()) + "_data_draw").c_str()));
  TH1* hMCDraw = dynamic_cast<TH1*>(hMC->Clone((std::string(hMC->GetName()) + "_mc_draw").c_str()));
  if (!hDataDraw || !hMCDraw) return;
  hDataDraw->SetDirectory(nullptr);
  hMCDraw->SetDirectory(nullptr);
  hDataDraw->Sumw2();
  hMCDraw->Sumw2();


  const double intData = hDataDraw->Integral();
  const double intMC = hMCDraw->Integral();
  if (intData > 0.0) hDataDraw->Scale(1.0 / intData);
  if (intMC > 0.0) hMCDraw->Scale(1.0 / intMC);
  if (!strcmp(hData->GetName(),"h_lead_spectra_record")) {
    hDataDraw->Rebin(10);
    hMCDraw->Rebin(10);
    hDataDraw->Scale(1.0 / hDataDraw->GetBinContent(hDataDraw->FindBin(21.5)));
    hMCDraw->Scale(1.0 / hMCDraw->GetBinContent(hMCDraw->FindBin(21.5)));
  }

  TH1* ratio = dynamic_cast<TH1*>(hDataDraw->Clone((std::string(hData->GetName()) + "_ratio").c_str()));
  ratio->SetDirectory(nullptr);
  ratio->Divide(hMCDraw);

  TCanvas* c = new TCanvas((outName + "_c").c_str(), "", 900, 900);
  c->cd();
  TPad* pTop = new TPad((outName + "_top").c_str(), "", 0.0, 0.30, 1.0, 1.0);
  TPad* pBot = new TPad((outName + "_bot").c_str(), "", 0.0, 0.0, 1.0, 0.30);
  pTop->SetBottomMargin(0.02);
  pBot->SetTopMargin(0.05);
  pBot->SetBottomMargin(0.28);
  pTop->Draw();
  pBot->Draw();

  pTop->cd();
  hDataDraw->SetMarkerStyle(20);
  hDataDraw->SetMarkerColor(kBlack);
  hDataDraw->SetLineColor(kBlack);
  hMCDraw->SetMarkerStyle(24);
  hMCDraw->SetMarkerColor(kRed + 1);
  hMCDraw->SetLineColor(kRed + 1);
  hDataDraw->SetTitle(hData->GetName());
  hDataDraw->GetYaxis()->SetTitle("Normalized entries");
  hDataDraw->GetXaxis()->SetLabelSize(0.0);
  hDataDraw->Draw("E1");
  hMCDraw->Draw("E1 SAME");

  TLegend* leg = new TLegend(0.62, 0.72, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hDataDraw, "Data", "lep");
  leg->AddEntry(hMCDraw, "MC", "lep");
  leg->Draw();
  if (!strcmp(hData->GetName(),"h_lead_spectra_record")) { pTop->SetLogy(1); }

  pBot->cd();
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(kBlue + 1);
  ratio->SetLineColor(kBlue + 1);
  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetTitleSize(0.09);
  ratio->GetYaxis()->SetTitleOffset(0.45);
  ratio->GetYaxis()->SetLabelSize(0.08);
  ratio->GetYaxis()->SetRangeUser(0.0, 2.0);
  ratio->GetXaxis()->SetTitle(hData->GetXaxis()->GetTitle());
  ratio->GetXaxis()->SetTitleSize(0.10);
  ratio->GetXaxis()->SetTitleOffset(1.05);
  ratio->GetXaxis()->SetLabelSize(0.085);
  ratio->Draw("E1");

  TLine* line1 = new TLine(ratio->GetXaxis()->GetXmin(), 1.0, ratio->GetXaxis()->GetXmax(), 1.0);
  line1->SetLineStyle(2);
  line1->Draw("SAME");

  c->SaveAs((outName + ".pdf").c_str());
  c->SaveAs((outName + ".png").c_str());
}

}  // namespace

void plot_mbd_vertex_eff_compare(const std::string& mcFileName = "analysis_syst_output/output_sim_zvtx_lt_60.000000_all_run28.root",
                                 const std::string& dataFileName = "analysis_syst_output/output_zvtx_lt_60.000000_all_51274.root",
                                 const std::string& outDir = "comparison_plots",
                                 const std::string& ratioOutFileName = "efficiency_ratio_scale_factors.root")
{
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  gSystem->mkdir(outDir.c_str(), true);

  TFile* fMC = TFile::Open(mcFileName.c_str(), "READ");
  TFile* fData = TFile::Open(dataFileName.c_str(), "READ");
  if (!fMC || fMC->IsZombie()) {
    std::cerr << "Failed to open MC file: " << mcFileName << std::endl;
    return;
  }
  if (!fData || fData->IsZombie()) {
    std::cerr << "Failed to open data file: " << dataFileName << std::endl;
    return;
  }

  TH2* hPassData = dynamic_cast<TH2*>(fData->Get("h_pass_cut_measure_mbdvertex"));
  TH2* hTotData = dynamic_cast<TH2*>(fData->Get("h_total_measure_mbdvertex"));
  TH2* hPassMC = dynamic_cast<TH2*>(fMC->Get("h_pass_cut_measure_mbdvertex"));
  TH2* hTotMC = dynamic_cast<TH2*>(fMC->Get("h_total_measure_mbdvertex"));

  TFile* fRatioOut = TFile::Open(ratioOutFileName.c_str(), "RECREATE");
  if (!fRatioOut || fRatioOut->IsZombie()) {
    std::cerr << "Failed to create ratio output file: " << ratioOutFileName << std::endl;
    return;
  }

  if (hPassData && hTotData && hPassMC && hTotMC) {
    TH1D* hPassDataX = hPassData->ProjectionX("hPassDataX");
    TH1D* hTotDataX = hTotData->ProjectionX("hTotDataX");
    TH1D* hPassMCX = hPassMC->ProjectionX("hPassMCX");
    TH1D* hTotMCX = hTotMC->ProjectionX("hTotMCX");

    TH1D* hEffDataX = BuildEfficiency(hPassDataX, hTotDataX, "hEffDataX");
    TH1D* hEffMCX = BuildEfficiency(hPassMCX, hTotMCX, "hEffMCX");
    TH1D* hEffRatioX = BuildRatio(hEffDataX, hEffMCX, "h_eff_ratio_jet_pt_data_over_mc");
    if (hEffRatioX) {
      hEffRatioX->SetTitle("Efficiency ratio Data/MC vs jet p_{T}");
      hEffRatioX->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
      hEffRatioX->GetYaxis()->SetTitle("Data / MC");
    }
    DrawEfficiencyWithRatio(hEffDataX, hEffMCX, outDir + "/efficiency_vs_jet_pt", "Jet p_{T} [GeV]");

    TH1D* hPassDataY = hPassData->ProjectionY("hPassDataY");
    TH1D* hTotDataY = hTotData->ProjectionY("hTotDataY");
    TH1D* hPassMCY = hPassMC->ProjectionY("hPassMCY");
    TH1D* hTotMCY = hTotMC->ProjectionY("hTotMCY");

    TH1D* hEffDataY = BuildEfficiency(hPassDataY, hTotDataY, "hEffDataY");
    TH1D* hEffMCY = BuildEfficiency(hPassMCY, hTotMCY, "hEffMCY");
    TH1D* hEffRatioY = BuildRatio(hEffDataY, hEffMCY, "h_eff_ratio_calo_et_data_over_mc");
    if (hEffRatioY) {
      hEffRatioY->SetTitle("Efficiency ratio Data/MC vs calo E_{T}");
      hEffRatioY->GetXaxis()->SetTitle("Calo E_{T} [GeV]");
      hEffRatioY->GetYaxis()->SetTitle("Data / MC");
    }
    DrawEfficiencyWithRatio(hEffDataY, hEffMCY, outDir + "/efficiency_vs_calo_et", "Calo E_{T} [GeV]");

    fRatioOut->cd();
    if (hEffDataX) hEffDataX->Write("h_eff_data_jet_pt");
    if (hEffMCX) hEffMCX->Write("h_eff_mc_jet_pt");
    if (hEffRatioX) hEffRatioX->Write();
    if (hEffDataY) hEffDataY->Write("h_eff_data_calo_et");
    if (hEffMCY) hEffMCY->Write("h_eff_mc_calo_et");
    if (hEffRatioY) hEffRatioY->Write();
  } else {
    std::cerr << "Could not find one or more required 2D histograms for efficiency plots." << std::endl;
  }

  TIter nextKey(fData->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = static_cast<TKey*>(nextKey()))) {
    TObject* objData = key->ReadObj();
    if (!objData || !objData->InheritsFrom(TH1::Class())) continue;

    const std::string histName = objData->GetName();
    TObject* objMC = fMC->Get(histName.c_str());
    if (!objMC || !objMC->InheritsFrom(TH1::Class())) continue;

    TH1* hData = dynamic_cast<TH1*>(objData);
    TH1* hMC = dynamic_cast<TH1*>(objMC);
    if (!hData || !hMC) continue;

    DrawSimpleComparison(hData, hMC, outDir + "/" + histName + "_data_vs_mc");
  }

  fRatioOut->Close();
  std::cout << "Saved plots in: " << outDir << std::endl;
  std::cout << "Saved efficiency ratio histograms to: " << ratioOutFileName << std::endl;
  fData->Close();
  fMC->Close();
}
