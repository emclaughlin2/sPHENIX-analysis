#include <TFile.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TF1.h>

void do_analysis(TH1 *hPass, TH1 *hTot, TFile *f_out_record, TFile *f_out, float fitmin, float fitmax, string prefix) {
  TEfficiency* eff = new TEfficiency(*hPass, *hTot);
  eff->SetName(("eff" + prefix).c_str());
  TBinomialEfficiencyFitter* bef = new TBinomialEfficiencyFitter(hPass, hTot);
  //TF1* fit = new TF1("fit", "[0] + [1]/pow(1+exp(-[3]*(x-[2])),[4])", 0, 100);
  //fit->SetParameters(0.0022, 0.965, 7.36, 0.467, 5.801);
  TF1* fit = new TF1("fit", "[0] + [1]/pow(1+exp(-[3]*(x-[2])),[4]) + [5]*pow(x,[6])", 0, 100);
  fit->SetParameters(0.0022, 0.965, 7.36, 0.467, 5.801, 0.0019, 0.8);
  fit->SetRange(fitmin, fitmax);
  TFitResultPtr fr = bef->Fit(fit, "RS");

  constexpr int N = 100;
  double xs[N], cis[N];
  double xmin = fitmin, xmax = fitmax;
  for(int i = 0; i < N; ++i) {
    xs[i] = xmin + (xmax - xmin) * i/(N - 1);
  }
  fr->GetConfidenceIntervals(N, 1, 1, xs, cis, 0.6827);  
  TGraphErrors* grCI = new TGraphErrors(N);
  for(int i = 0; i < N; ++i) {
    double x = xs[i];
    double y = fit->Eval(x);
    double dy = cis[i];
    grCI->SetPoint(i, x, y);
    grCI->SetPointError(i, 0, dy);
  }
  TGraph* gUpper = new TGraph(N);
  TGraph* gLower = new TGraph(N);
  for(int i = 0; i < N; ++i) {
    double x, y, dy;
    grCI->GetPoint(i, x, y);
    dy = grCI->GetErrorY(i);
    gUpper->SetPoint(i, x, y + dy);
    gLower->SetPoint(i, x, y - dy);
  }

  TF1* fit_up = (TF1*)fit->Clone(("jettrig" + prefix + "_up").c_str());
  TF1* fit_down = (TF1*)fit->Clone(("jettrig" + prefix + "_down").c_str());
  fit_up->SetName(("jettrig" + prefix + "_up").c_str());
  fit_down->SetName(("jettrig" + prefix + "_down").c_str());
  gUpper->Fit(fit_up, "R");
  gLower->Fit(fit_down, "R");
  fit_up->SetRange(0,100);
  fit_down->SetRange(0,100);
  fit->SetName(("jettrig" + prefix + "_nominal").c_str());
  fit->SetRange(0,100);

  f_out_record->cd();
  eff->Write();
  fit->Write();
  gUpper->Write(("g" + prefix + "_up").c_str());
  gLower->Write(("g" + prefix + "_down").c_str());
  fit_up->Write();
  fit_down->Write();
  f_out->cd();
  fit->Write();
  fit_up->Write();
  fit_down->Write();
}

void ppg09_jettrig_efficiency() {
  TFile *f_in = new TFile("test_analysis_output/trig12_zvtx_lt_10_output.root", "READ");
  TFile *f_out_record = new TFile("ppg09_trig12_zvtx_lt_10_output.root", "RECREATE");
  TFile *f_out = new TFile("ppg09_trig12_zvtx_lt_10_jetefficiency.root", "RECREATE");

  TH1F *h_leadingjet04pt_all_12 = (TH1F*)f_in->Get("h_leadingjet04pt_all_12"); h_leadingjet04pt_all_12->Rebin(10);
  TH1F *h_leadingjet04pt_all_34 = (TH1F*)f_in->Get("h_leadingjet04pt_all_34"); h_leadingjet04pt_all_34->Rebin(10);
  TH1F *h_leadingjet04pt_all_22 = (TH1F*)f_in->Get("h_leadingjet04pt_all_22"); h_leadingjet04pt_all_22->Rebin(10);

  do_analysis(h_leadingjet04pt_all_34, h_leadingjet04pt_all_12, f_out_record, f_out, 8, 30, "_04_pt_trig34");
  do_analysis(h_leadingjet04pt_all_22, h_leadingjet04pt_all_12, f_out_record, f_out, 8, 30, "_04_pt_trig22");

  f_in->Close();
  f_out_record->Close();
  f_out->Close();
}
