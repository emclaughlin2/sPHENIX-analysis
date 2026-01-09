#include "dlUtility.h"

double background(double* x, double* par)
{
  return par[0];
}

double gaussianpeak(double* x, double* par)
{
  double expt = -pow((x[0]-par[1])/par[2],2);
  return par[0]*exp(expt);
}

double fitf(double* x, double* par)
{
  double num = background(x,par) + gaussianpeak(x,&par[1]);
  if(num<1e-20) num = 1e-20;
  return num;
}

double bg2(double* val, double* par)
{
  return par[0];
}

double twogaus(double* val, double* par)
{
  double x = val[0];
  double y = val[1];
  double a  = par[0];
  double mx = par[1];
  double my = par[2];
  double sx = par[3];
  double sy = par[4];

  double expt = -(pow((x-mx)/sx,2)+pow((y-my)/sy,2));
  return a*exp(expt);
}

double fitf2(double* val, double*par)
{
  double num = bg2(val,par)+twogaus(val,&par[1]);
  if(num<1e-20) num = 1e-20;
  return num;
}

int timing_eff_test()
{
  const int nbincheck = 5;

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetEndErrorSize(0);

  TFile* inf = TFile::Open(
    "analysis_data_ana509_output/output_efrac_and_time_bkg_cut.root","READ"
  );

  TH3D* h3pttdt = (TH3D*)inf->Get("hpttdtdat");

  // ---- SINGLE X BIN DEFINITION (CRITICAL CHANGE)
  const int xfirst = 1;
  const int xlast  = h3pttdt->GetXaxis()->GetNbins();

  TCanvas* c = new TCanvas("c","",1000,1000);
  c->cd();
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);

  TH1D* nomcuts[2][nbincheck];
  TH2D* nomcuts2[nbincheck];

  TF1* fitsnom1[2][nbincheck];
  TF2* fitsnom2[nbincheck];

  TH1D* effs[3][3];

  float ltc[2] = {-8,4};
  float ltv[2] = {-9,5};
  float dtc[2] = {-5,1};
  float dtv[2] = {-6,2};

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      effs[i][j] = new TH1D(
        Form("eff_%d_%d",i,j),
        ";Uncalibrated Lead Jet p_{T} [GeV];Efficiency",
        9,10,100
      );

  // ===================== MAIN LOOP =====================
  for(int i=0;i<2;i++)
  {
    for(int j=0;j<nbincheck;j++)
    {
      // ----------- i == 0 : Y projection -------------
      if(i==0)
      {
        nomcuts[i][j] = h3pttdt->ProjectionY(
          Form("projY_%d_%d",i,j),
          xfirst, xlast,
          1, h3pttdt->GetZaxis()->GetNbins(),
          "e"
        );

        nomcuts2[j] = (TH2D*)h3pttdt->Project3D("zy");

        fitsnom1[i][j] = new TF1(Form("f1_%d_%d",i,j),fitf,-30,30,4);
        fitsnom1[i][j]->SetParameters(
          nomcuts[i][j]->GetBinContent(22),
          nomcuts[i][j]->GetMaximum(),0,
          nomcuts[i][j]->GetStdDev()
        );

        nomcuts[i][j]->Fit(fitsnom1[i][j],"LI");

        fitsnom2[j] = new TF2(Form("f2_%d",j),fitf2,-30,30,-30,30,6);
        fitsnom2[j]->SetParameters(
          nomcuts2[j]->GetMaximum(),0,0,0,
          nomcuts2[j]->GetStdDev(1),
          nomcuts2[j]->GetStdDev(2)
        );

        nomcuts2[j]->Fit(fitsnom2[j],"LI");

        c->Clear();
        nomcuts[i][j]->Draw("PE");
        fitsnom1[i][j]->Draw("SAME");

        c->SaveAs(Form("eff_plots/timeff_fit_t_%d_%d.pdf",i,j));
      }

      // ----------- i == 1 : Z projection -------------
      if(i==1)
      {
        nomcuts[i][j] = h3pttdt->ProjectionZ(
          Form("projZ_%d_%d",i,j),
          xfirst, xlast,
          1, h3pttdt->GetYaxis()->GetNbins(),
          "e"
        );

        fitsnom1[i][j] = new TF1(Form("f1_%d_%d",i,j),fitf,-30,30,4);
        fitsnom1[i][j]->SetParameters(
          nomcuts[i][j]->GetBinContent(27),
          nomcuts[i][j]->GetMaximum(),0,
          nomcuts[i][j]->GetRMS()
        );

        nomcuts[i][j]->Fit(fitsnom1[i][j],"LI");

        c->Clear();
        nomcuts[i][j]->Draw("PE");
        fitsnom1[i][j]->Draw("SAME");

        c->SaveAs(Form("eff_plots/timeff_fit_dt_%d_%d.pdf",i,j));
      }
    }
  }

  // ===================== OUTPUT =====================
  TFile* outf = TFile::Open("effs_timing.root","RECREATE");
  outf->cd();

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      effs[i][j]->Write();

  outf->Close();
  inf->Close();

  return 0;
}
