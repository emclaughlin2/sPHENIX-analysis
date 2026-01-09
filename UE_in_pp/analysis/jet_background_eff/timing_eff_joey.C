#include "../dlUtility.h"

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
  double a = par[0];
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

int timing_eff()
{
  const int nbincheck = 5;
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  cout << "test" << endl;
  TFile* inf = TFile::Open("../hists_mbdtimereq_out_datsam_slewed20251211.root","READ");
  TCanvas* c = new TCanvas("","",1000,1000);
  c->cd();
  gPad->SetTopMargin(0.05);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.15);
  gStyle->SetEndErrorSize(0);
  TH3D* h3pttdt = (TH3D*)inf->Get("hpttdtdat");

  TH1D* nomcuts[2][nbincheck];
  TF1* fitsnom1[2][nbincheck];
  TF2* fitsnom2[nbincheck];

  TH2D* nomcuts2[nbincheck];

  TH1D* effs[3][3];

  float ltc[2] = {-6,6};
  float ltv[2] = {-7,7};
  float dtc[2] = {-3,3};
  float dtv[2] = {-4,4};

  h3pttdt->Rebin3D(10,10,10);
  
  for(int i=0; i<3; ++i)
    {
      for(int j=0; j<3; ++j)
	{
	  effs[i][j] = new TH1D(("eff"+to_string(i)+"_"+to_string(j)).c_str(),";Uncalibrated Lead Jet p_{T} [GeV];Efficiency",9,10,100);

	}
    }

  for(int i=0; i<2; ++i)
    {
      for(int j=0; j<nbincheck; ++j)
	{
	  cout << i << " " << j << endl;
	  if(i==0)
	    {
	      nomcuts[i][j] = h3pttdt->ProjectionY((string(h3pttdt->GetName()+to_string(i)+"_"+to_string(j)+"_nomcuts")).c_str(),(j+1)*1+1,(j+2)*1,1,h3pttdt->GetZaxis()->GetNbins(),"e");
	      h3pttdt->GetXaxis()->SetRange((j+1)*1+1,(j+2)*1);
	      nomcuts2[j] = (TH2D*)h3pttdt->Project3D("zy");
	      fitsnom1[i][j] = new TF1(("f1nom"+to_string(i)+"_"+to_string(j)).c_str(),fitf,-30,30,4);
	      fitsnom1[i][j]->SetParLimits(3,0.00001,999999);
	      fitsnom1[i][j]->SetParLimits(0,1e-30,999999);
	      fitsnom1[i][j]->SetParLimits(1,0.00001,99999999);
	      fitsnom1[i][j]->SetParameters(nomcuts[i][j]->GetBinContent(22),nomcuts[i][j]->GetMaximum(),0,nomcuts[i][j]->GetStdDev());
	      nomcuts[i][j]->Fit(fitsnom1[i][j],"LI");
	      
	      fitsnom2[j] = new TF2(("f2nom"+to_string(j)).c_str(),fitf2,-30,30,-30,30,6);
	      fitsnom2[j]->SetParLimits(0,1e-30,999999);
	      fitsnom2[j]->SetParLimits(1,0.00001,99999999);
	      fitsnom2[j]->SetParLimits(4,0.00001,999999);
	      fitsnom2[j]->SetParLimits(5,0.00001,999999);
	      fitsnom2[j]->SetParameters(nomcuts2[j]->GetBinContent(22,27)>0?nomcuts2[j]->GetBinContent(22,27):1,nomcuts2[j]->GetMaximum(),0,0,nomcuts2[j]->GetStdDev(1),nomcuts2[j]->GetStdDev(2));
	      nomcuts2[j]->Fit(fitsnom2[j],"LI");

	      
	      TF1* usegaus = new TF1("usegaus",gaussianpeak,-30,30,3);
	      TF1* higaus = new TF1("higaus",gaussianpeak,-30,30,3);
	      TF1* logaus = new TF1("logaus",gaussianpeak,-30,30,3);
	      TF1* p0 = new TF1("p0","pol0",-30,30);
	      TF1* thefit = nomcuts[i][j]->GetFunction(("f1nom"+to_string(i)+"_"+to_string(j)).c_str());
	      usegaus->SetParameters(thefit->GetParameter(1),thefit->GetParameter(2),thefit->GetParameter(3));
	      higaus->SetParameters(thefit->GetParameter(1)+thefit->GetParError(1),thefit->GetParameter(2)+thefit->GetParError(2),thefit->GetParameter(3)+thefit->GetParError(3));
	      logaus->SetParameters(thefit->GetParameter(1)-thefit->GetParError(1),thefit->GetParameter(2)+thefit->GetParError(2),thefit->GetParameter(3)+thefit->GetParError(3));
	      double errshi[3] = {higaus->Integral(ltc[0],ltc[1])/higaus->Integral(-30,30)-usegaus->Integral(ltc[0],ltc[1])/usegaus->Integral(-30,30),higaus->Integral(ltv[0],ltv[1])/higaus->Integral(-30,30)-usegaus->Integral(ltv[0],ltv[1])/usegaus->Integral(-30,30),higaus->Integral(-7,3)/higaus->Integral(-30,30)-usegaus->Integral(-7,3)/usegaus->Integral(-30,30)};
	      double errslo[3] = {logaus->Integral(ltc[0],ltc[1])/logaus->Integral(-30,30)-usegaus->Integral(ltc[0],ltc[1])/usegaus->Integral(-30,30),logaus->Integral(ltv[0],ltv[1])/logaus->Integral(-30,30)-usegaus->Integral(ltv[0],ltv[1])/usegaus->Integral(-30,30),logaus->Integral(-7,3)/logaus->Integral(-30,30)-usegaus->Integral(-7,3)/usegaus->Integral(-30,30)};
	      effs[0][0]->SetBinContent(j+1,usegaus->Integral(ltc[0],ltc[1])/usegaus->Integral(-30,30));
	      effs[0][1]->SetBinContent(j+1,usegaus->Integral(ltv[0],ltv[1])/usegaus->Integral(-30,30));
	      effs[0][2]->SetBinContent(j+1,usegaus->Integral(-7,3)/usegaus->Integral(-30,30));
	      double errs[3] = {(abs(errshi[0])+abs(errslo[0]))/2,(abs(errshi[1])+abs(errslo[1]))/2,(abs(errshi[2])+abs(errslo[2]))/2};
	      for(int k=0; k<3; ++k)
		{
		  if(effs[0][k]->GetBinContent(j+1)+errs[k]>1)
		    {
		      errs[k] = 1-effs[0][k]->GetBinContent(j+1);
		    }
		}
	      effs[0][0]->SetBinError(j+1,errs[0]);
	      effs[0][1]->SetBinError(j+1,errs[1]);
	      effs[0][2]->SetBinError(j+1,errs[2]);

	      nomcuts[i][j]->SetMarkerColor(kBlack);
	      nomcuts[i][j]->SetMarkerSize(1);
	      nomcuts[i][j]->SetMarkerStyle(20);
	      nomcuts[i][j]->SetLineColor(kBlack);
	      nomcuts[i][j]->GetYaxis()->SetRangeUser(0,nomcuts[i][j]->GetMaximum()*1.1);
	      nomcuts[i][j]->GetYaxis()->SetTitle("Counts");
	      nomcuts[i][j]->Draw("PE");
	      thefit->SetLineColor(kViolet);
	      p0->SetParameters(thefit->GetParameter(0));
	      //p0->SetLineColor(kBlue);
	      usegaus->SetLineColor(kBlue);

	      TLegend* theleg = new TLegend(0.7,0.75,0.91,0.9);
	      theleg->AddEntry(thefit,"Fit signal+bg","l");
	      theleg->AddEntry(usegaus,"Signal","l");
	      thefit->SetNpx(600);
	      usegaus->SetNpx(600);
	      const int N = 1200;   // resolution
	      double xmin = thefit->GetXmin();
	      double xmax = thefit->GetXmax();

	      // Build the filled polygon
	      TGraph* g = new TGraph(2*N);
	      
	      for (int i = 0; i < N; i++) {
		double x = xmin + (xmax - xmin) * i / (N - 1);
		g->SetPoint(i, x, thefit->Eval(x));         // upper curve (left → right)
	      }
	      for (int i = 0; i < N; i++) {
		double x = xmax - (xmax - xmin) * i / (N - 1);
		g->SetPoint(N + i, x, p0->Eval(x));     // lower curve (right → left)
	      }

	      g->SetFillColorAlpha(kAzure+2, 0.35);
	      g->SetLineWidth(0);
	      g->Draw("F");

	      TBox* thebox = new TBox(-30,0,30,p0->GetParameter(0));

	      thebox->SetFillColorAlpha(kRed,0.35);
	      thebox->SetLineWidth(0);
	      thebox->Draw();
	      
	      thefit->Draw("SAME");
	      usegaus->Draw("SAME");
	
	      theleg->SetFillStyle(0);
	      theleg->SetFillColor(0);
	      theleg->SetBorderSize(0);
	      TLine* lines[2] = {new TLine(ltc[0],0,ltc[0],nomcuts[i][j]->GetMaximum()),new TLine(ltc[1],0,ltc[1],nomcuts[i][j]->GetMaximum())};
	      for(int k=0; k<2; ++k)
		{
		  lines[k]->SetLineColor(kGreen);
		  lines[k]->SetLineStyle(9);
		  lines[k]->SetLineWidth(2);
		  lines[k]->Draw();
		}
	      theleg->AddEntry(lines[0],"Cut Bounds","l");
	      theleg->Draw();

	      

	      maintexts(0.91,0.17,0,0.03,1,0);
	      drawText("Background included",0.17,0.8,0,kBlack,0.03);
	      drawText("No reconstructed",0.17,0.76,0,kBlack,0.03);
	      drawText("z_{vtx} requirement",0.17,0.72,0,kBlack,0.03);
	      drawText("Dijet cut applied",0.17,0.68,0,kBlack,0.03);

	      c->SaveAs(("../../images/effs/timeff_fit_t_"+to_string(i)+"_"+to_string(j)+".pdf").c_str());
	      
	      delete usegaus;
	      delete higaus;
	      delete logaus;
	      
	      TF2* usegaus2 = new TF2("usegaus2",twogaus,-30,30,-30,30,5);
	      TF2* higaus2 = new TF2("higaus2",twogaus,-30,30,-30,30,5);
	      TF2* logaus2 = new TF2("logaus2",twogaus,-30,30,-30,30,5);
	      
	      TF2* fit2gaus = (TF2*)nomcuts2[j]->GetFunction(("f2nom"+to_string(j)).c_str());
	      for(int k=0; k<5; ++k)
		{
		  usegaus2->SetParameter(k,fit2gaus->GetParameter(k+1));
		  higaus2->SetParameter(k,fit2gaus->GetParameter(k+1)+fit2gaus->GetParError(k+1));
		  logaus2->SetParameter(k,fit2gaus->GetParameter(k+1)-fit2gaus->GetParError(k+1));
		}
	      cout << usegaus2->Integral(ltc[0],ltc[1],dtc[0],dtc[1]) << endl;
	      cout << usegaus2->Integral(-30,30,-30,30) << endl;

	      double errshi2[3] = {higaus2->Integral(ltc[0],ltc[1],dtc[0],dtc[1])/higaus2->Integral(-30,30,-30,30)-usegaus2->Integral(ltc[0],ltc[1],dtc[0],dtc[1])/usegaus2->Integral(-30,30,-30,30),higaus2->Integral(ltv[0],ltv[1],dtv[0],4)/higaus2->Integral(-30,30,-30,30)-usegaus2->Integral(ltv[0],ltv[1],dtv[0],4)/usegaus2->Integral(-30,30,-30,30),higaus2->Integral(-7,3,-2,2)/higaus2->Integral(-30,30,-30,30)-usegaus2->Integral(-7,3,-2,2)/usegaus2->Integral(-30,30,-30,30)};
	      double errslo2[3] = {logaus2->Integral(ltc[0],ltc[1],dtc[0],dtc[1])/logaus2->Integral(-30,30,-30,30)-usegaus2->Integral(ltc[0],ltc[1],dtc[0],dtc[1])/usegaus2->Integral(-30,30,-30,30),logaus2->Integral(ltv[0],ltv[1],dtv[0],4)/logaus2->Integral(-30,30,-30,30)-usegaus2->Integral(ltv[0],ltv[1],dtv[0],4)/usegaus2->Integral(-30,30,-30,30),logaus2->Integral(-7,3,-2,2)/logaus2->Integral(-30,30,-30,30)-usegaus2->Integral(-7,3,-2,2)/usegaus2->Integral(-30,30,-30,30)};
	      
	      
	      effs[2][0]->SetBinContent(j+1,usegaus2->Integral(ltc[0],ltc[1],dtc[0],dtc[1])/usegaus2->Integral(-30,30,-30,30));
	      effs[2][1]->SetBinContent(j+1,usegaus2->Integral(ltv[0],ltv[1],dtv[0],dtv[1])/usegaus2->Integral(-30,30,-30,30));
	      effs[2][2]->SetBinContent(j+1,usegaus2->Integral(-7,3,-2,2)/usegaus2->Integral(-30,30,-30,30));

	      double errs2[3] = {(abs(errshi2[0])+abs(errslo2[0]))/2,(abs(errshi2[1])+abs(errslo2[1]))/2,(abs(errshi2[2])+abs(errslo2[2]))/2};
	      for(int k=0; k<3; ++k)
		{
		  if(effs[2][k]->GetBinContent(j+1)+errs2[k]>1)
		    {
		      errs2[k] = 1-effs[2][k]->GetBinContent(j+1);
		    }
		}
	      effs[2][0]->SetBinError(j+1,errs2[0]);
	      effs[2][1]->SetBinError(j+1,errs2[1]);
	      effs[2][2]->SetBinError(j+1,errs2[2]);



	      
	      delete usegaus2;
	      delete logaus2;
	      delete higaus2;
	    }
	  if(i==1)
	    {
	      nomcuts[i][j] = h3pttdt->ProjectionZ((string(h3pttdt->GetName()+to_string(i)+"_"+to_string(j)+"_nomcuts")).c_str(),(j+1)*1+1,(j+2)*1,1,h3pttdt->GetYaxis()->GetNbins(),"e");
	      fitsnom1[i][j] = new TF1(("f1nom"+to_string(i)+"_"+to_string(j)).c_str(),fitf,-30,30,4);
	      fitsnom1[i][j]->SetParLimits(3,0.00001,999999);
	      fitsnom1[i][j]->SetParLimits(0,1e-30,999999);
	      fitsnom1[i][j]->SetParLimits(1,0.00001,99999999);
	      fitsnom1[i][j]->SetParameters(nomcuts[i][j]->GetBinContent(27),nomcuts[i][j]->GetMaximum(),0,nomcuts[i][j]->GetRMS());
	      nomcuts[i][j]->Fit(fitsnom1[i][j],"LI");

	      TF1* p0 = new TF1("p0","pol0",-30,30);
	      TF1* usegaus = new TF1("usegaus",gaussianpeak,-30,30,3);
	      TF1* higaus = new TF1("higaus",gaussianpeak,-30,30,3);
	      TF1* logaus = new TF1("logaus",gaussianpeak,-30,30,3);
	      TF1* thefit = nomcuts[i][j]->GetFunction(("f1nom"+to_string(i)+"_"+to_string(j)).c_str());
	      usegaus->SetParameters(thefit->GetParameter(1),thefit->GetParameter(2),thefit->GetParameter(3));
	      higaus->SetParameters(thefit->GetParameter(1)+thefit->GetParError(1),thefit->GetParameter(2)+thefit->GetParError(2),thefit->GetParameter(3)+thefit->GetParError(3));
	      logaus->SetParameters(thefit->GetParameter(1)-thefit->GetParError(1),thefit->GetParameter(2)+thefit->GetParError(2),thefit->GetParameter(3)+thefit->GetParError(3));
	      double errshi[3] = {higaus->Integral(dtc[0],dtc[1])/higaus->Integral(-30,30)-usegaus->Integral(dtc[0],dtc[1])/usegaus->Integral(-30,30),higaus->Integral(dtv[0],dtv[1])/higaus->Integral(-30,30)-usegaus->Integral(dtv[0],dtv[1])/usegaus->Integral(-30,30),higaus->Integral(-7,3)/higaus->Integral(-30,30)-usegaus->Integral(-7,3)/usegaus->Integral(-30,30)};
	      double errslo[3] = {logaus->Integral(dtc[0],dtc[1])/logaus->Integral(-30,30)-usegaus->Integral(dtc[0],dtc[1])/usegaus->Integral(-30,30),logaus->Integral(dtv[0],dtv[1])/logaus->Integral(-30,30)-usegaus->Integral(dtv[0],dtv[1])/usegaus->Integral(-30,30),logaus->Integral(-7,3)/logaus->Integral(-30,30)-usegaus->Integral(-7,3)/usegaus->Integral(-30,30)};
	      effs[1][0]->SetBinContent(j+1,usegaus->Integral(dtc[0],dtc[1])/usegaus->Integral(-30,30));
	      effs[1][1]->SetBinContent(j+1,usegaus->Integral(dtv[0],dtv[1])/usegaus->Integral(-30,30));
	      effs[1][2]->SetBinContent(j+1,usegaus->Integral(-7,3)/usegaus->Integral(-30,30));
	      double errs[3] = {(abs(errshi[0])+abs(errslo[0]))/2,(abs(errshi[1])+abs(errslo[1]))/2,(abs(errshi[2])+abs(errslo[2]))/2};
	      for(int k=0; k<3; ++k)
		{
		  if(effs[1][k]->GetBinContent(j+1)+errs[k]>1)
		    {
		      errs[k] = 1-effs[1][k]->GetBinContent(j+1);
		    }
		}
	      effs[1][0]->SetBinError(j+1,errs[0]);
	      effs[1][1]->SetBinError(j+1,errs[1]);
	      effs[1][2]->SetBinError(j+1,errs[2]);


	      nomcuts[i][j]->SetMarkerColor(kBlack);
	      nomcuts[i][j]->SetMarkerSize(1);
	      nomcuts[i][j]->SetMarkerStyle(20);
	      nomcuts[i][j]->SetLineColor(kBlack);
	      nomcuts[i][j]->GetYaxis()->SetRangeUser(0,nomcuts[i][j]->GetMaximum()*1.1);
	      nomcuts[i][j]->Draw("PE");
	      thefit->SetLineColor(kViolet);
	      p0->SetParameters(thefit->GetParameter(0));
	      //p0->SetLineColor(kBlue);
	      usegaus->SetLineColor(kBlue);

	      TLegend* theleg = new TLegend(0.7,0.75,0.91,0.9);
	      theleg->AddEntry(thefit,"Fit signal+bg","l");
	      theleg->AddEntry(usegaus,"Signal","l");
	      usegaus->SetNpx(600);
	      thefit->SetNpx(600);

	      const int N = 1200;   // resolution
	      double xmin = thefit->GetXmin();
	      double xmax = thefit->GetXmax();

	      // Build the filled polygon
	      TGraph* g = new TGraph(2*N);
	      
	      for (int i = 0; i < N; i++) {
		double x = xmin + (xmax - xmin) * i / (N - 1);
		g->SetPoint(i, x, thefit->Eval(x));         // upper curve (left → right)
	      }
	      for (int i = 0; i < N; i++) {
		double x = xmax - (xmax - xmin) * i / (N - 1);
		g->SetPoint(N + i, x, p0->Eval(x));     // lower curve (right → left)
	      }

	      g->SetFillColorAlpha(kAzure+2, 0.35);
	      g->SetLineWidth(0);
	      g->Draw("F");

	      TBox* thebox = new TBox(-30,0,30,p0->GetParameter(0));

	      thebox->SetFillColorAlpha(kRed,0.35);
	      thebox->SetLineWidth(0);
	      thebox->Draw();
	      
	      thefit->Draw("SAME");
	      usegaus->Draw("SAME");
	      theleg->SetFillStyle(0);
	      theleg->SetFillColor(0);
	      theleg->SetBorderSize(0);
	      TLine* lines[2] = {new TLine(dtc[0],0,dtc[0],nomcuts[i][j]->GetMaximum()),new TLine(3,0,3,nomcuts[i][j]->GetMaximum())};
	      
	      for(int k=0; k<2; ++k)
		{
		  lines[k]->SetLineColor(kGreen);
		  lines[k]->SetLineStyle(9);
		  lines[k]->SetLineWidth(2);
		  lines[k]->Draw();
		}
	      theleg->AddEntry(lines[0],"Cut Bounds","l");
	      theleg->Draw();

	      maintexts(0.91,0.17,0,0.03,1,0);
	      drawText("Background included",0.17,0.8,0,kBlack,0.03);
	      drawText("No reconstructed",0.17,0.76,0,kBlack,0.03);
	      drawText("z_{vtx} requirement",0.17,0.72,0,kBlack,0.03);
	      drawText("Dijet cut applied",0.17,0.68,0,kBlack,0.03);

	      c->SaveAs(("../../images/effs/timeff_fit_dt_"+to_string(i)+"_"+to_string(j)+".pdf").c_str());
	      
	      delete usegaus;
	      delete logaus;
	      delete higaus;
	    }
	  
	}
    }

  TFile* file = TFile::Open("effs_timing.root","RECREATE");
  file->cd();
  int colors[3] = {kRed+2, kAzure+2, kSpring+2};
  string labels[3] = {"Nominal cut","Widened by 1 ns","Narrowed by 1 ns"};
  string type[3] = {"t_{lead} Cut Efficiency","#Delta-t Cut Efficiency","Both Timing Cuts Efficiency"};
  TLegend* leg = new TLegend(0.6,0.62,0.9,0.74);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  TH1D* effc[3][3];
  for(int i=0; i<3; ++i)
    {
      effs[i][1]->GetYaxis()->SetTitle(type[i].c_str());
      cout << type[i] << endl;
      for(int j=1; j>=0; --j)
	{
	  effc[i][j] = (TH1D*)effs[i][j]->Clone(("effc"+to_string(i)+to_string(j)).c_str());
	  effc[i][j]->SetLineColor(kBlack);
	  effc[i][j]->SetMarkerStyle(1);
	  cout <<labels[j] << endl;
	  effs[i][j]->Write();
	  effs[i][j]->GetYaxis()->SetRangeUser(0,1.1);
	  effs[i][j]->SetMarkerStyle(20);
	  effs[i][j]->SetMarkerSize(2);
	  effs[i][j]->SetMarkerColor(colors[j]);
	  effs[i][j]->SetLineColor(kBlack);
	  if(i==0) leg->AddEntry(effs[i][j],labels[j].c_str(),"p");
	  effs[i][j]->GetYaxis()->SetRangeUser(0.9,1);
	  if(j==1) effs[i][j]->Draw("P");
	  else effs[i][j]->Draw("SAME P");
	  effc[i][j]->Draw("SAME PE");
	  for(int k=0; k<effs[i][j]->GetNbinsX(); ++k)
	    {
	      cout << effs[i][j]->GetBinContent(k+1) << endl;
	    }
	  
	}
      leg->Draw();
      maintexts(0.9,0.6);
      //drawText(type[i].c_str(),0.6,0.78,0,kBlack,0.03);
      c->SaveAs(("../../images/effs/effstype"+to_string(i)+".pdf").c_str());
    }
  return 0;
}
