#include <TTree.h>
#include <TProfile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <vector>

using namespace std;

void clusterprofile() {
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
 

  string particle = "gamma";
  const char * filename = Form("/sphenix/tg/tg01/jets/samfred/SkimMaker/outtree_towerntuple_vtx150cm_Run24_%s_debug_merged.root",particle.c_str());
  TFile * f = TFile::Open(filename,"READ");
  TTree * t = (TTree*)f->Get("towerntup");

  t->SetBranchStatus("*",0);
  t->SetBranchStatus("truth_vz",1);
  t->SetBranchStatus("nClusters",1);
  t->SetBranchStatus("nClusters_mother",1);
  t->SetBranchStatus("nearbytowers_mother",1);
  t->SetBranchStatus("photon_4mom_mother",1);
  t->SetBranchStatus("truth_diphoton_4mom",1);
  t->SetBranchStatus("nPairs",1);
  t->SetBranchStatus("truth_isconverted",1);

  short nClusters;
  short nClusters_mother;
  float mother_rr;
  vector<vector<float>> * nearbytowers_mother = 0;
  TClonesArray * photon_4mom_mother = 0;
  TClonesArray * truth_diphoton_4mom = 0;
  float truth_eta;
  float truth_phi;
  short nPairs;
  float truth_vz;
  bool truth_isconverted;

  t->SetBranchAddress("nClusters",&nClusters);
  t->SetBranchAddress("nClusters_mother",&nClusters_mother);
  t->SetBranchAddress("nearbytowers_mother",&nearbytowers_mother);
  t->SetBranchAddress("photon_4mom_mother",&photon_4mom_mother);
  t->SetBranchAddress("truth_diphoton_4mom",&truth_diphoton_4mom);
  t->SetBranchAddress("nPairs",&nPairs);
  t->SetBranchAddress("truth_vz",&truth_vz);
  t->SetBranchAddress("truth_isconverted",&truth_isconverted);

  TH2D * energyprof = new TH2D("energyprof",";distance from center; total energy; energy fraction",1000,0,5,1000,0,1);
  TH2D * energyprof_threshcut = new TH2D("energyprof_threshcut",";distance from center; total energy; energy fraction",1000,0,5,1000,0,1);
  TH2D * clustershape[10];
  for (int i = 0; i < 10; i++) {
    clustershape[i] = new TH2D(Form("clustershape_eta%i",i),";etabin;phibin",100,-3.5,3.5,100,-3.5,3.5);
  }
  TProfile * energyprof1D = new TProfile("energyprof1D", ";distance from center;Average energy percentage",1000,0,5);
  TProfile * energyprofs1D[5][3];
  TH2D * energyprofs2D[5][3];
  int energybins[6] = {7,10,15,20,30,1000};
  float etabins[4] = {0,.3,.6,.9};
  for (int i = 0; i < 5; i ++) {
    for (int j = 0; j < 3; j++) {
      energyprofs1D[i][j] = new TProfile(Form("energyprof_%i_%iGeV_%i_%ieta",energybins[i],energybins[i+1],(int)(10*etabins[j]),(int)(10*etabins[j+1])), ";distance from center;Average energy percentage",1000,0,5);
      energyprofs2D[i][j] = new TH2D(Form("energyprof2D_%i_%iGeV_%i_%ieta",energybins[i],energybins[i+1],(int)(10*etabins[j]),(int)(10*etabins[j+1])), ";distance from center;Average energy percentage",1000,0,5,1000,0,1);
    }
  }
  TProfile * energyprof1D_threshcut = new TProfile("energyprof1D_threshcut", ";distance from center;Average energy percentage",1000,0,5);
  TH1D * rr = new TH1D("rr",";distance from center;counts",1000,0,1);

  gStyle->SetOptStat(0);
  int entries = t->GetEntries()/10;
  cout << "I see " << entries << " entries!" << endl;
  for (int ie = 0; ie < entries; ie++) {
    t->GetEntry(ie);
    if (ie % 5000 == 0) cout << "Event " << ie << " / " << entries << ": (" << static_cast<int>(100*static_cast<float>(ie)/static_cast<float>(entries)) << "%)" << endl;
    if (nClusters_mother != 1 || nClusters != 1) continue;
    if (truth_isconverted) continue;
    if (abs(truth_vz) > 30) continue;
    if (nearbytowers_mother->size() != 1) continue;
    TLorentzVector diphoton = *(TLorentzVector*) photon_4mom_mother->At(0);
    TLorentzVector truth_diphoton = *(TLorentzVector*) truth_diphoton_4mom->At(0);
    float diphotoneta = diphoton.Eta();
    float diphotonphi = diphoton.Phi();
    float diphotonpt = diphoton.Pt();
    if (diphotonpt < 1) continue;
    if (abs(diphotoneta) > 0.9) continue;

    float r_error = TMath::Sqrt(TMath::Power(diphotoneta-truth_diphoton.Eta(),2) + TMath::Power(diphotonphi-truth_diphoton.Phi(),2));
    rr->Fill(r_error);
    if (r_error > 0.05) continue;
    
    float eta_cm = 0;
    float phi_cm = 0;
    float e_tot = 0;
    for (int i = 0; i < 49; i++) {
      float energy = nearbytowers_mother->at(0).at(i);
      if (energy < 0.07) continue;
      eta_cm += energy * (i/7);
      phi_cm += energy * (i%7);
      e_tot += energy;
    }
    eta_cm /= e_tot;
    phi_cm /= e_tot;
    
    int h_etabin = (int)((diphotoneta+1.1)/2.2 * 10);
    if (h_etabin < 0) h_etabin = 0;
    if (h_etabin > 9) h_etabin = 9;
    for (int i = 0; i < 49; i++) {
      float energy = nearbytowers_mother->at(0).at(i);
      if (energy < 0.07) continue;
      clustershape[h_etabin]->Fill(i/7-eta_cm,i%7-phi_cm,energy);
    }

    //cout << ie << endl;
    for (int i = 0; i < 49; i++) {
      //cout <<i << " ";
      float dist = TMath::Sqrt(TMath::Power(eta_cm - (i/7),2) + TMath::Power(phi_cm - (i%7),2));
      float energy = nearbytowers_mother->at(0).at(i);
      int energybin;
      int etabin;
      for (int ienergy = 0; ienergy < 5; ienergy++) {
        if (energybins[ienergy] < diphoton.E() && diphoton.E() < energybins[ienergy + 1]) {energybin = ienergy; break;}
      }
      for (int ieta = 0; ieta < 3; ieta++) {
        if (etabins[ieta] < diphoton.Eta() && diphoton.Eta() < etabins[ieta + 1]) {etabin = ieta; break;}
      }
      energyprof->Fill(dist, energy/e_tot);
      energyprof1D->Fill(dist, energy/e_tot);
      energyprofs1D[energybin][etabin]->Fill(dist,energy/e_tot); 
      energyprofs2D[energybin][etabin]->Fill(dist,energy/e_tot);
      if (energy < 0.07) continue;
      energyprof_threshcut->Fill(dist, energy/e_tot);
      energyprof1D_threshcut->Fill(dist, energy/e_tot);

    }
    //cout << endl;
  }
  /*
  TCanvas * c = new TCanvas("c","");
  c->SaveAs("profiles.pdf[");
  energyprof->Draw("colz");
  c->SaveAs("profiles.pdf");
  energyprof1D->Draw();
  c->SaveAs("profiles.pdf");
  energyprof_threshcut->Draw("colz");
  c->SaveAs("profiles.pdf");
  energyprof1D_threshcut->Draw();
  c->SaveAs("profiles.pdf");
  rr->Draw();
  c->SaveAs("profiles.pdf");
  c->SaveAs("profiles.pdf]");
*/
  TFile * of = TFile::Open(Form("profilehists_%s.root",particle.c_str()),"RECREATE");
  energyprof->Write();
  energyprof1D->Write();
  energyprof_threshcut->Write();
  energyprof1D_threshcut->Write();
  rr->Write();
  for (int i = 0; i < 10; i++) {
    clustershape[i]->Write();
  }
  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 3; j++) {
      energyprofs1D[i][j]->Write();
      energyprofs2D[i][j]->Write();
    }
  }
  of->Close();
}
