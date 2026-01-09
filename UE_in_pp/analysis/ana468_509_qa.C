#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled);
void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_et);

////////////////////////////////////////// Main Function //////////////////////////////////////////
void ana468_509_qa(int runnumber = 47289)  {
  /////////////// General Set up ///////////////
  string outfilename = "ana468_509_qa_output/ana468_509_run" + to_string(runnumber) + "_jet10GeV_qa.root";
  TFile *f_out = new TFile(outfilename.c_str(), "RECREATE");

  /////////////// Read Files ///////////////
  const char* baseDirJet = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput";
  TChain chain("T");
  chain.Add(Form("%s/output_ana509_%d_*.root", baseDirJet, runnumber));
  chain.SetBranchStatus("*", 0);
  TChain chain2("T");
  chain2.Add(Form("%s/output_0mrad_ana468_%d_*.root", baseDirJet, runnumber));
  chain2.SetBranchStatus("*", 0);

  float zvertex; chain.SetBranchStatus("zvtx", 1); chain.SetBranchAddress("zvtx", &zvertex);
  std::vector<int>* gl1_trigger_vector_scaled = nullptr; chain.SetBranchStatus("triggerVector", 1); chain.SetBranchAddress("triggerVector", &gl1_trigger_vector_scaled);
  std::vector<float>* unsubjet_e = nullptr; chain.SetBranchStatus("e", 1); chain.SetBranchAddress("e", &unsubjet_e);
  std::vector<float>* unsubjet_pt = nullptr; chain.SetBranchStatus("pt", 1); chain.SetBranchAddress("pt", &unsubjet_pt);
  std::vector<float>* unsubjet_eta = nullptr; chain.SetBranchStatus("eta", 1); chain.SetBranchAddress("eta", &unsubjet_eta);
  std::vector<float>* unsubjet_phi = nullptr; chain.SetBranchStatus("phi", 1); chain.SetBranchAddress("phi", &unsubjet_phi);
  std::vector<float>* unsubjet_emcal_calo_e = nullptr; chain.SetBranchStatus("jetEmcalE", 1); chain.SetBranchAddress("jetEmcalE", &unsubjet_emcal_calo_e);
  std::vector<float>* unsubjet_ihcal_calo_e = nullptr; chain.SetBranchStatus("jetIhcalE", 1); chain.SetBranchAddress("jetIhcalE", &unsubjet_ihcal_calo_e);
  std::vector<float>* unsubjet_ohcal_calo_e = nullptr; chain.SetBranchStatus("jetOhcalE", 1); chain.SetBranchAddress("jetOhcalE", &unsubjet_ohcal_calo_e);
  int clsmult = 0; float cluster_e[10000] = {0.0}; float cluster_eta[10000] = {0.0}; float cluster_phi[10000] = {0.0};
  chain.SetBranchStatus("clsmult", 1); chain.SetBranchAddress("clsmult",&clsmult);
  chain.SetBranchStatus("cluster_e", 1); chain.SetBranchAddress("cluster_e",cluster_e);
  chain.SetBranchStatus("cluster_eta", 1); chain.SetBranchAddress("cluster_eta",cluster_eta);
  chain.SetBranchStatus("cluster_phi", 1); chain.SetBranchAddress("cluster_phi",cluster_phi);

  float zvertex2; chain2.SetBranchStatus("zvtx", 1); chain2.SetBranchAddress("zvtx", &zvertex2);
  int event2; chain2.SetBranchStatus("m_event",1); chain2.SetBranchAddress("m_event", &event2);
  std::vector<int>* gl1_trigger_vector_scaled2 = nullptr; chain2.SetBranchStatus("triggerVector", 1); chain2.SetBranchAddress("triggerVector", &gl1_trigger_vector_scaled2);
  std::vector<float>* unsubjet_e2 = nullptr; chain2.SetBranchStatus("e", 1); chain2.SetBranchAddress("e", &unsubjet_e2);
  std::vector<float>* unsubjet_pt2 = nullptr; chain2.SetBranchStatus("pt", 1); chain2.SetBranchAddress("pt", &unsubjet_pt2);
  std::vector<float>* unsubjet_eta2 = nullptr; chain2.SetBranchStatus("eta", 1); chain2.SetBranchAddress("eta", &unsubjet_eta2);
  std::vector<float>* unsubjet_phi2 = nullptr; chain2.SetBranchStatus("phi", 1); chain2.SetBranchAddress("phi", &unsubjet_phi2);
  std::vector<float>* unsubjet_emcal_calo_e2 = nullptr; chain2.SetBranchStatus("jetEmcalE", 1); chain2.SetBranchAddress("jetEmcalE", &unsubjet_emcal_calo_e2);
  std::vector<float>* unsubjet_ihcal_calo_e2 = nullptr; chain2.SetBranchStatus("jetIhcalE", 1); chain2.SetBranchAddress("jetIhcalE", &unsubjet_ihcal_calo_e2);
  std::vector<float>* unsubjet_ohcal_calo_e2 = nullptr; chain2.SetBranchStatus("jetOhcalE", 1); chain2.SetBranchAddress("jetOhcalE", &unsubjet_ohcal_calo_e2);
  int clsmult2 = 0; float cluster_e2[10000] = {0.0}; float cluster_eta2[10000] = {0.0}; float cluster_phi2[10000] = {0.0};
  chain2.SetBranchStatus("clsmult", 1); chain2.SetBranchAddress("clsmult",&clsmult2);
  chain2.SetBranchStatus("cluster_e", 1); chain2.SetBranchAddress("cluster_e",cluster_e2);
  chain2.SetBranchStatus("cluster_eta", 1); chain2.SetBranchAddress("cluster_eta",cluster_eta2);
  chain2.SetBranchStatus("cluster_phi", 1); chain2.SetBranchAddress("cluster_phi",cluster_phi2);

  /////////////// Histograms ///////////////
  TH2D *h_zvertex = new TH2D("h_zvertex", "", 200, -100, 100, 200, -100, 100);
  TH2D* h_trigger = new TH2D("h_trigger","",40,0,40,40,0,40);
  TH2D* h_lead_spectra_record = new TH2D("h_lead_spectra_record","", 1000, 0, 100,1000, 0, 100);
  TH2D* h_lead_emcal_e = new TH2D("h_lead_emcal_e","", 1000, 0, 100,1000, 0, 100);
  TH2D* h_lead_ihcal_e = new TH2D("h_lead_ihcal_e","", 1000, 0, 100,1000, 0, 100);
  TH2D* h_lead_ohcal_e = new TH2D("h_lead_ohcal_e","", 1000, 0, 100,1000, 0, 100);

  TH1D* h_lead_spectra2 = new TH1D("h_lead_spectra2","",1000,0,100);
  TH1D* h_trigger2 = new TH1D("h_trigger2","",40,0,40);
  TH1D* h_zvertex2 = new TH1D("h_zvertex2","",200,-100,100);
  TH1D* h_lead_emcal_e2 = new TH1D("h_lead_emcal_e2","", 1000, 0, 100);
  TH1D* h_lead_ihcal_e2 = new TH1D("h_lead_ihcal_e2","", 1000, 0, 100);
  TH1D* h_lead_ohcal_e2 = new TH1D("h_lead_ohcal_e2","", 1000, 0, 100);
  TH1D* h_lead_spectra1 = new TH1D("h_lead_spectra1","",1000,0,100);
  TH1D* h_trigger1 = new TH1D("h_trigger1","",40,0,40);
  TH1D* h_zvertex1 = new TH1D("h_zvertex1","",200,-100,100);
  TH1D* h_lead_emcal_e1 = new TH1D("h_lead_emcal_e1","", 1000, 0, 100);
  TH1D* h_lead_ihcal_e1 = new TH1D("h_lead_ihcal_e1","", 1000, 0, 100);
  TH1D* h_lead_ohcal_e1 = new TH1D("h_lead_ohcal_e1","", 1000, 0, 100);

  TH1D* h_run_ana509 = new TH1D("h_run_ana509","",4000,47250,51250);
  TH1D* h_run_ana468 = new TH1D("h_run_ana468","",4000,47250,51250);
  TH1D* h_run_ratio = new TH1D("h_run_ratio","",4000,47250,51250);

  Long64_t n_events = chain.GetEntries();
  Long64_t n_events2 = chain2.GetEntries();

  /////////////// Event Loop ///////////////
  std::cout << "Data analysis started." << std::endl;
  std::cout << "Total number of events: " << n_events << " " << n_events2 << std::endl;

  for (Long64_t entry = 0; entry < n_events; ++entry) {
    chain.GetEntry(entry);
    if (entry % 10000 == 0) { std::cout << entry << std::endl; }
    if (fabs(zvertex) > 30) { 
      continue; 
    }
    bool jet_trig = false;
    for (int i = 0; i < gl1_trigger_vector_scaled->size(); i++) {
      if (gl1_trigger_vector_scaled->at(i) == 22) {
        jet_trig = true; 
        break;
      }
    }
    if (!jet_trig) {
      continue;
    }
    std::vector<float> recoe_new, recopt_new, recoeta_new, recophi_new, recoemcal_new, recoihcal_new, recoohcal_new;
    for (size_t i = 0; i < unsubjet_eta->size(); ++i) {
      if (fabs(unsubjet_eta->at(i)) < 0.7) {
        recoe_new.push_back(unsubjet_e->at(i));
        recopt_new.push_back(unsubjet_pt->at(i));
        recoeta_new.push_back(unsubjet_eta->at(i));
        recophi_new.push_back(unsubjet_phi->at(i));
        recoemcal_new.push_back(unsubjet_emcal_calo_e->at(i));
        recoihcal_new.push_back(unsubjet_ihcal_calo_e->at(i));
        recoohcal_new.push_back(unsubjet_ohcal_calo_e->at(i));
      }
    }
    *unsubjet_e = std::move(recoe_new);
    *unsubjet_pt = std::move(recopt_new);
    *unsubjet_eta = std::move(recoeta_new);
    *unsubjet_phi = std::move(recophi_new);
    *unsubjet_emcal_calo_e = std::move(recoemcal_new);
    *unsubjet_ihcal_calo_e = std::move(recoihcal_new);
    *unsubjet_ohcal_calo_e = std::move(recoohcal_new);
    if (unsubjet_pt->size() < 1) {
      continue;
    }
    int ind_lead = -1;
    get_leading_jet(ind_lead, unsubjet_pt);
    if (unsubjet_pt->at(ind_lead) < 10.0) {
      continue;
    }

    h_zvertex1->Fill(zvertex);
    for (int i = 0; i < gl1_trigger_vector_scaled->size(); i++) {
      h_trigger1->Fill(gl1_trigger_vector_scaled->at(i));
    }

    h_lead_spectra1->Fill(unsubjet_pt->at(ind_lead));
    h_lead_emcal_e1->Fill(unsubjet_emcal_calo_e->at(ind_lead));
    h_lead_ihcal_e1->Fill(unsubjet_ihcal_calo_e->at(ind_lead));
    h_lead_ohcal_e1->Fill(unsubjet_ohcal_calo_e->at(ind_lead));

  }

  for (Long64_t entry = 0; entry < n_events2; ++entry) {
    chain2.GetEntry(entry);
    if (entry % 10000 == 0) { std::cout << entry << std::endl; }
    if (fabs(zvertex2) > 30) { 
      continue; 
    }
    bool jet_trig = false;
    for (int i = 0; i < gl1_trigger_vector_scaled2->size(); i++) {
      if (gl1_trigger_vector_scaled2->at(i) == 22) {
        jet_trig = true; 
        break;
      }
    }
    if (!jet_trig) {
      continue;
    }
    int ind_lead2 = -1;
    get_leading_jet(ind_lead2, unsubjet_pt2);
    if (unsubjet_pt2->at(ind_lead2) < 10.0) {
      continue;
    }

    h_zvertex2->Fill(zvertex2);
    for (int i = 0; i < gl1_trigger_vector_scaled2->size(); i++) {
      h_trigger2->Fill(gl1_trigger_vector_scaled2->at(i));
    }

    h_lead_spectra2->Fill(unsubjet_pt2->at(ind_lead2));
    h_lead_emcal_e2->Fill(unsubjet_emcal_calo_e2->at(ind_lead2));
    h_lead_ihcal_e2->Fill(unsubjet_ihcal_calo_e2->at(ind_lead2));
    h_lead_ohcal_e2->Fill(unsubjet_ohcal_calo_e2->at(ind_lead2));

  }

  h_run_ratio->Fill(runnumber, h_lead_spectra1->GetEntries()/h_lead_spectra2->GetEntries());
  h_run_ana509->Fill(runnumber, h_lead_spectra1->GetEntries());
  h_run_ana468->Fill(runnumber, h_lead_spectra2->GetEntries());


  /*
  int oldentry = 0;
  bool iterate = false;
  for (Long64_t entry = 0; entry < n_events; ++entry) {
  //for (Long64_t entry = 0; entry < 200; ++entry) {
    if (iterate) oldentry += 1;
    chain.GetEntry(entry);
    chain2.GetEntry(oldentry);
    std::cout << "Entry: " << entry << " " << oldentry << std::endl;
    std::cout << " Vertices: " << zvertex << " " << zvertex2 << std::endl;
    if (oldentry >= 3974) {
      break;
    }
    if (fabs(zvertex) > 30) { 
      iterate = false;
      continue; 
    }

    //std::cout << " Vertices: " << zvertex << " " << zvertex2 << std::endl;

    std::vector<float> recoe_new, recopt_new, recoeta_new, recophi_new, recoemcal_new, recoihcal_new, recoohcal_new;
    for (size_t i = 0; i < unsubjet_eta->size(); ++i) {
      if (fabs(unsubjet_eta->at(i)) < 0.7) {
        recoe_new.push_back(unsubjet_e->at(i));
        recopt_new.push_back(unsubjet_pt->at(i));
        recoeta_new.push_back(unsubjet_eta->at(i));
        recophi_new.push_back(unsubjet_phi->at(i));
        recoemcal_new.push_back(unsubjet_emcal_calo_e->at(i));
        recoihcal_new.push_back(unsubjet_ihcal_calo_e->at(i));
        recoohcal_new.push_back(unsubjet_ohcal_calo_e->at(i));
      }
    }
    *unsubjet_e = std::move(recoe_new);
    *unsubjet_pt = std::move(recopt_new);
    *unsubjet_eta = std::move(recoeta_new);
    *unsubjet_phi = std::move(recophi_new);
    *unsubjet_emcal_calo_e = std::move(recoemcal_new);
    *unsubjet_ihcal_calo_e = std::move(recoihcal_new);
    *unsubjet_ohcal_calo_e = std::move(recoohcal_new);
    if (unsubjet_pt->size() < 1) {
      iterate = false; 
      continue;
    }
    int ind_lead = -1;
    get_leading_jet(ind_lead, unsubjet_pt);
    if (unsubjet_pt->at(ind_lead) < 4.0) {
      iterate = false;
      continue;
    }
    std::cout << " Vertices: " << zvertex << " " << zvertex2 << std::endl;

    int ind_lead2 = -1;
    get_leading_jet(ind_lead2, unsubjet_pt2);

    std::cout << "Leading jets: " << unsubjet_pt->at(ind_lead) << " " << unsubjet_pt2->at(ind_lead2) << " IHCal: " << unsubjet_ihcal_calo_e->at(ind_lead) << " " << unsubjet_ihcal_calo_e2->at(ind_lead2) << " OHCal: " << unsubjet_ohcal_calo_e->at(ind_lead) << " " << unsubjet_ohcal_calo_e2->at(ind_lead2) << std::endl;

    if ((fabs(zvertex - zvertex2) > 0.05)) {
      iterate = false; 
      continue;
    }
    if (fabs(unsubjet_ihcal_calo_e->at(ind_lead) - unsubjet_ihcal_calo_e2->at(ind_lead2)) > 0.005 && fabs(unsubjet_ohcal_calo_e->at(ind_lead) - unsubjet_ohcal_calo_e2->at(ind_lead2)) > 0.05) {
      iterate = false;
      continue;
    }
    //std::cout << "Leading jets: " << unsubjet_pt->at(ind_lead) << " " << unsubjet_pt2->at(ind_lead2) << " IHCal: " << unsubjet_ihcal_calo_e->at(ind_lead) << " " << unsubjet_ihcal_calo_e2->at(ind_lead2) << " OHCal: " << unsubjet_ohcal_calo_e->at(ind_lead) << " " << unsubjet_ohcal_calo_e2->at(ind_lead2) << std::endl;

    //std::cout << "Leading jets: " << unsubjet_pt->at(ind_lead) << " " << unsubjet_pt2->at(ind_lead2) << std::endl;
    std::cout << "Event PASSED " << entry << " " << oldentry << std::endl;
    h_zvertex->Fill(zvertex, zvertex2);
    if (gl1_trigger_vector_scaled->size() == gl1_trigger_vector_scaled2->size()) {
      for (int i = 0; i < gl1_trigger_vector_scaled->size(); i++) {
        h_trigger->Fill(gl1_trigger_vector_scaled->at(i), gl1_trigger_vector_scaled2->at(i));
      }
    } else {
      std::cout << "Event " << entry << " triggers don't match old: ";
      for (int i = 0; i < gl1_trigger_vector_scaled->size(); i++) {
        std::cout << gl1_trigger_vector_scaled->at(i) << " ";
      }
      std::cout << " new: ";
      for (int i = 0; i < gl1_trigger_vector_scaled2->size(); i++) {
        std::cout << gl1_trigger_vector_scaled2->at(i) << " ";
      }
    }

    h_lead_spectra_record->Fill(unsubjet_pt->at(ind_lead),unsubjet_pt2->at(ind_lead2));
    h_lead_emcal_e->Fill(unsubjet_emcal_calo_e->at(ind_lead), unsubjet_emcal_calo_e2->at(ind_lead2));
    h_lead_ihcal_e->Fill(unsubjet_ihcal_calo_e->at(ind_lead), unsubjet_ihcal_calo_e2->at(ind_lead2));
    h_lead_ohcal_e->Fill(unsubjet_ohcal_calo_e->at(ind_lead), unsubjet_ohcal_calo_e2->at(ind_lead2));
    iterate = true;
  }
  */

  std::cout << h_lead_spectra1->GetEntries()/h_lead_spectra2->GetEntries() << std::endl;
  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();
  h_zvertex1->Write();
  h_lead_spectra1->Write();
  h_lead_emcal_e1->Write();
  h_lead_ihcal_e1->Write();
  h_lead_ohcal_e1->Write();
  h_zvertex2->Write();
  h_lead_spectra2->Write();
  h_lead_emcal_e2->Write();
  h_lead_ihcal_e2->Write();
  h_lead_ohcal_e2->Write();
  h_run_ratio->Write();
  h_run_ana468->Write();
  h_run_ana509->Write();
  f_out->Close();
  std::cout << "All done!" << std::endl;
}

////////////////////////////////////////// Functions //////////////////////////////////////////
bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled) {
  if (std::find(gl1_trigger_vector_scaled->begin(), gl1_trigger_vector_scaled->end(), 18) != gl1_trigger_vector_scaled->end()) { return true; }
  if (std::find(gl1_trigger_vector_scaled->begin(), gl1_trigger_vector_scaled->end(), 34) != gl1_trigger_vector_scaled->end()) { return true; }
  return false;
}

void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_et) {
  leadingjet_index = -1;
  float leadingjet_et = -9999;
  for (int ij = 0; ij < jet_et->size(); ++ij) {
    float jetet = jet_et->at(ij);
    if (jetet > leadingjet_et) {
      leadingjet_et = jetet;
      leadingjet_index = ij;
    }
  }
}