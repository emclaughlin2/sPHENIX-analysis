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
#include "../unfolding/unfold_Def.h"

bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled);
void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et);
void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_et);
bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi);

////////////////////////////////////////// Main Function //////////////////////////////////////////
void jet_background_data_eff(int runnumber = 51274, std::string bkg_cut = "dijet")  {
  /////////////// General Set up ///////////////
  string outfilename = "analysis_data_ana509_output/output_" + bkg_cut + "_and_time_bkg_cut_"+ to_string(runnumber) + ".root";
  TFile *f_out = new TFile(outfilename.c_str(), "RECREATE");

  /////////////// Read Files ///////////////
  const char* baseDirJet = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput";
  TChain chain("T");
  chain.Add(Form("%s/output_ana509_v2_%d_*.root", baseDirJet, runnumber));
  chain.SetBranchStatus("*", 0);

  float zvertex; chain.SetBranchStatus("zvtx", 1); chain.SetBranchAddress("zvtx", &zvertex);
  std::vector<int>* gl1_trigger_vector_scaled = nullptr; chain.SetBranchStatus("triggerVector", 1); chain.SetBranchAddress("triggerVector", &gl1_trigger_vector_scaled);
  std::vector<float>* unsubjet_e = nullptr; chain.SetBranchStatus("e", 1); chain.SetBranchAddress("e", &unsubjet_e);
  std::vector<float>* unsubjet_pt = nullptr; chain.SetBranchStatus("pt", 1); chain.SetBranchAddress("pt", &unsubjet_pt);
  std::vector<float>* unsubjet_eta = nullptr; chain.SetBranchStatus("eta", 1); chain.SetBranchAddress("eta", &unsubjet_eta);
  std::vector<float>* unsubjet_phi = nullptr; chain.SetBranchStatus("phi", 1); chain.SetBranchAddress("phi", &unsubjet_phi);
  std::vector<float>* unsubjet_emcal_calo_e = nullptr; chain.SetBranchStatus("jetEmcalE", 1); chain.SetBranchAddress("jetEmcalE", &unsubjet_emcal_calo_e);
  std::vector<float>* unsubjet_ihcal_calo_e = nullptr; chain.SetBranchStatus("jetIhcalE", 1); chain.SetBranchAddress("jetIhcalE", &unsubjet_ihcal_calo_e);
  std::vector<float>* unsubjet_ohcal_calo_e = nullptr; chain.SetBranchStatus("jetOhcalE", 1); chain.SetBranchAddress("jetOhcalE", &unsubjet_ohcal_calo_e);
  std::vector<float>* unsubjet_time = nullptr; chain.SetBranchStatus("jettime", 1); chain.SetBranchAddress("jettime", &unsubjet_time);
  std::vector<float>* unsubjet_unweighttime = nullptr; chain.SetBranchStatus("jetunweighttime",1); chain.SetBranchAddress("jetunweighttime", &unsubjet_unweighttime);
  float mbd_t0; chain.SetBranchStatus("mbd_t0", 1); chain.SetBranchAddress("mbd_t0", &mbd_t0);

  int emcaln = 0; float emcale[24576] = {0.0}; float emcaleta[24576] = {0.0}; float emcalphi[24576] = {0.0};
  int ihcaln = 0; float ihcale[1536] = {0.0}; float ihcaleta[1536] = {0.0}; float ihcalphi[1536] = {0.0};
  int ohcaln = 0; float ohcale[1536] = {0.0}; float ohcaleta[1536] = {0.0}; float ohcalphi[1536] = {0.0};
  int clsmult = 0; float cluster_e[10000] = {0.0}; float cluster_eta[10000] = {0.0}; float cluster_phi[10000] = {0.0};
  int cluster_ntowers[10000] = {0};
  chain.SetBranchStatus("clsmult", 1); chain.SetBranchAddress("clsmult",&clsmult);
  chain.SetBranchStatus("cluster_e", 1); chain.SetBranchAddress("cluster_e",cluster_e);
  chain.SetBranchStatus("cluster_eta", 1); chain.SetBranchAddress("cluster_eta",cluster_eta);
  chain.SetBranchStatus("cluster_phi", 1); chain.SetBranchAddress("cluster_phi",cluster_phi);
  chain.SetBranchStatus("cluster_ntowers", 1); chain.SetBranchAddress("cluster_ntowers",cluster_ntowers);

  /////////////// JES func ///////////////
  TFile *corrFile = new TFile("../unfolding/JES_Calib_Default.root", "READ");
  if (!corrFile) {
    std::cout << "Error: cannot open JES_Calib_Default.root" << std::endl;
    return;
  }
  TF1 *f_corr = (TF1*)corrFile->Get("JES_Calib_Default_Func");
  if (!f_corr) {
    std::cout << "Error: cannot open f_corr" << std::endl;
    return;
  }

  double tight_calibptbins[] = {22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 59, 63.5};
  int tight_calibnpt = sizeof(tight_calibptbins) / sizeof(tight_calibptbins[0]) - 1;

  /////////////// Histograms ///////////////
  TH1D *h_zvertex = new TH1D("h_zvertex", ";Z-vertex [cm]", 200, -100, 100);
  TH1D* h_deltaphi_record = new TH1D("h_deltaphi_record","",125,-2*M_PI,2*M_PI);
  TH1D* h_xj_record = new TH1D("h_xj_record","",20,0,1);
  TH1D* h_lead_spectra_record = new TH1D("h_lead_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
  TH1D* h_sub_spectra_record = new TH1D("h_sub_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_et_transverse_record = new TH1D("h_et_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
  TH2D* h_ue_pt_transverse_record = new TH2D("h_ue_pt_transverse_record","", calibnpt, 0, 1, calibnet, 0, 1);

  TH2D* h_total_measure = new TH2D("h_total_measure","", calibnpt, calibptbins, calibnet, calibetbins);
  TH2D* h_pass_time_cut_measure = new TH2D("h_pass_time_cut_measure","", calibnpt, calibptbins, calibnet, calibetbins);
  TH2D* h_total_measure_tight = new TH2D("h_total_measure_tight","", tight_calibnpt, tight_calibptbins, calibnet, calibetbins);
  TH2D* h_pass_time_cut_measure_tight = new TH2D("h_pass_time_cut_measure_tight","", tight_calibnpt, tight_calibptbins, calibnet, calibetbins);
  TH1D* h_total_measure_recojet = new TH1D("h_total_measure_recojet","", 100, 0, 100);
  TH1D* h_pass_time_cut_measure_recojet = new TH1D("h_pass_time_cut_measure_recojet","",100,0,100);
  TH1D* h_total_measure_calibjet = new TH1D("h_total_measure_calibjet","", tight_calibnpt, tight_calibptbins);
  TH1D* h_pass_time_cut_measure_calibjet = new TH1D("h_pass_time_cut_measure_calibjet","",tight_calibnpt, tight_calibptbins);
  TH2D* h_pass_syst_time_cut_measure = new TH2D("h_pass_syst_time_cut_measure","", calibnpt, calibptbins, calibnet, calibetbins);
  TH2D* h_pass_syst_time_cut_measure_tight = new TH2D("h_pass_syst_time_cut_measure_tight","", tight_calibnpt, tight_calibptbins, calibnet, calibetbins);
  TH1D* h_pass_syst_time_cut_measure_recojet = new TH1D("h_pass_syst_time_cut_measure_recojet","",100,0,100);
  TH1D* h_pass_syst_time_cut_measure_calibjet = new TH1D("h_pass_syst_time_cut_measure_calibjet","",tight_calibnpt, tight_calibptbins);

  TH1D* h_mbd_time = new TH1D("h_mbd_time","",200,-20,20);
  TH1D* h_jettime = new TH1D("h_jettime","",200,-20,20);
  TH1D* h_jet_mbd_timediff = new TH1D("h_jet_mbd_timediff","",200,-20,20);
  TH1D* h_jetunweighttime = new TH1D("h_jetunweighttime","",200,-20,20);
  TH2D* h_jet_mbd_time = new TH2D("h_jet_mbd_time","",200,-20,20,200,-20,20);
  TH1D* h_deltajettime = new TH1D("h_deltajettime","",200,-20,20);
  TH1D* h_deltajetunweighttime = new TH1D("h_deltajetunweighttime","",200,-20,20);
  TH2D* h_2D_deltajettime = new TH2D("h_2D_deltajettime","",200,-20,20,200,-20,20);
  TH2D* h_2D_deltajetunweighttime = new TH2D("h_2D_deltajetunweighttime","",200,-20,20,200,-20,20);
  TH2D* h_jettimecorr = new TH2D("h_jettimecorr","",200,-20,20,200,-20,20);
  TH2D* h_deltatcorr = new TH2D("h_deltatcorr","",200,-20,20,200,-20,20);
  TH1D* h_et_jet_mbd_timediff_gt_0 = new TH1D("h_et_jet_mbd_timediff_gt_0","",500,0,100);
  TH1D* h_et_jet_mbd_timediff_lt_0 = new TH1D("h_et_jet_mbd_timediff_lt_0","",500,0,100);
  TH1D* h_trans_et_jet_mbd_timediff_gt_0 = new TH1D("h_trans_et_jet_mbd_timediff_gt_0","",500,0,100);
  TH1D* h_trans_et_jet_mbd_timediff_lt_0 = new TH1D("h_trans_et_jet_mbd_timediff_lt_0","",500,0,100);
  TH2D* h_em_jet_frac_jettime = new TH2D("h_em_jet_frac_jettime","",5,0,1,200,-20,20);
  TH2D* h_em_jet_frac_jet_mbd_timediff = new TH2D("h_em_jet_frac_jet_mbd_timediff","",5,0,1,200,-20,20);
  TH2D* h_em_jet_frac_mbd_time = new TH2D("h_em_jet_frac_mbd_time","",5,0,1,200,-20,20);
  TH2D* h_jet_pt_jettime = new TH2D("h_jet_pt_jettime","",4,20,60,200,-20,20);
  TH2D* h_jet_pt_jet_mbd_timediff = new TH2D("h_jet_pt_jet_mbd_timediff","",4,20,60,200,-20,20);
  TH2D* h_jet_pt_mbd_time = new TH2D("h_jet_pt_mbd_time","",4,20,60,200,-20,20);
  TH3D* h_jet_pt_jet_time_jet_mbd_timediff = new TH3D("h_jet_pt_jet_time_jet_mbd_timediff","",2,20,40,200,-20,20,200,-20,20);

  TH3D* hpttdtdat = new TH3D("hpttdtdat","",10,0,100,600,-30,30,600,-30,30);
  
  /////////////// Event Loop ///////////////
  std::cout << "Data analysis started." << std::endl;
  int n_events = chain.GetEntries();
  std::cout << "Total number of events: " << n_events << std::endl;

  // Event variables setup.
  bool reco_bkg_cut = false; 
  string deltat_type = "jet_mbd";

  for (int ie = 0; ie < n_events; ++ie) { // event loop start
  //for (int ie = 0; ie < 20; ++ie) {
    // Load event.
    if (ie % 1000 == 0) {
      std::cout << "Processing event " << ie << "..." << std::endl;
    }
    chain.GetEntry(ie);

    // Trigger and Z-vertex cut.
    if (!check_bad_trigger(gl1_trigger_vector_scaled)) continue;
    if (isnan(zvertex)) { continue; }
    if (fabs(zvertex) > 60) continue;

    //////////////////////////// SETUP JET VARIABLES FOR UNFOLDING ////////////////////////////

    std::vector<float> recoe_new, recopt_new, recoeta_new, recophi_new, recoemcal_new, recoihcal_new, recoohcal_new, recotime_new, reconwtime_new;

    for (size_t i = 0; i < unsubjet_eta->size(); ++i) {
        if (!check_bad_jet_eta(unsubjet_eta->at(i), zvertex, jet_radius)) {
            recoe_new.push_back(unsubjet_e->at(i));
            recopt_new.push_back(unsubjet_pt->at(i));
            recoeta_new.push_back(unsubjet_eta->at(i));
            recophi_new.push_back(unsubjet_phi->at(i));
            recoemcal_new.push_back(unsubjet_emcal_calo_e->at(i));
            recoihcal_new.push_back(unsubjet_ihcal_calo_e->at(i));
            recoohcal_new.push_back(unsubjet_ohcal_calo_e->at(i));
            recotime_new.push_back(unsubjet_time->at(i));
            reconwtime_new.push_back(unsubjet_unweighttime->at(i));
        }
    }

    // Replace original contents
    *unsubjet_e = std::move(recoe_new);
    *unsubjet_pt = std::move(recopt_new);
    *unsubjet_eta = std::move(recoeta_new);
    *unsubjet_phi = std::move(recophi_new);
    *unsubjet_emcal_calo_e = std::move(recoemcal_new);
    *unsubjet_ihcal_calo_e = std::move(recoihcal_new);
    *unsubjet_ohcal_calo_e = std::move(recoohcal_new);
    *unsubjet_time = std::move(recotime_new);
    *unsubjet_unweighttime = std::move(reconwtime_new);

    // indices to find leading and subleading jets 
    int ind_lead = -1; int ind_sub = -1;
    float lead_e = 0;
    reco_bkg_cut = false;

    int nJetReq = 1;
    if (bkg_cut == "dijet") { nJetReq = 2; }
    if (unsubjet_pt->size() < nJetReq) { continue; }

    // find reco leading and subleading jets if nreco jets >= nJetReq
    TVector3 lead;
    if (bkg_cut == "dijet") {
        get_leading_subleading_jet(ind_lead, ind_sub, unsubjet_pt); 
        //std::cout << " sub/lead: " << unsubjet_e->at(ind_sub)/unsubjet_e->at(ind_lead) << " dPhi: " << get_dphi(unsubjet_phi->at(ind_lead), unsubjet_phi->at(ind_sub)) << std::endl; 
        if (unsubjet_e->at(ind_sub)/unsubjet_e->at(ind_lead) > 0.3 && match_leading_subleading_jet(unsubjet_phi->at(ind_lead), unsubjet_phi->at(ind_sub))) {
            reco_bkg_cut = true;
        } else {
            reco_bkg_cut = false;
        }
    } else if (bkg_cut == "efrac") {
        get_leading_jet(ind_lead, unsubjet_pt);
        float jete = unsubjet_e->at(ind_lead);
        //std::cout << " EMCal frac: " << unsubjet_emcal_calo_e->at(ind_lead)/jete << " IHCal frac: " << unsubjet_ihcal_calo_e->at(ind_lead)/jete << " OHCal frac: " << unsubjet_ohcal_calo_e->at(ind_lead)/jete << std::endl;
        if (unsubjet_emcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_emcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ihcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ohcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_ohcal_calo_e->at(ind_lead)/jete < 0.9) {
            reco_bkg_cut = true; 
        } else {
            reco_bkg_cut = false;
        }
    } else if (bkg_cut == "efrac_bkg_only") {
        get_leading_jet(ind_lead, unsubjet_pt);
        float jete = unsubjet_e->at(ind_lead);
        //std::cout << " EMCal frac: " << unsubjet_emcal_calo_e->at(ind_lead)/jete << " IHCal frac: " << unsubjet_ihcal_calo_e->at(ind_lead)/jete << " OHCal frac: " << unsubjet_ohcal_calo_e->at(ind_lead)/jete << std::endl;
        if (unsubjet_emcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_emcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ihcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ohcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_ohcal_calo_e->at(ind_lead)/jete < 0.9) {
            reco_bkg_cut = false; 
        } else {
            reco_bkg_cut = true;
        }
    } else {
        get_leading_jet(ind_lead, unsubjet_pt);
        float jete = unsubjet_e->at(ind_lead);
        reco_bkg_cut = true;
    }
    lead.SetPtEtaPhi(unsubjet_pt->at(ind_lead), unsubjet_eta->at(ind_lead), unsubjet_phi->at(ind_lead));
    lead_e = unsubjet_e->at(ind_lead);
    
    if (lead_e < 0) { continue; }

    if (lead.Pt() > 20.0 && lead.Pt() < 60.0 && bkg_cut == "efrac_bkg_only" && reco_bkg_cut && unsubjet_emcal_calo_e->at(ind_lead)/lead_e > 0.1 && unsubjet_emcal_calo_e->at(ind_lead)/lead_e < 0.9) {
        std::cout << "bkg_cut " << reco_bkg_cut << " EMCal frac: " << unsubjet_emcal_calo_e->at(ind_lead)/lead_e << " IHCal frac: " << unsubjet_ihcal_calo_e->at(ind_lead)/lead_e << " OHCal frac: " << unsubjet_ohcal_calo_e->at(ind_lead)/lead_e << " Jet Uncalib pT " << lead.Pt() << std::endl;
    }

    // Fill z-vertex histogram.
    h_zvertex->Fill(zvertex);

    float lead_time = unsubjet_time->at(ind_lead)*17.6;
    float lead_nwtime = unsubjet_unweighttime->at(ind_lead)*17.6;
    float lead_deltat, lead_nwdeltat;
    if (deltat_type == "dijet") {
      lead_deltat = (unsubjet_time->at(ind_lead) - unsubjet_time->at(ind_sub))*17.6;
      lead_nwdeltat = (unsubjet_unweighttime->at(ind_lead) - unsubjet_unweighttime->at(ind_sub))*17.6;
    } else {
      lead_deltat = unsubjet_time->at(ind_lead)*17.6 - mbd_t0;
      lead_nwdeltat = unsubjet_unweighttime->at(ind_lead)*17.6 - mbd_t0;
    }
    if (lead.Pt() > 20.0 && lead.Pt() < 60.0 && reco_bkg_cut) {
      h_jettime->Fill(lead_time);
      h_jetunweighttime->Fill(lead_nwtime);
      h_jettimecorr->Fill(lead_time, lead_nwtime);
      h_mbd_time->Fill(mbd_t0);
      h_jet_mbd_time->Fill(lead_time, mbd_t0);
      h_jet_mbd_timediff->Fill(lead_time - mbd_t0);
      h_deltajettime->Fill(lead_deltat);
      h_deltajetunweighttime->Fill(lead_nwdeltat);
      h_2D_deltajettime->Fill(lead_time, lead_deltat);
      h_2D_deltajetunweighttime->Fill(lead_nwtime, lead_nwdeltat);
      h_deltatcorr->Fill(lead_deltat, lead_nwdeltat);
      h_em_jet_frac_jettime->Fill(unsubjet_emcal_calo_e->at(ind_lead)/lead_e, lead_time);
      h_em_jet_frac_jet_mbd_timediff->Fill(unsubjet_emcal_calo_e->at(ind_lead)/lead_e, lead_time - mbd_t0);
      h_em_jet_frac_mbd_time->Fill(unsubjet_emcal_calo_e->at(ind_lead)/lead_e, mbd_t0);
      h_jet_pt_jettime->Fill(lead.Pt(), lead_time);
      h_jet_pt_jet_mbd_timediff->Fill(lead.Pt(), lead_time - mbd_t0);
      h_jet_pt_mbd_time->Fill(lead.Pt(), mbd_t0);
      h_jet_pt_jet_time_jet_mbd_timediff->Fill(lead.Pt(), lead_time, lead_time - mbd_t0);
      if (lead_time == 0) {
        std::cout << "Jet uncalib pT " << lead.Pt() << " eta " << lead.Eta() << " phi " << lead.Phi() << " bkg_cut " << reco_bkg_cut << " EMCal frac: " << unsubjet_emcal_calo_e->at(ind_lead)/lead_e << " IHCal frac: " << unsubjet_ihcal_calo_e->at(ind_lead)/lead_e << " OHCal frac: " << unsubjet_ohcal_calo_e->at(ind_lead)/lead_e << " Jet time " << lead_time << " Jet energy " << lead_e << std::endl;
        for (int i = 0; i < clsmult; i++) {
            TVector3 clus;
            clus.SetPtEtaPhi(cluster_e[i]/cosh(cluster_eta[i]),cluster_eta[i],cluster_phi[i]);
            double deltaR = lead.DeltaR(clus);
            //if (fabs(deltaR) < 0.4) {
                std::cout << "\t cluster e: " << cluster_e[i] << " eta: " << cluster_eta[i] << " phi: " << cluster_phi[i] << " ntowers: " << cluster_ntowers[i] << std::endl;
            //}
        }
      }
    }
    if (lead.Pt() > 5.0 && reco_bkg_cut) {
        hpttdtdat->Fill(lead.Pt(), lead_time, lead_time - mbd_t0);
    }

    TVector3 caliblead;
    caliblead.SetPtEtaPhi(f_corr->Eval(lead.Pt()), lead.Eta(), lead.Phi());

    //std::cout << "Reco lead: pt " << lead.Pt() << " e " << lead_e << " eta " << lead.Eta() << " phi " << lead.Phi()  << " bkg_cut " << reco_bkg_cut << std::endl;
    //std::cout << "Good reco lead: pt " << caliblead.Pt() << " eta " << caliblead.Eta() << " phi " << caliblead.Phi() << std::endl;
    
    //////////////////////////// SETUP UE VARIABLES //////////////////////////////

    // find reco ET information 
    float et_transverse = 0;
    float total_et = 0;
    for (int i = 0; i < clsmult; i++) {
        float dphi = get_dphi(lead.Phi(),cluster_phi[i]);
        if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += cluster_e[i]/cosh(cluster_eta[i]); }
        total_et += cluster_e[i]/cosh(cluster_eta[i]);
    }

    if (caliblead.Pt() > 21.0 && reco_bkg_cut) { 
      h_total_measure_recojet->Fill(lead.Pt());
      h_total_measure->Fill(caliblead.Pt(), et_transverse);
      h_total_measure_tight->Fill(caliblead.Pt(), et_transverse);
      h_total_measure_calibjet->Fill(caliblead.Pt());
      if (fabs(lead_time + 2.0) < 6.0 && lead_deltat >= -5.0 && lead_deltat <= 0.0) {
        h_pass_time_cut_measure_recojet->Fill(lead.Pt());
        h_pass_time_cut_measure->Fill(caliblead.Pt(), et_transverse);
        h_pass_time_cut_measure_tight->Fill(caliblead.Pt(), et_transverse);
        h_pass_time_cut_measure_calibjet->Fill(caliblead.Pt());
      }
      if (fabs(lead_time + 2.0) < 7.0 && lead_deltat >= -6.0 && lead_deltat <= 1.0) {
        h_pass_syst_time_cut_measure_recojet->Fill(lead.Pt());
        h_pass_syst_time_cut_measure->Fill(caliblead.Pt(), et_transverse);
        h_pass_syst_time_cut_measure_tight->Fill(caliblead.Pt(), et_transverse);
        h_pass_syst_time_cut_measure_calibjet->Fill(caliblead.Pt());
      }
    }

    if (lead.Pt() > 20.0 && lead.Pt() < 40.0 && reco_bkg_cut) {
      if (lead_time - mbd_t0 > 0) {
        h_et_jet_mbd_timediff_gt_0->Fill(total_et);
        h_trans_et_jet_mbd_timediff_gt_0->Fill(et_transverse);
      } else {
        h_et_jet_mbd_timediff_lt_0->Fill(total_et);
        h_trans_et_jet_mbd_timediff_lt_0->Fill(et_transverse);
      }
    }

    // Subleading jet E / Leading jet E > 0.3.
    // Delta phi between subleading jet and leading jet > 7*pi/8
    // | Leading jet T + 2 ns| < 5 ns (between -7 to 3 ns)
    // Delta T between subleading jet and leading jet < 3 ns (1 sample = 17.6 ns)

    //std::cout << "Reco ET " << et_transverse << std::endl;
    //std::cout << std::endl;

    //////////////////////////// RECORD QA PLOTS FOR ALL EVENTS IN UNFOLDING PROCEDURE ////////////////////////////
    h_lead_spectra_record->Fill(lead.Pt());
    h_et_transverse_record->Fill(et_transverse);
    h_ue_pt_transverse_record->Fill(lead.Pt(),et_transverse);

} // event loop end

  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();
  h_zvertex->Write();
  h_lead_spectra_record->Write();
  h_et_transverse_record->Write();
  h_ue_pt_transverse_record->Write();
  h_total_measure->Write();
  h_pass_time_cut_measure->Write();
  h_total_measure_tight->Write();
  h_pass_time_cut_measure_tight->Write();
  h_total_measure_recojet->Write();
  h_pass_time_cut_measure_recojet->Write();
  h_total_measure_calibjet->Write();
  h_pass_time_cut_measure_calibjet->Write();
  h_mbd_time->Write();
  h_jet_mbd_time->Write();
  h_jet_mbd_timediff->Write();
  h_jettime->Write();
  h_jetunweighttime->Write();
  h_deltajettime->Write();
  h_deltajetunweighttime->Write();
  h_2D_deltajettime->Write();
  h_2D_deltajetunweighttime->Write();
  h_jettimecorr->Write();
  h_deltatcorr->Write();
  h_et_jet_mbd_timediff_gt_0->Write();
  h_et_jet_mbd_timediff_lt_0->Write();
  h_trans_et_jet_mbd_timediff_gt_0->Write();
  h_trans_et_jet_mbd_timediff_lt_0->Write();
  h_em_jet_frac_jettime->Write();
  h_em_jet_frac_jet_mbd_timediff->Write();
  h_em_jet_frac_mbd_time->Write();
  h_jet_pt_jettime->Write();
  h_jet_pt_jet_mbd_timediff->Write();
  h_jet_pt_mbd_time->Write();
  h_pass_syst_time_cut_measure_recojet->Write();
  h_pass_syst_time_cut_measure->Write();
  h_pass_syst_time_cut_measure_tight->Write();
  h_pass_syst_time_cut_measure_calibjet->Write();
  h_jet_pt_jet_time_jet_mbd_timediff->Write();
  hpttdtdat->Write();
  f_out->Close();
  std::cout << "All done!" << std::endl;
}

////////////////////////////////////////// Functions //////////////////////////////////////////
bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled) {
  if (std::find(gl1_trigger_vector_scaled->begin(), gl1_trigger_vector_scaled->end(), 18) != gl1_trigger_vector_scaled->end()) { return true; }
  if (std::find(gl1_trigger_vector_scaled->begin(), gl1_trigger_vector_scaled->end(), 34) != gl1_trigger_vector_scaled->end()) { return true; }
  if (std::find(gl1_trigger_vector_scaled->begin(), gl1_trigger_vector_scaled->end(), 22) != gl1_trigger_vector_scaled->end()) { return true; }
  return false;
}

void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et) {
  leadingjet_index = -1;
  subleadingjet_index = -1;
  float leadingjet_et = -9999;
  float subleadingjet_et = -9999;
  for (int ij = 0; ij < jet_et->size(); ++ij) {
    float jetet = jet_et->at(ij);
    if (jetet > leadingjet_et) {
      subleadingjet_et = leadingjet_et;
      subleadingjet_index = leadingjet_index;
      leadingjet_et = jetet;
      leadingjet_index = ij;
    } else if (jetet > subleadingjet_et) {
      subleadingjet_et = jetet;
      subleadingjet_index = ij;
    }
  }
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

bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi) {
  float dijet_min_phi = 7*TMath::Pi()/8.;
  float dphi = get_dphi(leadingjet_phi, subleadingjet_phi);
  //std::cout << "delta phi: " << dphi << std::endl;
  return fabs(dphi) > dijet_min_phi;
}