#include <iostream>
#include <TH2D.h>
#include <TH1D.h>
#include <TChain.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <sstream> //std::ostringstsream
#include <fstream> //std::ifstream
#include <iostream> //std::cout, std::endl
#include <cmath>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <string>
#include <set>
#include <TVector3.h>
#include <map>
#include <vector>
#include <TDatabasePDG.h>
#include <tuple>
#include <TProfile.h>
#include <TProfile2D.h>
#include "TRandom.h"
#include "TRandom3.h"
#include "../unfolding/unfold_Def.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

using namespace std;

void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et);
void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_et, std::vector<float>* jet_eta);
bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi);
void get_calibjet(float& calibjet_pt, float& calibjet_eta, float& calibjet_phi, float jet_pt, float jet_eta, float jet_phi, bool reco_bkg_cut, TF1* f_corr, float jes_para, float jer_para);
void get_truthjet(float& goodtruthjet_pt, float& goodtruthjet_eta, float& goodtruthjet_phi, float jet_pt, float jet_eta, float jet_phi);
void match_meas_truth(float meas_eta, float meas_phi, bool& matched, float truth_eta, float truth_phi, float jet_radius);
void find_match_meas_truth(float meas_eta, float meas_phi, bool& matched, float& deltaR, float truth_eta, float truth_phi, float jet_radius);

TRandom3 randGen(1234);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                       //
//  Note: 10.14.25                                                                                                                       //
//  This macro creates efficiency plots for jet background cuts for 2D unfolding of leading jet pT and transverse region energy density  //
//                                                                                                                                       //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double MapToUniform(double et, const double* bin_edges, int nbins) {
    if (et < bin_edges[0] || et >= bin_edges[nbins]) return -1; // out of range
    for (int i = 0; i < nbins; ++i) {
        if (et >= bin_edges[i] && et < bin_edges[i + 1]) {
            double bin_center_fraction = (i + 0.5) / nbins;
            return bin_center_fraction;
        }
    }
    return -1;
}

void jet_background_eff_sim(std::string runtype = "mb", int start_seg = 0, int end_seg = 200, std::string bkg_cut = "dijet")  {
    
    ////////// General Set up //////////
    double weight_scale = 1.0, truthjet_pt_min = 0, truthjet_pt_max = 3000;
    string count_filename;
    if (runtype == "mb") {
        std::cout << "Not using MB dataset" << std::endl;
        //weight_scale = MB_scale;
        //truthjet_pt_min = 0;
        //truthjet_pt_max = 14;
    } else if (runtype == "jet10") {
        weight_scale = Jet10GeV_scale;
        truthjet_pt_min = 14;
        truthjet_pt_max = 17;
    } else if (runtype == "jet15") {
        weight_scale = Jet15GeV_scale;
        truthjet_pt_min = 17;
        truthjet_pt_max = 22;
    } else if (runtype == "jet20") {
        weight_scale = Jet20GeV_scale;
        truthjet_pt_min = 22;
        truthjet_pt_max = 35;
    } else if (runtype == "jet30") {
        weight_scale = Jet30GeV_scale;
        truthjet_pt_min = 35;
        truthjet_pt_max = 52;
    } else if (runtype == "jet50") {
        weight_scale = Jet50GeV_scale;
        truthjet_pt_min = 52;
        truthjet_pt_max = 71;
    } else if (runtype == "jet70") {
        weight_scale = Jet70GeV_scale;
        truthjet_pt_min = 71;
        truthjet_pt_max = 1000;
    } else if (runtype == "herwig_jet10") {
        weight_scale = Herwig_Jet10GeV_scale;
        truthjet_pt_min = 14;
        truthjet_pt_max = 35;
    } else if (runtype == "herwig_jet30") {
        weight_scale = Herwig_Jet30GeV_scale;
        truthjet_pt_min = 35;
        truthjet_pt_max = 3000;
    } else {
        std::cout << "Unknown runtype" << std::endl;
        return;
    }

    ////////// Files //////////
    string outfilename = "analysis_sim_run28_output/output_sim_" + runtype + "_" + bkg_cut + "_bkg_cut_" + to_string(start_seg) + "_" + to_string(end_seg) + ".root";
    TFile *f_out = new TFile(outfilename.c_str(), "RECREATE");
    TChain chain("T");
    for (int i = start_seg; i < end_seg; ++i) {
        chain.Add(Form("/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_%s_3sigma_output_%d.root", runtype.c_str(), i));
    }
    chain.SetBranchStatus("*", 0);

    float zvertex; chain.SetBranchStatus("zvtx", 1); chain.SetBranchAddress("zvtx", &zvertex);
    std::vector<float>* unsubjet_e = nullptr; chain.SetBranchStatus("e", 1); chain.SetBranchAddress("e", &unsubjet_e);
    std::vector<float>* unsubjet_pt = nullptr; chain.SetBranchStatus("pt", 1); chain.SetBranchAddress("pt", &unsubjet_pt);
    std::vector<float>* unsubjet_eta = nullptr; chain.SetBranchStatus("eta", 1); chain.SetBranchAddress("eta", &unsubjet_eta);
    std::vector<float>* unsubjet_phi = nullptr; chain.SetBranchStatus("phi", 1); chain.SetBranchAddress("phi", &unsubjet_phi);
    std::vector<float>* unsubjet_emcal_calo_e = nullptr; chain.SetBranchStatus("jetEmcalE", 1); chain.SetBranchAddress("jetEmcalE", &unsubjet_emcal_calo_e);
    std::vector<float>* unsubjet_ihcal_calo_e = nullptr; chain.SetBranchStatus("jetIhcalE", 1); chain.SetBranchAddress("jetIhcalE", &unsubjet_ihcal_calo_e);
    std::vector<float>* unsubjet_ohcal_calo_e = nullptr; chain.SetBranchStatus("jetOhcalE", 1); chain.SetBranchAddress("jetOhcalE", &unsubjet_ohcal_calo_e);

    std::vector<float>* truthjet_e = nullptr; chain.SetBranchStatus("truthE", 1); chain.SetBranchAddress("truthE", &truthjet_e);
    std::vector<float>* truthjet_pt = nullptr; chain.SetBranchStatus("truthPt", 1); chain.SetBranchAddress("truthPt", &truthjet_pt);
    std::vector<float>* truthjet_eta = nullptr; chain.SetBranchStatus("truthEta", 1); chain.SetBranchAddress("truthEta", &truthjet_eta);
    std::vector<float>* truthjet_phi = nullptr; chain.SetBranchStatus("truthPhi", 1); chain.SetBranchAddress("truthPhi", &truthjet_phi);

    int truthpar_n = 0; chain.SetBranchStatus("truthpar_n", 1); chain.SetBranchAddress("truthpar_n",&truthpar_n);
    float truthpar_e[100000] = {0.0}; chain.SetBranchStatus("truthpar_e", 1); chain.SetBranchAddress("truthpar_e",truthpar_e);
    float truthpar_eta[100000] = {0.0}; chain.SetBranchStatus("truthpar_eta", 1); chain.SetBranchAddress("truthpar_eta",truthpar_eta);
    float truthpar_phi[100000] = {0.0}; chain.SetBranchStatus("truthpar_phi", 1); chain.SetBranchAddress("truthpar_phi",truthpar_phi);
    int truthpar_pid[100000] = {0}; chain.SetBranchStatus("truthpar_pid", 1); chain.SetBranchAddress("truthpar_pid",truthpar_pid);

    int emcaln = 0; float emcale[24576] = {0.0}; float emcaleta[24576] = {0.0}; float emcalphi[24576] = {0.0};
    int ihcaln = 0; float ihcale[1536] = {0.0}; float ihcaleta[1536] = {0.0}; float ihcalphi[1536] = {0.0};
    int ohcaln = 0; float ohcale[1536] = {0.0}; float ohcaleta[1536] = {0.0}; float ohcalphi[1536] = {0.0};
    int clsmult = 0; float cluster_e[10000] = {0.0}; float cluster_eta[10000] = {0.0}; float cluster_phi[10000] = {0.0};
    int cluster_ntowers[2000]; int cluster_tower_calo[200][500]; int cluster_tower_ieta[200][500]; int cluster_tower_iphi[200][500]; float cluster_tower_e[200][500];
    
    chain.SetBranchStatus("clsmult", 1); chain.SetBranchAddress("clsmult",&clsmult);
    chain.SetBranchStatus("cluster_e", 1); chain.SetBranchAddress("cluster_e",cluster_e);
    chain.SetBranchStatus("cluster_eta", 1); chain.SetBranchAddress("cluster_eta",cluster_eta);
    chain.SetBranchStatus("cluster_phi", 1); chain.SetBranchAddress("cluster_phi",cluster_phi);
    chain.SetBranchStatus("cluster_ntowers", 1); chain.SetBranchAddress("cluster_ntowers",cluster_ntowers);
    chain.SetBranchStatus("cluster_tower_e", 1); chain.SetBranchAddress("cluster_tower_e",cluster_tower_e);
    chain.SetBranchStatus("cluster_tower_calo", 1); chain.SetBranchAddress("cluster_tower_calo",cluster_tower_calo);
    chain.SetBranchStatus("cluster_tower_ieta", 1); chain.SetBranchAddress("cluster_tower_ieta",cluster_tower_ieta);
    chain.SetBranchStatus("cluster_tower_iphi", 1); chain.SetBranchAddress("cluster_tower_iphi",cluster_tower_iphi);

    ////////// JES func //////////
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

    ////////// Histograms //////////
    TH1D* h_zvertex = new TH1D("h_zvertex", ";Z-vertex [cm]", 400, -200, 200);
    TH1D* h_lead_spectra_record = new TH1D("h_lead_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
    TH1D* h_lead_truth_spectra_record = new TH1D("h_lead_truth_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
    TH2D* h_jes_qa = new TH2D("h_jes_qa",";p^{truth}_{T} [GeV]; p^{reco}_{T}/p^{truth}_{T}", 50, 0, 100, 1000, 0, 10);
   
    double tight_calibptbins[] = {22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 53, 56, 59, 63.5};
    int tight_calibnpt = sizeof(tight_calibptbins) / sizeof(tight_calibptbins[0]) - 1;

    TH2D* h_efficiency_reco_bkg_cut = new TH2D("h_efficiency_reco_bkg_cut", "", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_match = new TH2D("h_efficiency_match", "", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_reco_pt_cut = new TH2D("h_efficiency_reco_pt_cut","", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_match_reco_bkg_cut = new TH2D("h_efficiency_match_reco_bkg_cut", "", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_reco_pt_cut_reco_bkg_cut = new TH2D("h_efficiency_reco_pt_cut_reco_bkg_cut","", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_match_reco_pt_cut_reco_bkg_cut = new TH2D("h_efficiency_match_reco_pt_cut_reco_bkg_cut", "", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_match_reco_pt_cut = new TH2D("h_efficiency_match_reco_pt_cut","", truthnpt, truthptbins, truthnet, truthetbins);

    TH2D* h_purity_truth_pt_cut = new TH2D("h_purity_truth_pt_cut","", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_purity_match = new TH2D("h_purity_match","", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_truth_qa = new TH2D("h_truth_qa","", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_qa = new TH2D("h_measure_qa","", calibnpt, calibptbins, calibnet, calibetbins);
    TH1D* efficiency_truth_match_test = new TH1D("efficiency_truth_match_test","", truthnpt, truthptbins);

    TH1D *h_et_transverse_record = new TH1D("h_et_transverse_record", ";#SigmaE_{T} [GeV]", calibnet, calibetbins);
    TH1D* h_et_truth_transverse_record = new TH1D("h_et_truth_transverse_record", ";#SigmaE_{T} [GeV]", truthnet, truthetbins);
    TH1D *h_nw_et_transverse_record = new TH1D("h_nw_et_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH1D* h_nw_et_truth_transverse_record = new TH1D("h_nw_et_truth_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH2D* h_ue_pt_transverse_record = new TH2D("h_ue_pt_transverse_record","", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_ue_pt_truth_transverse_record = new TH2D("h_ue_pt_truth_transverse_record","",truthnpt, truthptbins, truthnet, truthetbins);

    TH2D* h_total_measure = new TH2D("h_total_measure","", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_pass_cut_measure = new TH2D("h_pass_cut_measure","", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_total_measure_tight = new TH2D("h_total_measure_tight","", tight_calibnpt, tight_calibptbins, calibnet, calibetbins);
    TH2D* h_pass_cut_measure_tight = new TH2D("h_pass_cut_measure_tight","", tight_calibnpt, tight_calibptbins, calibnet, calibetbins);
    
    TH1D* h_total_measure_recojet = new TH1D("h_total_measure_recojet","", 100, 0, 100);
    TH1D* h_pass_cut_measure_recojet = new TH1D("h_pass_cut_measure_recojet","",100,0,100);
    TH1D* h_total_measure_recojet_matched = new TH1D("h_total_measure_recojet_matched","", 100, 0, 100);
    TH1D* h_pass_cut_measure_recojet_matched = new TH1D("h_pass_cut_measure_recojet_matched","",100,0,100);
    TH1D* h_total_measure_calibjet = new TH1D("h_total_measure_calibjet","", tight_calibnpt, tight_calibptbins);
    TH1D* h_pass_cut_measure_calibjet = new TH1D("h_pass_cut_measure_calibjet","",tight_calibnpt, tight_calibptbins);
    TH1D* h_total_alljet = new TH1D("h_total_alljet","", 100, 0, 100);
    TH1D* h_pass_cut_alljet = new TH1D("h_pass_cut_alljet","",100,0,100);

    ////////// Event Loop //////////
    std::cout << "Data analysis started." << std::endl;
    Long64_t nEntries = chain.GetEntries();
    std::cout << "Total number of events: " << nEntries << std::endl;
    
    // Event variables setup.
    bool reco_bkg_cut = false;
    float goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi;
    float calibjet_pt_dijet, calibjet_eta_dijet, calibjet_phi_dijet;
    bool calibjet_matched_dijet, qa_matched;
    
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 2000; ++entry) {
        if (entry % 1000 == 0) cout << "event " << entry << endl;
        chain.GetEntry(entry);

        // Z-vertex cut.
        if (isnan(zvertex)) { continue; }
        if (fabs(zvertex) > 60) { continue; }

        // Fill z-vertex histogram.
        h_zvertex->Fill(zvertex);

        //////////////////////////// SETUP JET VARIABLES FOR UNFOLDING ////////////////////////////

        // indices to find leading and subleading jets 
        int ind_truth_lead = -1; int ind_truth_sub = -1;
        int ind_lead = -1; int ind_sub = -1;
        float lead_e = 0; float truthlead_e = 0;
        reco_bkg_cut = false;

        //std::cout << "reco jets " << unsubjet_pt->size() << " truth jets " << truthjet_pt->size() << std::endl;
        
        std::vector<float> truthe_new, truthpt_new, trutheta_new, truthphi_new;
        std::vector<float> recoe_new, recopt_new, recoeta_new, recophi_new, recoemcal_new, recoihcal_new, recoohcal_new;
        
        for (size_t i = 0; i < truthjet_eta->size(); ++i) {
            if (!check_bad_jet_eta(truthjet_eta->at(i), zvertex, jet_radius)) {
                truthe_new.push_back(truthjet_e->at(i));
                truthpt_new.push_back(truthjet_pt->at(i));
                trutheta_new.push_back(truthjet_eta->at(i));
                truthphi_new.push_back(truthjet_phi->at(i));
            }
        }

        // Replace original contents
        *truthjet_e = std::move(truthe_new);
        *truthjet_pt = std::move(truthpt_new);
        *truthjet_eta = std::move(trutheta_new);
        *truthjet_phi = std::move(truthphi_new);
        /*
        for (size_t i = 0; i < unsubjet_eta->size(); ++i) {
            if (!check_bad_jet_eta(unsubjet_eta->at(i), zvertex, jet_radius)) {
                recoe_new.push_back(unsubjet_e->at(i));
                recopt_new.push_back(unsubjet_pt->at(i));
                recoeta_new.push_back(unsubjet_eta->at(i));
                recophi_new.push_back(unsubjet_phi->at(i));
                recoemcal_new.push_back(unsubjet_emcal_calo_e->at(i));
                recoihcal_new.push_back(unsubjet_ihcal_calo_e->at(i));
                recoohcal_new.push_back(unsubjet_ohcal_calo_e->at(i));
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
        */
        int nJetReq = 1;
        if (bkg_cut == "dijet") { nJetReq = 2; }
        if (unsubjet_pt->size() < nJetReq) {
            continue;
        }
        
        //std::cout << "reco jets " << unsubjet_pt->size() << " truth jets " << truthjet_pt->size() << std::endl;

        //find truth leading if ntruth jets >= 1
        TVector3 truthlead;
        if (truthjet_pt->size() >= 1) {
            get_leading_jet(ind_truth_lead, truthjet_pt, truthjet_eta);
            truthlead.SetPtEtaPhi(truthjet_pt->at(ind_truth_lead), truthjet_eta->at(ind_truth_lead), truthjet_phi->at(ind_truth_lead));
            truthlead_e = truthjet_e->at(ind_truth_lead);
            // if leading truth jet outside of trigger range, discard event
            if (truthlead.Pt() < truthjet_pt_min || truthlead.Pt() > truthjet_pt_max) { continue; } 
        } else {
            truthlead.SetPtEtaPhi(0,0,0);
            truthlead_e = 0;
        }

        // find reco leading and subleading jets if nreco jets >= nJetReq
        TVector3 lead;
        if (unsubjet_pt->size() >= nJetReq) {
            if (bkg_cut == "dijet") {
                get_leading_subleading_jet(ind_lead, ind_sub, unsubjet_pt); 
                //std::cout << " sub/lead: " << unsubjet_e->at(ind_sub)/unsubjet_e->at(ind_lead) << " dPhi: " << get_dphi(unsubjet_phi->at(ind_lead), unsubjet_phi->at(ind_sub)) << std::endl; 
                if (unsubjet_e->at(ind_sub)/unsubjet_e->at(ind_lead) > 0.3 && match_leading_subleading_jet(unsubjet_phi->at(ind_lead), unsubjet_phi->at(ind_sub))) {
                    reco_bkg_cut = true;
                } else {
                    reco_bkg_cut = false;
                }
            } else if (bkg_cut == "efrac") {
                get_leading_jet(ind_lead, unsubjet_pt, unsubjet_eta);
                float jete = unsubjet_e->at(ind_lead);
                //std::cout << " EMCal frac: " << unsubjet_emcal_calo_e->at(ind_lead)/jete << " IHCal frac: " << unsubjet_ihcal_calo_e->at(ind_lead)/jete << " OHCal frac: " << unsubjet_ohcal_calo_e->at(ind_lead)/jete << std::endl;
                if (unsubjet_emcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_emcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ihcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ohcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_ohcal_calo_e->at(ind_lead)/jete < 0.9) {
                    reco_bkg_cut = true; 
                } else {
                    reco_bkg_cut = false;
                }
            } else {
                reco_bkg_cut = true;
            }
            for (size_t i = 0; i < unsubjet_eta->size(); ++i) {
                if (!check_bad_jet_eta(unsubjet_eta->at(i), zvertex, jet_radius)) {
                    recoe_new.push_back(unsubjet_e->at(i));
                    recopt_new.push_back(unsubjet_pt->at(i));
                    recoeta_new.push_back(unsubjet_eta->at(i));
                    recophi_new.push_back(unsubjet_phi->at(i));
                    recoemcal_new.push_back(unsubjet_emcal_calo_e->at(i));
                    recoihcal_new.push_back(unsubjet_ihcal_calo_e->at(i));
                    recoohcal_new.push_back(unsubjet_ohcal_calo_e->at(i));
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
            if (unsubjet_pt->size() >= 1) {
                get_leading_jet(ind_lead, unsubjet_pt, unsubjet_eta);
                lead.SetPtEtaPhi(unsubjet_pt->at(ind_lead), unsubjet_eta->at(ind_lead), unsubjet_phi->at(ind_lead));
                lead_e = unsubjet_e->at(ind_lead);
            } else {
                lead.SetPtEtaPhi(0,0,0);
                lead_e = 0;
                reco_bkg_cut = false;
            }
        } else {
            lead.SetPtEtaPhi(0,0,0);
            lead_e = 0;
            reco_bkg_cut = false;
        }

        //std::cout << "Reco lead: pt " << lead.Pt() << " e " << lead_e << " eta " << lead.Eta() << " phi " << lead.Phi()  << " bkg_cut " << reco_bkg_cut << std::endl;
        //std::cout << "Truth lead: pt " << truthlead.Pt() << " truth e " << truthlead_e << " truth eta " << truthlead.Eta() << " truth phi " << truthlead.Phi() << std::endl;
        get_truthjet(goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi, truthlead.Pt(), truthlead.Eta(), truthlead.Phi());
        get_calibjet(calibjet_pt_dijet, calibjet_eta_dijet, calibjet_phi_dijet, lead.Pt(), lead.Eta(), lead.Phi(), reco_bkg_cut, f_corr, 1, 0.1);    
        //std::cout << "Good reco lead: pt " << calibjet_pt_dijet << " eta " << calibjet_eta_dijet << " phi " << calibjet_phi_dijet << " bkg_cut " << reco_bkg_cut << std::endl;
        //std::cout << "Good truth lead: pt " << goodtruthjet_pt << " truth eta " << goodtruthjet_eta << " truth phi " << goodtruthjet_phi << std::endl;
        match_meas_truth(calibjet_eta_dijet, calibjet_phi_dijet, calibjet_matched_dijet, goodtruthjet_eta, goodtruthjet_phi, jet_radius);     
        //std::cout << "Match " << calibjet_matched_dijet << std::endl;
        
        match_meas_truth(lead.Eta(),lead.Phi(),qa_matched,truthlead.Eta(),truthlead.Phi(),jet_radius);
        double calib_pt_qa = f_corr->Eval(lead.Pt()) * (1 + randGen.Gaus(0.0, 0.1)) * 1;
        //std::cout << "calib_pt_qa: " << calib_pt_qa << " calibjet_pt_dijet: " << calibjet_pt_dijet << std::endl;
        if (truthlead.Pt() >= truthptbins[0] && truthlead.Pt() < truthptbins[truthnpt] && qa_matched) {
            //std::cout << "matched: " << qa_matched << " respmatrix match: " << calibjet_matched_dijet << " truth pt: " << truthlead.Pt() << " ratio: " << f_corr->Eval(lead.Pt())/truthlead.Pt() << " weight: " << weight_scale << std::endl;
            h_jes_qa->Fill(truthlead.Pt(), f_corr->Eval(lead.Pt())/truthlead.Pt(), weight_scale);
        }

        if (lead.Pt() > 14) {
            h_total_measure_recojet->Fill(lead.Pt(), weight_scale);
            if (reco_bkg_cut) h_pass_cut_measure_recojet->Fill(lead.Pt(), weight_scale);
            if (truthlead.Pt() >= truthptbins[0] && truthlead.Pt() < truthptbins[truthnpt] && qa_matched) {
                h_total_measure_recojet_matched->Fill(lead.Pt(), weight_scale);
                if (reco_bkg_cut) {
                    h_pass_cut_measure_recojet_matched->Fill(lead.Pt(), weight_scale);
                }
            }
        }

        if (lead.Pt() > 14) {
            for (int i = 0; i < unsubjet_pt->size(); i++) {
                if (check_bad_jet_eta(unsubjet_eta->at(i), zvertex, jet_radius)) { continue; }
                h_total_alljet->Fill(unsubjet_pt->at(i), weight_scale);
                if (reco_bkg_cut) {
                    h_pass_cut_alljet->Fill(unsubjet_pt->at(i), weight_scale);
                }
            }
        }

        bool truth_match = false;
        float truth_deltaR0 = 999.0; float truth_deltaR = 999.0;
        for (int i = 0; i < unsubjet_pt->size(); i++) {
            find_match_meas_truth(unsubjet_eta->at(i),unsubjet_phi->at(i),truth_match,truth_deltaR,truthlead.Eta(),truthlead.Phi(),jet_radius);
            if (truth_deltaR < truth_deltaR0) { truth_deltaR0 = truth_deltaR; }
        }
        if (truthlead.Pt() >= truthptbins[0] && truthlead.Pt() < truthptbins[truthnpt] && truth_deltaR0 < 0.75*jet_radius) {
            efficiency_truth_match_test->Fill(truthlead.Pt(), weight_scale);
        }
    
        //////////////////////////// SETUP UE VARIABLES FOR UNFOLDING //////////////////////////////

        // find reco ET information 
        float et_transverse = 0;
        for (int i = 0; i < clsmult; i++) {
            float dphi = get_dphi(lead.Phi(),cluster_phi[i]);
            if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += cluster_e[i]/cosh(cluster_eta[i]); }
        }

        // find truth ET information
        float truth_et_transverse = 0;
        for (int i = 0; i < truthpar_n; i++) {
            if (fabs(truthpar_eta[i]) > 1.1) { continue; }
            if ((truthpar_pid[i] == 22 || truthpar_pid[i] == 111) && fabs(truthpar_e[i]) < 0.2) { continue; }
            else if (fabs(truthpar_e[i]) < 0.2) { continue; } // edited from 0.5 to 0.2
            float dphi = get_dphi(truthlead.Phi(),truthpar_phi[i]);
            if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { truth_et_transverse += truthpar_e[i]/cosh(truthpar_eta[i]); } 
        }

        // can toggle reco pT cut, reco bkg cut and match requirement for efficiency 
        if (truthlead.Pt() >= truthptbins[0] && truthlead.Pt() < truthptbins[truthnpt]) {
            h_truth_qa->Fill(truthlead.Pt(), truth_et_transverse, weight_scale);
            if (reco_bkg_cut) { h_efficiency_reco_bkg_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale); }
            if (qa_matched) { h_efficiency_match->Fill(truthlead.Pt(), truth_et_transverse, weight_scale); }
            if (calib_pt_qa >= calibptbins[0] && calib_pt_qa < calibptbins[calibnpt]) { h_efficiency_reco_pt_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale); }
            if (reco_bkg_cut && qa_matched) { h_efficiency_match_reco_bkg_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale); }
            if (reco_bkg_cut && calib_pt_qa >= calibptbins[0] && calib_pt_qa < calibptbins[calibnpt]) { h_efficiency_reco_pt_cut_reco_bkg_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale); }
            if (qa_matched && calib_pt_qa >= calibptbins[0] && calib_pt_qa < calibptbins[calibnpt]) { h_efficiency_match_reco_pt_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale); }
            if (reco_bkg_cut && qa_matched && calib_pt_qa >= calibptbins[0] && calib_pt_qa < calibptbins[calibnpt]) { h_efficiency_match_reco_pt_cut_reco_bkg_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale); }
        }

        // can toggle truth pT cut and match requirement for purity 
        if (calibjet_pt_dijet >= calibptbins[0] && calibjet_pt_dijet < calibptbins[calibnpt] && reco_bkg_cut) {
            h_measure_qa->Fill(calibjet_pt_dijet, et_transverse, weight_scale);
            if (truthlead.Pt() >= truthptbins[0] && truthlead.Pt() < truthptbins[truthnpt]) { h_purity_truth_pt_cut->Fill(calibjet_pt_dijet, et_transverse, weight_scale); }
            if (calibjet_matched_dijet) { h_purity_match->Fill(calibjet_pt_dijet, et_transverse, weight_scale); }
        }

        //////////////////////////// RECORD QA PLOTS FOR ALL EVENTS IN UNFOLDING PROCEDURE ////////////////////////////
        if (reco_bkg_cut) {
            h_lead_spectra_record->Fill(lead.Pt(), weight_scale);
            if (calibjet_pt_dijet > calibptbins[0]) {
                h_nw_et_transverse_record->Fill(et_transverse, weight_scale);
                h_et_transverse_record->Fill(et_transverse, weight_scale);
                h_ue_pt_transverse_record->Fill(lead.Pt(),et_transverse, weight_scale);
            }
        }
            
        h_lead_truth_spectra_record->Fill(truthlead.Pt(), weight_scale);
        if (goodtruthjet_pt > truthptbins[0]) {
            h_nw_et_truth_transverse_record->Fill(truth_et_transverse, weight_scale);
            h_et_truth_transverse_record->Fill(truth_et_transverse, weight_scale);
            h_ue_pt_truth_transverse_record->Fill(truthlead.Pt(),truth_et_transverse, weight_scale);
        }

        // Create efficiency plots here
        if (calibjet_pt_dijet >= calibptbins[0] && calibjet_pt_dijet < calibptbins[calibnpt]) {
            h_total_measure->Fill(calibjet_pt_dijet, et_transverse, weight_scale);
            h_total_measure_tight->Fill(calibjet_pt_dijet, et_transverse, weight_scale);
            h_total_measure_calibjet->Fill(calibjet_pt_dijet, weight_scale);
            if (reco_bkg_cut) {
                h_pass_cut_measure->Fill(calibjet_pt_dijet, et_transverse, weight_scale);
                h_pass_cut_measure_tight->Fill(calibjet_pt_dijet, et_transverse, weight_scale);
                h_pass_cut_measure_calibjet->Fill(calibjet_pt_dijet, weight_scale);
            }
        }
    }

    std::cout << "Writing histograms..." << std::endl;
    f_out->cd();

    h_zvertex->Write();  h_jes_qa->Write(); h_truth_qa->Write(); h_measure_qa->Write(); efficiency_truth_match_test->Write();
    h_efficiency_match_reco_bkg_cut->Write(); h_efficiency_reco_pt_cut_reco_bkg_cut->Write(); h_efficiency_match_reco_pt_cut_reco_bkg_cut->Write(); h_efficiency_match_reco_pt_cut->Write();
    h_purity_match->Write(); h_purity_truth_pt_cut->Write(); h_efficiency_match->Write(); h_efficiency_reco_pt_cut->Write(); h_efficiency_reco_bkg_cut->Write(); 
    h_lead_spectra_record->Write(); h_et_transverse_record->Write(); h_nw_et_transverse_record->Write(); h_ue_pt_transverse_record->Write(); 
    h_lead_truth_spectra_record->Write(); h_et_truth_transverse_record->Write(); h_nw_et_truth_transverse_record->Write(); h_ue_pt_truth_transverse_record->Write(); 
    h_total_measure->Write(); h_total_measure_tight->Write(); h_pass_cut_measure->Write(); h_pass_cut_measure_tight->Write();
    h_total_measure_recojet->Write(); h_total_measure_calibjet->Write(); h_pass_cut_measure_recojet->Write(); h_pass_cut_measure_calibjet->Write();
    h_total_measure_recojet_matched->Write(); h_pass_cut_measure_recojet_matched->Write(); h_total_alljet->Write(); h_pass_cut_alljet->Write();

    f_out->Close();
    std::cout << "All done!" << std::endl;

}

////////////////////////////////////////// Functions //////////////////////////////////////////
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

void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_et, std::vector<float>* jet_eta) {
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
  float dijet_min_phi = 3*TMath::Pi()/4.;
  float dphi = get_dphi(leadingjet_phi, subleadingjet_phi);
  return fabs(dphi) > dijet_min_phi;
}

void get_calibjet(float& calibjet_pt, float& calibjet_eta, float& calibjet_phi, float jet_pt, float jet_eta, float jet_phi, bool reco_bkg_cut, TF1* f_corr, float jes_para, float jer_para) {
    calibjet_pt = -9999;
    calibjet_eta = -9999;
    calibjet_phi = -9999;
    //if (!reco_bkg_cut) return;
    double calib_pt = f_corr->Eval(jet_pt) * (1 + randGen.Gaus(0.0, jer_para)) * jes_para;
    if (calib_pt < calibptbins[0] || calib_pt > calibptbins[calibnpt]) return;
    calibjet_pt = calib_pt;
    calibjet_eta = jet_eta;
    calibjet_phi = jet_phi;
}

void get_truthjet(float& goodtruthjet_pt, float& goodtruthjet_eta, float& goodtruthjet_phi, float jet_pt, float jet_eta, float jet_phi) {
    goodtruthjet_pt = -9999;
    goodtruthjet_eta = -9999;
    goodtruthjet_phi = -9999;
    if (jet_pt < truthptbins[0] || jet_pt > truthptbins[truthnpt]) return;
    goodtruthjet_pt = jet_pt;
    goodtruthjet_eta = jet_eta;
    goodtruthjet_phi = jet_phi;
}

void match_meas_truth(float meas_eta, float meas_phi, bool& matched, float truth_eta, float truth_phi, float jet_radius) {
    if (meas_eta == -9999 || meas_phi == -9999 || truth_eta == -9999 || truth_phi == -9999) {
        matched = false;
        return; 
    } else {
        float dR = get_dR(meas_eta, meas_phi, truth_eta, truth_phi);
        if (dR < jet_radius * 0.75) {
            //std::cout << "dR = " << dR << std::endl;
            matched = true;
            return;
        } else {
            //std::cout << "dR = " << dR << std::endl;
            matched = false;
            return;
        }
    }
}

void find_match_meas_truth(float meas_eta, float meas_phi, bool& matched, float& deltaR, float truth_eta, float truth_phi, float jet_radius) {
    if (meas_eta == -9999 || meas_phi == -9999 || truth_eta == -9999 || truth_phi == -9999) {
        matched = false;
        return; 
    } else {
        float dR = get_dR(meas_eta, meas_phi, truth_eta, truth_phi);
        deltaR = dR;
        if (dR < jet_radius * 0.75) {
            //std::cout << "dR = " << dR << std::endl;
            matched = true;
            return;
        } else {
            //std::cout << "dR = " << dR << std::endl;
            matched = false;
            return;
        }
    }
}