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
#include <TH2D.h>

using namespace std;

bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled);
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
//                                                                                                                                        //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mbd_vertex_eff_data(int runnumber = 51274, std::string bkg_cut = "dijet", float zvertex_cut = 20.0, float wide_zvertex_cut = 100.0)  {
    /////////////// General Set up ///////////////
    string outfilename = "analysis_data_ana509_output/output_zvtx_lt_" + to_string(floor(zvertex_cut)) + "_" + to_string(floor(wide_zvertex_cut)) + "_" + bkg_cut + "_"+ to_string(runnumber) + ".root";
    TFile *f_out = new TFile(outfilename.c_str(), "RECREATE");

    /////////////// Read Files ///////////////
    const char* baseDirJet = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput";
    TChain chain("T");
    chain.Add(Form("%s/output_mbd_eff_%d_*.root", baseDirJet, runnumber));
    chain.SetBranchStatus("*", 0);

    std::vector<int>* gl1_trigger_vector_scaled = nullptr; chain.SetBranchStatus("triggerVector", 1); chain.SetBranchAddress("triggerVector", &gl1_trigger_vector_scaled);
    float zvertex; chain.SetBranchStatus("zvtx", 1); chain.SetBranchAddress("zvtx", &zvertex);
    float trackvertex; chain.SetBranchStatus("zsvtx", 1); chain.SetBranchAddress("zsvtx", &trackvertex);
    float siliconvertex; chain.SetBranchStatus("zsiliconvtx", 1); chain.SetBranchAddress("zsiliconvtx", &siliconvertex);
    std::vector<float>* unsubjet_e = nullptr; chain.SetBranchStatus("e", 1); chain.SetBranchAddress("e", &unsubjet_e);
    std::vector<float>* unsubjet_pt = nullptr; chain.SetBranchStatus("pt", 1); chain.SetBranchAddress("pt", &unsubjet_pt);
    std::vector<float>* unsubjet_eta = nullptr; chain.SetBranchStatus("eta", 1); chain.SetBranchAddress("eta", &unsubjet_eta);
    std::vector<float>* unsubjet_phi = nullptr; chain.SetBranchStatus("phi", 1); chain.SetBranchAddress("phi", &unsubjet_phi);
    std::vector<float>* unsubjet_emcal_calo_e = nullptr; chain.SetBranchStatus("jetEmcalE", 1); chain.SetBranchAddress("jetEmcalE", &unsubjet_emcal_calo_e);
    std::vector<float>* unsubjet_ihcal_calo_e = nullptr; chain.SetBranchStatus("jetIhcalE", 1); chain.SetBranchAddress("jetIhcalE", &unsubjet_ihcal_calo_e);
    std::vector<float>* unsubjet_ohcal_calo_e = nullptr; chain.SetBranchStatus("jetOhcalE", 1); chain.SetBranchAddress("jetOhcalE", &unsubjet_ohcal_calo_e);
    std::vector<float>* svtxVector = nullptr; chain.SetBranchStatus("svtxVector", 1); chain.SetBranchAddress("svtxVector", &svtxVector);
    std::vector<int>* svtxBcoVector = nullptr; chain.SetBranchStatus("svtxBcoVector", 1); chain.SetBranchAddress("svtxBcoVector", &svtxBcoVector);
    std::vector<float>* siliconVector = nullptr; chain.SetBranchStatus("siliconVector", 1); chain.SetBranchAddress("siliconVector", &siliconVector);
    std::vector<int>* siliconBcoVector = nullptr; chain.SetBranchStatus("siliconBcoVector", 1); chain.SetBranchAddress("siliconBcoVector", &siliconBcoVector);

    int trkmult = 0; chain.SetBranchStatus("trkmult", 1); chain.SetBranchAddress("trkmult",&trkmult);
    int tr_nintt[2000] = {0}; chain.SetBranchStatus("tr_nintt", 1); chain.SetBranchAddress("tr_nintt", tr_nintt);
    int tr_ntpc[2000] = {0}; chain.SetBranchStatus("tr_ntpc", 1); chain.SetBranchAddress("tr_ntpc", tr_ntpc);
    int tr_crossing[2000] = {0}; chain.SetBranchStatus("tr_crossing", 1); chain.SetBranchAddress("tr_crossing", tr_crossing);
    float tr_pt[2000] = {0}; chain.SetBranchStatus("tr_pt", 1); chain.SetBranchAddress("tr_pt", tr_pt);

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

    double calibtotaletbins[441] = {0.0};
    for (int i = 0; i < 441; i++) {
        calibtotaletbins[i] = -20.0 + 0.5*i;
    }
    int calibtotalnet = sizeof(calibtotaletbins) / sizeof(calibtotaletbins[0]) - 1;

    ////////// Histograms //////////
    TH1D* h_lead_spectra_record = new TH1D("h_lead_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
    TH1D *h_et_transverse_record = new TH1D("h_et_transverse_record", ";#SigmaE_{T} [GeV]", calibnet, calibetbins);
    TH1D *h_nw_et_transverse_record = new TH1D("h_nw_et_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH2D* h_ue_pt_transverse_record = new TH2D("h_ue_pt_transverse_record","", calibnpt, calibptbins, calibnet, calibetbins);
    
    TH2D* h_measure_track_mbd_corr = new TH2D("h_measure_track_mbd_corr",";Track z-vertex [cm]; MBD z-vertex [cm]", 400, -200, 200, 400, -200, 200);
    TH2D* h_measure_mbd_track_corr = new TH2D("h_measure_mbd_track_corr",";Track z-vertex [cm]; MBD z-vertex [cm]", 400, -200, 200, 400, -200, 200);
    TH2D* h_measure_track_silicon_corr = new TH2D("h_measure_track_silicon_corr",";Silicon z-vertex [cm]; Track z-vertex [cm]", 400, -200, 200, 400, -200, 200);
    TH2D* h_measure_mbd_silicon_corr = new TH2D("h_measure_mbd_silicon_corr",";Silicon z-vertex [cm]; MBD z-vertex [cm]", 400, -200, 200, 400, -200, 200);
    TH2D* h_measure_silicon_track_corr = new TH2D("h_measure_silicon_track_corr",";Silicon z-vertex [cm]; Track z-vertex [cm]", 400, -200, 200, 400, -200, 200);
    TH2D* h_measure_silicon_mbd_corr = new TH2D("h_measure_silicon_mbd_corr",";Silicon z-vertex [cm]; MBD z-vertex [cm]", 400, -200, 200, 400, -200, 200);

    TH2D* h_total_measure_trackvertex = new TH2D("h_total_measure_trackvertex","", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_total_measure_mbdvertex = new TH2D("h_total_measure_mbdvertex","",calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_total_measure_siliconvertex = new TH2D("h_total_measure_siliconvertex","",calibnpt, calibptbins, calibnet, calibetbins);
    
    TH2D* h_pass_cut_measure_mbdvertex_trackvertex = new TH2D("h_pass_cut_measure_mbdvertex_trackvertex","", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_pass_cut_measure_trackvertex_mbdvertex = new TH2D("h_pass_cut_measure_trackvertex_mbdvertex","",calibnpt, calibptbins, calibnet, calibetbins);  
    TH2D* h_pass_cut_measure_mbdvertex_siliconvertex = new TH2D("h_pass_cut_measure_mbdvertex_siliconvertex","",calibnpt, calibptbins, calibnet, calibetbins); 
    TH2D* h_pass_cut_measure_trackvertex_siliconvertex = new TH2D("h_pass_cut_measure_trackvertex_siliconvertex","",calibnpt, calibptbins, calibnet, calibetbins); 
    TH2D* h_pass_cut_measure_siliconvertex_mbdvertex = new TH2D("h_pass_cut_measure_siliconvertex_mbdvertex","",calibnpt, calibptbins, calibnet, calibetbins); 
    TH2D* h_pass_cut_measure_siliconvertex_trackvertex = new TH2D("h_pass_cut_measure_siliconvertex_trackvertex","",calibnpt, calibptbins, calibnet, calibetbins);     
    
    TH1D* h_measure_track_mbd_diff = new TH1D("h_measure_track_mbd_diff",";MBD z-vertex - Track z-vertex [cm]", 400, -20.0, 20.0);
    TH1D* h_measure_mbd_track_diff = new TH1D("h_measure_mbd_track_diff",";MBD z-vertex - Track z-vertex [cm]", 400, -20.0, 20.0);
    TH1D* h_measure_track_silicon_diff = new TH1D("h_measure_track_silicon_diff",";Track z-vertex - Silicon z-vertex [cm]", 400, -20.0, 20.0);
    TH1D* h_measure_mbd_silicon_diff = new TH1D("h_measure_mbd_silicon_diff",";MBD z-vertex - Silicon z-vertex [cm]", 400, -20.0, 20.0);
    TH1D* h_measure_silicon_track_diff = new TH1D("h_measure_silicon_track_diff",";Track z-vertex - Silicon z-vertex [cm]", 400, -20.0, 20.0);
    TH1D* h_measure_silicon_mbd_diff = new TH1D("h_measure_silicon_mbd_diff",";MBD z-vertex - Silicon z-vertex [cm]", 400, -20.0, 20.0);

    TH1D* h_all_svtx_bco = new TH1D("h_all_svtx_bco","",600,-100,500);
    TH1D* h_all_silicon_bco = new TH1D("h_all_silicon_bco","",600,-100,500);
    
    TH1D* h_track_mult_no_pass_cut_measure_trackvertex_mbdvertex = new TH1D("h_track_mult_no_pass_cut_measure_trackvertex_mbdvertex","",100,0,100);
    TH1D* h_track_mult_pass_cut_measure_trackvertex_mbdvertex = new TH1D("h_track_mult_pass_cut_measure_trackvertex_mbdvertex","",100,0,100);
    TH1D* h_track_pt_spectra_no_pass_cut_measure_trackvertex_mbdvertex = new TH1D("h_track_pt_spectra_no_pass_cut_measure_trackvertex_mbdvertex","",100,0,100);
    TH1D* h_track_pt_spectra_pass_cut_measure_trackvertex_mbdvertex = new TH1D("h_track_pt_spectra_pass_cut_measure_trackvertex_mbdvertex","",100,0,100);
    TH1D* h_track_mult_no_pass_cut_measure_trackvertex_siliconvertex = new TH1D("h_track_mult_no_pass_cut_measure_trackvertex_siliconvertex","",100,0,100);
    TH1D* h_track_mult_pass_cut_measure_trackvertex_siliconvertex = new TH1D("h_track_mult_pass_cut_measure_trackvertex_siliconvertex","",100,0,100);
    TH1D* h_track_pt_spectra_no_pass_cut_measure_trackvertex_siliconvertex = new TH1D("h_track_pt_spectra_no_pass_cut_measure_trackvertex_siliconvertex","",100,0,100);
    TH1D* h_track_pt_spectra_pass_cut_measure_trackvertex_siliconvertex = new TH1D("h_track_pt_spectra_pass_cut_measure_trackvertex_siliconvertex","",100,0,100);
    TH1D* h_track_mult_no_pass_cut_measure_siliconvertex_mbdvertex = new TH1D("h_track_mult_no_pass_cut_measure_siliconvertex_mbdvertex","",100,0,100);
    TH1D* h_track_mult_pass_cut_measure_siliconvertex_mbdvertex = new TH1D("h_track_mult_pass_cut_measure_siliconvertex_mbdvertex","",100,0,100);
    TH1D* h_track_pt_spectra_no_pass_cut_measure_siliconvertex_mbdvertex = new TH1D("h_track_pt_spectra_no_pass_cut_measure_siliconvertex_mbdvertex","",100,0,100);
    TH1D* h_track_pt_spectra_pass_cut_measure_siliconvertex_mbdvertex = new TH1D("h_track_pt_spectra_pass_cut_measure_siliconvertex_mbdvertex","",100,0,100);

    TH1D* h_totalet_total_measure_trackvertex = new TH1D("h_totalet_total_measure_trackvertex","", calibtotalnet, calibtotaletbins);
    TH1D* h_totalet_total_measure_mbdvertex = new TH1D("h_totalet_total_measure_mbdvertex","", calibtotalnet, calibtotaletbins);
    TH1D* h_totalet_total_measure_siliconvertex = new TH1D("h_totalet_total_measure_siliconvertex","", calibtotalnet, calibtotaletbins);
    TH1D* h_totalet_pass_cut_measure_mbdvertex_trackvertex = new TH1D("h_totalet_pass_cut_measure_mbdvertex_trackvertex","", calibtotalnet, calibtotaletbins);
    TH1D* h_totalet_pass_cut_measure_trackvertex_mbdvertex = new TH1D("h_totalet_pass_cut_measure_trackvertex_mbdvertex","", calibtotalnet, calibtotaletbins);    
    TH1D* h_totalet_no_pass_cut_measure_mbdvertex_trackvertex = new TH1D("h_totalet_no_pass_cut_measure_mbdvertex_trackvertex","", calibtotalnet, calibtotaletbins);
    TH1D* h_totalet_no_pass_cut_measure_trackvertex_mbdvertex = new TH1D("h_totalet_no_pass_cut_measure_trackvertex_mbdvertex","", calibtotalnet, calibtotaletbins);    
    TH1D* h_totalet_pass_cut_measure_siliconvertex_trackvertex = new TH1D("h_totalet_pass_cut_measure_siliconvertex_trackvertex","", calibtotalnet, calibtotaletbins);
    TH1D* h_totalet_pass_cut_measure_trackvertex_siliconvertex = new TH1D("h_totalet_pass_cut_measure_trackvertex_siliconvertex","", calibtotalnet, calibtotaletbins);    
    TH1D* h_totalet_no_pass_cut_measure_siliconvertex_trackvertex = new TH1D("h_totalet_no_pass_cut_measure_siliconvertex_trackvertex","", calibtotalnet, calibtotaletbins);
    TH1D* h_totalet_no_pass_cut_measure_trackvertex_siliconvertex = new TH1D("h_totalet_no_pass_cut_measure_trackvertex_siliconvertex","", calibtotalnet, calibtotaletbins);    
    TH1D* h_totalet_pass_cut_measure_mbdvertex_siliconvertex = new TH1D("h_totalet_pass_cut_measure_mbdvertex_siliconvertex","", calibtotalnet, calibtotaletbins);
    TH1D* h_totalet_pass_cut_measure_siliconvertex_mbdvertex = new TH1D("h_totalet_pass_cut_measure_siliconvertex_mbdvertex","", calibtotalnet, calibtotaletbins);    
    TH1D* h_totalet_no_pass_cut_measure_mbdvertex_siliconvertex = new TH1D("h_totalet_no_pass_cut_measure_mbdvertex_siliconvertex","", calibtotalnet, calibtotaletbins);
    TH1D* h_totalet_no_pass_cut_measure_siliconvertex_mbdvertex = new TH1D("h_totalet_no_pass_cut_measure_siliconvertex_mbdvertex","", calibtotalnet, calibtotaletbins);    

    TH2D* h_totalet_calibjet_total_measure_trackvertex = new TH2D("h_totalet_calibjet_total_measure_trackvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);
    TH2D* h_totalet_calibjet_total_measure_mbdvertex = new TH2D("h_totalet_calibjet_total_measure_mbdvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);
    TH2D* h_totalet_calibjet_total_measure_siliconvertex = new TH2D("h_totalet_calibjet_total_measure_siliconvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);
    TH2D* h_totalet_calibjet_pass_cut_measure_mbdvertex_trackvertex = new TH2D("h_totalet_calibjet_pass_cut_measure_mbdvertex_trackvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);
    TH2D* h_totalet_calibjet_pass_cut_measure_trackvertex_mbdvertex = new TH2D("h_totalet_calibjet_pass_cut_measure_trackvertex_mbdvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);    
    TH2D* h_totalet_calibjet_no_pass_cut_measure_mbdvertex_trackvertex = new TH2D("h_totalet_calibjet_no_pass_cut_measure_mbdvertex_trackvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);
    TH2D* h_totalet_calibjet_no_pass_cut_measure_trackvertex_mbdvertex = new TH2D("h_totalet_calibjet_no_pass_cut_measure_trackvertex_mbdvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);    
    TH2D* h_totalet_calibjet_pass_cut_measure_siliconvertex_trackvertex = new TH2D("h_totalet_calibjet_pass_cut_measure_siliconvertex_trackvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);
    TH2D* h_totalet_calibjet_pass_cut_measure_trackvertex_siliconvertex = new TH2D("h_totalet_calibjet_pass_cut_measure_trackvertex_siliconvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);    
    TH2D* h_totalet_calibjet_no_pass_cut_measure_siliconvertex_trackvertex = new TH2D("h_totalet_calibjet_no_pass_cut_measure_siliconvertex_trackvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);
    TH2D* h_totalet_calibjet_no_pass_cut_measure_trackvertex_siliconvertex = new TH2D("h_totalet_calibjet_no_pass_cut_measure_trackvertex_siliconvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);    
    TH2D* h_totalet_calibjet_pass_cut_measure_mbdvertex_siliconvertex = new TH2D("h_totalet_calibjet_pass_cut_measure_mbdvertex_siliconvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);
    TH2D* h_totalet_calibjet_pass_cut_measure_siliconvertex_mbdvertex = new TH2D("h_totalet_calibjet_pass_cut_measure_siliconvertex_mbdvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);    
    TH2D* h_totalet_calibjet_no_pass_cut_measure_mbdvertex_siliconvertex = new TH2D("h_totalet_calibjet_no_pass_cut_measure_mbdvertex_siliconvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);
    TH2D* h_totalet_calibjet_no_pass_cut_measure_siliconvertex_mbdvertex = new TH2D("h_totalet_calibjet_no_pass_cut_measure_siliconvertex_mbdvertex","", calibnpt, calibptbins, calibtotalnet, calibtotaletbins);    


    ////////// Event Loop //////////
    std::cout << "Data analysis started." << std::endl;
    Long64_t nEntries = chain.GetEntries();
    std::cout << "Total number of events: " << nEntries << std::endl;
    
    // Event variables setup.
    bool reco_bkg_cut = false;
    float calibjet_pt_dijet, calibjet_eta_dijet, calibjet_phi_dijet;
    bool calibjet_matched_dijet, qa_matched;
    
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 2000; ++entry) {
        if (entry % 1000 == 0) cout << "event " << entry << endl;
        chain.GetEntry(entry);

        if (!check_bad_trigger(gl1_trigger_vector_scaled)) continue;
        if (fabs(trackvertex) > zvertex_cut && fabs(zvertex) > zvertex_cut) { continue;}

        //////////////////////////// SETUP JET VARIABLES FOR UNFOLDING ////////////////////////////

        // indices to find leading and subleading jets 
        int ind_lead = -1; int ind_sub = -1;
        float lead_e = 0;
        reco_bkg_cut = false;

        std::vector<float> recoe_new, recopt_new, recoeta_new, recophi_new, recoemcal_new, recoihcal_new, recoohcal_new;

        int nJetReq = 1;
        if (bkg_cut == "dijet") { nJetReq = 2; }
        if (unsubjet_pt->size() < nJetReq) { continue; }
        
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
                float vertex = fabs(zvertex) < 60.0 ? zvertex : 0.0;
                if (!check_bad_jet_eta(unsubjet_eta->at(i), vertex, jet_radius)) {
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

        get_calibjet(calibjet_pt_dijet, calibjet_eta_dijet, calibjet_phi_dijet, lead.Pt(), lead.Eta(), lead.Phi(), reco_bkg_cut, f_corr, 1, 0.1);    
        //std::cout << "Good reco lead: pt " << calibjet_pt_dijet << " eta " << calibjet_eta_dijet << " phi " << calibjet_phi_dijet << " bkg_cut " << reco_bkg_cut << std::endl;
    
        //////////////////////////// SETUP UE VARIABLES FOR UNFOLDING //////////////////////////////

        // find reco ET information 
        float et_transverse = 0;
        float total_et = 0;
        float radius_IH = 127.503;
        for (int i = 0; i < clsmult; i++) {
            float eta = cluster_eta[i];
            total_et += cluster_e[i]/cosh(cluster_eta[i]);
            float dphi = get_dphi(lead.Phi(),cluster_phi[i]);
            if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += cluster_e[i]/cosh(cluster_eta[i]); }
        }

        //////////////////////////// RECORD QA PLOTS FOR ALL EVENTS IN UNFOLDING PROCEDURE ////////////////////////////
        if (reco_bkg_cut) {
            h_lead_spectra_record->Fill(lead.Pt());
            if (calibjet_pt_dijet > calibptbins[0]) {
                h_nw_et_transverse_record->Fill(et_transverse);
                h_et_transverse_record->Fill(et_transverse);
                h_ue_pt_transverse_record->Fill(lead.Pt(),et_transverse);
            }
        }

        for (int i = 0; i < svtxBcoVector->size(); i++) {
            if (svtxBcoVector->at(i) > 500) {
                std::cout << "Large BCO value: " << svtxBcoVector->at(i) << " zvertex: " << svtxVector->at(i) << std::endl;
            }
            h_all_svtx_bco->Fill(svtxBcoVector->at(i));
            h_all_svtx_bco->Fill(svtxBcoVector->at(i), svtxVector->at(i));
        }

        for (int i = 0; i < siliconBcoVector->size(); i++) {
            //if (siliconBcoVector->at(i) > 500) {
            //    std::cout << "Large BCO value: " << siliconBcoVector->at(i) << " zvertex: " << siliconVector->at(i) << std::endl;
            //}
            h_all_silicon_bco->Fill(siliconBcoVector->at(i));
            h_all_silicon_bco->Fill(siliconBcoVector->at(i), siliconVector->at(i));
        }

        //std::cout << "Vertex info: track " << trackvertex << " mbd " << zvertex << std::endl;
        //std::cout << std::endl;
        // Create efficiency plots here
        if (calibjet_pt_dijet >= calibptbins[0] && calibjet_pt_dijet < calibptbins[calibnpt] && reco_bkg_cut) {
            if (fabs(trackvertex) < zvertex_cut) {
                h_total_measure_trackvertex->Fill(calibjet_pt_dijet, et_transverse);
                h_totalet_total_measure_trackvertex->Fill(total_et);
                h_totalet_calibjet_total_measure_trackvertex->Fill(calibjet_pt_dijet, total_et);
                if (fabs(zvertex) < wide_zvertex_cut) {
                    h_measure_track_mbd_corr->Fill(trackvertex, zvertex);
                    h_measure_track_mbd_diff->Fill(zvertex - trackvertex);
                    h_pass_cut_measure_mbdvertex_trackvertex->Fill(calibjet_pt_dijet, et_transverse);
                    h_totalet_pass_cut_measure_mbdvertex_trackvertex->Fill(total_et);
                    h_totalet_calibjet_pass_cut_measure_mbdvertex_trackvertex->Fill(calibjet_pt_dijet, total_et);
                } else {
                    //std::cout << "reco jet pt " << calibjet_pt_dijet << " eta " << calibjet_eta_dijet << " phi " << calibjet_phi_dijet << " vz_mbd " << zvertex << std::endl;
                    //for (int i = 0; i < clsmult; i++) {
                    //    std::cout << "\t cluster e " << cluster_e[i] << " eta " << cluster_eta[i] << " phi " << cluster_phi[i] << std::endl;
                    //}
                    h_totalet_no_pass_cut_measure_mbdvertex_trackvertex->Fill(total_et);
                    h_totalet_calibjet_no_pass_cut_measure_mbdvertex_trackvertex->Fill(calibjet_pt_dijet, total_et);
                }
                if (fabs(siliconvertex) < wide_zvertex_cut) {
                    h_measure_track_silicon_corr->Fill(siliconvertex, trackvertex);
                    h_measure_track_silicon_diff->Fill(trackvertex - siliconvertex);
                    h_pass_cut_measure_siliconvertex_trackvertex->Fill(calibjet_pt_dijet, et_transverse);
                    h_totalet_pass_cut_measure_siliconvertex_trackvertex->Fill(total_et);
                    h_totalet_calibjet_pass_cut_measure_siliconvertex_trackvertex->Fill(calibjet_pt_dijet, total_et);
                } else {
                    //std::cout << "reco jet pt " << calibjet_pt_dijet << " eta " << calibjet_eta_dijet << " phi " << calibjet_phi_dijet << " vz_silicon " << zvertex << std::endl;
                    //for (int i = 0; i < clsmult; i++) {
                    //    std::cout << "\t cluster e " << cluster_e[i] << " eta " << cluster_eta[i] << " phi " << cluster_phi[i] << std::endl;
                    //}
                    h_totalet_no_pass_cut_measure_siliconvertex_trackvertex->Fill(total_et);
                    h_totalet_calibjet_no_pass_cut_measure_siliconvertex_trackvertex->Fill(calibjet_pt_dijet, total_et);
                }   
            }
            if (fabs(siliconvertex) < zvertex_cut) {
                h_total_measure_siliconvertex->Fill(calibjet_pt_dijet, et_transverse);
                h_totalet_total_measure_siliconvertex->Fill(total_et);
                h_totalet_calibjet_total_measure_siliconvertex->Fill(calibjet_pt_dijet, total_et);
                if (fabs(zvertex) < wide_zvertex_cut) {
                    h_measure_silicon_mbd_corr->Fill(siliconvertex, zvertex);
                    h_measure_silicon_mbd_diff->Fill(zvertex - siliconvertex);
                    h_pass_cut_measure_mbdvertex_siliconvertex->Fill(calibjet_pt_dijet, et_transverse);
                    h_totalet_pass_cut_measure_mbdvertex_siliconvertex->Fill(total_et);
                    h_totalet_calibjet_pass_cut_measure_mbdvertex_siliconvertex->Fill(calibjet_pt_dijet, total_et);
                } else {
                    //std::cout << "reco jet pt " << calibjet_pt_dijet << " eta " << calibjet_eta_dijet << " phi " << calibjet_phi_dijet << " vz_mbd " << zvertex << std::endl;
                    //for (int i = 0; i < clsmult; i++) {
                    //    std::cout << "\t cluster e " << cluster_e[i] << " eta " << cluster_eta[i] << " phi " << cluster_phi[i] << std::endl;
                    //}
                    h_totalet_no_pass_cut_measure_mbdvertex_siliconvertex->Fill(total_et);
                    h_totalet_calibjet_no_pass_cut_measure_mbdvertex_siliconvertex->Fill(calibjet_pt_dijet, total_et);
                }
                if (fabs(trackvertex) < wide_zvertex_cut) {
                    h_measure_silicon_track_corr->Fill(siliconvertex, trackvertex);
                    h_measure_silicon_track_diff->Fill(trackvertex - siliconvertex);
                    h_pass_cut_measure_trackvertex_siliconvertex->Fill(calibjet_pt_dijet, et_transverse);
                    h_totalet_pass_cut_measure_trackvertex_siliconvertex->Fill(total_et);
                    h_totalet_calibjet_pass_cut_measure_trackvertex_siliconvertex->Fill(calibjet_pt_dijet, total_et);
                    int ntracks = 0;
                    for (int i = 0; i < trkmult; i++) {
                        if (tr_nintt[i] > 1 && tr_ntpc[i] > 25 && tr_crossing[i] == 0) {
                            h_track_pt_spectra_pass_cut_measure_trackvertex_siliconvertex->Fill(tr_pt[i]);
                            ntracks++;
                        }
                    }
                    h_track_mult_pass_cut_measure_trackvertex_siliconvertex->Fill(ntracks);
                } else {
                    int ntracks = 0;
                    for (int i = 0; i < trkmult; i++) {
                        if (tr_nintt[i] > 1 && tr_ntpc[i] > 25 && tr_crossing[i] == 0) {
                            h_track_pt_spectra_no_pass_cut_measure_trackvertex_siliconvertex->Fill(tr_pt[i]);
                            ntracks++;
                        }
                    }
                    h_track_mult_no_pass_cut_measure_trackvertex_siliconvertex->Fill(ntracks);
                    h_totalet_no_pass_cut_measure_trackvertex_siliconvertex->Fill(total_et);
                    h_totalet_calibjet_no_pass_cut_measure_trackvertex_siliconvertex->Fill(calibjet_pt_dijet, total_et);
                }
            }
            if (fabs(zvertex) < zvertex_cut) {
                h_total_measure_mbdvertex->Fill(calibjet_pt_dijet, et_transverse);
                h_totalet_total_measure_mbdvertex->Fill(total_et);
                h_totalet_calibjet_total_measure_mbdvertex->Fill(calibjet_pt_dijet, total_et);
                if (fabs(trackvertex) < wide_zvertex_cut) {
                    h_measure_mbd_track_corr->Fill(trackvertex, zvertex);
                    h_measure_mbd_track_diff->Fill(zvertex - trackvertex);
                    h_pass_cut_measure_trackvertex_mbdvertex->Fill(calibjet_pt_dijet, et_transverse);
                    h_totalet_pass_cut_measure_trackvertex_mbdvertex->Fill(total_et);
                    h_totalet_calibjet_pass_cut_measure_trackvertex_mbdvertex->Fill(calibjet_pt_dijet, total_et);
                    int ntracks = 0;
                    for (int i = 0; i < trkmult; i++) {
                        if (tr_nintt[i] > 1 && tr_ntpc[i] > 25 && tr_crossing[i] == 0) {
                            h_track_pt_spectra_pass_cut_measure_trackvertex_mbdvertex->Fill(tr_pt[i]);
                            ntracks++;
                        }
                    }
                    h_track_mult_pass_cut_measure_trackvertex_mbdvertex->Fill(ntracks);
                } else {
                    int ntracks = 0;
                    for (int i = 0; i < trkmult; i++) {
                        if (tr_nintt[i] > 1 && tr_ntpc[i] > 25 && tr_crossing[i] == 0) {
                            h_track_pt_spectra_no_pass_cut_measure_trackvertex_mbdvertex->Fill(tr_pt[i]);
                            ntracks++;
                        }
                    }
                    h_track_mult_no_pass_cut_measure_trackvertex_mbdvertex->Fill(ntracks);
                    h_totalet_no_pass_cut_measure_trackvertex_mbdvertex->Fill(total_et);
                    h_totalet_calibjet_no_pass_cut_measure_trackvertex_mbdvertex->Fill(calibjet_pt_dijet, total_et);
                }
                if (fabs(siliconvertex) < wide_zvertex_cut) {
                    h_measure_mbd_silicon_corr->Fill(siliconvertex, zvertex);
                    h_measure_mbd_silicon_diff->Fill(zvertex - siliconvertex);
                    h_pass_cut_measure_siliconvertex_mbdvertex->Fill(calibjet_pt_dijet, et_transverse);
                    h_totalet_pass_cut_measure_siliconvertex_mbdvertex->Fill(total_et);
                    h_totalet_calibjet_pass_cut_measure_siliconvertex_mbdvertex->Fill(calibjet_pt_dijet, total_et);
                    int ntracks = 0;
                    for (int i = 0; i < trkmult; i++) {
                        if (tr_nintt[i] > 1 && tr_ntpc[i] > 25 && tr_crossing[i] == 0) {
                            h_track_pt_spectra_pass_cut_measure_siliconvertex_mbdvertex->Fill(tr_pt[i]);
                            ntracks++;
                        }
                    }
                    h_track_mult_pass_cut_measure_siliconvertex_mbdvertex->Fill(ntracks);
                } else {
                    int ntracks = 0;
                    for (int i = 0; i < trkmult; i++) {
                        if (tr_nintt[i] > 1 && tr_ntpc[i] > 25 && tr_crossing[i] == 0) {
                            h_track_pt_spectra_no_pass_cut_measure_siliconvertex_mbdvertex->Fill(tr_pt[i]);
                            ntracks++;
                        }
                    }
                    h_track_mult_no_pass_cut_measure_siliconvertex_mbdvertex->Fill(ntracks);
                    h_totalet_no_pass_cut_measure_siliconvertex_mbdvertex->Fill(total_et);
                    h_totalet_calibjet_no_pass_cut_measure_siliconvertex_mbdvertex->Fill(calibjet_pt_dijet, total_et);
                }
            }
        }
    }

    std::cout << "Writing histograms..." << std::endl;
    f_out->cd();

    h_lead_spectra_record->Write(); h_et_transverse_record->Write(); h_nw_et_transverse_record->Write(); h_ue_pt_transverse_record->Write();
    h_measure_track_mbd_corr->Write(); h_measure_mbd_track_corr->Write(); h_measure_track_silicon_corr->Write(); h_measure_mbd_silicon_corr->Write(); h_measure_silicon_track_corr->Write(); h_measure_silicon_mbd_corr->Write();
    h_total_measure_trackvertex->Write(); h_total_measure_mbdvertex->Write(); h_total_measure_siliconvertex->Write();
    h_pass_cut_measure_mbdvertex_trackvertex->Write(); h_pass_cut_measure_trackvertex_mbdvertex->Write(); h_pass_cut_measure_mbdvertex_siliconvertex->Write(); h_pass_cut_measure_trackvertex_siliconvertex->Write(); h_pass_cut_measure_siliconvertex_mbdvertex->Write(); h_pass_cut_measure_siliconvertex_trackvertex->Write();      
    h_measure_track_mbd_diff->Write(); h_measure_mbd_track_diff->Write(); h_measure_track_silicon_diff->Write(); h_measure_mbd_silicon_diff->Write(); h_measure_silicon_track_diff->Write(); h_measure_silicon_mbd_diff->Write();
    h_all_svtx_bco->Write(); h_all_silicon_bco->Write();  
    h_track_mult_no_pass_cut_measure_trackvertex_mbdvertex->Write(); h_track_mult_pass_cut_measure_trackvertex_mbdvertex->Write(); h_track_pt_spectra_no_pass_cut_measure_trackvertex_mbdvertex->Write(); h_track_pt_spectra_pass_cut_measure_trackvertex_mbdvertex->Write();  
    h_track_mult_no_pass_cut_measure_trackvertex_siliconvertex->Write(); h_track_mult_pass_cut_measure_trackvertex_siliconvertex->Write(); h_track_pt_spectra_no_pass_cut_measure_trackvertex_siliconvertex->Write(); h_track_pt_spectra_pass_cut_measure_trackvertex_siliconvertex->Write();  
    h_track_mult_no_pass_cut_measure_siliconvertex_mbdvertex->Write(); h_track_mult_pass_cut_measure_siliconvertex_mbdvertex->Write(); h_track_pt_spectra_no_pass_cut_measure_siliconvertex_mbdvertex->Write(); h_track_pt_spectra_pass_cut_measure_siliconvertex_mbdvertex->Write();  
    h_totalet_total_measure_trackvertex->Write(); h_totalet_total_measure_mbdvertex->Write(); h_totalet_total_measure_siliconvertex->Write(); 
    h_totalet_pass_cut_measure_mbdvertex_trackvertex->Write(); h_totalet_pass_cut_measure_trackvertex_mbdvertex->Write(); h_totalet_no_pass_cut_measure_mbdvertex_trackvertex->Write(); h_totalet_no_pass_cut_measure_trackvertex_mbdvertex->Write();     
    h_totalet_pass_cut_measure_siliconvertex_trackvertex->Write(); h_totalet_pass_cut_measure_trackvertex_siliconvertex->Write(); h_totalet_no_pass_cut_measure_siliconvertex_trackvertex->Write(); h_totalet_no_pass_cut_measure_trackvertex_siliconvertex->Write();     
    h_totalet_pass_cut_measure_mbdvertex_siliconvertex->Write(); h_totalet_pass_cut_measure_siliconvertex_mbdvertex->Write(); h_totalet_no_pass_cut_measure_mbdvertex_siliconvertex->Write(); h_totalet_no_pass_cut_measure_siliconvertex_mbdvertex->Write();     
    h_totalet_calibjet_total_measure_trackvertex->Write(); h_totalet_calibjet_total_measure_mbdvertex->Write(); h_totalet_calibjet_total_measure_siliconvertex->Write(); 
    h_totalet_calibjet_pass_cut_measure_mbdvertex_trackvertex->Write(); h_totalet_calibjet_pass_cut_measure_trackvertex_mbdvertex->Write(); h_totalet_calibjet_no_pass_cut_measure_mbdvertex_trackvertex->Write(); h_totalet_calibjet_no_pass_cut_measure_trackvertex_mbdvertex->Write();     
    h_totalet_calibjet_pass_cut_measure_siliconvertex_trackvertex->Write(); h_totalet_calibjet_pass_cut_measure_trackvertex_siliconvertex->Write(); h_totalet_calibjet_no_pass_cut_measure_siliconvertex_trackvertex->Write(); h_totalet_calibjet_no_pass_cut_measure_trackvertex_siliconvertex->Write();     
    h_totalet_calibjet_pass_cut_measure_mbdvertex_siliconvertex->Write(); h_totalet_calibjet_pass_cut_measure_siliconvertex_mbdvertex->Write(); h_totalet_calibjet_no_pass_cut_measure_mbdvertex_siliconvertex->Write(); h_totalet_calibjet_no_pass_cut_measure_siliconvertex_mbdvertex->Write();     
    f_out->Close();
    std::cout << "All done!" << std::endl;

}

////////////////////////////////////////// Functions //////////////////////////////////////////
bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled) {
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
  float dijet_min_phi = 7*TMath::Pi()/8.;
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