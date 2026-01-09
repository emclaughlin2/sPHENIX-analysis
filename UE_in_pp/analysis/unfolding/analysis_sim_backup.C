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
#include "unfold_Def.h"
#include "RooUnfold.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

using namespace std;

void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_pt);
void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et);
bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi);
void filter_jet(bool& jet_filter, float jet_e, float jet_pt, float jet_eta, float zvertex, float jet_radius);
void get_jet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius);
void get_calibjet(float& calibjet_pt, float& calibjet_eta, float& calibjet_phi, bool jet_filter, float jet_pt, float jet_eta, float jet_phi, bool jet_background, TF1* f_corr, float jes_para, float jer_para);
void get_truthjet(float& goodtruthjet_pt, float& goodtruthjet_eta, float& goodtruthjet_phi, bool truthjet_filter, float jet_pt, float jet_eta, float jet_phi);
void match_meas_truth(float meas_eta, float meas_phi, bool& matched, float truth_eta, float truth_phi, float jet_radius);
void fill_response_matrix(TH2D*& h_truth, TH2D*& h_meas, RooUnfoldResponse*& h_resp, TH2D*& h_fake, TH2D*& h_miss, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale);

TF1* f_jer = new TF1("f_jer", "sqrt( 0.095077098*0.095077098 + (0.63134847*0.63134847/x) + (2.1664610*2.1664610/(x*x)) )", 0, 100);
TRandom3 randGen(1234);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                               //
//  Note: 5.13.25                                                                                                                //
//  This macro creates response matrices for 2D unfolding of leading jet pT and transverse region energy density                 //
//                                                                                                                               //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void analysis_sim(std::string runtype = "mb", int start_seg = 0, int end_seg = 200, bool clusters = true, bool emcal_clusters = false)  {
    
    ////////// General Set up //////////
    double weight_scale = 1.0, truthjet_pt_min = 0, truthjet_pt_max = 3000;
    if (runtype == "mb") {
        weight_scale = MB_scale;
        truthjet_pt_min = 0;
        truthjet_pt_max = 14;
    } else if (runtype == "jet10") {
        weight_scale = Jet10GeV_scale;
        truthjet_pt_min = 14;
        truthjet_pt_max = 22;
    } else if (runtype == "jet20") {
        weight_scale = Jet20GeV_scale;
        truthjet_pt_min = 22;
        truthjet_pt_max = 35;
    } else if (runtype == "jet30") {
        weight_scale = Jet30GeV_scale;
        truthjet_pt_min = 35;
        truthjet_pt_max = 3000;
    } else {
        std::cout << "Unknown runtype" << std::endl;
        return;
    }

    ////////// Files //////////
    TFile *f_out = new TFile(Form("analysis_sim_output/output_sim_%s_%d_%d.root", runtype.c_str(), start_seg, end_seg), "RECREATE");
    TChain chain("T");
    for (int i = start_seg; i < end_seg; ++i) {
        chain.Add(Form("/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run22_%s_3sigma_output_%d.root", runtype.c_str(), i));
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
    if (!clusters) {
        chain.SetBranchStatus("emcaln", 1); chain.SetBranchAddress("emcaln",&emcaln);
        chain.SetBranchStatus("emcale", 1); chain.SetBranchAddress("emcale",emcale);
        chain.SetBranchStatus("emcaleta", 1); chain.SetBranchAddress("emcaleta",emcaleta);
        chain.SetBranchStatus("emcalphi", 1); chain.SetBranchAddress("emcalphi",emcalphi);
        chain.SetBranchStatus("ihcaln", 1); chain.SetBranchAddress("ihcaln",&ihcaln);
        chain.SetBranchStatus("ihcale", 1); chain.SetBranchAddress("ihcale",ihcale);
        chain.SetBranchStatus("ihcaleta", 1); chain.SetBranchAddress("ihcaleta",ihcaleta);
        chain.SetBranchStatus("ihcalphi", 1); chain.SetBranchAddress("ihcalphi",ihcalphi);
        chain.SetBranchStatus("ohcaln", 1); chain.SetBranchAddress("ohcaln",&ohcaln);
        chain.SetBranchStatus("ohcale", 1); chain.SetBranchAddress("ohcale",ohcale);
        chain.SetBranchStatus("ohcaleta", 1); chain.SetBranchAddress("ohcaleta",ohcaleta);
        chain.SetBranchStatus("ohcalphi", 1); chain.SetBranchAddress("ohcalphi",ohcalphi);
    } else if (clusters && !emcal_clusters) {
        chain.SetBranchStatus("clsmult", 1); chain.SetBranchAddress("clsmult",&clsmult);
        chain.SetBranchStatus("cluster_e", 1); chain.SetBranchAddress("cluster_e",cluster_e);
        chain.SetBranchStatus("cluster_eta", 1); chain.SetBranchAddress("cluster_eta",cluster_eta);
        chain.SetBranchStatus("cluster_phi", 1); chain.SetBranchAddress("cluster_phi",cluster_phi);
    } else {
        chain.SetBranchStatus("emcal_clsmult", 1); chain.SetBranchAddress("emcal_clsmult",&clsmult);
        chain.SetBranchStatus("emcal_cluster_e", 1); chain.SetBranchAddress("emcal_cluster_e",cluster_e);
        chain.SetBranchStatus("emcal_cluster_eta", 1); chain.SetBranchAddress("emcal_cluster_eta",cluster_eta);
        chain.SetBranchStatus("emcal_cluster_phi", 1); chain.SetBranchAddress("emcal_cluster_phi",cluster_phi);
    }

    ////////// JES func //////////
    TFile *corrFile = new TFile("JES_Calib_Default.root", "READ");
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
    TH1D* h_zvertex = new TH1D("h_zvertex", ";Z-vertex [cm]", 4000, -2000, 2000);
    TH1D *h_recojet_pt_record = new TH1D("h_recojet_pt_record", ";p_{T} [GeV]", 1000, 0, 100);
    TH1D *h_recojet_pt_record_dijet = new TH1D("h_recojet_pt_record_dijet", ";p_{T} [GeV]", 1000, 0, 100);
    TH1D *h_recojet_pt_record_frac = new TH1D("h_recojet_pt_record_frac", ";p_{T} [GeV]", 1000, 0, 100);
    TH1D *h_reco_et_transverse_record_dijet = new TH1D("h_reco_et_transverse_record_dijet", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH1D* h_truth_et_transverse_record_dijet = new TH1D("h_truth_et_transverse_record_dijet", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH1D *h_reco_et_transverse_record_dijet_weighted = new TH1D("h_reco_et_transverse_record_dijet_weighted", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH1D* h_truth_et_transverse_record_dijet_weighted = new TH1D("h_truth_et_transverse_record_dijet_weighted", ";#SigmaE_{T} [GeV]", 7000, -20, 50);

    TH2D* h_truth_calib_dijet = new TH2D("h_truth_calib_dijet", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_dijet = new TH2D("h_measure_calib_dijet", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_dijet = new TH2D("h_fake_calib_dijet", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_dijet = new TH2D("h_miss_calib_dijet", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_dijet = new RooUnfoldResponse("h_respmatrix_calib_dijet",""); h_respmatrix_calib_dijet->Setup(h_measure_calib_dijet, h_truth_calib_dijet);

    TH2D* h_truth_calib_dijet_jesdown = new TH2D("h_truth_calib_dijet_jesdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_dijet_jesdown = new TH2D("h_measure_calib_dijet_jesdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_dijet_jesdown = new TH2D("h_fake_calib_dijet_jesdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_dijet_jesdown = new TH2D("h_miss_calib_dijet_jesdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_dijet_jesdown = new RooUnfoldResponse("h_respmatrix_calib_dijet_jesdown",""); h_respmatrix_calib_dijet_jesdown->Setup(h_measure_calib_dijet_jesdown, h_truth_calib_dijet_jesdown);

    TH2D* h_truth_calib_dijet_jesup = new TH2D("h_truth_calib_dijet_jesup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_dijet_jesup = new TH2D("h_measure_calib_dijet_jesup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_dijet_jesup = new TH2D("h_fake_calib_dijet_jesup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_dijet_jesup = new TH2D("h_miss_calib_dijet_jesup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_dijet_jesup = new RooUnfoldResponse("h_respmatrix_calib_dijet_jesup",""); h_respmatrix_calib_dijet_jesup->Setup(h_measure_calib_dijet_jesup, h_truth_calib_dijet_jesup);

    TH2D* h_truth_calib_dijet_jerdown = new TH2D("h_truth_calib_dijet_jerdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_dijet_jerdown = new TH2D("h_measure_calib_dijet_jerdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_dijet_jerdown = new TH2D("h_fake_calib_dijet_jerdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_dijet_jerdown = new TH2D("h_miss_calib_dijet_jerdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_dijet_jerdown = new RooUnfoldResponse("h_respmatrix_calib_dijet_jerdown",""); h_respmatrix_calib_dijet_jerdown->Setup(h_measure_calib_dijet_jerdown, h_truth_calib_dijet_jerdown);

    TH2D* h_truth_calib_dijet_jerup = new TH2D("h_truth_calib_dijet_jerup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_dijet_jerup = new TH2D("h_measure_calib_dijet_jerup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_dijet_jerup = new TH2D("h_fake_calib_dijet_jerup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_dijet_jerup = new TH2D("h_miss_calib_dijet_jerup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_dijet_jerup = new RooUnfoldResponse("h_respmatrix_calib_dijet_jerup",""); h_respmatrix_calib_dijet_jerup->Setup(h_measure_calib_dijet_jerup, h_truth_calib_dijet_jerup);

    TH2D* h_truth_calib_frac = new TH2D("h_truth_calib_frac", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_frac = new TH2D("h_measure_calib_frac", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_frac = new TH2D("h_fake_calib_frac", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_frac = new TH2D("h_miss_calib_frac", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_frac = new RooUnfoldResponse("h_respmatrix_calib_frac",""); h_respmatrix_calib_frac->Setup(h_measure_calib_frac, h_truth_calib_frac);

    TH2D* h_truth_calib_frac_jesdown = new TH2D("h_truth_calib_frac_jesdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_frac_jesdown = new TH2D("h_measure_calib_frac_jesdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_frac_jesdown = new TH2D("h_fake_calib_frac_jesdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_frac_jesdown = new TH2D("h_miss_calib_frac_jesdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_frac_jesdown = new RooUnfoldResponse("h_respmatrix_calib_frac_jesdown",""); h_respmatrix_calib_frac_jesdown->Setup(h_measure_calib_frac_jesdown, h_truth_calib_frac_jesdown);

    TH2D* h_truth_calib_frac_jesup = new TH2D("h_truth_calib_frac_jesup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_frac_jesup = new TH2D("h_measure_calib_frac_jesup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_frac_jesup = new TH2D("h_fake_calib_frac_jesup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_frac_jesup = new TH2D("h_miss_calib_frac_jesup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_frac_jesup = new RooUnfoldResponse("h_respmatrix_calib_frac_jesup",""); h_respmatrix_calib_frac_jesup->Setup(h_measure_calib_frac_jesup, h_truth_calib_frac_jesup);

    TH2D* h_truth_calib_frac_jerdown = new TH2D("h_truth_calib_frac_jerdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_frac_jerdown = new TH2D("h_measure_calib_frac_jerdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_frac_jerdown = new TH2D("h_fake_calib_frac_jerdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_frac_jerdown = new TH2D("h_miss_calib_frac_jerdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_frac_jerdown = new RooUnfoldResponse("h_respmatrix_calib_frac_jerdown",""); h_respmatrix_calib_frac_jerdown->Setup(h_measure_calib_frac_jerdown, h_truth_calib_frac_jerdown);

    TH2D* h_truth_calib_frac_jerup = new TH2D("h_truth_calib_frac_jerup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_frac_jerup = new TH2D("h_measure_calib_frac_jerup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_frac_jerup = new TH2D("h_fake_calib_frac_jerup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_frac_jerup = new TH2D("h_miss_calib_frac_jerup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_frac_jerup = new RooUnfoldResponse("h_respmatrix_calib_frac_jerup",""); h_respmatrix_calib_frac_jerup->Setup(h_measure_calib_frac_jerup, h_truth_calib_frac_jerup);

    TH2D* h_truth_calib_frac_half1 = new TH2D("h_truth_calib_frac_half1", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_frac_half1 = new TH2D("h_measure_calib_frac_half1", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_frac_half1 = new TH2D("h_fake_calib_frac_half1", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_frac_half1 = new TH2D("h_miss_calib_frac_half1", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_frac_half1 = new RooUnfoldResponse("h_respmatrix_calib_frac_half1",""); h_respmatrix_calib_frac_half1->Setup(h_measure_calib_frac_half1, h_truth_calib_frac_half1);

    TH2D* h_truth_calib_frac_half2 = new TH2D("h_truth_calib_frac_half2", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_calib_frac_half2 = new TH2D("h_measure_calib_frac_half2", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_fake_calib_frac_half2 = new TH2D("h_fake_calib_frac_half2", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_miss_calib_frac_half2 = new TH2D("h_miss_calib_frac_half2", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, truthptbins, calibnet, calibetbins);
    RooUnfoldResponse* h_respmatrix_calib_frac_half2 = new RooUnfoldResponse("h_respmatrix_calib_frac_half2",""); h_respmatrix_calib_frac_half2->Setup(h_measure_calib_frac_half2, h_truth_calib_frac_half2);

    ////////// Event Loop //////////
    std::cout << "Data analysis started." << std::endl;
    Long64_t nEntries = chain.GetEntries();
    std::cout << "Total number of events: " << nEntries << std::endl;
    
    // Event variables setup.
    bool background_dijet = false; bool background_frac = false;
    bool jet_filter = false; bool truthjet_filter = false;
    float goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi;
    float calibjet_pt_dijet, calibjet_eta_dijet, calibjet_phi_dijet;
    float calibjet_pt_dijet_jesdown, calibjet_eta_dijet_jesdown, calibjet_phi_dijet_jesdown;
    float calibjet_pt_dijet_jesup, calibjet_eta_dijet_jesup, calibjet_phi_dijet_jesup;
    float calibjet_pt_dijet_jerdown, calibjet_eta_dijet_jerdown, calibjet_phi_dijet_jerdown;
    float calibjet_pt_dijet_jerup, calibjet_eta_dijet_jerup, calibjet_phi_dijet_jerup;
    float calibjet_pt_frac, calibjet_eta_frac, calibjet_phi_frac;
    float calibjet_pt_frac_jesdown, calibjet_eta_frac_jesdown, calibjet_phi_frac_jesdown;
    float calibjet_pt_frac_jesup, calibjet_eta_frac_jesup, calibjet_phi_frac_jesup;
    float calibjet_pt_frac_jerdown, calibjet_eta_frac_jerdown, calibjet_phi_frac_jerdown;
    float calibjet_pt_frac_jerup, calibjet_eta_frac_jerup, calibjet_phi_frac_jerup;
    bool calibjet_matched_dijet, calibjet_matched_dijet_jesdown, calibjet_matched_dijet_jesup, calibjet_matched_dijet_jerdown, calibjet_matched_dijet_jerup;
    bool calibjet_matched_frac, calibjet_matched_frac_jesdown, calibjet_matched_frac_jesup, calibjet_matched_frac_jerdown, calibjet_matched_frac_jerup;

    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 100; ++entry) {
        if (entry % 1000 == 0) cout << "event " << entry << endl;
        chain.GetEntry(entry);

        // Z-vertex cut.
        if (isnan(zvertex)) { continue; }
        if (fabs(zvertex) > 30) { continue; }

        // Fill z-vertex histogram.
        h_zvertex->Fill(zvertex);

        //////////////////////////// SETUP JET VARIABLES FOR UNFOLDING ////////////////////////////

        // Get truth leading jet
        int leadingtruthjet_index = -9999;
        get_leading_jet(leadingtruthjet_index, truthjet_pt);
        if (leadingtruthjet_index < 0) { continue; }
        float truthlead_pt = truthjet_pt->at(leadingtruthjet_index);
        float truthlead_e = truthjet_e->at(leadingtruthjet_index);
        float truthlead_eta = truthjet_eta->at(leadingtruthjet_index);
        float truthlead_phi = truthjet_phi->at(leadingtruthjet_index);

        //std::cout << "Before truth cut: truth pt " << truthlead_pt << " truth e " << truthlead_e << " truth eta " << truthlead_eta << " truth phi " << truthlead_phi << std::endl;
        if (truthlead_pt < truthjet_pt_min || truthlead_pt > truthjet_pt_max) { continue; }
        //std::cout << "After truth cut: truth pt " << truthlead_pt << " truth e " << truthlead_e << " truth eta " << truthlead_eta << " truth phi " << truthlead_phi << std::endl;

        // Get reco leading jet and subleading jet. Do basic jet cuts.
        int leadingunsubjet_index = -1;
        int subleadingunsubjet_index = -1;
        get_leading_subleading_jet(leadingunsubjet_index, subleadingunsubjet_index, unsubjet_pt);
        if (leadingunsubjet_index < 0) { continue; }

        // Beam background cut. Event level recording.
        background_dijet = true;
        background_frac = true;
        
        // Dijet cut.
        bool match_dijet = false;
        if (subleadingunsubjet_index >= 0 && unsubjet_e->at(subleadingunsubjet_index) / (float) unsubjet_e->at(leadingunsubjet_index) > 0.3) {
          match_dijet = match_leading_subleading_jet(unsubjet_phi->at(leadingunsubjet_index), unsubjet_phi->at(subleadingunsubjet_index));
        }
        if (match_dijet) { background_dijet = false; }
        
        // Fraction cut.
        double emfrac = unsubjet_emcal_calo_e->at(leadingunsubjet_index) / (double)(unsubjet_e->at(leadingunsubjet_index));
        double ihfrac = unsubjet_ihcal_calo_e->at(leadingunsubjet_index) / (double)(unsubjet_e->at(leadingunsubjet_index));
        double ohfrac = unsubjet_ohcal_calo_e->at(leadingunsubjet_index) / (double)(unsubjet_e->at(leadingunsubjet_index));
        if (emfrac > 0.1 || emfrac < 0.9 || ohfrac > 0.1 || ohfrac < 0.9 || ihfrac < 0.9) { background_frac = false; }

        float lead_pt = unsubjet_pt->at(leadingunsubjet_index);
        float lead_e = unsubjet_e->at(leadingunsubjet_index);
        float lead_eta = unsubjet_eta->at(leadingunsubjet_index);
        float lead_phi = unsubjet_phi->at(leadingunsubjet_index);

        //std::cout << "Reco lead: pt " << lead_pt << " e " << lead_e << " eta " << lead_eta << " phi " << lead_phi << std::endl;

        jet_filter = false; truthjet_filter = false;
        filter_jet(jet_filter, lead_e, lead_pt, lead_eta, zvertex, jet_radius);
        filter_jet(truthjet_filter, truthlead_e, truthlead_pt, truthlead_eta, zvertex, jet_radius);
        if (!jet_filter) h_recojet_pt_record->Fill(lead_pt);
        if (!jet_filter && !background_dijet) h_recojet_pt_record_dijet->Fill(lead_pt);
        if (!jet_filter && !background_frac) h_recojet_pt_record_frac->Fill(lead_pt);

        //std::cout << "Reco lead: pt " << lead_pt << " e " << lead_e << " eta " << lead_eta << " phi " << lead_phi << " jet filter " << jet_filter << " dijet bkg " << background_dijet << " frac bkg " << background_frac << std::endl;
        //std::cout << "Truth lead: pt " << truthlead_pt << " truth e " << truthlead_e << " truth eta " << truthlead_eta << " truth phi " << truthlead_phi << " jet filter " << truthjet_filter << std::endl;

        get_truthjet(goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi, truthjet_filter, truthlead_pt, truthlead_eta, truthlead_phi);
        get_calibjet(calibjet_pt_dijet, calibjet_eta_dijet, calibjet_phi_dijet, jet_filter, lead_pt, lead_eta, lead_phi, background_dijet, f_corr, 1, 0.1);
        get_calibjet(calibjet_pt_dijet_jesdown, calibjet_eta_dijet_jesdown, calibjet_phi_dijet_jesdown, jet_filter, lead_pt, lead_eta, lead_phi, background_dijet, f_corr, 0.94, 0.1);
        get_calibjet(calibjet_pt_dijet_jesup, calibjet_eta_dijet_jesup, calibjet_phi_dijet_jesup, jet_filter, lead_pt, lead_eta, lead_phi, background_dijet, f_corr, 1.06, 0.1);
        get_calibjet(calibjet_pt_dijet_jerdown, calibjet_eta_dijet_jerdown, calibjet_phi_dijet_jerdown, jet_filter, lead_pt, lead_eta, lead_phi, background_dijet, f_corr, 1, 0.05);
        get_calibjet(calibjet_pt_dijet_jerup, calibjet_eta_dijet_jerup, calibjet_phi_dijet_jerup, jet_filter, lead_pt, lead_eta, lead_phi, background_dijet, f_corr, 1, 0.15);
        get_calibjet(calibjet_pt_frac, calibjet_eta_frac, calibjet_phi_frac, jet_filter, lead_pt, lead_eta, lead_phi, background_frac, f_corr, 1, 0.1);
        get_calibjet(calibjet_pt_frac_jesdown, calibjet_eta_frac_jesdown, calibjet_phi_frac_jesdown, jet_filter, lead_pt, lead_eta, lead_phi, background_frac, f_corr, 0.94, 0.1);
        get_calibjet(calibjet_pt_frac_jesup, calibjet_eta_frac_jesup, calibjet_phi_frac_jesup, jet_filter, lead_pt, lead_eta, lead_phi, background_frac, f_corr, 1.06, 0.1);
        get_calibjet(calibjet_pt_frac_jerdown, calibjet_eta_frac_jerdown, calibjet_phi_frac_jerdown, jet_filter, lead_pt, lead_eta, lead_phi, background_frac, f_corr, 1, 0.05);
        get_calibjet(calibjet_pt_frac_jerup, calibjet_eta_frac_jerup, calibjet_phi_frac_jerup, jet_filter, lead_pt, lead_eta, lead_phi, background_frac, f_corr, 1, 0.15);

        //std::cout << "Good reco lead: pt " << calibjet_pt_frac << " eta " << calibjet_eta_frac << " phi " << calibjet_phi_frac << std::endl;
        //std::cout << "Good truth lead: pt " << goodtruthjet_pt << " truth eta " << goodtruthjet_eta << " truth phi " << goodtruthjet_phi << std::endl;
        
        match_meas_truth(calibjet_eta_dijet, calibjet_phi_dijet, calibjet_matched_dijet, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_dijet_jesdown, calibjet_phi_dijet_jesdown, calibjet_matched_dijet_jesdown, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_dijet_jesup, calibjet_phi_dijet_jesup, calibjet_matched_dijet_jesup, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_dijet_jerdown, calibjet_phi_dijet_jerdown, calibjet_matched_dijet_jerdown, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_dijet_jerup, calibjet_phi_dijet_jerup, calibjet_matched_dijet_jerup, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_frac, calibjet_phi_frac, calibjet_matched_frac, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_frac_jesdown, calibjet_phi_frac_jesdown, calibjet_matched_frac_jesdown, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_frac_jesup, calibjet_phi_frac_jesup, calibjet_matched_frac_jesup, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_frac_jerdown, calibjet_phi_frac_jerdown, calibjet_matched_frac_jerdown, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_frac_jerup, calibjet_phi_frac_jerup, calibjet_matched_frac_jerup, goodtruthjet_eta, goodtruthjet_phi, jet_radius);

        //std::cout << "Match " << calibjet_matched_frac << std::endl;
    
        //////////////////////////// SETUP UE VARIABLES FOR UNFOLDING //////////////////////////////

        // find reco ET information 
        float et_transverse = 0;
        if (!clusters) {
            for (int i = 0; i < emcaln; i++) {
                float dphi = get_dphi(lead_phi,emcalphi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += emcale[i]/cosh(emcaleta[i]); } 
            }
            for (int i = 0; i < ihcaln; i++) {
                float dphi = get_dphi(lead_phi,ihcalphi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += ihcale[i]/cosh(ihcaleta[i]); } 
            }
            for (int i = 0; i < ohcaln; i++) {
                float dphi = get_dphi(lead_phi,ohcalphi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += ohcale[i]/cosh(ohcaleta[i]); } 
            }
        } else {
            for (int i = 0; i < clsmult; i++) {
                float dphi = get_dphi(lead_phi,cluster_phi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += cluster_e[i]/cosh(cluster_eta[i]); }
            }
        }

        // find truth ET information
        float truth_et_transverse = 0;
        for (int i = 0; i < truthpar_n; i++) {
            if (fabs(truthpar_eta[i]) > 1.1) { continue; }
            if ((truthpar_pid[i] == 22 || truthpar_pid[i] == 111) && fabs(truthpar_e[i]) < 0.2) { continue; }
            else if (fabs(truthpar_e[i]) < 0.2) { continue; } // edited from 0.5 to 0.2
            float dphi = get_dphi(truthlead_phi,truthpar_phi[i]);
            if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { truth_et_transverse += truthpar_e[i]/cosh(truthpar_eta[i]); } 
        }

        //std::cout << "Truth ET " << truth_et_transverse << " Reco ET " << et_transverse << std::endl;
        //std::cout << std::endl;

        if (!jet_filter && !background_dijet) {
            h_reco_et_transverse_record_dijet->Fill(et_transverse);
            h_reco_et_transverse_record_dijet_weighted->Fill(et_transverse, weight_scale);
        }
        if (!truthjet_filter) {
            h_truth_et_transverse_record_dijet->Fill(truth_et_transverse);
            h_truth_et_transverse_record_dijet_weighted->Fill(truth_et_transverse, weight_scale);
        }


        //////////////////////////// FILL RESPONSE MATRICES ////////////////////////////
        fill_response_matrix(h_truth_calib_dijet, h_measure_calib_dijet, h_respmatrix_calib_dijet, h_fake_calib_dijet, h_miss_calib_dijet, calibjet_pt_dijet, calibjet_matched_dijet, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        fill_response_matrix(h_truth_calib_dijet_jesdown, h_measure_calib_dijet_jesdown, h_respmatrix_calib_dijet_jesdown, h_fake_calib_dijet_jesdown, h_miss_calib_dijet_jesdown, calibjet_pt_dijet_jesdown, calibjet_matched_dijet_jesdown, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        fill_response_matrix(h_truth_calib_dijet_jesup, h_measure_calib_dijet_jesup, h_respmatrix_calib_dijet_jesup, h_fake_calib_dijet_jesup, h_miss_calib_dijet_jesup, calibjet_pt_dijet_jesup, calibjet_matched_dijet_jesup, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        fill_response_matrix(h_truth_calib_dijet_jerdown, h_measure_calib_dijet_jerdown, h_respmatrix_calib_dijet_jerdown, h_fake_calib_dijet_jerdown, h_miss_calib_dijet_jerdown, calibjet_pt_dijet_jerdown, calibjet_matched_dijet_jerdown, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        fill_response_matrix(h_truth_calib_dijet_jerup, h_measure_calib_dijet_jerup, h_respmatrix_calib_dijet_jerup, h_fake_calib_dijet_jerup, h_miss_calib_dijet_jerup, calibjet_pt_dijet_jerup, calibjet_matched_dijet_jerup, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        fill_response_matrix(h_truth_calib_frac, h_measure_calib_frac, h_respmatrix_calib_frac, h_fake_calib_frac, h_miss_calib_frac, calibjet_pt_frac, calibjet_matched_frac, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        fill_response_matrix(h_truth_calib_frac_jesdown, h_measure_calib_frac_jesdown, h_respmatrix_calib_frac_jesdown, h_fake_calib_frac_jesdown, h_miss_calib_frac_jesdown, calibjet_pt_frac_jesdown, calibjet_matched_frac_jesdown, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        fill_response_matrix(h_truth_calib_frac_jesup, h_measure_calib_frac_jesup, h_respmatrix_calib_frac_jesup, h_fake_calib_frac_jesup, h_miss_calib_frac_jesup, calibjet_pt_frac_jesup, calibjet_matched_frac_jesup, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        fill_response_matrix(h_truth_calib_frac_jerdown, h_measure_calib_frac_jerdown, h_respmatrix_calib_frac_jerdown, h_fake_calib_frac_jerdown, h_miss_calib_frac_jerdown, calibjet_pt_frac_jerdown, calibjet_matched_frac_jerdown, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        fill_response_matrix(h_truth_calib_frac_jerup, h_measure_calib_frac_jerup, h_respmatrix_calib_frac_jerup, h_fake_calib_frac_jerup, h_miss_calib_frac_jerup, calibjet_pt_frac_jerup, calibjet_matched_frac_jerup, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        if (entry % 2 == 0){
            fill_response_matrix(h_truth_calib_frac_half1, h_measure_calib_frac_half1, h_respmatrix_calib_frac_half1, h_fake_calib_frac_half1, h_miss_calib_frac_half1, calibjet_pt_frac, calibjet_matched_frac, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        } else {
            fill_response_matrix(h_truth_calib_frac_half2, h_measure_calib_frac_half2, h_respmatrix_calib_frac_half2, h_fake_calib_frac_half2, h_miss_calib_frac_half2, calibjet_pt_frac, calibjet_matched_frac, goodtruthjet_pt, et_transverse, truth_et_transverse, weight_scale);
        }
    }

    std::cout << "Writing histograms..." << std::endl;
    f_out->cd();

    h_zvertex->Write();
    h_recojet_pt_record->Write();
    h_recojet_pt_record_dijet->Write();
    h_recojet_pt_record_frac->Write();
    h_reco_et_transverse_record_dijet->Write();
    h_truth_et_transverse_record_dijet->Write();
    h_reco_et_transverse_record_dijet_weighted->Write();
    h_truth_et_transverse_record_dijet_weighted->Write();

    h_truth_calib_dijet->Write();
    h_measure_calib_dijet->Write();
    h_respmatrix_calib_dijet->Write();
    h_fake_calib_dijet->Write();
    h_miss_calib_dijet->Write();

    h_truth_calib_dijet_jesdown->Write();
    h_measure_calib_dijet_jesdown->Write();
    h_respmatrix_calib_dijet_jesdown->Write();
    h_fake_calib_dijet_jesdown->Write();
    h_miss_calib_dijet_jesdown->Write();

    h_truth_calib_dijet_jesup->Write();
    h_measure_calib_dijet_jesup->Write();
    h_respmatrix_calib_dijet_jesup->Write();
    h_fake_calib_dijet_jesup->Write();
    h_miss_calib_dijet_jesup->Write();

    h_truth_calib_dijet_jerdown->Write();
    h_measure_calib_dijet_jerdown->Write();
    h_respmatrix_calib_dijet_jerdown->Write();
    h_fake_calib_dijet_jerdown->Write();
    h_miss_calib_dijet_jerdown->Write();

    h_truth_calib_dijet_jerup->Write();
    h_measure_calib_dijet_jerup->Write();
    h_respmatrix_calib_dijet_jerup->Write();
    h_fake_calib_dijet_jerup->Write();
    h_miss_calib_dijet_jerup->Write();

    h_truth_calib_frac->Write();
    h_measure_calib_frac->Write();
    h_respmatrix_calib_frac->Write();
    h_fake_calib_frac->Write();
    h_miss_calib_frac->Write();

    h_truth_calib_frac_jesdown->Write();
    h_measure_calib_frac_jesdown->Write();
    h_respmatrix_calib_frac_jesdown->Write();
    h_fake_calib_frac_jesdown->Write();
    h_miss_calib_frac_jesdown->Write();

    h_truth_calib_frac_jesup->Write();
    h_measure_calib_frac_jesup->Write();
    h_respmatrix_calib_frac_jesup->Write();
    h_fake_calib_frac_jesup->Write();
    h_miss_calib_frac_jesup->Write();

    h_truth_calib_frac_jerdown->Write();
    h_measure_calib_frac_jerdown->Write();
    h_respmatrix_calib_frac_jerdown->Write();
    h_fake_calib_frac_jerdown->Write();
    h_miss_calib_frac_jerdown->Write();

    h_truth_calib_frac_jerup->Write();
    h_measure_calib_frac_jerup->Write();
    h_respmatrix_calib_frac_jerup->Write();
    h_fake_calib_frac_jerup->Write();
    h_miss_calib_frac_jerup->Write();

    h_truth_calib_frac_half1->Write();
    h_measure_calib_frac_half1->Write();
    h_respmatrix_calib_frac_half1->Write();
    h_fake_calib_frac_half1->Write();
    h_miss_calib_frac_half1->Write();

    h_truth_calib_frac_half2->Write();
    h_measure_calib_frac_half2->Write();
    h_respmatrix_calib_frac_half2->Write();
    h_fake_calib_frac_half2->Write();
    h_miss_calib_frac_half2->Write();

    f_out->Close();
    std::cout << "All done!" << std::endl;

}

////////////////////////////////////////// Functions //////////////////////////////////////////
void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_pt) {
  leadingjet_index = -1;
  float leadingjet_pt = -9999;
  for (int ij = 0; ij < jet_pt->size(); ++ij) {
    float jetpt = jet_pt->at(ij);
    if (jetpt > leadingjet_pt) {
      leadingjet_pt = jetpt;
      leadingjet_index = ij;
    }
  }
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

bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi) {
  float dijet_min_phi = 3*TMath::Pi()/4.;
  float dphi = get_dphi(leadingjet_phi, subleadingjet_phi);
  return dphi > dijet_min_phi;
}

void filter_jet(bool& jet_filter, float jet_e, float jet_pt, float jet_eta, float zvertex, float jet_radius) {
    jet_filter = (jet_e < 0 || check_bad_jet_eta(jet_eta, zvertex, jet_radius));
}

void get_jet_filter(std::vector<bool>& jet_filter, std::vector<float>* jet_e, std::vector<float>* jet_pt, std::vector<float>* jet_eta, float zvertex, float jet_radius) {
  jet_filter.clear();
  int njet = jet_e->size();
  for (int ij = 0; ij < njet; ++ij) {
    jet_filter.push_back(jet_e->at(ij) < 0 || check_bad_jet_eta(jet_eta->at(ij), zvertex, jet_radius));
  }
}

void get_calibjet(float& calibjet_pt, float& calibjet_eta, float& calibjet_phi, bool jet_filter, float jet_pt, float jet_eta, float jet_phi, bool jet_background, TF1* f_corr, float jes_para, float jer_para) {
    calibjet_pt = -9999;
    calibjet_eta = -9999;
    calibjet_phi = -9999;
    if (jet_filter || jet_background) return;
    double calib_pt = f_corr->Eval(jet_pt) * (1 + randGen.Gaus(0.0, jer_para)) * jes_para;
    if (calib_pt < calibptbins[0] || calib_pt > calibptbins[calibnpt]) return;
    calibjet_pt = calib_pt;
    calibjet_eta = jet_eta;
    calibjet_phi = jet_phi;
}

void get_truthjet(float& goodtruthjet_pt, float& goodtruthjet_eta, float& goodtruthjet_phi, bool truthjet_filter, float jet_pt, float jet_eta, float jet_phi) {
    goodtruthjet_pt = -9999;
    goodtruthjet_eta = -9999;
    goodtruthjet_phi = -9999;
    if (truthjet_filter) return;
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

void fill_response_matrix(TH2D*& h_truth, TH2D*& h_meas, RooUnfoldResponse*& h_resp, TH2D*& h_fake, TH2D*& h_miss, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale) {
    if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet] && truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
        if (matched) {
            h_meas->Fill(meas_pt, meas_et, weight_scale);
            h_truth->Fill(truth_pt, truth_et, weight_scale);
            h_resp->Fill(meas_pt, meas_et, truth_pt, truth_et, weight_scale);
        } else {
            h_meas->Fill(meas_pt, meas_et, weight_scale);
            h_fake->Fill(meas_pt, meas_et, weight_scale);
            h_resp->Fake(meas_pt, meas_et, weight_scale);
            h_truth->Fill(truth_pt, truth_et, weight_scale);
            h_miss->Fill(truth_pt, truth_et, weight_scale);
            h_resp->Miss(truth_pt, truth_et, weight_scale);
        }
    } else if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet]) {
        h_meas->Fill(meas_pt, meas_et, weight_scale);
        h_fake->Fill(meas_pt, meas_et, weight_scale);
        h_resp->Fake(meas_pt, meas_et, weight_scale);

    } else if (truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
        h_truth->Fill(truth_pt, truth_et, weight_scale);
        h_miss->Fill(truth_pt, truth_et, weight_scale);
        h_resp->Miss(truth_pt, truth_et, weight_scale);
    }
}
