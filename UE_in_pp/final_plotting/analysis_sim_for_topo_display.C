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

void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et);
void get_leading_subleading_subsubleading_jet(int& leadingjet_index, int& subleadingjet_index, int& subsubleadingjet_index, std::vector<float>* jet_et);
void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_et);
bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi);
void get_calibjet(float& calibjet_pt, float& calibjet_eta, float& calibjet_phi, float jet_pt, float jet_eta, float jet_phi, bool reco_cut, TF1* f_corr, float jes_para, float jer_para);
void get_truthjet(float& goodtruthjet_pt, float& goodtruthjet_eta, float& goodtruthjet_phi, float jet_pt, float jet_eta, float jet_phi, bool truth_cut);
void match_meas_truth(float meas_eta, float meas_phi, bool& matched, float truth_eta, float truth_phi, float jet_radius);
void find_match_meas_truth(float meas_eta, float meas_phi, bool& matched, float& deltaR, float truth_eta, float truth_phi, float jet_radius);
void fill_response_matrix(TH2D*& h_truth, TH2D*& h_meas, RooUnfoldResponse*& h_resp, TH2D*& h_fake, TH2D*& h_miss, RooUnfoldResponse*& h_count, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale);
void fill_response_matrix_full(TH2D*& h_truth, TH2D*& h_meas, RooUnfoldResponse*& h_resp, TH2D*& h_fake, TH2D*& h_miss, RooUnfoldResponse*& h_count, TH2D*& h_count_fake, TH2D*& h_count_miss, RooUnfoldResponse*& h_jetpt, RooUnfoldResponse*& h_caloet, RooUnfoldResponse*& h_count_jetpt, RooUnfoldResponse*& h_count_caloet, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale);
void fill_trim_response_matrix(RooUnfoldResponse*& h_resp, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale, const TMatrixD& counts_matrix, TH2D* h_counts_measured, TH2D* h_counts_truth, float trim_value);
void fill_reweighted_trim_response_matrix(RooUnfoldResponse*& h_resp, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale, const TMatrixD& counts_matrix, TH2D* h_counts_measured, TH2D* h_counts_truth, float trim_value, TH2D* h_prior_weights);
void fill_trim_1D_response_matrices(RooUnfoldResponse*& h_jetpt_resp, RooUnfoldResponse*& h_caloet_resp, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale, const TMatrixD& counts_matrix, TH2D* h_counts_measured, TH2D* h_counts_truth, float trim_value);
void fill_reweighted_trim_1D_response_matrices(RooUnfoldResponse*& h_jetpt_resp, RooUnfoldResponse*& h_caloet_resp, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale, const TMatrixD& counts_matrix, TH2D* h_counts_measured, TH2D* h_counts_truth, float trim_value, TH1D* h_jetpt_prior_weights, TH1D* h_caloet_prior_weights);

TRandom3 randGen(1234);
TRandom3 randGen2(123);
TRandom3 clusERandGen(34);
TRandom3 randGen3(67);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                               //
//  Note: 5.13.25                                                                                                                //
//  This macro creates response matrices for 2D unfolding of leading jet pT and transverse region energy density                 //
//                                                                                                                               //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

void OutputVarBinMapping(TH1D*& h_binning, const double* bins, int nbins) {
  if (h_binning->GetNbinsX() + 1 != nbins) {
    std::cout << "Error in uniform to variable binning output, sizes don't match: ";
    std::cout << h_binning->GetNbinsX() + 1 << " " << nbins << std::endl;
  }
  for (int i = 0; i < h_binning->GetNbinsX() + 1; i++) {
    h_binning->SetBinContent(i,bins[i]);
  }
}

float ohcal_scale_eta[] = {0.925,0.943,0.9920,1.025,1.0443,1.0700,1.069,1.0540,1.0005,1.006,1.0410,1.0620,1.07,1.0247,1.0068,0.9759,1.0084,1.050,1.040,1.019,1.0012,1.0326,0.952,0.997};

void analysis_sim(std::string runtype = "mb", int start_seg = 0, int end_seg = 200, int iter = 1, std::string bkg_cut = "dijet", bool clusters = true, bool emcal_clusters = false)  {
    
    ////////// General Set up //////////
    double weight_scale = 1.0, truthjet_pt_min = 0, truthjet_pt_max = 3000;
    string count_filename;
    if (runtype == "jet5") {
        weight_scale = Jet5GeV_scale;
        truthjet_pt_min = 7;
        truthjet_pt_max = 14;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_jet5.root";
    } else if (runtype == "jet12") {
        weight_scale = Jet12GeV_scale;
        truthjet_pt_min = 14;
        truthjet_pt_max = 21;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_jet12.root";
    } else if (runtype == "jet20") {
        weight_scale = Jet20GeV_scale;
        truthjet_pt_min = 21;
        truthjet_pt_max = 32;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_jet20.root";
    } else if (runtype == "jet30") {
        weight_scale = Jet30GeV_scale;
        truthjet_pt_min = 32;
        truthjet_pt_max = 42;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_jet30.root";
    } else if (runtype == "jet40") {
        weight_scale = Jet40GeV_scale;
        truthjet_pt_min = 42;
        truthjet_pt_max = 52;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_jet40.root";
    } else if (runtype == "jet50") {
        weight_scale = Jet50GeV_scale;
        truthjet_pt_min = 52;
        truthjet_pt_max = 62;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_jet50.root";
    } else if (runtype == "jet60") {
        weight_scale = Jet60GeV_scale;
        truthjet_pt_min = 62;
        truthjet_pt_max = 1000;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_jet60.root";
    } else if (runtype == "herwig_jet12") {
        weight_scale = Herwig_Jet12GeV_scale;
        truthjet_pt_min = 14;
        truthjet_pt_max = 21;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_herwig_jet12.root";
    } else if (runtype == "herwig_jet20") {
        weight_scale = Herwig_Jet20GeV_scale;
        truthjet_pt_min = 21;
        truthjet_pt_max = 32;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_herwig_jet20.root";
    } else if (runtype == "herwig_jet30") {
        weight_scale = Herwig_Jet30GeV_scale;
        truthjet_pt_min = 32;
        truthjet_pt_max = 42;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_herwig_jet30.root";
    } else if (runtype == "herwig_jet40") {
        weight_scale = Herwig_Jet40GeV_scale;
        truthjet_pt_min = 42;
        truthjet_pt_max = 52;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_herwig_jet40.root";
    } else if (runtype == "herwig_jet50") {
        weight_scale = Herwig_Jet50GeV_scale;
        truthjet_pt_min = 52;
        truthjet_pt_max = 1000;
        if (iter > 1) count_filename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_1_herwig_jet50.root";
    } else {
        std::cout << "Unknown runtype" << std::endl;
        return;
    }

    TFile *f_zvertex;
    TH1D* h_vertex_weight;
    string zvertex_filename = "output_vertex_reweight_run28_" + bkg_cut + "_bkg_cut_iter_1.root";
    f_zvertex = new TFile(zvertex_filename.c_str(),"READ");
    TH1D* htemp = (TH1D*)f_zvertex->Get("vertex_weight");
    if (!htemp) {
        std::cerr << "ERROR: Could not find histogram 'vertex_weight'!" << std::endl;
        exit(1);
    }
    h_vertex_weight = (TH1D*)htemp->Clone("vertex_weight_clone");
    h_vertex_weight->SetDirectory(0);  // detach from file
    f_zvertex->Close();

    TFile *f_count;
    RooUnfoldResponse* h_count[13];
    std::vector<std::string> count_matrix_names = {"h_respmatrix_calib_dijet_counts","h_respmatrix_calib_dijet_jesdown_counts",
    "h_respmatrix_calib_dijet_jesup_counts","h_respmatrix_calib_dijet_jerdown_counts","h_respmatrix_calib_dijet_jerup_counts",
    "h_respmatrix_calib_dijet_clus_smear_counts","h_respmatrix_calib_dijet_ohcal_mc_data_var_counts","h_respmatrix_calib_dijet_2sigma_noise_counts","h_respmatrix_calib_dijet_4sigma_noise_counts",
    "h_respmatrix_calib_dijet_half1_counts","h_respmatrix_calib_dijet_half2_counts",
    "h_jetpt_respmatrix_counts","h_caloet_respmatrix_counts"};
    std::vector<TMatrixD> counts_matrix;
    //counts_matrix.resize(7);
    TH2D* counts_measured[11]; 
    TH2D* counts_truth[11]; 
    TH1D* counts1D_measured[2];
    TH1D* counts1D_truth[2];

    if (iter > 1) {
        f_count = new TFile(count_filename.c_str(), "READ");
        for (int i = 0; i < count_matrix_names.size(); i++) {
            RooUnfoldResponse* resp = dynamic_cast<RooUnfoldResponse*>(f_count->Get(count_matrix_names[i].c_str()));
            if (!resp) {
                std::cerr << "ERROR: missing RooUnfoldResponse '" << count_matrix_names[i] << "'" << std::endl;
                exit(1);
            }
            const TMatrixD& mat_ref = resp->Mresponse(false);
            counts_matrix.push_back(mat_ref);
            std::cout << "Matrix " << i << " size: " << counts_matrix.back().GetNrows() << " x " << counts_matrix.back().GetNcols() << std::endl;
            TH1* h_meas = resp->Hmeasured();
            TH1* h_true = resp->Htruth();
            if (!h_meas || !h_true) {
                std::cerr << "ERROR: Measured or truth histogram missing for response " << i << std::endl;
                exit(1);
            }
            if (i < 11) {
                counts_measured[i] = dynamic_cast<TH2D*>(h_meas->Clone());
                counts_truth[i]    = dynamic_cast<TH2D*>(h_true->Clone());
                counts_measured[i]->SetDirectory(0);
                counts_truth[i]->SetDirectory(0);
            }
            else {
                counts1D_measured[i-11] = dynamic_cast<TH1D*>(h_meas->Clone());
                counts1D_truth[i-11]    = dynamic_cast<TH1D*>(h_true->Clone());
                counts1D_measured[i-11]->SetDirectory(0);
                counts1D_truth[i-11]->SetDirectory(0);
            }
            h_count[i] = dynamic_cast<RooUnfoldResponse*>(resp->Clone());
            h_count[i]->SetName((count_matrix_names[i] + "_clone").c_str());
        }
        f_count->Close();
        std::cout << "counts_measured: " << counts_measured[0]->GetEntries() << " counts_truth: " << counts_truth[0]->GetEntries() << std::endl;
    }
    /*
    if (iter > 1) {
        f_count = new TFile(count_filename.c_str(),"READ");
        for (int i = 0; i < count_matrix_names.size(); i++) {
            h_count[i] = dynamic_cast<RooUnfoldResponse*>(f_count->Get(count_matrix_names[i].c_str()));
            const TMatrixD& mat_ref = h_count[i]->Mresponse(false);
            TMatrixD mat_copy(mat_ref);
            counts_matrix.push_back(mat_copy);
            std::cout << "Matrix " << i << " size: " << counts_matrix[i].GetNrows() << " x " << counts_matrix[i].GetNcols() << std::endl;
            if (i < 7) {
                counts_measured[i] = dynamic_cast<TH2D*>(h_count[i]->Hmeasured());
                counts_truth[i] = dynamic_cast<TH2D*>(h_count[i]->Htruth());
            } else {
                counts1D_measured[i-7] = dynamic_cast<TH1D*>(h_count[i]->Hmeasured());
                counts1D_truth[i-7] = dynamic_cast<TH1D*>(h_count[i]->Htruth());
            }
        }
        f_count->Close();
        std::cout << "counts_measured: " << counts_measured[0]->GetEntries() << " counts_truth: " << counts_truth[0]->GetEntries() << std::endl;
    }
    */
    TFile *f_reweight; 
    TH2D* weights[11][3];
    TH1D* weights1D[2][3];
    std::vector<std::string> pw_syst = {"weights_calib_dijet","weights_calib_dijet_jesdown","weights_calib_dijet_jesup","weights_calib_dijet_jerdown","weights_calib_dijet_jerup",
    "weights_calib_dijet_clus_smear","weights_calib_dijet_ohcal_mc_data_var","weights_calib_dijet_2sigma_noise","weights_calib_dijet_4sigma_noise","weights_calib_dijet_half1","weights_calib_dijet_half2",
    "weights_jetpt","weights_caloet"};
    std::vector<std::string> pw_trim = {"","_trim_5","_trim_10"};
    if (iter > 2) {
        if (runtype == "herwig_jet12" || runtype == "herwig_jet20" || runtype == "herwig_jet30" || runtype == "herwig_jet40" || runtype == "herwig_jet50") {
            std::string rwfilename = "sphenix_primary_run28_output_files/output_herwig_reweighted_respmatrix_run28_iter_2_"+bkg_cut+"_bkg_cut.root";
            f_reweight = new TFile(rwfilename.c_str(), "READ");
        } else {
            std::string rwfilename = "sphenix_primary_run28_output_files/output_reweighted_respmatrix_run28_iter_2_"+bkg_cut+"_bkg_cut.root";
            f_reweight = new TFile(rwfilename.c_str(), "READ");
        }
        for (int i = 0; i < pw_syst.size(); i++) {
            for (int j = 0; j < pw_trim.size(); j++) {
                if (i < 11) { weights[i][j] = dynamic_cast<TH2D*>(f_reweight->Get((pw_syst[i]+pw_trim[j]).c_str())); }
                else { weights1D[i-11][j] = dynamic_cast<TH1D*>(f_reweight->Get((pw_syst[i]+pw_trim[j]).c_str())); }
            }
        }
    }

    TFile *f_deadmap = new TFile("hists_EMCalHotMap_and_NoCalibMap_ana509_2024p022_v001_53244cdb.root","READ");
    TH2F* EMCal_deadmap = dynamic_cast<TH2F*>(f_deadmap->Get("h_temp"));
    EMCal_deadmap->SetDirectory(0);
    std::cout << "deadmap loaded: " << EMCal_deadmap->GetEntries() << std::endl;

    ////////// Files //////////
    //TFile *f_out = new TFile(Form("analysis_sim_output/output_dijet_sim_iter_%d_%s_%d_%d.root", iter, runtype.c_str(), start_seg, end_seg), "RECREATE");
    string outfilename = "sphenix_primary_analysis_sim_run28_output/output_" + bkg_cut + "_bkg_cut_sim_iter_" + to_string(iter) + "_" + runtype + "_" + to_string(start_seg) + "_" + to_string(end_seg) + ".root";
    TFile *f_out = new TFile(outfilename.c_str(), "RECREATE");
    TChain chain("T");
    for (int i = start_seg; i < end_seg; ++i) {
        //chain.Add(Form("/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run21_%s_3sigma_output_%d.root", runtype.c_str(), i)); 
        //chain.Add(Form("run21_test_files/sim_run21_%s_3sigma_output_%d.root", runtype.c_str(), i)); 
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
    std::vector<float>* calibjet_pt = nullptr; chain.SetBranchStatus("calibpt", 1); chain.SetBranchAddress("calibpt", &calibjet_pt);
    std::vector<float>* truthjet_e = nullptr; chain.SetBranchStatus("truthE", 1); chain.SetBranchAddress("truthE", &truthjet_e);
    std::vector<float>* truthjet_pt = nullptr; chain.SetBranchStatus("truthPt", 1); chain.SetBranchAddress("truthPt", &truthjet_pt);
    std::vector<float>* truthjet_eta = nullptr; chain.SetBranchStatus("truthEta", 1); chain.SetBranchAddress("truthEta", &truthjet_eta);
    std::vector<float>* truthjet_phi = nullptr; chain.SetBranchStatus("truthPhi", 1); chain.SetBranchAddress("truthPhi", &truthjet_phi);

    int truthpar_n = 0; chain.SetBranchStatus("truthpar_n", 1); chain.SetBranchAddress("truthpar_n",&truthpar_n);
    float truthpar_e[10000] = {0.0}; chain.SetBranchStatus("truthpar_e", 1); chain.SetBranchAddress("truthpar_e",truthpar_e);
    float truthpar_eta[10000] = {0.0}; chain.SetBranchStatus("truthpar_eta", 1); chain.SetBranchAddress("truthpar_eta",truthpar_eta);
    float truthpar_phi[10000] = {0.0}; chain.SetBranchStatus("truthpar_phi", 1); chain.SetBranchAddress("truthpar_phi",truthpar_phi);
    int truthpar_pid[10000] = {0}; chain.SetBranchStatus("truthpar_pid", 1); chain.SetBranchAddress("truthpar_pid",truthpar_pid);
    float truthpar_pt[10000] = {0.0}; chain.SetBranchStatus("truthpar_pt", 1); chain.SetBranchAddress("truthpar_pt",truthpar_pt);

    //int emcaln = 0; float emcale[24576] = {0.0}; float emcaleta[24576] = {0.0}; float emcalphi[24576] = {0.0};
    //int ihcaln = 0; float ihcale[1536] = {0.0}; float ihcaleta[1536] = {0.0}; float ihcalphi[1536] = {0.0};
    //int ohcaln = 0; float ohcale[1536] = {0.0}; float ohcaleta[1536] = {0.0}; float ohcalphi[1536] = {0.0};
    int clsmult = 0; float cluster_e[1000] = {0.0}; float cluster_eta[1000] = {0.0}; float cluster_phi[1000] = {0.0};
    int clsmult2 = 0; float cluster2_e[1000] = {0.0}; float cluster2_eta[1000] = {0.0}; float cluster2_phi[1000] = {0.0};
    int clsmult4 = 0; float cluster4_e[1000] = {0.0}; float cluster4_eta[1000] = {0.0}; float cluster4_phi[1000] = {0.0};
    int cluster_ntowers[1000]; int cluster_tower_calo[1000][500]; int cluster_tower_ieta[1000][500]; int cluster_tower_iphi[1000][500]; float cluster_tower_e[1000][500];
    //int cluster2_ntowers[1000]; int cluster2_tower_calo[1000][500]; int cluster2_tower_ieta[1000][500]; int cluster2_tower_iphi[1000][500]; float cluster2_tower_e[1000][500];
    //int cluster4_ntowers[1000]; int cluster4_tower_calo[1000][500]; int cluster4_tower_ieta[1000][500]; int cluster4_tower_iphi[1000][500]; float cluster4_tower_e[1000][500];
    if (!clusters) {
        /*
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
        */
    } else if (clusters && !emcal_clusters) {
        chain.SetBranchStatus("clsmult", 1); chain.SetBranchAddress("clsmult",&clsmult);
        chain.SetBranchStatus("cluster_e", 1); chain.SetBranchAddress("cluster_e",cluster_e);
        chain.SetBranchStatus("cluster_eta", 1); chain.SetBranchAddress("cluster_eta",cluster_eta);
        chain.SetBranchStatus("cluster_phi", 1); chain.SetBranchAddress("cluster_phi",cluster_phi);
        chain.SetBranchStatus("cluster_ntowers", 1); chain.SetBranchAddress("cluster_ntowers",cluster_ntowers);
        chain.SetBranchStatus("cluster_tower_e", 1); chain.SetBranchAddress("cluster_tower_e",cluster_tower_e);
        chain.SetBranchStatus("cluster_tower_calo", 1); chain.SetBranchAddress("cluster_tower_calo",cluster_tower_calo);
        chain.SetBranchStatus("cluster_tower_ieta", 1); chain.SetBranchAddress("cluster_tower_ieta",cluster_tower_ieta);
        chain.SetBranchStatus("cluster_tower_iphi", 1); chain.SetBranchAddress("cluster_tower_iphi",cluster_tower_iphi);
        chain.SetBranchStatus("clsmult2", 1); chain.SetBranchAddress("clsmult2",&clsmult2);
        chain.SetBranchStatus("cluster2_e", 1); chain.SetBranchAddress("cluster2_e",cluster2_e);
        chain.SetBranchStatus("cluster2_eta", 1); chain.SetBranchAddress("cluster2_eta",cluster2_eta);
        chain.SetBranchStatus("cluster2_phi", 1); chain.SetBranchAddress("cluster2_phi",cluster2_phi);
        //chain.SetBranchStatus("cluster2_ntowers", 1); chain.SetBranchAddress("cluster2_ntowers",cluster2_ntowers);
        //chain.SetBranchStatus("cluster2_tower_e", 1); chain.SetBranchAddress("cluster2_tower_e",cluster2_tower_e);
        //chain.SetBranchStatus("cluster2_tower_calo", 1); chain.SetBranchAddress("cluster2_tower_calo",cluster2_tower_calo);
        //chain.SetBranchStatus("cluster2_tower_ieta", 1); chain.SetBranchAddress("cluster2_tower_ieta",cluster2_tower_ieta);
        //chain.SetBranchStatus("cluster2_tower_iphi", 1); chain.SetBranchAddress("cluster2_tower_iphi",cluster2_tower_iphi);
        chain.SetBranchStatus("clsmult4", 1); chain.SetBranchAddress("clsmult4",&clsmult4);
        chain.SetBranchStatus("cluster4_e", 1); chain.SetBranchAddress("cluster4_e",cluster4_e);
        chain.SetBranchStatus("cluster4_eta", 1); chain.SetBranchAddress("cluster4_eta",cluster4_eta);
        chain.SetBranchStatus("cluster4_phi", 1); chain.SetBranchAddress("cluster4_phi",cluster4_phi);
        //chain.SetBranchStatus("cluster4_ntowers", 1); chain.SetBranchAddress("cluster4_ntowers",cluster4_ntowers);
        //chain.SetBranchStatus("cluster4_tower_e", 1); chain.SetBranchAddress("cluster4_tower_e",cluster4_tower_e);
        //chain.SetBranchStatus("cluster4_tower_calo", 1); chain.SetBranchAddress("cluster4_tower_calo",cluster4_tower_calo);
        //chain.SetBranchStatus("cluster4_tower_ieta", 1); chain.SetBranchAddress("cluster4_tower_ieta",cluster4_tower_ieta);
        //chain.SetBranchStatus("cluster4_tower_iphi", 1); chain.SetBranchAddress("cluster4_tower_iphi",cluster4_tower_iphi);
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

    ///////// JES and JER variations ////////////////
    double JER_smear = 0.101;
    double JER_var = 0.018;
    double JES_var = 0.025;

    double efracbins[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};

    ////////// Histograms //////////
    TH1D* h_zvertex_noreweight = new TH1D("h_zvertex_noreweight", ";Z-vertex [cm]", 400, -200, 200);
    TH1D* h_zvertex = new TH1D("h_zvertex", ";Z-vertex [cm]", 400, -200, 200);
    TH1D* h_deltaphi_record = new TH1D("h_deltaphi_record","",125,-2*M_PI,2*M_PI);
    TH1D* h_xj_record = new TH1D("h_xj_record","",20,0,1);
    TH1D* h_uncalib_lead_spectra_record = new TH1D("h_uncalib_lead_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
    TH1D* h_lead_spectra_record = new TH1D("h_lead_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
    TH1D* h_sub_spectra_record = new TH1D("h_sub_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
    TH1D* h_truth_deltaphi_record = new TH1D("h_truth_deltaphi_record","",125,-2*M_PI,2*M_PI);
    TH1D* h_truth_xj_record = new TH1D("h_truth_xj_record","",20,0,1);
    TH1D* h_lead_truth_spectra_record = new TH1D("h_lead_truth_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
    TH1D* h_sub_truth_spectra_record = new TH1D("h_sub_truth_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
    TH2D* h_jes_qa = new TH2D("h_jes_qa",";p^{truth}_{T} [GeV]; p^{reco}_{T}/p^{truth}_{T}", 50, 0, 100, 1000, 0, 10);
    TH2D* h_jes_qa_reco = new TH2D("h_jes_qa_reco",";p^{reco}_{T} [GeV]; p^{reco}_{T}/p^{truth}_{T}", 50, 0, 100, 1000, 0, 10);
    TH2D* h_efficiency_reco_bkg_cut = new TH2D("h_efficiency_reco_bkg_cut", "", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_match = new TH2D("h_efficiency_match", "", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_reco_pt_cut = new TH2D("h_efficiency_reco_pt_cut","", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_match_reco_bkg_cut = new TH2D("h_efficiency_match_reco_bkg_cut", "", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_reco_pt_cut_reco_bkg_cut = new TH2D("h_efficiency_reco_pt_cut_reco_bkg_cut","", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_match_reco_pt_cut_reco_bkg_cut = new TH2D("h_efficiency_match_reco_pt_cut_reco_bkg_cut", "", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efficiency_match_reco_pt_cut = new TH2D("h_efficiency_match_reco_pt_cut","", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_efrac_lead_pt_dijet = new TH2D("h_efrac_lead_pt_dijet","", calibnpt, calibptbins, 15, efracbins);
    TProfile* h_truth_avg_pt_vs_jet_pt = new TProfile("h_truth_avg_pt_vs_jet_pt","", truthnpt, truthptbins);

    TH2D* h_purity_truth_pt_cut = new TH2D("h_purity_truth_pt_cut","", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_purity_match = new TH2D("h_purity_match","", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_truth_qa = new TH2D("h_truth_qa","", truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_measure_qa = new TH2D("h_measure_qa","", calibnpt, calibptbins, calibnet, calibetbins);

    TH1D* efficiency_truth_match_test = new TH1D("efficiency_truth_match_test","", truthnpt, truthptbins);
    TH1D *h_et_transverse_record = new TH1D("h_et_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH1D* h_et_truth_transverse_record = new TH1D("h_et_truth_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH1D *h_nw_et_transverse_record = new TH1D("h_nw_et_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH1D* h_nw_et_truth_transverse_record = new TH1D("h_nw_et_truth_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH2D* h_ue_pt_transverse_record = new TH2D("h_ue_pt_transverse_record","", calibnpt, calibptbins, calibnet, calibetbins);
    TH2D* h_ue_pt_truth_transverse_record = new TH2D("h_ue_pt_truth_transverse_record","",truthnpt, truthptbins, truthnet, truthetbins);
    TH2D* h_ue_pt_transverse_record_full = new TH2D("h_ue_pt_transverse_record_full","", 100,0,100,700,-20,50);
    TH2D* h_ue_pt_truth_transverse_record_full = new TH2D("h_ue_pt_truth_transverse_record_full","",100,0,100,500,0,50);
    TH1D* h_thres_et_truth_transverse_record = new TH1D("h_thres_et_truth_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH1D* h_thres_nw_et_truth_transverse_record = new TH1D("h_thres_nw_et_truth_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH2D* h_thres_ue_pt_truth_transverse_record = new TH2D("h_thres_ue_pt_truth_transverse_record","",truthnpt, truthptbins, truthnet, truthetbins);
    TH1D* h_reco_et_truth_transverse_record = new TH1D("h_reco_et_truth_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH1D* h_reco_nw_et_truth_transverse_record = new TH1D("h_reco_nw_et_truth_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH2D* h_reco_ue_pt_truth_transverse_record = new TH2D("h_reco_ue_pt_truth_transverse_record","",truthnpt, truthptbins, truthnet, truthetbins);
    TH1D* h_reco_thres_et_truth_transverse_record = new TH1D("h_reco_thres_et_truth_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH1D* h_reco_thres_nw_et_truth_transverse_record = new TH1D("h_reco_thres_nw_et_truth_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
    TH2D* h_reco_thres_ue_pt_truth_transverse_record = new TH2D("h_reco_thres_ue_pt_truth_transverse_record","",truthnpt, truthptbins, truthnet, truthetbins);
    
    TH2D* h_truth_calib_dijet = new TH2D("h_truth_calib_dijet", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_measure_calib_dijet = new TH2D("h_measure_calib_dijet", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_fake_calib_dijet = new TH2D("h_fake_calib_dijet", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_miss_calib_dijet = new TH2D("h_miss_calib_dijet", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_counts_fake_calib_dijet = new TH2D("h_counts_fake_calib_dijet", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_counts_miss_calib_dijet = new TH2D("h_counts_miss_calib_dijet", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", truthnpt, 0, 1, calibnet, 0, 1);
    TH1D* h_jetpt_truth = new TH1D("h_jetpt_truth","", truthnpt, 0, 1); 
    TH1D* h_jetpt_measure = new TH1D("h_jetpt_measure","", calibnpt, 0, 1);
    TH1D* h_caloet_truth = new TH1D("h_caloet_truth","", truthnet, 0, 1); 
    TH1D* h_caloet_measure = new TH1D("h_caloet_measure","", calibnet, 0, 1);
    
    TH2D* h_truth_calib_dijet_jesdown = new TH2D("h_truth_calib_dijet_jesdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_measure_calib_dijet_jesdown = new TH2D("h_measure_calib_dijet_jesdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_fake_calib_dijet_jesdown = new TH2D("h_fake_calib_dijet_jesdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_miss_calib_dijet_jesdown = new TH2D("h_miss_calib_dijet_jesdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_truth_calib_dijet_jesup = new TH2D("h_truth_calib_dijet_jesup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_measure_calib_dijet_jesup = new TH2D("h_measure_calib_dijet_jesup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_fake_calib_dijet_jesup = new TH2D("h_fake_calib_dijet_jesup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_miss_calib_dijet_jesup = new TH2D("h_miss_calib_dijet_jesup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    
    TH2D* h_truth_calib_dijet_jerdown = new TH2D("h_truth_calib_dijet_jerdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_measure_calib_dijet_jerdown = new TH2D("h_measure_calib_dijet_jerdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_fake_calib_dijet_jerdown = new TH2D("h_fake_calib_dijet_jerdown", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_miss_calib_dijet_jerdown = new TH2D("h_miss_calib_dijet_jerdown", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_truth_calib_dijet_jerup = new TH2D("h_truth_calib_dijet_jerup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_measure_calib_dijet_jerup = new TH2D("h_measure_calib_dijet_jerup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_fake_calib_dijet_jerup = new TH2D("h_fake_calib_dijet_jerup", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_miss_calib_dijet_jerup = new TH2D("h_miss_calib_dijet_jerup", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    
    TH2D* h_truth_calib_dijet_half1 = new TH2D("h_truth_calib_dijet_half1", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_measure_calib_dijet_half1 = new TH2D("h_measure_calib_dijet_half1", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_fake_calib_dijet_half1 = new TH2D("h_fake_calib_dijet_half1", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_miss_calib_dijet_half1 = new TH2D("h_miss_calib_dijet_half1", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_truth_calib_dijet_half2 = new TH2D("h_truth_calib_dijet_half2", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_measure_calib_dijet_half2 = new TH2D("h_measure_calib_dijet_half2", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_fake_calib_dijet_half2 = new TH2D("h_fake_calib_dijet_half2", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_miss_calib_dijet_half2 = new TH2D("h_miss_calib_dijet_half2", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);

    TH2D* h_truth_calib_dijet_clus_smear = new TH2D("h_truth_calib_dijet_clus_smear", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_measure_calib_dijet_clus_smear = new TH2D("h_measure_calib_dijet_clus_smear", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_fake_calib_dijet_clus_smear = new TH2D("h_fake_calib_dijet_clus_smear", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_miss_calib_dijet_clus_smear = new TH2D("h_miss_calib_dijet_clus_smear", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_truth_calib_dijet_ohcal_mc_data_var = new TH2D("h_truth_calib_dijet_ohcal_mc_data_var", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_measure_calib_dijet_ohcal_mc_data_var = new TH2D("h_measure_calib_dijet_ohcal_mc_data_var", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_fake_calib_dijet_ohcal_mc_data_var = new TH2D("h_fake_calib_dijet_ohcal_mc_data_var", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_miss_calib_dijet_ohcal_mc_data_var = new TH2D("h_miss_calib_dijet_ohcal_mc_data_var", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_truth_calib_dijet_2sigma_noise = new TH2D("h_truth_calib_dijet_2sigma_noise", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_measure_calib_dijet_2sigma_noise = new TH2D("h_measure_calib_dijet_2sigma_noise", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_fake_calib_dijet_2sigma_noise = new TH2D("h_fake_calib_dijet_2sigma_noise", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_miss_calib_dijet_2sigma_noise = new TH2D("h_miss_calib_dijet_2sigma_noise", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_truth_calib_dijet_4sigma_noise = new TH2D("h_truth_calib_dijet_4sigma_noise", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);
    TH2D* h_measure_calib_dijet_4sigma_noise = new TH2D("h_measure_calib_dijet_4sigma_noise", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_fake_calib_dijet_4sigma_noise = new TH2D("h_fake_calib_dijet_4sigma_noise", ";p_{T}^{Calib jet} [GeV];#SigmaE_{T}^{Reco} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
    TH2D* h_miss_calib_dijet_4sigma_noise = new TH2D("h_miss_calib_dijet_4sigma_noise", ";p_{T}^{Truth jet} [GeV];#SigmaE_{T}^{Truth} [GeV]", truthnpt, 0, 1, truthnet, 0, 1);

    TH1D* h_exclusive_dijet_passcut = new TH1D("h_exclusive_dijet_passcut","",truthnpt, truthptbins);
    TH1D* h_exclusive_dijet_total = new TH1D("h_exclusive_dijet_total","",truthnpt, truthptbins);

    TH2D* h_total_topo_event_display = new TH2D("h_total_topo_event_display","",24,0,24,64,0,64);
    TH2D* h_total_topo_emcal_event_display = new TH2D("h_total_topo_emcal_event_display","",96,0,96,256,0,256);
    TH2D* h_total_topo_ihcal_event_display = new TH2D("h_total_topo_ihcal_event_display","",24,0,24,64,0,64);
    TH2D* h_total_topo_ohcal_event_display = new TH2D("h_total_topo_ohcal_event_display","",24,0,24,64,0,64);   
    TH2D* h_topo_event_display[30];
    TH2D* h_topo_emcal_event_display[30];
    TH2D* h_topo_ihcal_event_display[30];
    TH2D* h_topo_ohcal_event_display[30];
    for (int i = 0; i < 30; i++) {
        h_topo_event_display[i] = new TH2D(Form("h_topo_event_display_%d",i),"",24,0,24,64,0,64);
        h_topo_emcal_event_display[i] = new TH2D(Form("h_topo_emcal_event_display_%d",i),"",96,0,96,256,0,256);
        h_topo_ihcal_event_display[i] = new TH2D(Form("h_topo_ihcal_event_display_%d",i),"",24,0,24,64,0,64);
        h_topo_ohcal_event_display[i] = new TH2D(Form("h_topo_ohcal_event_display_%d",i),"",24,0,24,64,0,64);   
    }

    std::vector<int> trim_val = {0, 5, 10};
    std::vector<std::string> respmatrix_tags = {"", "_counts", "_trim_5", "_trim_10", "_reweight", "_reweight_trim_5", "_reweight_trim_10"};
    RooUnfoldResponse* h_respmatrix_calib_dijet[7];
    RooUnfoldResponse* h_respmatrix_calib_dijet_jesdown[7];
    RooUnfoldResponse* h_respmatrix_calib_dijet_jesup[7];
    RooUnfoldResponse* h_respmatrix_calib_dijet_jerdown[7];
    RooUnfoldResponse* h_respmatrix_calib_dijet_jerup[7];
    RooUnfoldResponse* h_respmatrix_calib_dijet_half1[7];
    RooUnfoldResponse* h_respmatrix_calib_dijet_half2[7];
    RooUnfoldResponse* h_respmatrix_calib_dijet_clus_smear[7];
    RooUnfoldResponse* h_respmatrix_calib_dijet_ohcal_mc_data_var[7];
    RooUnfoldResponse* h_respmatrix_calib_dijet_2sigma_noise[7];
    RooUnfoldResponse* h_respmatrix_calib_dijet_4sigma_noise[7];
    RooUnfoldResponse* h_jetpt_respmatrix[7];
    RooUnfoldResponse* h_caloet_respmatrix[7];

    for (int i = 0; i < 7; i++) {
        h_respmatrix_calib_dijet[i] = new RooUnfoldResponse(("h_respmatrix_calib_dijet" + respmatrix_tags[i]).c_str(),""); h_respmatrix_calib_dijet[i]->Setup(h_measure_calib_dijet, h_truth_calib_dijet);
        h_respmatrix_calib_dijet_jesdown[i] = new RooUnfoldResponse(("h_respmatrix_calib_dijet_jesdown" + respmatrix_tags[i]).c_str(),""); h_respmatrix_calib_dijet_jesdown[i]->Setup(h_measure_calib_dijet_jesdown, h_truth_calib_dijet_jesdown);
        h_respmatrix_calib_dijet_jesup[i] = new RooUnfoldResponse(("h_respmatrix_calib_dijet_jesup" + respmatrix_tags[i]).c_str(),""); h_respmatrix_calib_dijet_jesup[i]->Setup(h_measure_calib_dijet_jesup, h_truth_calib_dijet_jesup);
        h_respmatrix_calib_dijet_jerdown[i] = new RooUnfoldResponse(("h_respmatrix_calib_dijet_jerdown" + respmatrix_tags[i]).c_str(),""); h_respmatrix_calib_dijet_jerdown[i]->Setup(h_measure_calib_dijet_jerdown, h_truth_calib_dijet_jerdown);
        h_respmatrix_calib_dijet_jerup[i] = new RooUnfoldResponse(("h_respmatrix_calib_dijet_jerup" + respmatrix_tags[i]).c_str(),""); h_respmatrix_calib_dijet_jerup[i]->Setup(h_measure_calib_dijet_jerup, h_truth_calib_dijet_jerup);
        h_respmatrix_calib_dijet_half1[i] = new RooUnfoldResponse(("h_respmatrix_calib_dijet_half1" + respmatrix_tags[i]).c_str(),""); h_respmatrix_calib_dijet_half1[i]->Setup(h_measure_calib_dijet_half1, h_truth_calib_dijet_half1);
        h_respmatrix_calib_dijet_half2[i] = new RooUnfoldResponse(("h_respmatrix_calib_dijet_half2" + respmatrix_tags[i]).c_str(),""); h_respmatrix_calib_dijet_half2[i]->Setup(h_measure_calib_dijet_half2, h_truth_calib_dijet_half2);
        h_respmatrix_calib_dijet_clus_smear[i] = new RooUnfoldResponse(("h_respmatrix_calib_dijet_clus_smear" + respmatrix_tags[i]).c_str(),""); h_respmatrix_calib_dijet_clus_smear[i]->Setup(h_measure_calib_dijet_clus_smear, h_truth_calib_dijet_clus_smear);
        h_respmatrix_calib_dijet_ohcal_mc_data_var[i] = new RooUnfoldResponse(("h_respmatrix_calib_dijet_ohcal_mc_data_var" + respmatrix_tags[i]).c_str(),""); h_respmatrix_calib_dijet_ohcal_mc_data_var[i]->Setup(h_measure_calib_dijet_ohcal_mc_data_var, h_truth_calib_dijet_ohcal_mc_data_var);
        h_respmatrix_calib_dijet_2sigma_noise[i] = new RooUnfoldResponse(("h_respmatrix_calib_dijet_2sigma_noise" + respmatrix_tags[i]).c_str(),""); h_respmatrix_calib_dijet_2sigma_noise[i]->Setup(h_measure_calib_dijet_2sigma_noise, h_truth_calib_dijet_2sigma_noise);
        h_respmatrix_calib_dijet_4sigma_noise[i] = new RooUnfoldResponse(("h_respmatrix_calib_dijet_4sigma_noise" + respmatrix_tags[i]).c_str(),""); h_respmatrix_calib_dijet_4sigma_noise[i]->Setup(h_measure_calib_dijet_4sigma_noise, h_truth_calib_dijet_4sigma_noise);
        //h_jetpt_respmatrix[i] = new RooUnfoldResponse(h_jetpt_measure, h_jetpt_truth,("h_jetpt_respmatrix" + respmatrix_tags[i]).c_str(),"");
        //h_jetpt_respmatrix[i]->UseDensityStatus();
        h_jetpt_respmatrix[i] = new RooUnfoldResponse(("h_jetpt_respmatrix" + respmatrix_tags[i]).c_str(),""); h_jetpt_respmatrix[i]->Setup(h_jetpt_measure, h_jetpt_truth);
        h_caloet_respmatrix[i] = new RooUnfoldResponse(("h_caloet_respmatrix" + respmatrix_tags[i]).c_str(),""); h_caloet_respmatrix[i]->Setup(h_caloet_measure, h_caloet_truth);
    }

    TH1D* h_calib_jet_pt_tight = new TH1D("h_calib_jet_pt_tight","",50,15,65);
    TH1D* h_truth_jet_pt_tight = new TH1D("h_truth_jet_pt_tight","",70,15,85);
    TH1D* h_calib_calo_et_tight = new TH1D("h_calib_calo_et_tight","",110,-2,20);
    TH1D* h_truth_calo_et_tight = new TH1D("h_truth_calo_et_tight","",100,0,20);

    TH1D* h_calib_jet_pt_uni_tight = new TH1D("h_calib_jet_pt_uni_tight","",50,0,1);
    TH1D* h_truth_jet_pt_uni_tight = new TH1D("h_truth_jet_pt_uni_tight","",70,0,1);
    TH1D* h_calib_calo_et_uni_tight = new TH1D("h_calib_calo_et_uni_tight","",110,0,1);
    TH1D* h_truth_calo_et_uni_tight = new TH1D("h_truth_calo_et_uni_tight","",100,0,1);

    TH1D* h_truthptbins = new TH1D("h_truthptbins","",truthnpt, 0, 1);
    TH1D* h_truthetbins = new TH1D("h_truthetbins","",truthnet, 0, 1);
    TH1D* h_calibptbins = new TH1D("h_calibptbins","",calibnpt, 0, 1);
    TH1D* h_calibetbins = new TH1D("h_calibetbins","",calibnet, 0, 1);
    
    ////////// Event Loop //////////
    std::cout << "Data analysis started." << std::endl;
    Long64_t nEntries = chain.GetEntries();
    std::cout << "Total number of events: " << nEntries << std::endl;
    
    // Event variables setup.
    bool reco_cut = true; bool truth_cut = true;
    float goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi;
    float calibjet_pt_dijet, calibjet_eta_dijet, calibjet_phi_dijet;
    float calibjet_pt_dijet_jesdown, calibjet_eta_dijet_jesdown, calibjet_phi_dijet_jesdown, calibjet_pt_dijet_jesup, calibjet_eta_dijet_jesup, calibjet_phi_dijet_jesup;
    float calibjet_pt_dijet_jerdown, calibjet_eta_dijet_jerdown, calibjet_phi_dijet_jerdown, calibjet_pt_dijet_jerup, calibjet_eta_dijet_jerup, calibjet_phi_dijet_jerup;
    bool calibjet_matched_dijet, calibjet_matched_dijet_jesdown, calibjet_matched_dijet_jesup, calibjet_matched_dijet_jerdown, calibjet_matched_dijet_jerup, qa_matched;
    double emfrac = 0;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 20; ++entry) {
        if (entry % 1000 == 0) cout << "event " << entry << endl;
        chain.GetEntry(entry);

        // Z-vertex cut.
        if (isnan(zvertex)) { continue; }
        if (fabs(zvertex) > 60) { continue; } 

        // Fill z-vertex histogram.
        h_zvertex_noreweight->Fill(zvertex);
        int zvertex_bin = h_vertex_weight->FindBin(zvertex);
        float vertex_weight = h_vertex_weight->GetBinContent(zvertex_bin);
        h_zvertex->Fill(zvertex, vertex_weight);

        //std::cout << "zvertex: " << zvertex << " vertex_bin: " << zvertex_bin << " vertex_weight: " << vertex_weight << std::endl;

        //////////////////////////// SETUP JET VARIABLES FOR UNFOLDING ////////////////////////////

        // indices to find leading and subleading jets 
        int ind_truth_lead = -1; int ind_truth_sub = -1; int ind_truth_subsub = -1; int ind_lead = -1; int ind_sub = -1; int ind_subsub = -1;
        float lead_e = 0; float truthlead_e = 0; float sub_e = 0; float truthsub_e = 0; float subsub_e = 0; float truthsubsub_e = 0;
        reco_cut = false; truth_cut = false;
        
        if (truthjet_pt->size() < 1) { continue; }
        get_leading_jet(ind_truth_lead, truthjet_pt); 
        if (truthjet_pt->at(ind_truth_lead) < truthjet_pt_min || truthjet_pt->at(ind_truth_lead) > truthjet_pt_max) { continue; } 

        //std::cout << "reco jets " << unsubjet_pt->size() << " truth jets " << truthjet_pt->size() << std::endl;

        // check number of jets above 5 GeV
        int Njet = 0;
        for (size_t i = 0; i < unsubjet_pt->size(); i++) {
          if (unsubjet_pt->at(i) >= 5.0) {
            Njet++;
          }
        }
        if (Njet >= 9) { 
            std::cout << "EVENT WITH 9 OR MORE JETS" << std::endl;
            continue; 
        }

        std::vector<float> truthe_new, truthpt_new, trutheta_new, truthphi_new;
        std::vector<float> recoe_new, recopt_new, recoeta_new, recophi_new, recoemcal_new, recoihcal_new, recoohcal_new, recocalibpt_new;

        for (size_t i = 0; i < truthjet_eta->size(); ++i) {
            if (std::fabs(truthjet_eta->at(i)) <= 0.7) {
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

        for (size_t i = 0; i < unsubjet_eta->size(); ++i) {
            if (!check_bad_jet_eta(unsubjet_eta->at(i), zvertex, jet_radius) && fabs(unsubjet_eta->at(i)) <= 0.7 && unsubjet_e->at(i) > 0.0) {
            //if (fabs(unsubjet_eta->at(i)) <= 0.7) {
                recoe_new.push_back(unsubjet_e->at(i));
                recopt_new.push_back(unsubjet_pt->at(i));
                recoeta_new.push_back(unsubjet_eta->at(i));
                recophi_new.push_back(unsubjet_phi->at(i));
                recoemcal_new.push_back(unsubjet_emcal_calo_e->at(i));
                recoihcal_new.push_back(unsubjet_ihcal_calo_e->at(i));
                recoohcal_new.push_back(unsubjet_ohcal_calo_e->at(i));
                recocalibpt_new.push_back(calibjet_pt->at(i));
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
        *calibjet_pt = std::move(recocalibpt_new);

        int nJetReq = 1;
        if (bkg_cut == "dijet") { nJetReq = 2; }
        // if both nreco jets < 1 and nTruthJet jets < 1, discard event
        if (unsubjet_pt->size() < nJetReq && truthjet_pt->size() < nJetReq) {
            continue;
        }

        emfrac = 0;
        TVector3 truthlead, truthsub, lead, sub;
        if (bkg_cut == "dijet") {
            if (truthjet_pt->size() >= nJetReq) { // if ntruth jets >= 2
                get_leading_subleading_jet(ind_truth_lead, ind_truth_sub, truthjet_pt); 
                if (truthjet_e->at(ind_truth_sub)/truthjet_e->at(ind_truth_lead) > 0.3 && match_leading_subleading_jet(truthjet_phi->at(ind_truth_lead), truthjet_phi->at(ind_truth_sub))) { // if leading and subleading jets match
                    h_exclusive_dijet_total->Fill(truthjet_pt->at(ind_truth_lead));
                    if (truthjet_pt->size() == 2) { // if ntruth jets == 2
                        truth_cut = true;
                        h_exclusive_dijet_passcut->Fill(truthjet_pt->at(ind_truth_lead));
                    } else { // if ntruth jets > 2
                        get_leading_subleading_subsubleading_jet(ind_truth_lead, ind_truth_sub, ind_truth_subsub, truthjet_pt); 
                        if (truthjet_e->at(ind_truth_subsub)/truthjet_e->at(ind_truth_lead) < 0.5) {
                            truth_cut = true;
                            h_exclusive_dijet_passcut->Fill(truthjet_pt->at(ind_truth_lead));
                        } else {
                            truth_cut = false;
                        }
                        //if (truthjet_pt->at(ind_truth_lead) > 17.0) { std::cout << "truthlead: " << truthjet_pt->at(ind_truth_lead) << " truthsub: " << truthjet_pt->at(ind_truth_sub) << " truthsubsub: " << truthjet_pt->at(ind_truth_subsub) << " truthsubsub/truthlead: " << truthjet_e->at(ind_truth_subsub)/truthjet_e->at(ind_truth_lead) << " truthsub/truthlead: " << truthjet_e->at(ind_truth_sub)/truthjet_e->at(ind_truth_lead) << " truth_cut: " << truth_cut << std::endl; }
                    }
                } else { // if leading and subleading jets do not match
                    truth_cut = false;
                }
                truthlead.SetPtEtaPhi(truthjet_pt->at(ind_truth_lead), truthjet_eta->at(ind_truth_lead), truthjet_phi->at(ind_truth_lead));
                truthlead_e = truthjet_e->at(ind_truth_lead);
                truthsub.SetPtEtaPhi(truthjet_pt->at(ind_truth_sub), truthjet_eta->at(ind_truth_sub), truthjet_phi->at(ind_truth_sub));
                truthsub_e = truthjet_e->at(ind_truth_sub);
            } else { // if ntruth jets < 2
                truthlead.SetPtEtaPhi(0,0,0);
                truthlead_e = 0;
                truthsub.SetPtEtaPhi(0,0,0);
                truthsub_e = 0;
                truth_cut = true;
            }
            if (unsubjet_pt->size() >= nJetReq) { // if nreco jets >= 2
                get_leading_subleading_jet(ind_lead, ind_sub, unsubjet_pt); 
                //std::cout << " sub/lead: " << unsubjet_e->at(ind_sub)/unsubjet_e->at(ind_lead) << " dPhi: " << get_dphi(unsubjet_phi->at(ind_lead), unsubjet_phi->at(ind_sub)) << std::endl; 
                //std::cout << "lead e: " << unsubjet_e->at(ind_lead) << " lead uncalib pt: " << unsubjet_pt->at(ind_lead) << " lead calib pt: " << calibjet_pt->at(ind_lead) << " sub e: " << unsubjet_e->at(ind_sub) << " sub uncalib pt: " << unsubjet_pt->at(ind_sub) << " sub calib pt: " << calibjet_pt->at(ind_sub) << std::endl;
                if (unsubjet_e->at(ind_sub)/unsubjet_e->at(ind_lead) > 0.3 && match_leading_subleading_jet(unsubjet_phi->at(ind_lead), unsubjet_phi->at(ind_sub))) { // if leading and subleading jets match
                    reco_cut = true;
                    if (unsubjet_pt->size() == 2) { // if nreco jets == 2
                        reco_cut = true;
                    } else { // if nreco jets > 2
                        get_leading_subleading_subsubleading_jet(ind_lead, ind_sub, ind_subsub, unsubjet_pt); 
                        if (unsubjet_e->at(ind_subsub)/unsubjet_e->at(ind_lead) < 0.5) {
                            reco_cut = true;
                        } else {
                            reco_cut = false;
                        }
                    }
                } else {
                    reco_cut = false;
                }
                //lead.SetPtEtaPhi(unsubjet_pt->at(ind_lead), unsubjet_eta->at(ind_lead), unsubjet_phi->at(ind_lead));
                lead.SetPtEtaPhi(calibjet_pt->at(ind_lead), unsubjet_eta->at(ind_lead), unsubjet_phi->at(ind_lead));
                lead_e = unsubjet_e->at(ind_lead);
                //sub.SetPtEtaPhi(unsubjet_pt->at(ind_sub), unsubjet_eta->at(ind_sub), unsubjet_phi->at(ind_sub));
                sub.SetPtEtaPhi(calibjet_pt->at(ind_sub), unsubjet_eta->at(ind_sub), unsubjet_phi->at(ind_sub));
                sub_e = unsubjet_e->at(ind_sub);
                emfrac = unsubjet_emcal_calo_e->at(ind_lead)/unsubjet_e->at(ind_lead);
            } else { // if nreco jets < 2
                lead.SetPtEtaPhi(0,0,0);
                lead_e = 0;
                sub.SetPtEtaPhi(0,0,0);
                sub_e = 0;
                reco_cut = true;
            }
        } else {
            //find truth leading if ntruth jets >= 1
            if (truthjet_pt->size() >= nJetReq) {
                get_leading_jet(ind_truth_lead, truthjet_pt);
                truthlead.SetPtEtaPhi(truthjet_pt->at(ind_truth_lead), truthjet_eta->at(ind_truth_lead), truthjet_phi->at(ind_truth_lead));
                truthlead_e = truthjet_e->at(ind_truth_lead);
            } else {
                truthlead.SetPtEtaPhi(0,0,0);
                truthlead_e = 0;
            }

            // find reco leading and subleading jets if nreco jets >= nJetReq
            if (unsubjet_pt->size() >= nJetReq) {
                get_leading_jet(ind_lead, unsubjet_pt);
                //lead.SetPtEtaPhi(unsubjet_pt->at(ind_lead), unsubjet_eta->at(ind_lead), unsubjet_phi->at(ind_lead));
                lead.SetPtEtaPhi(calibjet_pt->at(ind_lead), unsubjet_eta->at(ind_lead), unsubjet_phi->at(ind_lead));
                lead_e = unsubjet_e->at(ind_lead);
                emfrac = unsubjet_emcal_calo_e->at(ind_lead)/unsubjet_e->at(ind_lead);
                // edited to test unfolding procedure: apply efraction cut to leading jet
                  //std::cout << " EMCal frac: " << unsubjet_emcal_calo_e->at(ind_lead)/lead_e << " IHCal frac: " << unsubjet_ihcal_calo_e->at(ind_lead)/lead_e << " OHCal frac: " << unsubjet_ohcal_calo_e->at(ind_lead)/lead_e << std::endl;
                  if (unsubjet_emcal_calo_e->at(ind_lead)/lead_e > 0.1 && unsubjet_emcal_calo_e->at(ind_lead)/lead_e < 0.9 && unsubjet_ihcal_calo_e->at(ind_lead)/lead_e < 0.9 && unsubjet_ohcal_calo_e->at(ind_lead)/lead_e > 0.1 && unsubjet_ohcal_calo_e->at(ind_lead)/lead_e < 0.9) {
                    reco_cut = true; 
                  } else {
                    reco_cut = false;
                  }
            } else {
                lead.SetPtEtaPhi(0,0,0);
                lead_e = 0;
                reco_cut = true; // edited to test unfolding procedure
            }
            truth_cut = true;
            //reco_cut = true; edited to test unfolding procedure 
        }

        //std::cout << "Reco lead: pt " << lead.Pt() << " e " << lead_e << " eta " << lead.Eta() << " phi " << lead.Phi()  << " reco_cut " << reco_cut << std::endl;
        //std::cout << "Truth lead: pt " << truthlead.Pt() << " truth e " << truthlead_e << " truth eta " << truthlead.Eta() << " truth phi " << truthlead.Phi() << " truth_cut " << truth_cut << std::endl;

        get_truthjet(goodtruthjet_pt, goodtruthjet_eta, goodtruthjet_phi, truthlead.Pt(), truthlead.Eta(), truthlead.Phi(), truth_cut);
        get_calibjet(calibjet_pt_dijet, calibjet_eta_dijet, calibjet_phi_dijet, lead.Pt(), lead.Eta(), lead.Phi(), reco_cut, f_corr, 1, JER_smear);
        get_calibjet(calibjet_pt_dijet_jesdown, calibjet_eta_dijet_jesdown, calibjet_phi_dijet_jesdown, lead.Pt(), lead.Eta(), lead.Phi(), reco_cut, f_corr, 1.0-JES_var, JER_smear);
        get_calibjet(calibjet_pt_dijet_jesup, calibjet_eta_dijet_jesup, calibjet_phi_dijet_jesup, lead.Pt(), lead.Eta(), lead.Phi(), reco_cut, f_corr, 1.0+JES_var, JER_smear);
        get_calibjet(calibjet_pt_dijet_jerdown, calibjet_eta_dijet_jerdown, calibjet_phi_dijet_jerdown, lead.Pt(), lead.Eta(), lead.Phi(), reco_cut, f_corr, 1.0, JER_smear-JER_var);
        get_calibjet(calibjet_pt_dijet_jerup, calibjet_eta_dijet_jerup, calibjet_phi_dijet_jerup, lead.Pt(), lead.Eta(), lead.Phi(), reco_cut, f_corr, 1, JER_smear+JER_var);
        
        //std::cout << "Good reco lead: pt " << calibjet_pt_dijet << " eta " << calibjet_eta_dijet << " phi " << calibjet_phi_dijet << std::endl;
        //std::cout << "Good truth lead: pt " << goodtruthjet_pt << " truth eta " << goodtruthjet_eta << " truth phi " << goodtruthjet_phi << std::endl;
        
        match_meas_truth(calibjet_eta_dijet, calibjet_phi_dijet, calibjet_matched_dijet, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_dijet_jesdown, calibjet_phi_dijet_jesdown, calibjet_matched_dijet_jesdown, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_dijet_jesup, calibjet_phi_dijet_jesup, calibjet_matched_dijet_jesup, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_dijet_jerdown, calibjet_phi_dijet_jerdown, calibjet_matched_dijet_jerdown, goodtruthjet_eta, goodtruthjet_phi, jet_radius);
        match_meas_truth(calibjet_eta_dijet_jerup, calibjet_phi_dijet_jerup, calibjet_matched_dijet_jerup, goodtruthjet_eta, goodtruthjet_phi, jet_radius);

        //std::cout << "Match " << calibjet_matched_dijet << std::endl;
        /*
        if (bkg_cut == "dijet") std::cout << "Reco sub: pt " << sub.Pt() << " e " << sub_e << " eta " << sub.Eta() << " phi " << sub.Phi()  << " reco_cut " << reco_cut << std::endl;
        if (bkg_cut == "dijet") std::cout << "Truth sub: pt " << truthsub.Pt() << " truth e " << truthsub_e << " truth eta " << truthsub.Eta() << " truth phi " << truthsub.Phi() << " truth_cut " << truth_cut << std::endl;

        if (bkg_cut == "dijet" && calibjet_matched_dijet) { match_meas_truth(sub.Eta(), sub.Phi(), calibjet_matched_dijet, truthsub.Eta(), truthsub.Phi(), jet_radius); }  
        if (bkg_cut == "dijet" && calibjet_matched_dijet_jesdown) { match_meas_truth(sub.Eta(), sub.Phi(), calibjet_matched_dijet_jesdown, truthsub.Eta(), truthsub.Phi(), jet_radius); }       
        if (bkg_cut == "dijet" && calibjet_matched_dijet_jesup) { match_meas_truth(sub.Eta(), sub.Phi(), calibjet_matched_dijet_jesup, truthsub.Eta(), truthsub.Phi(), jet_radius); }       
        if (bkg_cut == "dijet" && calibjet_matched_dijet_jerdown) { match_meas_truth(sub.Eta(), sub.Phi(), calibjet_matched_dijet_jerdown, truthsub.Eta(), truthsub.Phi(), jet_radius); }       
        if (bkg_cut == "dijet" && calibjet_matched_dijet_jerup) { match_meas_truth(sub.Eta(), sub.Phi(), calibjet_matched_dijet_jerup, truthsub.Eta(), truthsub.Phi(), jet_radius); }       
        */
        //std::cout << "Match after subleading match " << calibjet_matched_dijet << std::endl;

        match_meas_truth(lead.Eta(),lead.Phi(),qa_matched,truthlead.Eta(),truthlead.Phi(),jet_radius);
        double calib_pt_qa = lead.Pt() * (1 + randGen2.Gaus(0.0, JER_smear)) * 1;
        //std::cout << "calib_pt_qa: " << calib_pt_qa << " calibjet_pt_dijet: " << calibjet_pt_dijet << std::endl;
        if (truthlead.Pt() >= truthptbins[0] && truthlead.Pt() < truthptbins[truthnpt] && qa_matched) {
            //std::cout << "matched: " << qa_matched << " respmatrix match: " << calibjet_matched_dijet << " truth pt: " << truthlead.Pt() << " ratio: " << f_corr->Eval(lead.Pt())/truthlead.Pt() << " weight: " << weight_scale << std::endl;
            h_jes_qa->Fill(truthlead.Pt(), lead.Pt()/truthlead.Pt(), weight_scale*vertex_weight);
            h_jes_qa_reco->Fill(lead.Pt(), lead.Pt()/truthlead.Pt(), weight_scale*vertex_weight);
        }
        /*
        bool truth_match = false;
        float truth_deltaR0 = 999.0; float truth_deltaR = 999.0;
        for (int i = 0; i < unsubjet_pt->size(); i++) {
            find_match_meas_truth(unsubjet_eta->at(i),unsubjet_phi->at(i),truth_match,truth_deltaR,truthlead.Eta(),truthlead.Phi(),jet_radius);
            if (truth_deltaR < truth_deltaR0) { truth_deltaR0 = truth_deltaR; }
        }
        if (truthlead.Pt() >= truthptbins[0] && truthlead.Pt() < truthptbins[truthnpt] && truth_deltaR0 < 0.75*jet_radius) {
            efficiency_truth_match_test->Fill(truthlead.Pt(), weight_scale);
        }
        */
        //////////////////////////// SETUP UE VARIABLES FOR UNFOLDING //////////////////////////////

        // find reco ET information 
        float et_transverse = 0;
        float et_transverse_clus_smear = 0;
        float et_transverse_ohcal_mc_data_var = 0;
        float et_transverse_2sigma = 0;
        float et_transverse_4sigma = 0;
        if (!clusters) {
            /*
            for (int i = 0; i < emcaln; i++) {
                float dphi = get_dphi(lead.Phi(),emcalphi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += emcale[i]/cosh(emcaleta[i]); } 
            }
            for (int i = 0; i < ihcaln; i++) {
                float dphi = get_dphi(lead.Phi(),ihcalphi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += ihcale[i]/cosh(ihcaleta[i]); } 
            }
            for (int i = 0; i < ohcaln; i++) {
                float dphi = get_dphi(lead.Phi(),ohcalphi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += ohcale[i]/cosh(ohcaleta[i]); } 
            }
            */
        } else {
            
            for (int i = 0; i < clsmult; i++) {
                float dphi = get_dphi(lead.Phi(),cluster_phi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0 && cluster_e[i]/cluster_eta[i] > -1.0) { 
                    et_transverse += cluster_e[i]/cosh(cluster_eta[i]); 
                    et_transverse_clus_smear += (cluster_e[i] * (1.0 + clusERandGen.Gaus(0,0.08)))/(cosh(cluster_eta[i]));
                }
            }
            
            for (int i = 0; i < clsmult2; i++) {
                float dphi = get_dphi(lead.Phi(),cluster2_phi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0 && cluster2_e[i]/cluster2_eta[i] > -1.0) { 
                    et_transverse_2sigma += cluster2_e[i]/cosh(cluster2_eta[i]); 
                }
            }

            for (int i = 0; i < clsmult4; i++) {
                float dphi = get_dphi(lead.Phi(),cluster4_phi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0 && cluster4_e[i]/cluster4_eta[i] > -1.0) { 
                    et_transverse_4sigma += cluster4_e[i]/cosh(cluster4_eta[i]); 
                }
            }
            
            et_transverse_ohcal_mc_data_var = et_transverse;
            for (int i = 0; i < clsmult; i++) {
                float dphi = get_dphi(lead.Phi(),cluster_phi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0 && cluster_e[i]/cluster_eta[i] > -1.0) {
                  for (int j = 0; j < cluster_ntowers[i]; j++) {
                    if (cluster_tower_calo[i][j] == 2) {
                        et_transverse_ohcal_mc_data_var += (ohcal_scale_eta[cluster_tower_ieta[i][j]]-1.0)*cluster_tower_e[i][j]/(cosh(cluster_eta[i]));
                    }
                  }
                }
            }
        }

        // find truth ET information
        float truth_et_transverse = 0;
        for (int i = 0; i < truthpar_n; i++) {
            if (fabs(truthpar_eta[i]) > 1.1) { continue; }
            if ((truthpar_pid[i] == 22 || truthpar_pid[i] == 111) && fabs(truthpar_e[i]) < 0.2) { continue; } // edited back to 0.2 for sPHENIX primary particle list
            else if (fabs(truthpar_e[i]) < 0.5) { continue; } // edited back to 0.5
            float dphi = get_dphi(truthlead.Phi(),truthpar_phi[i]);
            if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { 
                truth_et_transverse += truthpar_e[i]/cosh(truthpar_eta[i]); 
                if (truthlead.Pt() >= truthptbins[0] && truthlead.Pt() < truthptbins[truthnpt]) {
                    h_truth_avg_pt_vs_jet_pt->Fill(truthlead.Pt(), truthpar_pt[i]);
                }
            } 
        }
        float thres_truth_et_transverse = 0;
        for (int i = 0; i < truthpar_n; i++) {
            if (fabs(truthpar_eta[i]) > 1.1) { continue; }
            if (fabs(truthpar_e[i]) < 0.4352) { continue; }
            float dphi = get_dphi(truthlead.Phi(),truthpar_phi[i]);
            if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { thres_truth_et_transverse += truthpar_e[i]/cosh(truthpar_eta[i]); } 
        }
        float reco_truth_et_transverse = 0;
        for (int i = 0; i < truthpar_n; i++) {
            if (fabs(truthpar_eta[i]) > 1.1) { continue; }
            if ((truthpar_pid[i] == 22 || truthpar_pid[i] == 111) && fabs(truthpar_e[i]) > 0.2) {
                float dphi = get_dphi(truthlead.Phi(),truthpar_phi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { reco_truth_et_transverse += truthpar_e[i]/cosh(truthpar_eta[i]); } 
            } else if (fabs(truthpar_e[i] > 0.2)) {
                float dphi = get_dphi(truthlead.Phi(),truthpar_phi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { reco_truth_et_transverse += 0.36*truthpar_e[i]/cosh(truthpar_eta[i]); } 
            }
        }
        float reco_thres_truth_et_transverse = 0;
        for (int i = 0; i < truthpar_n; i++) {
            if (fabs(truthpar_eta[i]) > 1.1) { continue; }
            if ((truthpar_pid[i] == 22 || truthpar_pid[i] == 111) && fabs(truthpar_e[i]) > 0.4352) {
                float dphi = get_dphi(truthlead.Phi(),truthpar_phi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { reco_thres_truth_et_transverse += truthpar_e[i]/cosh(truthpar_eta[i]); } 
            } else if (fabs(truthpar_e[i] > 0.4352)) {
                float dphi = get_dphi(truthlead.Phi(),truthpar_phi[i]);
                if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { reco_thres_truth_et_transverse += 0.36*truthpar_e[i]/cosh(truthpar_eta[i]); } 
            }
        }

        // can toggle reco pT cut, reco bkg cut and match requirement for efficiency 
        if (truthlead.Pt() >= truthptbins[0] && truthlead.Pt() < truthptbins[truthnpt]) {
            h_truth_qa->Fill(truthlead.Pt(), truth_et_transverse, weight_scale*vertex_weight);
            if (reco_cut) { h_efficiency_reco_bkg_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale*vertex_weight); }
            if (qa_matched) { h_efficiency_match->Fill(truthlead.Pt(), truth_et_transverse, weight_scale*vertex_weight); }
            if (calib_pt_qa >= calibptbins[0] && calib_pt_qa < calibptbins[calibnpt]) { h_efficiency_reco_pt_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale*vertex_weight); }
            if (reco_cut && qa_matched) { h_efficiency_match_reco_bkg_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale*vertex_weight); }
            if (reco_cut && calib_pt_qa >= calibptbins[0] && calib_pt_qa < calibptbins[calibnpt]) { h_efficiency_reco_pt_cut_reco_bkg_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale*vertex_weight); }
            if (qa_matched && calib_pt_qa >= calibptbins[0] && calib_pt_qa < calibptbins[calibnpt]) { h_efficiency_match_reco_pt_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale*vertex_weight); }
            if (reco_cut && qa_matched && calib_pt_qa >= calibptbins[0] && calib_pt_qa < calibptbins[calibnpt]) { h_efficiency_match_reco_pt_cut_reco_bkg_cut->Fill(truthlead.Pt(), truth_et_transverse, weight_scale*vertex_weight); }
        }

        // can toggle truth pT cut and match requirement for purity 
        if (calibjet_pt_dijet >= calibptbins[0] && calibjet_pt_dijet < calibptbins[calibnpt] && reco_cut) {
            h_measure_qa->Fill(calibjet_pt_dijet, et_transverse, weight_scale*vertex_weight);
            if (truthlead.Pt() >= truthptbins[0] && truthlead.Pt() < truthptbins[truthnpt]) { h_purity_truth_pt_cut->Fill(calibjet_pt_dijet, et_transverse, weight_scale*vertex_weight); }
            if (calibjet_matched_dijet) { h_purity_match->Fill(calibjet_pt_dijet, et_transverse, weight_scale*vertex_weight); }
        }

        //////////////////////////// RECORD QA PLOTS FOR ALL EVENTS IN UNFOLDING PROCEDURE ////////////////////////////
        if ((goodtruthjet_pt > truthptbins[0] && truth_cut) || (calibjet_pt_dijet > calibptbins[0] && reco_cut)) {
            h_lead_spectra_record->Fill(calibjet_pt_dijet, weight_scale*vertex_weight);
            h_lead_truth_spectra_record->Fill(truthlead.Pt(), weight_scale*vertex_weight);
        }
        if (reco_cut) {
            if (calibjet_pt_dijet > calibptbins[0]) {
                //h_uncalib_lead_spectra_record->Fill(unsubjet_pt->at(ind_lead), weight_scale*vertex_weight);
                h_lead_spectra_record->Fill(calib_pt_qa, weight_scale*vertex_weight);
                h_nw_et_transverse_record->Fill(et_transverse, weight_scale*vertex_weight);
                h_et_transverse_record->Fill(et_transverse, weight_scale*vertex_weight);
                h_ue_pt_transverse_record->Fill(lead.Pt(),et_transverse, weight_scale*vertex_weight);
                h_ue_pt_transverse_record_full->Fill(lead.Pt(),et_transverse, weight_scale*vertex_weight);
                h_efrac_lead_pt_dijet->Fill(calibjet_pt_dijet, emfrac, weight_scale*vertex_weight);
            }
        }
            
        if (goodtruthjet_pt > truthptbins[0]) {
            //h_lead_truth_spectra_record->Fill(truthlead.Pt(), weight_scale*vertex_weight);
            h_nw_et_truth_transverse_record->Fill(truth_et_transverse, weight_scale*vertex_weight);
            h_et_truth_transverse_record->Fill(truth_et_transverse, weight_scale*vertex_weight);
            h_ue_pt_truth_transverse_record->Fill(truthlead.Pt(),truth_et_transverse, weight_scale*vertex_weight);
            h_ue_pt_truth_transverse_record_full->Fill(truthlead.Pt(),truth_et_transverse, weight_scale*vertex_weight);
            h_thres_nw_et_truth_transverse_record->Fill(thres_truth_et_transverse, weight_scale*vertex_weight);
            h_thres_et_truth_transverse_record->Fill(thres_truth_et_transverse, weight_scale*vertex_weight);
            h_thres_ue_pt_truth_transverse_record->Fill(truthlead.Pt(),thres_truth_et_transverse, weight_scale*vertex_weight);
            h_reco_nw_et_truth_transverse_record->Fill(reco_truth_et_transverse, weight_scale*vertex_weight);
            h_reco_et_truth_transverse_record->Fill(reco_truth_et_transverse, weight_scale*vertex_weight);
            h_reco_ue_pt_truth_transverse_record->Fill(truthlead.Pt(),reco_truth_et_transverse, weight_scale*vertex_weight);
            h_reco_thres_nw_et_truth_transverse_record->Fill(reco_thres_truth_et_transverse, weight_scale*vertex_weight);
            h_reco_thres_et_truth_transverse_record->Fill(reco_thres_truth_et_transverse, weight_scale*vertex_weight);
            h_reco_thres_ue_pt_truth_transverse_record->Fill(truthlead.Pt(),reco_thres_truth_et_transverse, weight_scale*vertex_weight);
        }

        double uni_meas_et = MapToUniform(et_transverse, calibetbins, calibnet);
        double uni_truth_et = MapToUniform(truth_et_transverse, truthetbins, truthnet);
        double uni_meas_pt = MapToUniform(calibjet_pt_dijet, calibptbins, calibnpt);
        double uni_truth_pt = MapToUniform(goodtruthjet_pt, truthptbins, truthnpt);
        double uni_meas_pt_jesdown = MapToUniform(calibjet_pt_dijet_jesdown, calibptbins, calibnpt);
        double uni_meas_pt_jesup = MapToUniform(calibjet_pt_dijet_jesup, calibptbins, calibnpt);
        double uni_meas_pt_jerdown = MapToUniform(calibjet_pt_dijet_jerdown, calibptbins, calibnpt);
        double uni_meas_pt_jerup = MapToUniform(calibjet_pt_dijet_jerup, calibptbins, calibnpt);
        double uni_meas_et_clus_smear = MapToUniform(et_transverse_clus_smear, calibetbins, calibnet);
        double uni_meas_et_ohcal_mc_data_var = MapToUniform(et_transverse_ohcal_mc_data_var, calibetbins, calibnet);
        double uni_meas_et_2sigma_noise = MapToUniform(et_transverse_2sigma, calibetbins, calibnet);
        double uni_meas_et_4sigma_noise = MapToUniform(et_transverse_4sigma, calibetbins, calibnet);

        //std::cout << "Uniform values: meas pt: " << uni_meas_pt << " meas et: " << uni_meas_et << " truth pt: " << uni_truth_pt << " truth et: " << uni_truth_et << std::endl;

        if (uni_meas_pt >= 0 && uni_meas_et >= 0 && uni_meas_et <= 1) {
            h_calib_jet_pt_tight->Fill(calibjet_pt_dijet, weight_scale*vertex_weight);
            h_calib_calo_et_tight->Fill(et_transverse, weight_scale*vertex_weight);
            h_calib_jet_pt_uni_tight->Fill(uni_meas_pt, weight_scale*vertex_weight);
            h_calib_calo_et_uni_tight->Fill(uni_meas_et, weight_scale*vertex_weight);
        }
        if (uni_truth_pt >= 0 && uni_truth_et >= 0 && uni_truth_et <= 1) {
            h_truth_jet_pt_tight->Fill(goodtruthjet_pt, weight_scale*vertex_weight);
            h_truth_calo_et_tight->Fill(truth_et_transverse, weight_scale*vertex_weight);
            h_truth_jet_pt_uni_tight->Fill(uni_truth_pt, weight_scale*vertex_weight);
            h_truth_calo_et_uni_tight->Fill(uni_truth_et, weight_scale*vertex_weight);
        }

        //////////////////////////// FILL RESPONSE MATRICES ////////////////////////////
        fill_response_matrix_full(h_truth_calib_dijet, h_measure_calib_dijet, h_respmatrix_calib_dijet[0], h_fake_calib_dijet, h_miss_calib_dijet, h_respmatrix_calib_dijet[1], h_counts_fake_calib_dijet, h_counts_miss_calib_dijet, h_jetpt_respmatrix[0], h_caloet_respmatrix[0], h_jetpt_respmatrix[1], h_caloet_respmatrix[1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight);
        fill_response_matrix(h_truth_calib_dijet_jesdown, h_measure_calib_dijet_jesdown, h_respmatrix_calib_dijet_jesdown[0], h_fake_calib_dijet_jesdown, h_miss_calib_dijet_jesdown, h_respmatrix_calib_dijet_jesdown[1], uni_meas_pt_jesdown, calibjet_matched_dijet_jesdown, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight);
        fill_response_matrix(h_truth_calib_dijet_jesup, h_measure_calib_dijet_jesup, h_respmatrix_calib_dijet_jesup[0], h_fake_calib_dijet_jesup, h_miss_calib_dijet_jesup, h_respmatrix_calib_dijet_jesup[1], uni_meas_pt_jesup, calibjet_matched_dijet_jesup, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight);
        fill_response_matrix(h_truth_calib_dijet_jerdown, h_measure_calib_dijet_jerdown, h_respmatrix_calib_dijet_jerdown[0], h_fake_calib_dijet_jerdown, h_miss_calib_dijet_jerdown, h_respmatrix_calib_dijet_jerdown[1], uni_meas_pt_jerdown, calibjet_matched_dijet_jerdown, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight);
        fill_response_matrix(h_truth_calib_dijet_jerup, h_measure_calib_dijet_jerup, h_respmatrix_calib_dijet_jerup[0], h_fake_calib_dijet_jerup, h_miss_calib_dijet_jerup, h_respmatrix_calib_dijet_jerup[1], uni_meas_pt_jerup, calibjet_matched_dijet_jerup, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight);
        fill_response_matrix(h_truth_calib_dijet_clus_smear, h_measure_calib_dijet_clus_smear, h_respmatrix_calib_dijet_clus_smear[0], h_fake_calib_dijet_clus_smear, h_miss_calib_dijet_clus_smear, h_respmatrix_calib_dijet_clus_smear[1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_clus_smear, uni_truth_et, weight_scale*vertex_weight);
        fill_response_matrix(h_truth_calib_dijet_ohcal_mc_data_var, h_measure_calib_dijet_ohcal_mc_data_var, h_respmatrix_calib_dijet_ohcal_mc_data_var[0], h_fake_calib_dijet_ohcal_mc_data_var, h_miss_calib_dijet_ohcal_mc_data_var, h_respmatrix_calib_dijet_ohcal_mc_data_var[1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_ohcal_mc_data_var, uni_truth_et, weight_scale*vertex_weight);
        fill_response_matrix(h_truth_calib_dijet_2sigma_noise, h_measure_calib_dijet_2sigma_noise, h_respmatrix_calib_dijet_2sigma_noise[0], h_fake_calib_dijet_2sigma_noise, h_miss_calib_dijet_2sigma_noise, h_respmatrix_calib_dijet_2sigma_noise[1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_2sigma_noise, uni_truth_et, weight_scale*vertex_weight);
        fill_response_matrix(h_truth_calib_dijet_4sigma_noise, h_measure_calib_dijet_4sigma_noise, h_respmatrix_calib_dijet_4sigma_noise[0], h_fake_calib_dijet_4sigma_noise, h_miss_calib_dijet_4sigma_noise, h_respmatrix_calib_dijet_4sigma_noise[1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_4sigma_noise, uni_truth_et, weight_scale*vertex_weight);
        if (entry % 2 == 0){ fill_response_matrix(h_truth_calib_dijet_half1, h_measure_calib_dijet_half1, h_respmatrix_calib_dijet_half1[0], h_fake_calib_dijet_half1, h_miss_calib_dijet_half1, h_respmatrix_calib_dijet_half1[1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight); }   
        else { fill_response_matrix(h_truth_calib_dijet_half2, h_measure_calib_dijet_half2, h_respmatrix_calib_dijet_half2[0], h_fake_calib_dijet_half2, h_miss_calib_dijet_half2, h_respmatrix_calib_dijet_half2[1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight); }

        if (iter > 1) {
            for (int i = 1; i < 3; i++) {
                fill_trim_response_matrix(h_respmatrix_calib_dijet[i+1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[0], counts_measured[0], counts_truth[0], trim_val[i]);
                fill_trim_response_matrix(h_respmatrix_calib_dijet_jesdown[i+1], uni_meas_pt_jesdown, calibjet_matched_dijet_jesdown, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[1], counts_measured[1], counts_truth[1], trim_val[i]);
                fill_trim_response_matrix(h_respmatrix_calib_dijet_jesup[i+1], uni_meas_pt_jesup, calibjet_matched_dijet_jesup, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[2], counts_measured[2], counts_truth[2], trim_val[i]);
                fill_trim_response_matrix(h_respmatrix_calib_dijet_jerdown[i+1], uni_meas_pt_jerdown, calibjet_matched_dijet_jerdown, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[3], counts_measured[3], counts_truth[3], trim_val[i]);
                fill_trim_response_matrix(h_respmatrix_calib_dijet_jerup[i+1], uni_meas_pt_jerup, calibjet_matched_dijet_jerup, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[4], counts_measured[4], counts_truth[4], trim_val[i]);
                fill_trim_response_matrix(h_respmatrix_calib_dijet_clus_smear[i+1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_clus_smear, uni_truth_et, weight_scale*vertex_weight, counts_matrix[5], counts_measured[5], counts_truth[5], trim_val[i]);
                fill_trim_response_matrix(h_respmatrix_calib_dijet_ohcal_mc_data_var[i+1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_ohcal_mc_data_var, uni_truth_et, weight_scale*vertex_weight, counts_matrix[6], counts_measured[6], counts_truth[6], trim_val[i]);
                fill_trim_response_matrix(h_respmatrix_calib_dijet_2sigma_noise[i+1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_2sigma_noise, uni_truth_et, weight_scale*vertex_weight, counts_matrix[7], counts_measured[7], counts_truth[7], trim_val[i]);
                fill_trim_response_matrix(h_respmatrix_calib_dijet_4sigma_noise[i+1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_4sigma_noise, uni_truth_et, weight_scale*vertex_weight, counts_matrix[8], counts_measured[8], counts_truth[8], trim_val[i]);
                if (entry % 2 == 0) { fill_trim_response_matrix(h_respmatrix_calib_dijet_half1[i+1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[9], counts_measured[9], counts_truth[9], trim_val[i]); }
                else { fill_trim_response_matrix(h_respmatrix_calib_dijet_half2[i+1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[10], counts_measured[10], counts_truth[10], trim_val[i]); }
                fill_trim_1D_response_matrices(h_jetpt_respmatrix[i+1], h_caloet_respmatrix[i+1], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[0], counts_measured[0], counts_truth[0], trim_val[i]);
            }
        }

        if (iter > 2) {
            for (int i = 0; i < 3; i++) {
                fill_reweighted_trim_response_matrix(h_respmatrix_calib_dijet[i+4], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[0], counts_measured[0], counts_truth[0], trim_val[i], weights[0][i]);
                fill_reweighted_trim_response_matrix(h_respmatrix_calib_dijet_jesdown[i+4], uni_meas_pt_jesdown, calibjet_matched_dijet_jesdown, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[1], counts_measured[1], counts_truth[1], trim_val[i], weights[1][i]);
                fill_reweighted_trim_response_matrix(h_respmatrix_calib_dijet_jesup[i+4], uni_meas_pt_jesup, calibjet_matched_dijet_jesup, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[2], counts_measured[2], counts_truth[2], trim_val[i], weights[2][i]);
                fill_reweighted_trim_response_matrix(h_respmatrix_calib_dijet_jerdown[i+4], uni_meas_pt_jerdown, calibjet_matched_dijet_jerdown, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[3], counts_measured[3], counts_truth[3], trim_val[i], weights[3][i]);
                fill_reweighted_trim_response_matrix(h_respmatrix_calib_dijet_jerup[i+4], uni_meas_pt_jerup, calibjet_matched_dijet_jerup, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[4], counts_measured[4], counts_truth[4], trim_val[i], weights[4][i]);
                fill_reweighted_trim_response_matrix(h_respmatrix_calib_dijet_clus_smear[i+4], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_clus_smear, uni_truth_et, weight_scale*vertex_weight, counts_matrix[5], counts_measured[5], counts_truth[5], trim_val[i], weights[5][i]);
                fill_reweighted_trim_response_matrix(h_respmatrix_calib_dijet_ohcal_mc_data_var[i+4], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_ohcal_mc_data_var, uni_truth_et, weight_scale*vertex_weight, counts_matrix[6], counts_measured[6], counts_truth[6], trim_val[i], weights[6][i]);
                fill_reweighted_trim_response_matrix(h_respmatrix_calib_dijet_2sigma_noise[i+4], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_2sigma_noise, uni_truth_et, weight_scale*vertex_weight, counts_matrix[7], counts_measured[7], counts_truth[7], trim_val[i], weights[7][i]);
                fill_reweighted_trim_response_matrix(h_respmatrix_calib_dijet_4sigma_noise[i+4], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et_4sigma_noise, uni_truth_et, weight_scale*vertex_weight, counts_matrix[8], counts_measured[8], counts_truth[8], trim_val[i], weights[8][i]);
                if (entry % 2 == 0) { fill_reweighted_trim_response_matrix(h_respmatrix_calib_dijet_half1[i+4], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[9], counts_measured[9], counts_truth[9], trim_val[i], weights[9][i]); }
                else { fill_reweighted_trim_response_matrix(h_respmatrix_calib_dijet_half2[i+4], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[10], counts_measured[10], counts_truth[10], trim_val[i], weights[10][i]); }
                fill_reweighted_trim_1D_response_matrices(h_jetpt_respmatrix[i+4], h_caloet_respmatrix[i+4], uni_meas_pt, calibjet_matched_dijet, uni_truth_pt, uni_meas_et, uni_truth_et, weight_scale*vertex_weight, counts_matrix[0], counts_measured[0], counts_truth[0], trim_val[i], weights1D[0][i], weights1D[1][i]);
            } 
        }

        if (entry == 10113 && reco_cut && calibjet_pt_dijet) {
            //std::cout << "Event: " << entry << " lead PT " << lead.Pt() << " lead phi " << lead.Phi() << " lead eta " << lead.Eta() << " et_transverse " << et_transverse << std::endl;
            for (int i = 0; i < clsmult; i++) {
                std::cout << "event " << entry << " cluster id " << i << " e " << cluster_e[i] << " ntowers " << cluster_ntowers[i] << " eta " << cluster_eta[i] << " phi " << cluster_phi[i] << std::endl;
                for (int n = 0; n < cluster_ntowers[i]; n++) {
                    if (cluster_tower_calo[i][n] == 1) {
                        h_total_topo_emcal_event_display->Fill(cluster_tower_ieta[i][n],cluster_tower_iphi[i][n],cluster_tower_e[i][n]);
                        h_total_topo_event_display->Fill(cluster_tower_ieta[i][n]/4,cluster_tower_iphi[i][n]/4,cluster_tower_e[i][n]);
                        h_topo_emcal_event_display[i]->Fill(cluster_tower_ieta[i][n],cluster_tower_iphi[i][n]);
                        h_topo_event_display[i]->Fill(cluster_tower_ieta[i][n]/4,cluster_tower_iphi[i][n]/4);
                    } else if (cluster_tower_calo[i][n] == 2) {
                        h_total_topo_ohcal_event_display->Fill(cluster_tower_ieta[i][n],cluster_tower_iphi[i][n],cluster_tower_e[i][n]);
                        h_total_topo_event_display->Fill(cluster_tower_ieta[i][n],cluster_tower_iphi[i][n],cluster_tower_e[i][n]);
                        h_topo_ohcal_event_display[i]->Fill(cluster_tower_ieta[i][n],cluster_tower_iphi[i][n]);
                        h_topo_event_display[i]->Fill(cluster_tower_ieta[i][n],cluster_tower_iphi[i][n]);
                    } else if (cluster_tower_calo[i][n] == 3) {
                        h_total_topo_ihcal_event_display->Fill(cluster_tower_ieta[i][n],cluster_tower_iphi[i][n],cluster_tower_e[i][n]);
                        h_total_topo_event_display->Fill(cluster_tower_ieta[i][n],cluster_tower_iphi[i][n],cluster_tower_e[i][n]);
                        h_topo_ihcal_event_display[i]->Fill(cluster_tower_ieta[i][n],cluster_tower_iphi[i][n]);
                        h_topo_event_display[i]->Fill(cluster_tower_ieta[i][n],cluster_tower_iphi[i][n]);
                    }
                    std::cout << "cluster tower: calo " << cluster_tower_calo[i][n] << " ieta " << cluster_tower_ieta[i][n] << " iphi " << cluster_tower_iphi[i][n] << " energy " << cluster_tower_e[i][n] << std::endl;
                }
                std::cout << std::endl;

            }
        }

    }

    OutputVarBinMapping(h_truthptbins, truthptbins, int(sizeof(truthptbins)/sizeof(truthptbins[0])));
    OutputVarBinMapping(h_truthetbins, truthetbins, int(sizeof(truthetbins)/sizeof(truthetbins[0])));
    OutputVarBinMapping(h_calibptbins, calibptbins, int(sizeof(calibptbins)/sizeof(calibptbins[0])));
    OutputVarBinMapping(h_calibetbins, calibetbins, int(sizeof(calibetbins)/sizeof(calibetbins[0])));

    std::cout << "Writing histograms..." << std::endl;
    f_out->cd();

    h_zvertex->Write(); h_zvertex_noreweight->Write(); h_jes_qa->Write(); h_jes_qa_reco->Write(); h_truth_qa->Write(); h_measure_qa->Write(); h_uncalib_lead_spectra_record->Write();
    efficiency_truth_match_test->Write();
    h_calib_jet_pt_tight->Write(); h_truth_jet_pt_tight->Write(); h_calib_calo_et_tight->Write(); h_truth_calo_et_tight->Write();
    h_calib_jet_pt_uni_tight->Write(); h_truth_jet_pt_uni_tight->Write(); h_calib_calo_et_uni_tight->Write(); h_truth_calo_et_uni_tight->Write();
    h_efficiency_match_reco_bkg_cut->Write(); h_efficiency_reco_pt_cut_reco_bkg_cut->Write(); h_efficiency_match_reco_pt_cut_reco_bkg_cut->Write(); h_efficiency_match_reco_pt_cut->Write();
    h_purity_match->Write(); h_purity_truth_pt_cut->Write(); h_efficiency_match->Write(); h_efficiency_reco_pt_cut->Write(); h_efficiency_reco_bkg_cut->Write(); 
    h_lead_spectra_record->Write(); h_et_transverse_record->Write(); h_nw_et_transverse_record->Write(); h_ue_pt_transverse_record->Write(); 
    h_lead_truth_spectra_record->Write(); h_et_truth_transverse_record->Write(); h_nw_et_truth_transverse_record->Write(); h_ue_pt_truth_transverse_record->Write(); 
    h_thres_et_truth_transverse_record->Write(); h_thres_nw_et_truth_transverse_record->Write(); h_thres_ue_pt_truth_transverse_record->Write(); 
    h_reco_et_truth_transverse_record->Write(); h_reco_nw_et_truth_transverse_record->Write(); h_reco_ue_pt_truth_transverse_record->Write(); 
    h_reco_thres_et_truth_transverse_record->Write(); h_reco_thres_nw_et_truth_transverse_record->Write(); h_reco_thres_ue_pt_truth_transverse_record->Write(); 
    h_truth_calib_dijet->Write(); h_measure_calib_dijet->Write(); h_fake_calib_dijet->Write(); h_miss_calib_dijet->Write(); h_counts_fake_calib_dijet->Write(); h_counts_miss_calib_dijet->Write(); 
    h_truth_calib_dijet_jesdown->Write(); h_measure_calib_dijet_jesdown->Write(); h_fake_calib_dijet_jesdown->Write(); h_miss_calib_dijet_jesdown->Write();
    h_truth_calib_dijet_jesup->Write(); h_measure_calib_dijet_jesup->Write(); h_fake_calib_dijet_jesup->Write(); h_miss_calib_dijet_jesup->Write();
    h_truth_calib_dijet_jerdown->Write(); h_measure_calib_dijet_jerdown->Write();  h_fake_calib_dijet_jerdown->Write(); h_miss_calib_dijet_jerdown->Write();
    h_truth_calib_dijet_jerup->Write(); h_measure_calib_dijet_jerup->Write(); h_fake_calib_dijet_jerup->Write(); h_miss_calib_dijet_jerup->Write();
    h_truth_calib_dijet_clus_smear->Write(); h_measure_calib_dijet_clus_smear->Write();  h_fake_calib_dijet_clus_smear->Write(); h_miss_calib_dijet_clus_smear->Write();
    h_truth_calib_dijet_ohcal_mc_data_var->Write(); h_measure_calib_dijet_ohcal_mc_data_var->Write();  h_fake_calib_dijet_ohcal_mc_data_var->Write(); h_miss_calib_dijet_ohcal_mc_data_var->Write();
    h_truth_calib_dijet_2sigma_noise->Write(); h_measure_calib_dijet_2sigma_noise->Write();  h_fake_calib_dijet_2sigma_noise->Write(); h_miss_calib_dijet_2sigma_noise->Write();
    h_truth_calib_dijet_4sigma_noise->Write(); h_measure_calib_dijet_4sigma_noise->Write();  h_fake_calib_dijet_4sigma_noise->Write(); h_miss_calib_dijet_4sigma_noise->Write();
    h_truth_calib_dijet_half1->Write(); h_measure_calib_dijet_half1->Write(); h_fake_calib_dijet_half1->Write(); h_miss_calib_dijet_half1->Write();
    h_truth_calib_dijet_half2->Write(); h_measure_calib_dijet_half2->Write(); h_fake_calib_dijet_half2->Write(); h_miss_calib_dijet_half2->Write();
    h_efrac_lead_pt_dijet->Write();
    h_truth_avg_pt_vs_jet_pt->Write();
    h_ue_pt_transverse_record_full->Write(); h_ue_pt_truth_transverse_record_full->Write();
    h_exclusive_dijet_passcut->Write(); h_exclusive_dijet_total->Write();

    h_total_topo_event_display->Write(); h_total_topo_emcal_event_display->Write(); h_total_topo_ihcal_event_display->Write(); h_total_topo_ohcal_event_display->Write();
    for (int i = 0; i < 30; i++) {
        h_topo_event_display[i]->Write(); h_topo_emcal_event_display[i]->Write(); h_topo_ihcal_event_display[i]->Write(); h_topo_ohcal_event_display[i]->Write();
    }

    for (int i = 0; i < 7; i++) {
        h_respmatrix_calib_dijet[i]->Write();
        h_respmatrix_calib_dijet_jesdown[i]->Write();
        h_respmatrix_calib_dijet_jesup[i]->Write();
        h_respmatrix_calib_dijet_jerdown[i]->Write();
        h_respmatrix_calib_dijet_jerup[i]->Write();
        h_respmatrix_calib_dijet_clus_smear[i]->Write();
        h_respmatrix_calib_dijet_ohcal_mc_data_var[i]->Write();
        h_respmatrix_calib_dijet_2sigma_noise[i]->Write();
        h_respmatrix_calib_dijet_4sigma_noise[i]->Write();
        h_respmatrix_calib_dijet_half1[i]->Write();
        h_respmatrix_calib_dijet_half2[i]->Write();
        h_jetpt_respmatrix[i]->Write();
        h_caloet_respmatrix[i]->Write();
    }
    h_truthptbins->Write(); h_truthetbins->Write(); h_calibptbins->Write(); h_calibetbins->Write();
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

void get_leading_subleading_subsubleading_jet(int& leadingjet_index, int& subleadingjet_index, int& subsubleadingjet_index, std::vector<float>* jet_et) {
  leadingjet_index = -1;
  subleadingjet_index = -1;
  subsubleadingjet_index = -1;
  float leadingjet_et = -9999;
  float subleadingjet_et = -9999;
  float subsubleadingjet_et = -9999;
  for (int ij = 0; ij < jet_et->size(); ++ij) {
    float jetet = jet_et->at(ij);
      if (jetet > leadingjet_et) {
        subsubleadingjet_et = subleadingjet_et;
        subsubleadingjet_index = subleadingjet_index;
        subleadingjet_et = leadingjet_et;
        subleadingjet_index = leadingjet_index;
        leadingjet_et = jetet;
        leadingjet_index = ij;
    } else if (jetet > subleadingjet_et) {
        subsubleadingjet_et = subleadingjet_et;
        subsubleadingjet_index = subleadingjet_index;
        subleadingjet_et = jetet;
        subleadingjet_index = ij;
    } else if (jetet > subsubleadingjet_et) {
        subsubleadingjet_et = jetet;
        subsubleadingjet_index = ij;
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
  float dphi = get_dphi(leadingjet_phi, subleadingjet_phi);
  return fabs(dphi) > dijet_min_phi;
}

void get_calibjet(float& calibjet_pt, float& calibjet_eta, float& calibjet_phi, float jet_pt, float jet_eta, float jet_phi, bool reco_cut, TF1* f_corr, float jes_para, float jer_para) {
    calibjet_pt = -9999;
    calibjet_eta = -9999;
    calibjet_phi = -9999;
    if (!reco_cut) return;
    //double calib_pt = f_corr->Eval(jet_pt) * (1 + randGen.Gaus(0.0, jer_para)) * jes_para;
    double calib_pt = jet_pt * (1 + randGen.Gaus(0.0, jer_para)) * jes_para;
    if (calib_pt < calibptbins[0] || calib_pt > calibptbins[calibnpt]) return;
    calibjet_pt = calib_pt;
    calibjet_eta = jet_eta;
    calibjet_phi = jet_phi;
}

void get_truthjet(float& goodtruthjet_pt, float& goodtruthjet_eta, float& goodtruthjet_phi, float jet_pt, float jet_eta, float jet_phi, bool truth_cut) {
    goodtruthjet_pt = -9999;
    goodtruthjet_eta = -9999;
    goodtruthjet_phi = -9999;
    if (!truth_cut) return;
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

void fill_response_matrix(TH2D*& h_truth, TH2D*& h_meas, RooUnfoldResponse*& h_resp, TH2D*& h_fake, TH2D*& h_miss, RooUnfoldResponse*& h_count, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale) {
    //if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet] && truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
    if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1 && truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        if (matched) {
            h_meas->Fill(meas_pt, meas_et, weight_scale);
            h_truth->Fill(truth_pt, truth_et, weight_scale);
            h_resp->Fill(meas_pt, meas_et, truth_pt, truth_et, weight_scale);
            h_count->Fill(meas_pt, meas_et, truth_pt, truth_et);
        } else {
            h_meas->Fill(meas_pt, meas_et, weight_scale);
            h_fake->Fill(meas_pt, meas_et, weight_scale);
            h_resp->Fake(meas_pt, meas_et, weight_scale);
            h_truth->Fill(truth_pt, truth_et, weight_scale);
            h_miss->Fill(truth_pt, truth_et, weight_scale);
            h_resp->Miss(truth_pt, truth_et, weight_scale);
            h_count->Fake(meas_pt, meas_et);
            h_count->Miss(truth_pt, truth_et);
        }
    //} else if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet]) {
    } else if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1) {
        h_meas->Fill(meas_pt, meas_et, weight_scale);
        h_fake->Fill(meas_pt, meas_et, weight_scale);
        h_resp->Fake(meas_pt, meas_et, weight_scale);
        h_count->Fake(meas_pt, meas_et);

    //} else if (truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
    } else if (truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        h_truth->Fill(truth_pt, truth_et, weight_scale);
        h_miss->Fill(truth_pt, truth_et, weight_scale);
        h_resp->Miss(truth_pt, truth_et, weight_scale);
        h_count->Miss(truth_pt, truth_et);
    }
    
}        

void fill_response_matrix_full(TH2D*& h_truth, TH2D*& h_meas, RooUnfoldResponse*& h_resp, TH2D*& h_fake, TH2D*& h_miss, RooUnfoldResponse*& h_count, TH2D*& h_count_fake, TH2D*& h_count_miss, RooUnfoldResponse*& h_jetpt, RooUnfoldResponse*& h_caloet, RooUnfoldResponse*& h_count_jetpt, RooUnfoldResponse*& h_count_caloet, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale) {
    //if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet] && truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
    if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1 && truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        if (matched) {
            // measured, truth and response matrix 
            h_meas->Fill(meas_pt, meas_et, weight_scale);
            h_truth->Fill(truth_pt, truth_et, weight_scale);
            h_resp->Fill(meas_pt, meas_et, truth_pt, truth_et, weight_scale);
                        
            // QA on response matrix 
            h_count->Fill(meas_pt, meas_et, truth_pt, truth_et);
            h_jetpt->Fill(meas_pt, truth_pt, (double)weight_scale);
            h_count_jetpt->Fill(meas_pt, truth_pt);
            h_caloet->Fill(meas_et, truth_et, (double)weight_scale);
            h_count_caloet->Fill(meas_et, truth_et);

            /*
            std::cout << "meas_pt = " << meas_pt << ", truth_pt = " << truth_pt << ", weight_scale = " << weight_scale << std::endl;
            std::cout << "meas_hist range: " << h_jetpt->Hmeasured()->GetXaxis()->GetXmin() << " to " << h_jetpt->Hmeasured()->GetXaxis()->GetXmax() << std::endl;
            std::cout << "truth_hist range: " << h_jetpt->Htruth()->GetXaxis()->GetXmin() << " to " << h_jetpt->Htruth()->GetXaxis()->GetXmax() << std::endl;
            std::cout << "Filled entries in h_jetpt Hresponse: " << h_jetpt->Hresponse()->GetEntries() << ", integral: " << h_jetpt->Hresponse()->Integral() << std::endl;
            
            int xbin = h_jetpt->Hresponse()->GetXaxis()->FindBin(meas_pt);
            int ybin = h_jetpt->Hresponse()->GetYaxis()->FindBin(truth_pt);
            double bincontent = h_jetpt->Hresponse()->GetBinContent(xbin, ybin);
            std::cout << "h_jetpt bin content at (meas_pt, truth_pt): " << bincontent << std::endl;

            std::cout << "meas_et = " << meas_et << ", truth_et = " << truth_et << ", weight_scale = " << weight_scale << std::endl;
            std::cout << "meas_hist range: " << h_caloet->Hmeasured()->GetXaxis()->GetXmin() << " to " << h_caloet->Hmeasured()->GetXaxis()->GetXmax() << std::endl;
            std::cout << "truth_hist range: " << h_caloet->Htruth()->GetXaxis()->GetXmin() << " to " << h_caloet->Htruth()->GetXaxis()->GetXmax() << std::endl;
            std::cout << "Filled entries in h_caloet Hresponse: " << h_caloet->Hresponse()->GetEntries() << ", integral: " << h_caloet->Hresponse()->Integral() << std::endl;

            xbin = h_caloet->Hresponse()->GetXaxis()->FindBin(meas_et);
            ybin = h_caloet->Hresponse()->GetYaxis()->FindBin(truth_et);
            bincontent = h_caloet->Hresponse()->GetBinContent(xbin, ybin);
            std::cout << "h_caloet bin content at (meas_et, truth_et): " << bincontent << std::endl;

            int xbin_meas = h_resp->Hmeasured()->GetXaxis()->FindBin(meas_pt);
            int ybin_meas = h_resp->Hmeasured()->GetYaxis()->FindBin(meas_et);
            int xbin_truth = h_resp->Htruth()->GetXaxis()->FindBin(truth_pt);
            int ybin_truth = h_resp->Htruth()->GetYaxis()->FindBin(truth_et);
            int nbins_meas_pt = h_resp->Hmeasured()->GetXaxis()->GetNbins();
            int nbins_truth_pt = h_resp->Htruth()->GetXaxis()->GetNbins();
            int meas_flat_bin = (ybin_meas - 1) * nbins_meas_pt + xbin_meas;  // 1-based ROOT indexing
            int truth_flat_bin = (ybin_truth - 1) * nbins_truth_pt + xbin_truth;
            double bin_content = h_resp->Hresponse()->GetBinContent(meas_flat_bin, truth_flat_bin);

            std::cout << "meas_pt = " << meas_pt << ", meas_et = " << meas_et  << ", truth_pt = " << truth_pt << ", truth_et = " << truth_et  << ", weight_scale = " << weight_scale << std::endl;
            std::cout << "measured pt axis range: " << h_resp->Hmeasured()->GetXaxis()->GetXmin() << " to " << h_resp->Hmeasured()->GetXaxis()->GetXmax() << std::endl;
            std::cout << "measured et axis range: " << h_resp->Hmeasured()->GetYaxis()->GetXmin() << " to " << h_resp->Hmeasured()->GetYaxis()->GetXmax() << std::endl;
            std::cout << "truth pt axis range: " << h_resp->Htruth()->GetXaxis()->GetXmin() << " to " << h_resp->Htruth()->GetXaxis()->GetXmax() << std::endl;
            std::cout << "truth et axis range: " << h_resp->Htruth()->GetYaxis()->GetXmin() << " to " << h_resp->Htruth()->GetYaxis()->GetXmax() << std::endl;
            std::cout << "Filled entries in h_resp Hresponse: " << h_resp->Hresponse()->GetEntries() << ", integral: " << h_resp->Hresponse()->Integral() << std::endl;
            std::cout << "h_resp bin content at ((meas_pt, meas_et), (truth_pt, truth_et)): " << bin_content << std::endl;
            */
        } else {
            h_meas->Fill(meas_pt, meas_et, weight_scale);
            h_fake->Fill(meas_pt, meas_et, weight_scale);
            h_resp->Fake(meas_pt, meas_et, weight_scale);
            h_truth->Fill(truth_pt, truth_et, weight_scale);
            h_miss->Fill(truth_pt, truth_et, weight_scale);
            h_resp->Miss(truth_pt, truth_et, weight_scale);

            // QA on response matrix 
            h_count->Fake(meas_pt, meas_et);
            h_count->Miss(truth_pt, truth_et);
            h_count_fake->Fill(meas_pt, meas_et);
            h_count_miss->Fill(truth_pt, truth_et);
            h_jetpt->Fake(meas_pt, weight_scale);
            h_jetpt->Miss(truth_pt, weight_scale);
            h_count_jetpt->Fake(meas_pt);
            h_count_jetpt->Miss(truth_pt);
            h_caloet->Fake(meas_et, weight_scale);
            h_caloet->Miss(truth_et, weight_scale);
            h_count_caloet->Fake(meas_et);
            h_count_caloet->Miss(truth_et);

        }
    //} else if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet]) {
    } else if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1) {
        h_meas->Fill(meas_pt, meas_et, weight_scale);
        h_fake->Fill(meas_pt, meas_et, weight_scale);
        h_resp->Fake(meas_pt, meas_et, weight_scale);

        // QA on response matrix 
        h_count->Fake(meas_pt, meas_et);
        h_count_fake->Fill(meas_pt, meas_et);
        h_jetpt->Fake(meas_pt, weight_scale);
        h_count_jetpt->Fake(meas_pt);
        h_caloet->Fake(meas_et, weight_scale);
        h_count_caloet->Fake(meas_et);
    //} else if (truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[calibnet]) {
    } else if (truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        h_truth->Fill(truth_pt, truth_et, weight_scale);
        h_miss->Fill(truth_pt, truth_et, weight_scale);
        h_resp->Miss(truth_pt, truth_et, weight_scale);

        // QA on response matrix 
        h_count->Miss(truth_pt, truth_et);
        h_count_miss->Fill(truth_pt, truth_et);
        h_jetpt->Miss(truth_pt, weight_scale);
        h_count_jetpt->Miss(truth_pt);
        h_caloet->Miss(truth_et, weight_scale);
        h_count_caloet->Miss(truth_et);
    }

    /*
        // 1D unfolding syst
    if (meas_pt >= calibptbins[0] && truth_pt >= truthptbins[0]) {
        if (matched) {
            h_jetpt->Fill(meas_pt, truth_pt, weight_scale);
            h_count_jetpt->Fill(meas_pt, truth_pt);
        } else {
            h_jetpt->Fake(meas_pt, weight_scale);
            h_jetpt->Miss(truth_pt, weight_scale);
            h_count_jetpt->Fake(meas_pt);
            h_count_jetpt->Miss(truth_pt);
        }
    } else if (meas_pt >= calibptbins[0]) {
        h_jetpt->Fake(meas_pt, weight_scale);
        h_count_jetpt->Fake(meas_pt);
    } else if (truth_pt >= truthptbins[0]) {
        h_jetpt->Miss(truth_pt, weight_scale);
        h_count_jetpt->Miss(truth_pt);
    }

    if (meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
        h_caloet->Fill(meas_et, truth_et, weight_scale);
        h_count_caloet->Fill(meas_et, truth_et);
    } else if (meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet]) {
        h_caloet->Fake(meas_et, weight_scale);
        h_count_caloet->Fake(meas_et);
    } else if (truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
        h_caloet->Miss(truth_et, weight_scale);
        h_count_caloet->Miss(truth_et);
    }
    */
}

void fill_trim_response_matrix(RooUnfoldResponse*& h_resp, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale, const TMatrixD& counts_matrix, TH2D* h_counts_measured, TH2D* h_counts_truth, float trim_value) { 
    //if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet] && truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
    if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1 && truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        // only fill trimmed bins 
        int meas_bin = h_counts_measured->GetNbinsX() * (h_counts_measured->GetYaxis()->FindBin(meas_et) - 1) + (h_counts_measured->GetXaxis()->FindBin(meas_pt) - 1);
        int truth_bin = h_counts_truth->GetNbinsX() * (h_counts_truth->GetYaxis()->FindBin(truth_et) - 1) + (h_counts_truth->GetXaxis()->FindBin(truth_pt) - 1);
        //std::cout << "meas (pt,et) " << meas_pt << " " << meas_et << " meas_bin " << meas_bin << " truth (pt,et) " << truth_pt << " " << truth_et << " truth_bin " << truth_bin << std::endl;
        if (matched && counts_matrix(meas_bin, truth_bin) >= trim_value) {
             h_resp->Fill(meas_pt, meas_et, truth_pt, truth_et, weight_scale);
        } else {
            h_resp->Fake(meas_pt, meas_et, weight_scale);
            h_resp->Miss(truth_pt, truth_et, weight_scale);
        }
    //} else if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet]) {
    } else if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1) {
        h_resp->Fake(meas_pt, meas_et, weight_scale);

    //} else if (truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
    } else if (truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        h_resp->Miss(truth_pt, truth_et, weight_scale);
    }
}

void fill_reweighted_trim_response_matrix(RooUnfoldResponse*& h_resp, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale, const TMatrixD& counts_matrix, TH2D* h_counts_measured, TH2D* h_counts_truth, float trim_value, TH2D* h_prior_weights) { 
    double prior_weight = h_prior_weights->GetBinContent(h_prior_weights->FindBin(truth_pt,truth_et));
    //if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet] && truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
    if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1 && truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        // only fill trimmed bins 
        int meas_bin = h_counts_measured->GetNbinsX() * (h_counts_measured->GetYaxis()->FindBin(meas_et) - 1) + (h_counts_measured->GetXaxis()->FindBin(meas_pt) - 1);
        int truth_bin = h_counts_truth->GetNbinsX() * (h_counts_truth->GetYaxis()->FindBin(truth_et) - 1) + (h_counts_truth->GetXaxis()->FindBin(truth_pt) - 1);
        if (matched && counts_matrix(meas_bin, truth_bin) >= trim_value) {
             h_resp->Fill(meas_pt, meas_et, truth_pt, truth_et, weight_scale*prior_weight);
        } else {
            h_resp->Fake(meas_pt, meas_et, weight_scale);
            h_resp->Miss(truth_pt, truth_et, weight_scale*prior_weight);
        }
    //} else if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet]) {
    } else if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1) {
        h_resp->Fake(meas_pt, meas_et, weight_scale);

    //} else if (truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
    } else if (truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        h_resp->Miss(truth_pt, truth_et, weight_scale*prior_weight);
    }
}

void fill_trim_1D_response_matrices(RooUnfoldResponse*& h_jetpt_resp, RooUnfoldResponse*& h_caloet_resp, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale, const TMatrixD& counts_matrix, TH2D* h_counts_measured, TH2D* h_counts_truth, float trim_value) { 
    // JET PT AND CALO ET 1D UNFOLDING
    //if (meas_pt >= calibptbins[0] && truth_pt >= truthptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
    if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1 && truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        // only fill trimmed bins 
        int meas_bin = h_counts_measured->GetNbinsX() * (h_counts_measured->GetYaxis()->FindBin(meas_et) - 1) + (h_counts_measured->GetXaxis()->FindBin(meas_pt) - 1);
        int truth_bin = h_counts_truth->GetNbinsX() * (h_counts_truth->GetYaxis()->FindBin(truth_et) - 1) + (h_counts_truth->GetXaxis()->FindBin(truth_pt) - 1);
        if (matched && counts_matrix(meas_bin, truth_bin) >= trim_value) {
            h_jetpt_resp->Fill(meas_pt, truth_pt, weight_scale);
            h_caloet_resp->Fill(meas_et, truth_et, weight_scale);
        } else {
            h_jetpt_resp->Fake(meas_pt, weight_scale);
            h_jetpt_resp->Miss(truth_pt, weight_scale);
            h_caloet_resp->Fake(meas_et, weight_scale);
            h_caloet_resp->Miss(truth_et, weight_scale);
        }
    //} else if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet]) {
    } else if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1) {
        h_jetpt_resp->Fake(meas_pt, weight_scale);
        h_caloet_resp->Fake(meas_et, weight_scale);

    //} else if (truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
    } else if (truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        h_jetpt_resp->Miss(truth_pt, weight_scale);
        h_caloet_resp->Miss(truth_et, weight_scale);
    }
}

void fill_reweighted_trim_1D_response_matrices(RooUnfoldResponse*& h_jetpt_resp, RooUnfoldResponse*& h_caloet_resp, float meas_pt, float matched, float truth_pt, float meas_et, float truth_et, float weight_scale, const TMatrixD& counts_matrix, TH2D* h_counts_measured, TH2D* h_counts_truth, float trim_value, TH1D* h_jetpt_prior_weights, TH1D* h_caloet_prior_weights) { 
    // 1D UNFOLDING
    double jetpt_prior_weight = h_jetpt_prior_weights->GetBinContent(h_jetpt_prior_weights->FindBin(truth_pt));
    double caloet_prior_weight = h_caloet_prior_weights->GetBinContent(h_caloet_prior_weights->FindBin(truth_et));
    //if (meas_pt >= calibptbins[0] && truth_pt >= truthptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
    if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1 && truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        // only fill trimmed bins 
        int meas_bin = h_counts_measured->GetNbinsX() * (h_counts_measured->GetYaxis()->FindBin(meas_et) - 1) + (h_counts_measured->GetXaxis()->FindBin(meas_pt) - 1);
        int truth_bin = h_counts_truth->GetNbinsX() * (h_counts_truth->GetYaxis()->FindBin(truth_et) - 1) + (h_counts_truth->GetXaxis()->FindBin(truth_pt) - 1);
        if (matched && counts_matrix(meas_bin, truth_bin) >= trim_value) {
             h_jetpt_resp->Fill(meas_pt, truth_pt, weight_scale*jetpt_prior_weight);
             h_caloet_resp->Fill(meas_et, truth_et, weight_scale*caloet_prior_weight);
        } else {
            h_jetpt_resp->Fake(meas_pt, weight_scale);
            h_jetpt_resp->Miss(truth_pt, weight_scale*jetpt_prior_weight);
            h_caloet_resp->Fake(meas_et, weight_scale);
            h_caloet_resp->Miss(truth_et, weight_scale*caloet_prior_weight);
        }
    //} else if (meas_pt >= calibptbins[0] && meas_et >= calibetbins[0] && meas_et <= calibetbins[calibnet]) {
    } else if (meas_pt >= 0 && meas_et >= 0 && meas_et <= 1) {
        h_jetpt_resp->Fake(meas_pt, weight_scale);
        h_caloet_resp->Fake(meas_et, weight_scale);

    //} else if (truth_pt >= truthptbins[0] && truth_et >= truthetbins[0] && truth_et <= truthetbins[truthnet]) {
    } else if (truth_pt >= 0 && truth_et >= 0 && truth_et <= 1) {
        h_jetpt_resp->Miss(truth_pt, weight_scale*jetpt_prior_weight);
        h_caloet_resp->Miss(truth_et, weight_scale*caloet_prior_weight);
    }
}
