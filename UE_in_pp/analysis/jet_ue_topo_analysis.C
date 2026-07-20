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
#include "unfolding/unfold_Def.h"
#include <TProfile.h>

bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled);
void get_leading_subleading_jet(int& leadingjet_index, int& subleadingjet_index, std::vector<float>* jet_et);
void get_leading_jet(int& leadingjet_index, std::vector<float>* jet_et);
bool match_leading_subleading_jet(float leadingjet_phi, float subleadingjet_phi);

std::vector<int> run_numbers;
std::vector<double> collision_rates;
std::vector<double> pileup_rates;
std::vector<double> pileup_bins = {0.0,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2};
std::vector<double> pileup_correction = {0.02,0.018,0.021}; // udpated correction with MBD cut, old correction:{0.0353,0.0341,0.0379}, correction determined from both, dijet and efrac cuts

TRandom3 randGen(1234);

// get pile up rates for correction
void get_pileup_rates(const char* infilename) {
  std::ifstream infile(infilename);
  if (!infile) {
      std::cerr << "Error: Unable to open input file!" << std::endl;
      return;
  }
  std::string line;
  while (std::getline(infile, line)) {
      std::istringstream iss(line);
      int run_number;
      double collision_rate;  
      if (!(iss >> run_number >> collision_rate)) {
          std::cerr << "Error: Malformed line in input file!: " << line << std::endl;
          continue;
      }
      run_numbers.push_back(run_number);
      collision_rates.push_back(collision_rate);

      // Calculate lambda
      double denominator = 78000.0 * 111.0;
      double lambda = collision_rate / denominator;
      double exp_neg_lambda = exp(-lambda);
      
      // Calculate pk
      double pk = 1 - lambda * exp_neg_lambda - exp_neg_lambda;
      
      // Calculate pileup rate
      double pileup_rate = pk / (pk + lambda * exp_neg_lambda);
      pileup_rates.push_back(pileup_rate);
  }
  infile.close();
}

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

////////////////////////////////////////// Main Function //////////////////////////////////////////
void jet_ue_topo_analysis(int sim = 0, int runnumber = 51274, std::string bkg_cut = "dijet", bool clusters = true, std::string runtype = "jet12")  {
  /////////////// General Set up ///////////////
  string outfilename;
  if (sim) {
    outfilename = "jet_ue_topo_analysis/output_sim_" + runtype + "_pu_correct_" + bkg_cut + "_bkg_cut_" + to_string(runnumber) + ".root";
  } else {
    outfilename = "jet_ue_topo_analysis/output_pu_correct_ana509_v2_" + bkg_cut + "_bkg_cut_"+ to_string(runnumber) + ".root";
  }
  TFile *f_out = new TFile(outfilename.c_str(), "RECREATE");

  /////////////// Read Files ///////////////
  double weight_scale = 1.0, truthjet_pt_min = 0, truthjet_pt_max = 3000;
  if (sim == 0) {
    get_pileup_rates("/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/macro/mbdliverate_grl100625.txt");
  } else {
    if (runtype == "jet12") {
        weight_scale = Jet12GeV_scale;
        truthjet_pt_min = 14;
        truthjet_pt_max = 21;
    } else if (runtype == "jet20") {
        weight_scale = Jet20GeV_scale;
        truthjet_pt_min = 21;
        truthjet_pt_max = 32;
    } else if (runtype == "jet30") {
        weight_scale = Jet30GeV_scale;
        truthjet_pt_min = 32;
        truthjet_pt_max = 42;
    } else if (runtype == "jet40") {
        weight_scale = Jet40GeV_scale;
        truthjet_pt_min = 42;
        truthjet_pt_max = 52;
    } else if (runtype == "jet50") {
        weight_scale = Jet50GeV_scale;
        truthjet_pt_min = 52;
        truthjet_pt_max = 62;
    } else if (runtype == "jet60") {
        weight_scale = Jet60GeV_scale;
        truthjet_pt_min = 62;
        truthjet_pt_max = 1000;
    }
  }
  //get_pileup_rates("/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/analysis/pileup/0mrad_collision_rates.txt");
  //get_pileup_rates("/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/analysis/pileup/updated_1.5mrad_collision_rates.txt");
  const char* baseDirJet = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput";
  TChain chain("T");
  if (sim) {
    chain.Add(Form("%s/sim_run28_%s_3sigma_output_0.root", baseDirJet, runtype.c_str()));
  } else {
    chain.Add(Form("%s/output_ana509_v2_%d_*.root", baseDirJet, runnumber));
  }
  /*
  if (runnumber < 51274) {
    chain.Add(Form("%s/output_0mrad_ana468_%d_*.root", baseDirJet, runnumber));
  } else {
    chain.Add(Form("%s/output1.5mrad_ana468_%d_*.root", baseDirJet, runnumber));
  }
  */
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
  float mbd_t0; chain.SetBranchStatus("mbd_t0", 1); chain.SetBranchAddress("mbd_t0", &mbd_t0);
  std::vector<float>* truthjet_pt = nullptr; chain.SetBranchStatus("truthPt", 1); chain.SetBranchAddress("truthPt", &truthjet_pt);
  std::vector<float>* truthjet_eta = nullptr; chain.SetBranchStatus("truthEta", 1); chain.SetBranchAddress("truthEta", &truthjet_eta);
  std::vector<float>* truthjet_phi = nullptr; chain.SetBranchStatus("truthPhi", 1); chain.SetBranchAddress("truthPhi", &truthjet_phi);

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
  } else {
      chain.SetBranchStatus("clsmult", 1); chain.SetBranchAddress("clsmult",&clsmult); // edited to use lower topoclsuter threshold
      chain.SetBranchStatus("cluster_e", 1); chain.SetBranchAddress("cluster_e",cluster_e);
      chain.SetBranchStatus("cluster_eta", 1); chain.SetBranchAddress("cluster_eta",cluster_eta);
      chain.SetBranchStatus("cluster_phi", 1); chain.SetBranchAddress("cluster_phi",cluster_phi);
  } 

  /////////////// JES func ///////////////
  TFile *corrFile = new TFile("unfolding/JES_Calib_Default.root", "READ");
  if (!corrFile) {
    std::cout << "Error: cannot open JES_Calib_Default.root" << std::endl;
    return;
  }
  TF1 *f_corr = (TF1*)corrFile->Get("JES_Calib_Default_Func");
  if (!f_corr) {
    std::cout << "Error: cannot open f_corr" << std::endl;
    return;
  }

  /////////////// Trigger Efficiency ///////////////
  TFile *f_trigger = new TFile("unfolding/output_jettrigeff.root", "READ");
  if (!f_trigger) {
    std::cout << "Error: cannot open output_jettrigeff.root" << std::endl;
    return;
  }
  TF1* f_turnon = (TF1*)f_trigger->Get("f_fit");
  TF1* f_turnon_down = (TF1*)f_trigger->Get("f_fit_shiftdown"); // for trigger efficiency uncertainty
  TF1* f_turnon_up = (TF1*)f_trigger->Get("f_fit_shiftup"); // for trigger efficiency uncertainty
  if (!f_turnon || !f_turnon_down || !f_turnon_up) {
    std::cout << "Error: cannot open trigger efficiency functions" << std::endl;
    return;
  }

  TFile* f_zvertex;
  TH1D* h_vertex_weight;
  string zvertex_filename = "output_vertex_reweight_run28_dijet_bkg_cut_iter_1.root";
  f_zvertex = new TFile(zvertex_filename.c_str(), "READ");
  TH1D* htemp = (TH1D*) f_zvertex->Get("vertex_weight");
  if (!htemp)
  {
    std::cerr << "ERROR: Could not find histogram 'vertex_weight'!" << std::endl;
    exit(1);
  }
  h_vertex_weight = (TH1D*) htemp->Clone("vertex_weight_clone");
  h_vertex_weight->SetDirectory(0);  // detach from file
  f_zvertex->Close();

  ////////////// Jet Background Cut Efficiency ////////////
  TFile *f_bkgcut = new TFile("jet_background_eff/tefficiency_fits.root");
  TF1* f_bkgcut_eff = (TF1*)f_bkgcut->Get("fit_nominal");
  TSpline3* f_bkgcut_eff_down = (TSpline3*)f_bkgcut->Get("fit_minus1");
  TSpline3* f_bkgcut_eff_up = (TSpline3*)f_bkgcut->Get("fit_plus1");

  ////////////// Jet Background Timing Cut Efficiency ////////////
  double f_timingcut_eff = 0.95;
  double f_timingcut_eff_up = 0.99;
  double timingcut_scale = 1.0 / f_timingcut_eff;
  double timingcut_scale_up = 1.0 / f_timingcut_eff_up;
  double lead_time_cut[2] = {-8.0,4.0}; // leading time cut
  double lead_time_var[2] = {-9.0,5.0}; // leading time variation
  double deltat_mbd_cut[2] = {-5.0,1.0}; // delta time cut
  double deltat_dijet_cut[2] = {-3.0,3.0}; // delta time cut
  double deltat_mbd_var[2] = {-6.0,2.0}; // delta time variation
  double deltat_dijet_var[2] = {-4.0,4.0}; // delta timme variation

  float ptbins[] = {5,6,7,8,9,10,12,14,16,18,20,22,24,27,30,35,40,50,60,80};
  int nptbins = sizeof(ptbins) / sizeof(ptbins[0]) - 1;
  float etmin = -20;
  float etmax = 80;
  int netbins = int((10*(etmax - etmin)));
  float etbins[netbins+1];
  for (int i = 0; i <= netbins; i++) {
    etbins[i] = etmin + i*(etmax - etmin) / netbins;
  }

  /////////////// Histograms ///////////////
  TH1D *h_zvertex = new TH1D("h_zvertex", ";Z-vertex [cm]", 120, -60, 60);
  TH1D* h_deltaphi_record = new TH1D("h_deltaphi_record","",125,-2*M_PI,2*M_PI);
  TH1D* h_xj_record = new TH1D("h_xj_record","",20,0,1);
  TH1D* h_lead_spectra_record = new TH1D("h_lead_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
  TH1D* h_sub_spectra_record = new TH1D("h_sub_spectra_record",";p_{T} [GeV]", 1000, 0, 100);
  TH1D *h_et_transverse_record = new TH1D("h_et_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
  TH2D* h_ue_pt_transverse_record = new TH2D("h_ue_pt_transverse_record","", 100, 0, 100, 700, -20, 50);
  TH2D *h_calibjet_pt_dijet_eff = new TH2D("h_calibjet_pt_dijet_eff", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, 0, 1, calibnet, 0, 1); // with trigger efficiency + beam background efficiency + timing cut efficiency correction applied
  TH2D *h_calibjet_pt_dijet_effdown = new TH2D("h_calibjet_pt_dijet_effdown", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, 0, 1, calibnet, 0, 1); // for trigger efficiency uncertainty
  TH2D *h_calibjet_pt_dijet_effup = new TH2D("h_calibjet_pt_dijet_effup", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, 0, 1, calibnet, 0, 1); // for trigger efficiency uncertainty
  TH2D *h_calibjet_pt_dijet_pu_correct_et = new TH2D("h_calibjet_pt_dijet_pu_correct_et", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, 0, 1, calibnet, 0, 1); // with pileup correction
  TH2D *h_calibjet_pt_pu_up_et = new TH2D("h_calibjet_pt_pu_up_et", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
  TH2D *h_calibjet_pt_pu_down_et = new TH2D("h_calibjet_pt_pu_down_et", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, 0, 1, calibnet, 0, 1);
  TH2D *h_calibjet_pt_timingeffup = new TH2D("h_calibjet_pt_timingeffup", ";p_{T} [GeV];#SigmaE_{T} [GeV]", calibnpt, 0, 1, calibnet, 0, 1); // for timing cut efficiency uncertainty
  TH1D* h_jetpt = new TH1D("h_jetpt","", calibnpt, 0, 1); TH1D* h_caloet = new TH1D("h_caloet","", calibnet, 0, 1);

  TH1D* h_calib_jet_pt_tight = new TH1D("h_calib_jet_pt_tight","",50,15,65);
  TH1D* h_calib_calo_et_tight = new TH1D("h_calib_calo_et_tight","",110,-2,20);
  TH1D* h_calib_jet_pt_uni_tight = new TH1D("h_calib_jet_pt_uni_tight","",50,0,1);
  TH1D* h_calib_calo_et_uni_tight = new TH1D("h_calib_calo_et_uni_tight","",110,0,1);
  
  // QA histograms //
  TProfile* h_unscale_average_et = new TProfile("h_unscale_average_et","", 7000, 47000, 54000);
  TProfile* h_scale_average_et = new TProfile("h_scale_average_et","", 7000, 47000, 54000);
  TProfile* h_scale_njet = new TProfile("h_scale_njet","", 7000, 47000, 54000);
  TProfile* h_unscale_ejet = new TProfile("h_unscale_ejet","",7000,47000,54000);

  TH1D *h_pu_correct_et_transverse_record = new TH1D("h_pu_correct_et_transverse_record", ";#SigmaE_{T} [GeV]", 7000, -20, 50);
  TH2D* h_pu_correct_ue_pt_transverse_record = new TH2D("h_pu_correct_ue_pt_transverse_record","", 100, 0, 100, 700, -20, 50);
  TH2D* h_et_transverse_pu_correct_et_transverse_record = new TH2D("h_et_transverse_pu_correct_et_transverse_record","", 700, -20, 50, 700, -20, 50);
  TH3D* h_et_transverse_pu_correction_pu_rate_record = new TH3D("h_et_transverse_pu_correction_pu_rate_record","", 700, -20, 50, 100, 0, 1, 100, 0, 0.2);

  TH1F* h_ue_towards = new TH1F("h_ue_towards","",netbins, etbins);
  TH1F* h_ue_transverse = new TH1F("h_ue_transverse","",netbins, etbins);
  TH1F* h_ue_away = new TH1F("h_ue_away","",netbins, etbins);

  TH1F* h_et_towards = new TH1F("h_et_towards","",netbins, etbins);
  TH1F* h_et_transverse = new TH1F("h_et_transverse","",netbins, etbins);
  TH1F* h_et_away = new TH1F("h_et_away","",netbins, etbins);

  // create topocluster histograms 
  int topo_thresholds[] = {-9999,0,100,200,300,500,1000,2000};

  TH1F* h_ntopo_towards[8];
  TH1F* h_ntopo_transverse[8];
  TH1F* h_ntopo_away[8];

  TH1F* h_topo_towards[8];
  TH1F* h_topo_transverse[8];
  TH1F* h_topo_away[8];

  TH2F* h_2D_topo_towards[8];
  TH2F* h_2D_topo_transverse[8];
  TH2F* h_2D_topo_away[8];

  TH1F* h_sume_topo_towards[8];
  TH1F* h_sume_topo_transverse[8];
  TH1F* h_sume_topo_away[8];

  for (int i = 0; i < 8; i++) {
    h_ntopo_towards[i] = new TH1F(Form("h_ntopo%d_towards",topo_thresholds[i]),"",200,0,200);
    h_ntopo_transverse[i] = new TH1F(Form("h_ntopo%d_transverse",topo_thresholds[i]),"",200,0,200);
    h_ntopo_away[i] = new TH1F(Form("h_ntopo%d_away",topo_thresholds[i]),"",200,0,200);

    h_topo_towards[i] = new TH1F(Form("h_topo%d_towards",topo_thresholds[i]),"",netbins, etbins);
    h_topo_transverse[i] = new TH1F(Form("h_topo%d_transverse",topo_thresholds[i]),"",netbins, etbins);
    h_topo_away[i] = new TH1F(Form("h_topo%d_away",topo_thresholds[i]),"",netbins, etbins);

    h_2D_topo_towards[i] = new TH2F(Form("h_2D_topo%d_towards",topo_thresholds[i]),"",24,-1.1,1.1,32,0,M_PI);
    h_2D_topo_transverse[i] = new TH2F(Form("h_2D_topo%d_transverse",topo_thresholds[i]),"",24,-1.1,1.1,32,0,M_PI);
    h_2D_topo_away[i] = new TH2F(Form("h_2D_topo%d_away",topo_thresholds[i]),"",24,-1.1,1.1,32,0,M_PI);

    h_sume_topo_towards[i] = new TH1F(Form("h_sume_topo%d_towards",topo_thresholds[i]),"",netbins, etbins);
    h_sume_topo_transverse[i] = new TH1F(Form("h_sume_topo%d_transverse",topo_thresholds[i]),"",netbins, etbins);
    h_sume_topo_away[i] = new TH1F(Form("h_sume_topo%d_away",topo_thresholds[i]),"",netbins, etbins);
  }

  TProfile* h_ue_pt_towards = new TProfile("h_ue_pt_towards","",nptbins, ptbins);
  TProfile* h_ue_pt_transverse = new TProfile("h_ue_pt_transverse","",nptbins, ptbins);
  TProfile* h_ue_pt_away = new TProfile("h_ue_pt_away","",nptbins, ptbins);

  TH1D* h_calibptbins = new TH1D("h_calibptbins","",calibnpt, 0, 1);
  TH1D* h_calibetbins = new TH1D("h_calibetbins","",calibnet, 0, 1);

  /////////////// Event Loop ///////////////
  std::cout << "Data analysis started." << std::endl;
  int n_events = chain.GetEntries();
  std::cout << "Total number of events: " << n_events << std::endl;

  // Event variables setup.
  bool reco_bkg_cut = false; 
  float pileup = 1.0;

  for (int ie = 0; ie < n_events; ++ie) { // event loop start
  //for (int ie = 0; ie < 20; ++ie) {
    // Load event.
    if (ie % 1000 == 0) {
      std::cout << "Processing event " << ie << "..." << std::endl;
    }
    int njets = 0; float ejets = 0.0;
    chain.GetEntry(ie);

    // Trigger and Z-vertex cut.
    if (!sim && !check_bad_trigger(gl1_trigger_vector_scaled)) continue;
    if (isnan(zvertex)) { continue; }
    if (fabs(zvertex) > 60) continue;  // edited from 60 cm to 30 cm
    double vertex_weight = 1.0;
    int ind_truth_lead = -1;
    if (sim) {
      int vertex_bin = h_vertex_weight->FindBin(zvertex);
      vertex_weight = h_vertex_weight->GetBinContent(vertex_bin);
      h_zvertex->Fill(zvertex, vertex_weight * weight_scale);
      if (truthjet_pt->size() < 1) {
        continue;
      }
      get_leading_jet(ind_truth_lead, truthjet_pt);
      if (truthjet_pt->at(ind_truth_lead) < truthjet_pt_min || truthjet_pt->at(ind_truth_lead) > truthjet_pt_max) {
        continue;
      }
    }
    else {
      h_zvertex->Fill(zvertex, weight_scale);
    }

    //////////////////////////// SETUP JET VARIABLES FOR UNFOLDING ////////////////////////////

    // indices to find leading and subleading jets 
    int ind_lead = -1; int ind_sub = -1;
    float lead_e = 0;
    reco_bkg_cut = false;

    //std::cout << "reco jets " << unsubjet_pt->size() << std::endl;

    // require 1 jet for inclusive, 2 jets for dijet
    int nJetReq = 1;
    if (bkg_cut == "dijet") { nJetReq = 2; }
    if (unsubjet_pt->size() < nJetReq) { continue; }

    // check number of jets above 5 GeV
    int Njet = 0;
    for (size_t i = 0; i < unsubjet_pt->size(); i++) {
      if (unsubjet_pt->at(i) >= 5.0) {
        Njet++;
      }
    }

    // get reco jets with eta in calorimeter acceptance
    std::vector<float> recoe_new, recopt_new, recoeta_new, recophi_new, recoemcal_new, recoihcal_new, recoohcal_new, recotime_new;
    for (size_t i = 0; i < unsubjet_eta->size(); ++i) {
      if (!check_bad_jet_eta(unsubjet_eta->at(i), zvertex, jet_radius) && fabs(unsubjet_eta->at(i)) < 0.7 && unsubjet_e->at(i) > 0.0) { 
        recoe_new.push_back(unsubjet_e->at(i));
        recopt_new.push_back(unsubjet_pt->at(i));
        recoeta_new.push_back(unsubjet_eta->at(i));
        recophi_new.push_back(unsubjet_phi->at(i));
        recoemcal_new.push_back(unsubjet_emcal_calo_e->at(i));
        recoihcal_new.push_back(unsubjet_ihcal_calo_e->at(i));
        recoohcal_new.push_back(unsubjet_ohcal_calo_e->at(i));
        recotime_new.push_back(unsubjet_time->at(i));
      }
    }

    *unsubjet_e = std::move(recoe_new);
    *unsubjet_pt = std::move(recopt_new);
    *unsubjet_eta = std::move(recoeta_new);
    *unsubjet_phi = std::move(recophi_new);
    *unsubjet_emcal_calo_e = std::move(recoemcal_new);
    *unsubjet_ihcal_calo_e = std::move(recoihcal_new);
    *unsubjet_ohcal_calo_e = std::move(recoohcal_new);
    *unsubjet_time = std::move(recotime_new);

    //std::cout << "reco jets " << unsubjet_pt->size() << std::endl;

    if (unsubjet_pt->size() < nJetReq) { continue; }

    // find reco leading and subleading jets if nreco jets >= nJetReq
    TVector3 lead;
    if (bkg_cut == "dijet") {
      // jet leading and subleading jets in event and apply dijet cut
      get_leading_subleading_jet(ind_lead, ind_sub, unsubjet_pt); 
      //std::cout << " sub/lead: " << unsubjet_e->at(ind_sub)/unsubjet_e->at(ind_lead) << " dPhi: " << get_dphi(unsubjet_phi->at(ind_lead), unsubjet_phi->at(ind_sub)) << std::endl; 
      if (unsubjet_e->at(ind_sub)/unsubjet_e->at(ind_lead) > 0.3 && match_leading_subleading_jet(unsubjet_phi->at(ind_lead), unsubjet_phi->at(ind_sub))) {
        reco_bkg_cut = true;
      } else {
        reco_bkg_cut = false;
      }
    } else {
      // apply efraction cut to leading jet
      get_leading_jet(ind_lead, unsubjet_pt);
      float jete = unsubjet_e->at(ind_lead);
      //std::cout << " EMCal frac: " << unsubjet_emcal_calo_e->at(ind_lead)/jete << " IHCal frac: " << unsubjet_ihcal_calo_e->at(ind_lead)/jete << " OHCal frac: " << unsubjet_ohcal_calo_e->at(ind_lead)/jete << std::endl;
      if (unsubjet_emcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_emcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ihcal_calo_e->at(ind_lead)/jete < 0.9 && unsubjet_ohcal_calo_e->at(ind_lead)/jete > 0.1 && unsubjet_ohcal_calo_e->at(ind_lead)/jete < 0.9) {
        reco_bkg_cut = true; 
      } else {
        reco_bkg_cut = false;
      }
    }
    if (!reco_bkg_cut) { continue; }

    // Njet and timing cuts
    if (Njet >= 9) {
      reco_bkg_cut = false; 
    }
    float lead_time = unsubjet_time->at(ind_lead)*17.6;
    if (lead_time < lead_time_cut[0] || lead_time > lead_time_cut[1]) {
      reco_bkg_cut = false;
    }
    float lead_deltat;
    if (bkg_cut == "dijet") {
      lead_deltat = (unsubjet_time->at(ind_lead) - unsubjet_time->at(ind_sub))*17.6;
      if (lead_deltat < deltat_dijet_cut[0] || lead_deltat > deltat_dijet_cut[1]) {
        reco_bkg_cut = false;
      }
    } 
    //else {
    lead_deltat = unsubjet_time->at(ind_lead)*17.6 - mbd_t0;
    if (lead_deltat < deltat_mbd_cut[0] || lead_deltat > deltat_mbd_cut[1]) {
      reco_bkg_cut = false;
    }
    //}

    if (!reco_bkg_cut) { continue; }
  
    // find leading jet after bkg and calo acceptance cuts
    if (unsubjet_pt->size() < nJetReq) { continue; }
    get_leading_jet(ind_lead, unsubjet_pt);
    lead.SetPtEtaPhi(unsubjet_pt->at(ind_lead), unsubjet_eta->at(ind_lead), unsubjet_phi->at(ind_lead));
    lead_e = unsubjet_e->at(ind_lead);

    // Trigger efficiency.
    double jettrigeff = f_turnon->Eval(lead.Pt());
    if (jettrigeff < 0.01) jettrigeff = 0.01;
    double jettrig_scale = 1.0 / jettrigeff;

    double jettrigeff_down = f_turnon_down->Eval(lead.Pt());
    if (jettrigeff_down < 0.01) jettrigeff_down = 0.01;
    double jettrig_scale_down = 1.0 / jettrigeff_down;

    double jettrigeff_up = f_turnon_up->Eval(lead.Pt());
    if (jettrigeff_up < 0.01) jettrigeff_up = 0.01;
    double jettrig_scale_up = 1.0 / jettrigeff_up;

    //std::cout << "Reco lead: pt " << lead.Pt() << " e " << lead_e << " eta " << lead.Eta() << " phi " << lead.Phi()  << " bkg_cut " << reco_bkg_cut << std::endl;

    TVector3 caliblead;
    if (!sim) {
      caliblead.SetPtEtaPhi(f_corr->Eval(lead.Pt()), lead.Eta(), lead.Phi());
    } else {
      caliblead.SetPtEtaPhi(f_corr->Eval(lead.Pt()) * (1 + randGen.Gaus(0.0, 0.1)) * 1, lead.Eta(), lead.Phi());
    }

    //std::cout << "Good reco lead: pt " << caliblead.Pt() << " eta " << caliblead.Eta() << " phi " << caliblead.Phi() << std::endl;
    //if (bkg_cut == "dijet") { std::cout << "reco sub: pt " << unsubjet_pt->at(ind_sub) << " eta " << unsubjet_eta->at(ind_sub) << " phi " << unsubjet_phi->at(ind_sub) << std::endl; }
    
    //////////////////////////// SETUP UE VARIABLES //////////////////////////////

    // define energy variables in the towards, transverse and away regions
    float et_towards = 0;
    float et_transverse = 0;
    float et_away = 0;
    int ntopo_towards[] = {0, 0, 0, 0, 0, 0, 0, 0};
    int ntopo_transverse[] = {0, 0, 0, 0, 0, 0, 0, 0};
    int ntopo_away[] = {0, 0, 0, 0, 0, 0, 0, 0};
    float sume_topo_towards[] = {0, 0, 0, 0, 0, 0, 0, 0};
    float sume_topo_transverse[] = {0, 0, 0, 0, 0, 0, 0, 0};
    float sume_topo_away[] = {0, 0, 0, 0, 0, 0, 0, 0};
    if (!clusters) {
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
    } else {
      for (int i = 0; i < clsmult; i++) {
        TVector3 cls;
        cls.SetPtEtaPhi(cluster_e[i] / cosh(cluster_eta[i]), cluster_eta[i], cluster_phi[i]);  // define cluster vector
        float dphi = lead.DeltaPhi(cls);                                                       // find the deltaphi between leading jet and cluster
        if (fabs(dphi) < M_PI / 3.0) {  // towards region
          et_towards += cluster_e[i] / cosh(cluster_eta[i]);
          for (int j = 0; j < 8; j++) {
            if (cluster_e[i] > float(topo_thresholds[j] / 1000.0)) {
              ntopo_towards[j] += 1;
              sume_topo_towards[j] += cluster_e[i] / cosh(cluster_eta[i]);
              h_topo_towards[j]->Fill(cluster_e[i] / cosh(cluster_eta[i]), weight_scale*vertex_weight);
            }
          }
        }
        else if (fabs(dphi) > M_PI / 3.0 && fabs(dphi) < (2.0 * M_PI) / 3.0) {  // transverse region
          et_transverse += cluster_e[i] / cosh(cluster_eta[i]);
          for (int j = 0; j < 8; j++) {
            if (cluster_e[i] > float(topo_thresholds[j] / 1000.0)) {
              ntopo_transverse[j] += 1;
              sume_topo_transverse[j] += cluster_e[i] / cosh(cluster_eta[i]);
              h_topo_transverse[j]->Fill(cluster_e[i] / cosh(cluster_eta[i]), weight_scale*vertex_weight);
            }
          }
        }
        else if (fabs(dphi) > (2.0 * M_PI) / 3.0) {  // away region
          et_away += cluster_e[i] / cosh(cluster_eta[i]);
          for (int j = 0; j < 8; j++) {
            if (cluster_e[i] > float(topo_thresholds[j] / 1000.0)) {
              ntopo_away[j] += 1;
              sume_topo_away[j] += cluster_e[i] / cosh(cluster_eta[i]);
              h_topo_away[j]->Fill(cluster_e[i] / cosh(cluster_eta[i]), weight_scale*vertex_weight);
            }
          }
        }
      }
    }

    float pu_correct_et_transverse = et_transverse;
    float pu_up_et_transverse = et_transverse;
    float pu_down_et_transverse = et_transverse;
    if (!sim) {
      for (int i = 0; i < run_numbers.size(); i++) {
        if (run_numbers[i] == runnumber) {
          pu_correct_et_transverse -= pileup_correction[0]*100.0*pileup_rates[i];
          pu_up_et_transverse -= pileup_correction[1]*100.0*pileup_rates[i];
          pu_down_et_transverse -= pileup_correction[2]*100*pileup_rates[i];
          //std::cout << "pile up rate: " << pileup_rates[i] << " correction: " << pileup_correction[0]*pileup_rates[i]*100.0 << " " << et_transverse << " " << pu_correct_et_transverse << std::endl;
          break;
        }
      }
    }

    //std::cout << "Reco ET " << et_transverse << std::endl;
    //std::cout << std::endl;

    //////////////////////////// RECORD QA PLOTS FOR ALL EVENTS IN UNFOLDING PROCEDURE ////////////////////////////
    h_lead_spectra_record->Fill(caliblead.Pt(), weight_scale*vertex_weight);
    h_et_transverse_record->Fill(et_transverse, weight_scale*vertex_weight);
    h_ue_pt_transverse_record->Fill(caliblead.Pt(),et_transverse, weight_scale*vertex_weight);

    h_pu_correct_et_transverse_record->Fill(pu_correct_et_transverse, weight_scale*vertex_weight);
    h_pu_correct_ue_pt_transverse_record->Fill(caliblead.Pt(), pu_correct_et_transverse, weight_scale*vertex_weight);
    h_et_transverse_pu_correct_et_transverse_record->Fill(et_transverse, pu_correct_et_transverse, weight_scale*vertex_weight);

    pileup = 1.0;
    if (!sim) {
      for (int i = 0; i < run_numbers.size(); i++) {
        if (run_numbers[i] == runnumber) { pileup += pileup_rates[i]; }
      }
    }

    h_et_transverse_pu_correction_pu_rate_record->Fill(et_transverse, et_transverse - pu_correct_et_transverse, pileup-1.0 * weight_scale*vertex_weight);

    double uni_meas_et = MapToUniform(et_transverse, calibetbins, calibnet);
    double uni_meas_pt = MapToUniform(caliblead.Pt(), calibptbins, calibnpt);
    double uni_pu_correct_meas_et = MapToUniform(pu_correct_et_transverse, calibetbins, calibnet);
    double uni_pu_up_meas_et = MapToUniform(pu_up_et_transverse, calibetbins, calibnet);
    double uni_pu_down_meas_et = MapToUniform(pu_down_et_transverse, calibetbins, calibnet);

    if (uni_meas_pt >= 0 && uni_meas_pt < 1 && uni_meas_et >= 0 && uni_meas_et < 1) {
      h_calib_jet_pt_tight->Fill(caliblead.Pt(), weight_scale*vertex_weight);
      h_calib_calo_et_tight->Fill(et_transverse, weight_scale*vertex_weight);
      h_calib_jet_pt_uni_tight->Fill(uni_meas_pt, weight_scale*vertex_weight);
      h_calib_calo_et_uni_tight->Fill(uni_meas_et, weight_scale*vertex_weight);
      if (!sim) {
        h_unscale_average_et->Fill(runnumber, et_transverse, weight_scale*vertex_weight);
        h_scale_average_et->Fill(runnumber, et_transverse/pileup, weight_scale*vertex_weight);
        for (int i = 0; i < unsubjet_pt->size(); i++) {
          if (unsubjet_pt->at(i) > 15.0) {
            njets++;
            ejets += unsubjet_pt->at(i);
          }
        }
        h_scale_njet->Fill(runnumber, njets/pileup, weight_scale*vertex_weight);
        h_unscale_ejet->Fill(runnumber, ejets/njets, weight_scale*vertex_weight);
      }
      h_et_towards->Fill(et_towards, weight_scale*vertex_weight);
      h_et_transverse->Fill(et_transverse, weight_scale*vertex_weight);
      h_et_away->Fill(et_away, weight_scale*vertex_weight);
      h_ue_towards->Fill(et_towards / (secteta * sectphi), weight_scale*vertex_weight);
      h_ue_transverse->Fill(et_transverse / (secteta * sectphi), weight_scale*vertex_weight);
      h_ue_away->Fill(et_away / (secteta * sectphi), weight_scale*vertex_weight);

      for (int i = 0; i < 8; i++)
      {
        h_ntopo_towards[i]->Fill(ntopo_towards[i], weight_scale*vertex_weight);
        h_ntopo_transverse[i]->Fill(ntopo_transverse[i], weight_scale*vertex_weight);
        h_ntopo_away[i]->Fill(ntopo_away[i], weight_scale*vertex_weight);
        h_sume_topo_towards[i]->Fill(sume_topo_towards[i], weight_scale*vertex_weight);
        h_sume_topo_transverse[i]->Fill(sume_topo_transverse[i], weight_scale*vertex_weight);
        h_sume_topo_away[i]->Fill(sume_topo_away[i], weight_scale*vertex_weight);
      }

      h_ue_pt_towards->Fill(lead.Pt(), et_towards, weight_scale*vertex_weight);
      h_ue_pt_transverse->Fill(lead.Pt(), et_transverse, weight_scale*vertex_weight);
      h_ue_pt_away->Fill(lead.Pt(), et_away, weight_scale*vertex_weight);
    }
    
    //////////////////////////// FILL HISTOGRAMS FOR UNFOLDING ////////////////////////////
    //if (caliblead.Pt() >= calibptbins[0] && caliblead.Pt() <= calibptbins[calibnpt] && et_transverse >= calibetbins[0] && et_transverse <= calibetbins[calibnet]) {
    if (uni_meas_pt >= 0 && uni_meas_pt < 1 && uni_meas_et >= 0 && uni_meas_et < 1) {
      if (!sim) {
        h_calibjet_pt_dijet_eff->Fill(uni_meas_pt, uni_meas_et, jettrig_scale*timingcut_scale);
        h_calibjet_pt_dijet_effdown->Fill(uni_meas_pt, uni_meas_et, jettrig_scale_down*timingcut_scale);
        h_calibjet_pt_dijet_effup->Fill(uni_meas_pt, uni_meas_et, jettrig_scale_up*timingcut_scale);
        h_calibjet_pt_timingeffup->Fill(uni_meas_pt, uni_meas_et, jettrig_scale*timingcut_scale_up);
        h_jetpt->Fill(uni_meas_pt, jettrig_scale*timingcut_scale);
        h_caloet->Fill(uni_meas_et, jettrig_scale*timingcut_scale);
      } else {
        h_calibjet_pt_dijet_eff->Fill(uni_meas_pt, uni_meas_et, weight_scale*vertex_weight);
        h_calibjet_pt_dijet_effdown->Fill(uni_meas_pt, uni_meas_et, weight_scale*vertex_weight);
        h_calibjet_pt_dijet_effup->Fill(uni_meas_pt, uni_meas_et, weight_scale*vertex_weight);
        h_calibjet_pt_timingeffup->Fill(uni_meas_pt, uni_meas_et, weight_scale*vertex_weight);
        h_jetpt->Fill(uni_meas_pt, weight_scale*vertex_weight);
        h_caloet->Fill(uni_meas_et, weight_scale*vertex_weight);
      }
    }
    //if (caliblead.Pt() >= calibptbins[0] && caliblead.Pt() <= calibptbins[calibnpt] && pu_correct_et_transverse >= calibetbins[0] && pu_correct_et_transverse <= calibetbins[calibnet]) {
    if (!sim) {
     if (uni_meas_pt >= 0 && uni_meas_pt < 1 && uni_pu_correct_meas_et >= 0 && uni_pu_correct_meas_et < 1) { h_calibjet_pt_dijet_pu_correct_et->Fill(uni_meas_pt, uni_pu_correct_meas_et, jettrig_scale*timingcut_scale); }
     if (uni_meas_pt >= 0 && uni_meas_pt < 1 && uni_pu_up_meas_et >= 0 && uni_pu_up_meas_et < 1) { h_calibjet_pt_pu_up_et->Fill(uni_meas_pt, uni_pu_up_meas_et, jettrig_scale*timingcut_scale); }
     if (uni_meas_pt >= 0 && uni_meas_pt < 1 && uni_pu_down_meas_et >= 0 && uni_pu_down_meas_et < 1) { h_calibjet_pt_pu_down_et->Fill(uni_meas_pt, uni_pu_down_meas_et, jettrig_scale*timingcut_scale); }
    }

  } // event loop end

  OutputVarBinMapping(h_calibptbins, calibptbins, int(sizeof(calibptbins)/sizeof(calibptbins[0])));
  OutputVarBinMapping(h_calibetbins, calibetbins, int(sizeof(calibetbins)/sizeof(calibetbins[0])));

  // Write histograms.
  std::cout << "Writing histograms..." << std::endl;
  f_out->cd();
  h_zvertex->Write();
  h_lead_spectra_record->Write();
  h_et_transverse_record->Write();
  h_ue_pt_transverse_record->Write();
  h_calibjet_pt_dijet_eff->Write();
  h_calibjet_pt_dijet_effdown->Write();
  h_calibjet_pt_dijet_effup->Write();
  h_calibjet_pt_dijet_pu_correct_et->Write();
  h_calibjet_pt_pu_up_et->Write();
  h_calibjet_pt_pu_down_et->Write();
  h_calibjet_pt_timingeffup->Write();
  h_jetpt->Write();
  h_caloet->Write();
  h_unscale_average_et->Write();
  h_scale_average_et->Write();
  h_scale_njet->Write();
  h_unscale_ejet->Write();
  h_calib_jet_pt_tight->Write();
  h_calib_calo_et_tight->Write();
  h_calib_jet_pt_uni_tight->Write();
  h_calib_calo_et_uni_tight->Write();
  h_calibptbins->Write();
  h_calibetbins->Write();
  h_pu_correct_et_transverse_record->Write();
  h_pu_correct_ue_pt_transverse_record->Write();
  h_et_transverse_pu_correct_et_transverse_record->Write();
  h_et_transverse_pu_correction_pu_rate_record->Write();
  h_et_towards->Write();
  h_et_transverse->Write();
  h_et_away->Write();
  h_ue_towards->Write();
  h_ue_transverse->Write();
  h_ue_away->Write();
  for (int i = 0; i < 8; i++)
  {
    h_ntopo_towards[i]->Write();
    h_ntopo_transverse[i]->Write();
    h_ntopo_away[i]->Write();
    h_sume_topo_towards[i]->Write();
    h_sume_topo_transverse[i]->Write();
    h_sume_topo_away[i]->Write();
  }
  h_ue_pt_towards->Write();
  h_ue_pt_transverse->Write();
  h_ue_pt_away->Write();
  f_out->Close();
  std::cout << "All done!" << std::endl;
}

////////////////////////////////////////// Functions //////////////////////////////////////////
bool check_bad_trigger(std::vector<int>* gl1_trigger_vector_scaled) {
  if (std::find(gl1_trigger_vector_scaled->begin(), gl1_trigger_vector_scaled->end(), 22) != gl1_trigger_vector_scaled->end()) { return true; }
  //if (std::find(gl1_trigger_vector_scaled->begin(), gl1_trigger_vector_scaled->end(), 34) != gl1_trigger_vector_scaled->end()) { return true; }
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
  float dphi = get_dphi(leadingjet_phi, subleadingjet_phi);
  //std::cout << "delta phi: " << dphi << std::endl;
  return fabs(dphi) > dijet_min_phi;
}