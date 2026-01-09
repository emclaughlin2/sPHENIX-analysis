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
#include "TH1D.h"
#include <RooUnfold.h>
#include <RooUnfoldResponse.h>
#include <RooUnfoldBayes.h>

using namespace std;

//R__LOAD_LIBRARY(/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/analysis/roounfold/libRooUnfold.so)

bool isInRange(float truthJetPt, float mcWeight);
void getLeadSubleadJet(std::vector<float> *pt, std::vector<float> *eta, int &ind_lead, int &ind_sub);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                              //
//  Note: 11.26.24                                                                                                              //
//  This macro creates dijet QA plots and does 1D unfolding of leading jet pT for back-to-back dijet events                     //
//                                                                                                                              //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ue_unfolding(string filename = "ue_unfolding.root") {

    //double ptbins[] = {5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,16,17,18,20,22,24,26,28,31,35,40,50};
    //float dRMax = 0.3;
    //float leadptmin = 10;
    //float subptmin = 5;

    std::cout << "RooUnfold library path: " << gSystem->Which(gSystem->GetDynamicPath(), "libRooUnfold.so", kReadPermission) << std::endl;

    double ptbins[] = {20,22,24,26,28,31,35,40,50};
    int nptbins = sizeof(ptbins) / sizeof(ptbins[0]) - 1;
    //double etbins[] = {0.5,0.6,0.7,0.85,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.9,3.1,3.5,3.9,4.4,5,5.7,6.5,7.5,9,25};
    double etbins[] = {0.5,1.5,2.5,3.5,4.5,5.5,7,9,11,15,25,50};
    int netbins = sizeof(etbins) / sizeof(etbins[0]) - 1;


    // truth level bounds should be bigger on both low and high side 

    float dRMax = 0.3;
    float leadptmin = 20;
    float subptmin = 15;
    float ETmin = 0.5;

    bool doUnfolding = false;
    bool sim = false;
    bool clusters = true;
    bool emcal_clusters = false;
    bool applyCorr = true;
    
    // event and jet histograms 
    TH1F* h_pass_deltaphi = new TH1F("h_pass_deltaphi","",125,-2*M_PI,2*M_PI);
    TH1F* h_pass_xj = new TH1F("h_pass_xj","",20,0,1);
    TH1F* h_pass_spectra = new TH1F("h_pass_spectra","",50,0,50);
    TH2F* h_pass_aj_ptavg = new TH2F("h_pass_aj_ptavg","",100,0,100,100,0,1);
    TH1F* h_pass_truth_deltaphi = new TH1F("h_pass_truth_deltaphi","",125,-2*M_PI,2*M_PI);
    TH1F* h_pass_truth_xj = new TH1F("h_pass_truth_xj","",20,0,1);
    TH1F* h_pass_truth_spectra = new TH1F("h_pass_truth_spectra","",50,0,50);
    TH2F* h_pass_truth_aj_ptavg = new TH2F("h_pass_truth_aj_ptavg","",100,0,100,100,0,1);

    // calorimeter tower/cluster histograms 
    TH1F* h_et_towards = new TH1F("h_et_towards","",700,-20,50);
    TH1F* h_et_transverse = new TH1F("h_et_transverse","",700,-20,50);
    TH1F* h_et_away = new TH1F("h_et_away","",700,-20,50);
    TH1F* h_ue_towards = new TH1F("h_ue_towards","",700,-20,50);
    TH1F* h_ue_transverse = new TH1F("h_ue_transverse","",700,-20,50);
    TH1F* h_ue_away = new TH1F("h_ue_away","",700,-20,50);
    TH1F* h_et_truth_towards = new TH1F("h_et_truth_towards","",700,-20,50);
    TH1F* h_et_truth_transverse = new TH1F("h_et_truth_transverse","",700,-20,50);
    TH1F* h_et_truth_away = new TH1F("h_et_truth_away","",700,-20,50);
    TH1F* h_ue_truth_towards = new TH1F("h_ue_truth_towards","",700,-20,50);
    TH1F* h_ue_truth_transverse = new TH1F("h_ue_truth_transverse","",700,-20,50);
    TH1F* h_ue_truth_away = new TH1F("h_ue_truth_away","",700,-20,50);
    int topo_thresholds[] = {-9999,0,100,200,300,500};
    TH1F* h_ntopo_towards[6];
    TH1F* h_ntopo_transverse[6];
    TH1F* h_ntopo_away[6];
    TH1F* h_sume_topo_towards[6];
    TH1F* h_sume_topo_transverse[6];
    TH1F* h_sume_topo_away[6];
    for (int i = 0; i < 6; i++) {
        h_ntopo_towards[i] = new TH1F(Form("h_ntopo%d_towards",topo_thresholds[i]),"",200,0,200);
        h_ntopo_transverse[i] = new TH1F(Form("h_ntopo%d_transverse",topo_thresholds[i]),"",200,0,200);
        h_ntopo_away[i] = new TH1F(Form("h_ntopo%d_away",topo_thresholds[i]),"",200,0,200);
        h_sume_topo_towards[i] = new TH1F(Form("h_sume_topo%d_towards",topo_thresholds[i]),"",550,-10,100);
        h_sume_topo_transverse[i] = new TH1F(Form("h_sume_topo%d_transverse",topo_thresholds[i]),"",550,-10,100);
        h_sume_topo_away[i] = new TH1F(Form("h_sume_topo%d_away",topo_thresholds[i]),"",550,-10,100);
    }

    // jet vs calorimeter tower/cluster histograms 
    TProfile* h_ue_xj_towards = new TProfile("h_ue_xj_towards","",20,0,1);
    TProfile* h_ue_xj_transverse = new TProfile("h_ue_xj_transverse","",20,0,1);
    TProfile* h_ue_xj_away = new TProfile("h_ue_xj_away","",20,0,1);
    TProfile* h_ue_pt_towards = new TProfile("h_ue_pt_towards","",nptbins, ptbins);
    TProfile* h_ue_pt_transverse = new TProfile("h_ue_pt_transverse","",nptbins, ptbins);
    TProfile* h_ue_pt_away = new TProfile("h_ue_pt_away","",nptbins, ptbins);
    TProfile* h_ue_xj_truth_towards = new TProfile("h_ue_xj_truth_towards","",20,0,1);
    TProfile* h_ue_xj_truth_transverse = new TProfile("h_ue_xj_truth_transverse","",20,0,1);
    TProfile* h_ue_xj_truth_away = new TProfile("h_ue_xj_truth_away","",20,0,1);
    TProfile* h_ue_pt_truth_towards = new TProfile("h_ue_pt_truth_towards","",nptbins, ptbins);
    TProfile* h_ue_pt_truth_transverse = new TProfile("h_ue_pt_truth_transverse","",nptbins, ptbins);
    TProfile* h_ue_pt_truth_away = new TProfile("h_ue_pt_truth_away","",nptbins, ptbins);

    //Needed for gaus function
    TRandom3 obj;

    //random number in each event
    TRandom3 Random;

    // unfolding histograms 
    TProfile* jes_ratio = new TProfile("jes_ratio","",nptbins,ptbins);

    //defining Meas and Truth Histograms
    TH1D* hMeas = new TH1D("hMeas","",nptbins,ptbins);
    TH1D* hTruth = new TH1D("hTruth","",nptbins,ptbins);
    //TH2D* hMeas = new TH2D("hMeas","",nptbins,ptbins,netbins,etbins);
    //TH2D* hTruth = new TH2D("hTruth","",nptbins,ptbins,netbins,etbins);

    // closure test histograms 
    TH2D* hMeasHalf = new TH2D("hMeasHalf","",nptbins,ptbins,netbins,etbins);
    TH2D* hTruthHalf = new TH2D("hTruthHalf","",nptbins,ptbins,netbins,etbins);

    //RooUnfoldResponse* resp_full = new RooUnfoldResponse("resp_full",""); resp_full->Setup(hMeas, hTruth);
    //RooUnfoldResponse* resp_half = new RooUnfoldResponse("resp_half",""); resp_half->Setup(hMeasHalf, hTruthHalf);

    //making response matrices
    RooUnfoldResponse *resp_full = new RooUnfoldResponse(hMeas,hTruth,"resp_full","",true);
    RooUnfoldResponse *resp_half = new RooUnfoldResponse(hMeasHalf,hTruthHalf,"resp_half","",true);

    //histograms for errors
    //TH2D* hResponseTruthMeasFull = new TH2D("hResponseTruthMeasFull","",etbins,etmin,etmax,etbins,etmin,etmax);
    //TH2D* hResponseTruthMeasHalf = new TH2D("hResponseTruthMeasHalf","",etbins,etmin,etmax,etbins,etmin,etmax);

    TChain chain("T");
    const char* inputDirectory = "/sphenix/tg/tg01/jets/egm2153/";
    //TString wildcardPath = TString::Format("%sJetValOutput/sim_truth_jet_output.root", inputDirectory); // run 15 dataset
    //chain.Add(wildcardPath);
    
    //for (int i = 1; i < 3334; i++) {
    //    TString wildcardPath = TString::Format("%sUEinppOutput/sim_run22_jet10_3sigma_output_%d.root", inputDirectory, i); // run 22 dataset
    //    chain.Add(wildcardPath);
    //}

    for (int i = 1; i < 3; i++) {
        TString wildcardPath = TString::Format("%sUEinppOutput/sim_run22_jet10_3sigma_output_%d.root", inputDirectory, i); // run 22 dataset
        chain.Add(wildcardPath);
    }

    int m_event;
    int nJet;
    int nTruthJet;
    float zvtx;
    float deltaeta = 0.0916667;
    float deltaphi = 0.0981748;
    float secteta = 2.2;
    float sectphi = (2.0*M_PI)/3.0;

    vector<int> *triggerVector = nullptr;
    vector<float> *eta = nullptr;
    vector<float> *phi = nullptr;
    vector<float> *e = nullptr;
    vector<float> *pt = nullptr;
    vector<float> *truthEta = nullptr;
    vector<float> *truthPhi = nullptr;
    vector<float> *truthE = nullptr;
    vector<float> *truthPt = nullptr;

    int emcaln = 0;
    float emcale[24576] = {0.0};
    float emcaleta[24576] = {0.0};
    float emcalphi[24576] = {0.0};
    int ihcaln = 0;
    float ihcale[1536] = {0.0};
    float ihcaleta[1536] = {0.0};
    float ihcalphi[1536] = {0.0};
    int ohcaln = 0;
    float ohcale[1536] = {0.0};
    float ohcaleta[1536] = {0.0};
    float ohcalphi[1536] = {0.0};
    int clsmult = 0;
    float cluster_e[10000] = {0.0};
    float cluster_eta[10000] = {0.0};
    float cluster_phi[10000] = {0.0};

    int truthpar_n = 0;
    float truthpar_e[100000] = {0.0};
    float truthpar_eta[100000] = {0.0};
    float truthpar_phi[100000] = {0.0};
    int truthpar_pid[100000] = {0};

    chain.SetBranchStatus("*", 0);

    chain.SetBranchStatus("m_event",1);
    chain.SetBranchStatus("nJet",1);
    chain.SetBranchStatus("nTruthJet", 1);
    chain.SetBranchStatus("zvtx",1);
    chain.SetBranchStatus("triggerVector",1);
    chain.SetBranchStatus("eta",1);
    chain.SetBranchStatus("phi",1);
    chain.SetBranchStatus("e",1);
    chain.SetBranchStatus("pt",1);
    chain.SetBranchStatus("truthEta",1);
    chain.SetBranchStatus("truthPhi",1);
    chain.SetBranchStatus("truthE",1);
    chain.SetBranchStatus("truthPt",1);

    chain.SetBranchAddress("m_event",&m_event);
    chain.SetBranchAddress("nJet",&nJet);
    chain.SetBranchAddress("nTruthJet",&nTruthJet);
    chain.SetBranchAddress("zvtx",&zvtx);
    chain.SetBranchAddress("triggerVector",&triggerVector);
    chain.SetBranchAddress("eta",&eta);
    chain.SetBranchAddress("phi",&phi);
    chain.SetBranchAddress("e",&e);
    chain.SetBranchAddress("pt",&pt);
    chain.SetBranchAddress("truthEta",&truthEta);
    chain.SetBranchAddress("truthPhi",&truthPhi);
    chain.SetBranchAddress("truthE",&truthE);
    chain.SetBranchAddress("truthPt",&truthPt);

    if (!clusters) {
        chain.SetBranchStatus("emcaln", 1);
        chain.SetBranchStatus("emcale",1);
        chain.SetBranchStatus("emcaleta", 1);
        chain.SetBranchStatus("emcalphi", 1);
        chain.SetBranchStatus("ihcaln", 1);
        chain.SetBranchStatus("ihcale", 1);
        chain.SetBranchStatus("ihcaleta", 1);
        chain.SetBranchStatus("ihcalphi", 1);
        chain.SetBranchStatus("ohcaln", 1);
        chain.SetBranchStatus("ohcale", 1);
        chain.SetBranchStatus("ohcaleta", 1);
        chain.SetBranchStatus("ohcalphi", 1);

        chain.SetBranchAddress("emcaln",&emcaln);
        chain.SetBranchAddress("emcale",emcale);
        chain.SetBranchAddress("emcaleta",emcaleta);
        chain.SetBranchAddress("emcalphi",emcalphi);
        chain.SetBranchAddress("ihcaln",&ihcaln);
        chain.SetBranchAddress("ihcale",ihcale);
        chain.SetBranchAddress("ihcaleta",ihcaleta);
        chain.SetBranchAddress("ihcalphi",ihcalphi);
        chain.SetBranchAddress("ohcaln",&ohcaln);
        chain.SetBranchAddress("ohcale",ohcale);
        chain.SetBranchAddress("ohcaleta",ohcaleta);
        chain.SetBranchAddress("ohcalphi",ohcalphi);
    }

    if (!emcal_clusters) {
        chain.SetBranchStatus("clsmult", 1);
        chain.SetBranchStatus("cluster_e", 1);
        chain.SetBranchStatus("cluster_eta", 1);
        chain.SetBranchStatus("cluster_phi", 1);
        chain.SetBranchAddress("clsmult",&clsmult);
        chain.SetBranchAddress("cluster_e",cluster_e);
        chain.SetBranchAddress("cluster_eta",cluster_eta);
        chain.SetBranchAddress("cluster_phi",cluster_phi);
    } else {
        chain.SetBranchStatus("emcal_clsmult", 1);
        chain.SetBranchStatus("emcal_cluster_e", 1);
        chain.SetBranchStatus("emcal_cluster_eta", 1);
        chain.SetBranchStatus("emcal_cluster_phi", 1);
        chain.SetBranchAddress("emcal_clsmult",&clsmult);
        chain.SetBranchAddress("emcal_cluster_e",cluster_e);
        chain.SetBranchAddress("emcal_cluster_eta",cluster_eta);
        chain.SetBranchAddress("emcal_cluster_phi",cluster_phi);
    }

    chain.SetBranchStatus("truthpar_n", 1);
    chain.SetBranchStatus("truthpar_e", 1);
    chain.SetBranchStatus("truthpar_eta", 1);
    chain.SetBranchStatus("truthpar_phi", 1);
    chain.SetBranchStatus("truthpar_pid", 1);
    chain.SetBranchAddress("truthpar_n",&truthpar_n);
    chain.SetBranchAddress("truthpar_e",truthpar_e);
    chain.SetBranchAddress("truthpar_eta",truthpar_eta);
    chain.SetBranchAddress("truthpar_phi",truthpar_phi);
    chain.SetBranchAddress("truthpar_pid",truthpar_pid);

    TFile *corrFile = new TFile("JES_Calib_Default.root", "READ");
    if (!corrFile) {
        std::cout << "Error: cannot open JES_Calib_Default.root" << std::endl;
        return;
    }
    TF1 *correction = (TF1*)corrFile->Get("JES_Calib_Default_Func");
    if (!correction) {
        std::cout << "Error: cannot open corection" << std::endl;
        return;
    }

    int eventnumber = 0;
    int events = 0;
    Long64_t nEntries = chain.GetEntries();
    std::cout << nEntries << std::endl;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
    //for (Long64_t entry = 0; entry < 2000; ++entry) {
        chain.GetEntry(entry);
        if (eventnumber % 10000 == 0) cout << "event " << eventnumber << endl;
        eventnumber++;

        bool negJet = false;
        for (int i = 0; i < nJet; i++) {
            if (e->at(i) < 0) {
                negJet = true;
            }
        }

        // require event with |zvtx| < 30 cm 
        if (isnan(zvtx)) { continue; }
        if (zvtx < -30 || zvtx > 30) { continue; }
        if (negJet) { continue; }
        /*
        for (int i = 0; i < truthEta->size();) {
            if (fabs(truthEta->at(i)) > 0.7) {
                truthEta->erase(truthEta->begin() + i);
                truthPt->erase(truthPt->begin() + i);
                truthE->erase(truthE->begin() + i);
                truthPhi->erase(truthPhi->begin() + i);
            } else {
                ++i;
            }
        }
        int nTruthJet = truthPt->size();
        */

        //////////////////////////// SETUP JET VARIABLES FOR UNFOLDING ////////////////////////////

        // indices to find leading and subleading jets 
        int ind_truth_lead = -1;
        int ind_truth_sub = -1;
        int ind_lead = -1;
        int ind_sub = -1;

        // if both nreco jets < 2 and nTruthJet jets < 2, discard event
        if (nJet < 2 && nTruthJet < 2) {
            continue;
        }

        // find reco leading and subleading jets if nreco jets >= 2
        TVector3 lead, sub;
        if (nJet >= 2) {
            getLeadSubleadJet(pt, eta, ind_lead, ind_sub);  
            lead.SetPtEtaPhi(correction->Eval(pt->at(ind_lead)), eta->at(ind_lead), phi->at(ind_lead));
            sub.SetPtEtaPhi(correction->Eval(pt->at(ind_sub)), eta->at(ind_sub), phi->at(ind_sub));
        } else {
            lead.SetPtEtaPhi(0,0,0);
            sub.SetPtEtaPhi(0,0,0);
        }

        //find truth leading and subleading jets if nTruthJet jets >= 2
        TVector3 truthlead, truthsub;
        if (nTruthJet >= 2) {
            getLeadSubleadJet(truthPt, truthEta, ind_truth_lead, ind_truth_sub);
            truthlead.SetPtEtaPhi(truthPt->at(ind_truth_lead), truthEta->at(ind_truth_lead), truthPhi->at(ind_truth_lead));
            truthsub.SetPtEtaPhi(truthPt->at(ind_truth_sub), truthEta->at(ind_truth_sub), truthPhi->at(ind_truth_sub));
        } else {
            truthlead.SetPtEtaPhi(0,0,0);
            truthsub.SetPtEtaPhi(0,0,0);
        }

        //////////////////////////// SETUP UE VARIABLES FOR UNFOLDING //////////////////////////////

        // find reco ET information 
        float et_towards = 0;
        float et_transverse = 0;
        float et_away = 0;
        int ntopo_towards[] = {0,0,0,0,0,0};
        int ntopo_transverse[] = {0,0,0,0,0,0};
        int ntopo_away[] = {0,0,0,0,0,0};
        float sume_topo_towards[] = {0,0,0,0,0,0};
        float sume_topo_transverse[] = {0,0,0,0,0,0};
        float sume_topo_away[] = {0,0,0,0,0,0};

        if (!clusters) {
            for (int i = 0; i < emcaln; i++) {
                TVector3 em;
                em.SetPtEtaPhi(emcale[i]/cosh(emcaleta[i]),emcaleta[i],emcalphi[i]);
                float dphi = lead.DeltaPhi(em);
                if (fabs(dphi) < M_PI/3.0) { et_towards += emcale[i]/cosh(emcaleta[i]); } 
                else if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += emcale[i]/cosh(emcaleta[i]); } 
                else if (fabs(dphi) > (2.0*M_PI)/3.0) { et_away += emcale[i]/cosh(emcaleta[i]); }
            }
            for (int i = 0; i < ihcaln; i++) {
                TVector3 ih;
                ih.SetPtEtaPhi(ihcale[i]/cosh(ihcaleta[i]),ihcaleta[i],ihcalphi[i]);
                float dphi = lead.DeltaPhi(ih);
                if (fabs(dphi) < M_PI/3.0) { et_towards += ihcale[i]/cosh(ihcaleta[i]); } 
                else if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += ihcale[i]/cosh(ihcaleta[i]); } 
                else if (fabs(dphi) > (2.0*M_PI)/3.0) { et_away += ihcale[i]/cosh(ihcaleta[i]); }
            }
            for (int i = 0; i < ohcaln; i++) {
                TVector3 oh;
                oh.SetPtEtaPhi(ohcale[i]/cosh(ohcaleta[i]),ohcaleta[i],ohcalphi[i]);
                float dphi = lead.DeltaPhi(oh);
                if (fabs(dphi) < M_PI/3.0) { et_towards += ohcale[i]/cosh(ohcaleta[i]); } 
                else if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { et_transverse += ohcale[i]/cosh(ohcaleta[i]); } 
                else if (fabs(dphi) > (2.0*M_PI)/3.0) { et_away += ohcale[i]/cosh(ohcaleta[i]); }
            }
        } else {
            for (int i = 0; i < clsmult; i++) {
                TVector3 cls;
                cls.SetPtEtaPhi(cluster_e[i]/cosh(cluster_eta[i]),cluster_eta[i],cluster_phi[i]);
                float dphi = lead.DeltaPhi(cls);
                if (fabs(dphi) < M_PI/3.0) {
                    et_towards += cluster_e[i]/cosh(cluster_eta[i]);
                    for (int j = 0; j < 6; j++) {
                        if (cluster_e[i] > float(topo_thresholds[j]/1000.0)) {
                            ntopo_towards[j] += 1;
                            sume_topo_towards[j] += cluster_e[i]/cosh(cluster_eta[i]);
                        }
                    }
                } else if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) {
                    et_transverse += cluster_e[i]/cosh(cluster_eta[i]);
                    for (int j = 0; j < 6; j++) {
                        if (cluster_e[i] > float(topo_thresholds[j]/1000.0)) {
                            ntopo_transverse[j] += 1;
                            sume_topo_transverse[j] += cluster_e[i]/cosh(cluster_eta[i]);
                        }
                    }
                } else if (fabs(dphi) > (2.0*M_PI)/3.0) {
                    et_away += cluster_e[i]/cosh(cluster_eta[i]);
                    for (int j = 0; j < 6; j++) {
                        if (cluster_e[i] > float(topo_thresholds[j]/1000.0)) {
                            ntopo_away[j] += 1;
                            sume_topo_away[j] += cluster_e[i]/cosh(cluster_eta[i]);
                        }
                    }
                }
            }
        }

        // find truth ET information
        float truth_et_towards = 0;
        float truth_et_transverse = 0;
        float truth_et_away = 0;

        for (int i = 0; i < truthpar_n; i++) {
            if (fabs(truthpar_eta[i]) > 1.1) { continue; }
            if ((truthpar_pid[i] == 22 || truthpar_pid[i] == 111) && fabs(truthpar_e[i]) < 0.2) { continue; }
            else if (fabs(truthpar_e[i]) < 0.5) { continue; }
            TVector3 truth;
            truth.SetPtEtaPhi(truthpar_e[i]/cosh(truthpar_eta[i]),truthpar_eta[i],truthpar_phi[i]);
            float dphi = truthlead.DeltaPhi(truth);
            if (fabs(dphi) < M_PI/3.0) { truth_et_towards += truthpar_e[i]/cosh(truthpar_eta[i]); } 
            else if (fabs(dphi) > M_PI/3.0 && fabs(dphi) < (2.0*M_PI)/3.0) { truth_et_transverse += truthpar_e[i]/cosh(truthpar_eta[i]); } 
            else if (fabs(dphi) > (2.0*M_PI)/3.0) { truth_et_away += truthpar_e[i]/cosh(truthpar_eta[i]); }
        }

        //////////////////////////// RECORD QA PLOTS FOR ALL EVENTS IN UNFOLDING PROCEDURE ////////////////////////////

        // if reco or truth info passes dijet pT and back to back criteria, record event 
        if ((nJet >= 2 && lead.Pt() > leadptmin && sub.Pt() > subptmin && fabs(lead.DeltaPhi(sub)) > 2.75 && et_transverse > ETmin) || (nTruthJet >= 2 && truthlead.Pt() > leadptmin && truthsub.Pt() > subptmin && fabs(truthlead.DeltaPhi(truthsub)) > 2.75 && truth_et_transverse > ETmin)) {
            double  choice = Random.Rndm();
            
            if (nJet >= 2 && lead.Pt() > leadptmin && sub.Pt() > subptmin && fabs(lead.DeltaPhi(sub)) > 2.75 && et_transverse > ETmin) {
                h_pass_deltaphi->Fill(lead.DeltaPhi(sub));
                h_pass_xj->Fill(sub.Pt()/lead.Pt());
                h_pass_spectra->Fill(lead.Pt());
                h_pass_spectra->Fill(sub.Pt());
                h_pass_aj_ptavg->Fill((lead.Pt()+sub.Pt())/2.0, (lead.Pt()-sub.Pt())/(lead.Pt()+sub.Pt()));

                h_et_towards->Fill(et_towards);
                h_et_transverse->Fill(et_transverse);
                h_et_away->Fill(et_away);
                h_ue_towards->Fill(et_towards/(secteta*sectphi));
                h_ue_transverse->Fill(et_transverse/(secteta*sectphi));
                h_ue_away->Fill(et_away/(secteta*sectphi));

                for (int i = 0; i < 6; i++) {
                    h_ntopo_towards[i]->Fill(ntopo_towards[i]);
                    h_ntopo_transverse[i]->Fill(ntopo_transverse[i]);
                    h_ntopo_away[i]->Fill(ntopo_away[i]);
                    h_sume_topo_towards[i]->Fill(sume_topo_towards[i]);
                    h_sume_topo_transverse[i]->Fill(sume_topo_transverse[i]);
                    h_sume_topo_away[i]->Fill(sume_topo_away[i]);
                }

                h_ue_xj_towards->Fill(sub.Pt()/lead.Pt(),et_towards);
                h_ue_xj_transverse->Fill(sub.Pt()/lead.Pt(),et_transverse);
                h_ue_xj_away->Fill(sub.Pt()/lead.Pt(),et_away);
                h_ue_pt_towards->Fill(lead.Pt(),et_towards);
                h_ue_pt_transverse->Fill(lead.Pt(),et_transverse);
                h_ue_pt_away->Fill(lead.Pt(),et_away);
            }
            
            if (nTruthJet >= 2 && truthlead.Pt() > leadptmin && truthsub.Pt() > subptmin && fabs(truthlead.DeltaPhi(truthsub)) > 2.75 && truth_et_transverse > ETmin) {
                h_pass_truth_deltaphi->Fill(truthlead.DeltaPhi(truthsub));
                h_pass_truth_xj->Fill(truthsub.Pt()/truthlead.Pt());
                h_pass_truth_spectra->Fill(truthlead.Pt());
                h_pass_truth_spectra->Fill(truthsub.Pt());
                h_pass_truth_aj_ptavg->Fill((truthlead.Pt()+truthsub.Pt())/2.0, (truthlead.Pt()-truthsub.Pt())/(truthlead.Pt()+truthsub.Pt()));

                h_et_truth_towards->Fill(truth_et_towards);
                h_et_truth_transverse->Fill(truth_et_transverse);
                h_et_truth_away->Fill(truth_et_away);
                h_ue_truth_towards->Fill(truth_et_towards/(secteta*sectphi));
                h_ue_truth_transverse->Fill(truth_et_transverse/(secteta*sectphi));
                h_ue_truth_away->Fill(truth_et_away/(secteta*sectphi));

                h_ue_xj_truth_towards->Fill(truthsub.Pt()/truthlead.Pt(),truth_et_towards);
                h_ue_xj_truth_transverse->Fill(truthsub.Pt()/truthlead.Pt(),truth_et_transverse);
                h_ue_xj_truth_away->Fill(truthsub.Pt()/truthlead.Pt(),truth_et_away);
                h_ue_pt_truth_towards->Fill(truthlead.Pt(),truth_et_towards);
                h_ue_pt_truth_transverse->Fill(truthlead.Pt(),truth_et_transverse);
                h_ue_pt_truth_away->Fill(truthlead.Pt(),truth_et_away);
            }
           
            //////////////////////////// START UNFOLDING PROCEDURE ////////////////////////////

            if ((nJet >= 2 && lead.Pt() > leadptmin && sub.Pt() > subptmin && fabs(lead.DeltaPhi(sub)) > 2.75 && et_transverse > ETmin) && (nTruthJet >= 2 && truthlead.Pt() > leadptmin && truthsub.Pt() > subptmin && fabs(truthlead.DeltaPhi(truthsub)) > 2.75) && truth_et_transverse > ETmin) {
                if (truthlead.DeltaR(lead) < dRMax && truthsub.DeltaR(sub) < dRMax) { // should this match be both the leading and subleading? 
                    // MATCH 
                    hMeas->Fill(lead.Pt());
                    hTruth->Fill(truthlead.Pt());
                    resp_full->Fill(lead.Pt(),truthlead.Pt());
                    jes_ratio->Fill(truthlead.Pt(),lead.Pt()/truthlead.Pt());

                    //hMeas->Fill(lead.Pt(),et_transverse);
                    //hTruth->Fill(truthlead.Pt(),truth_et_transverse);
                    //resp_full->Fill(lead.Pt(),et_transverse,truthlead.Pt(),truth_et_transverse);
                    //jes_ratio->Fill(truthlead.Pt(),lead.Pt()/truthlead.Pt());

                    if (choice > 0.5) {
                        hTruthHalf->Fill(truthlead.Pt(),truth_et_transverse);
                        resp_half->Fill(lead.Pt(),et_transverse,truthlead.Pt(),truth_et_transverse);
                    } else {
                        hMeasHalf->Fill(lead.Pt(),et_transverse);
                    }
                } else {
                    // FAKE AND MISS
                    hTruth->Fill(truthlead.Pt());
                    resp_full->Miss(truthlead.Pt());
                    //hTruth->Fill(truthlead.Pt(),truth_et_transverse);
                    //resp_full->Miss(truthlead.Pt(),truth_et_transverse);
                    

                    if (choice > 0.5) {
                        hTruthHalf->Fill(truthlead.Pt(),truth_et_transverse);
                        resp_half->Miss(truthlead.Pt(),truth_et_transverse);
                    }
                    hMeas->Fill(lead.Pt());
                    resp_full->Fake(lead.Pt());
                    //hMeas->Fill(lead.Pt(),et_transverse);
                    //resp_full->Fake(lead.Pt(),et_transverse);
                    

                    if (choice > 0.5) {
                        resp_half->Fake(lead.Pt(),et_transverse);
                    } else {
                        hMeasHalf->Fill(lead.Pt(),et_transverse);
                    }
                }
            } else if (nJet >= 2 && lead.Pt() > leadptmin && sub.Pt() > subptmin && fabs(lead.DeltaPhi(sub)) > 2.75 && et_transverse > ETmin) {
                // FAKE
                hMeas->Fill(lead.Pt());
                resp_full->Fake(lead.Pt());
                //hMeas->Fill(lead.Pt(),et_transverse);
                //resp_full->Fake(lead.Pt(),et_transverse);
                if (choice > 0.5) {
                    resp_half->Fake(lead.Pt(),et_transverse);
                } else {
                    hMeasHalf->Fill(lead.Pt(),et_transverse);
                }
            } else if (nTruthJet >= 2 && truthlead.Pt() > leadptmin && truthsub.Pt() > subptmin && fabs(truthlead.DeltaPhi(truthsub)) > 2.75 && truth_et_transverse > ETmin) {
                // MISS
                hTruth->Fill(truthlead.Pt());
                resp_full->Miss(truthlead.Pt());
                //hTruth->Fill(truthlead.Pt(),truth_et_transverse);
                //resp_full->Miss(truthlead.Pt(),truth_et_transverse);
                if (choice > 0.5) {
                    hTruthHalf->Fill(truthlead.Pt(),truth_et_transverse);
                    resp_half->Miss(truthlead.Pt(),truth_et_transverse);
                }
            }

            events++;
        }

    }

    hMeas->Print("all");
    resp_full->Hmeasured()->Print("all");
    resp_full->Hfakes()->Print("all");
    hTruth->Print("all");
    resp_full->Htruth()->Print("all");
    //resp_full->Hresponse()->Print("all");

    // unfolding 
    std::cout << resp_full->HasFakes() << std::endl;
    RooUnfoldBayes *full_unfold = new RooUnfoldBayes(resp_full,hMeas,1,false,true);
    full_unfold->SetVerbose(3);
    TH1D*  hReco = (TH1D*)full_unfold->Hunfold();
    hReco->SetName("hReco");
    /*
    TMatrixD full_cov_matrix = full_unfold->Ereco(RooUnfold::kCovToy);
    for (int i = 1; i <= hReco->GetNbinsX(); ++i) {
        for (int j = 1; j <= hReco->GetNbinsY(); ++j) {
            int globalBin = hReco->GetBin(i, j);  // Get global bin index
            double error = sqrt(full_cov_matrix(globalBin - 1, globalBin - 1)); // Convert to zero-based index
            hReco->SetBinError(i, j, error);
        }
    }
    */

    //RooUnfoldBayes *half_unfold1 = new RooUnfoldBayes(resp_half,hMeasHalf,1);
    //TH2D* hRecoHalf1 = (TH2D*)half_unfold1->Hunfold();
    //hRecoHalf1->SetName("hRecoHalf1");
    /*
    TMatrixD half_cov_matrix1 = half_unfold1->Ereco(RooUnfold::kCovToy);
    for (int i = 1; i <= hRecoHalf1->GetNbinsX(); ++i) {
        for (int j = 1; j <= hRecoHalf1->GetNbinsY(); ++j) {
            int globalBin = hRecoHalf1->GetBin(i, j);  // Get global bin index
            double error = sqrt(half_cov_matrix1(globalBin - 1, globalBin - 1)); // Convert to zero-based index
            hRecoHalf1->SetBinError(i, j, error);
        }
    }
    */

    //RooUnfoldBayes *half_unfold2 = new RooUnfoldBayes(resp_half,hMeasHalf,2);
    //TH2D* hRecoHalf2 = (TH2D*)half_unfold2->Hunfold();
    //hRecoHalf2->SetName("hRecoHalf2");
    /*TMatrixD half_cov_matrix2 = half_unfold2->Ereco(RooUnfold::kCovToy);
    for (int i = 1; i <= hRecoHalf2->GetNbinsX(); ++i) {
        for (int j = 1; j <= hRecoHalf2->GetNbinsY(); ++j) {
            int globalBin = hRecoHalf2->GetBin(i, j);  // Get global bin index
            double error = sqrt(half_cov_matrix2(globalBin - 1, globalBin - 1)); // Convert to zero-based index
            hRecoHalf2->SetBinError(i, j, error);
        }
    }*/

    //RooUnfoldBayes *half_unfold3 = new RooUnfoldBayes(resp_half,hMeasHalf,3);
    //TH2D* hRecoHalf3 = (TH2D*)half_unfold3->Hunfold();
    //hRecoHalf3->SetName("hRecoHalf3");
    /*TMatrixD half_cov_matrix3 = half_unfold3->Ereco(RooUnfold::kCovToy);
    for (int i = 1; i <= hRecoHalf3->GetNbinsX(); ++i) {
        for (int j = 1; j <= hRecoHalf3->GetNbinsY(); ++j) {
            int globalBin = hRecoHalf3->GetBin(i, j);  // Get global bin index
            double error = sqrt(half_cov_matrix3(globalBin - 1, globalBin - 1)); // Convert to zero-based index
            hRecoHalf3->SetBinError(i, j, error);
        }
    }*/

    //RooUnfoldBayes *half_unfold4 = new RooUnfoldBayes(resp_half,hMeasHalf,4);
    //TH2D* hRecoHalf4 = (TH2D*)half_unfold4->Hunfold();
    //hRecoHalf4->SetName("hRecoHalf4");
    /*TMatrixD half_cov_matrix4 = half_unfold4->Ereco(RooUnfold::kCovToy);
    for (int i = 1; i <= hRecoHalf4->GetNbinsX(); ++i) {
        for (int j = 1; j <= hRecoHalf4->GetNbinsY(); ++j) {
            int globalBin = hRecoHalf4->GetBin(i, j);  // Get global bin index
            double error = sqrt(half_cov_matrix4(globalBin - 1, globalBin - 1)); // Convert to zero-based index
            hRecoHalf4->SetBinError(i, j, error);
        }
    }*/

    //RooUnfoldBayes *half_unfold5 = new RooUnfoldBayes(resp_half,hMeasHalf,5);
    //TH2D* hRecoHalf5 = (TH2D*)half_unfold5->Hunfold();
    //hRecoHalf5->SetName("hRecoHalf5");
    /*TMatrixD half_cov_matrix5 = half_unfold5->Ereco(RooUnfold::kCovToy);
    for (int i = 1; i <= hRecoHalf5->GetNbinsX(); ++i) {
        for (int j = 1; j <= hRecoHalf5->GetNbinsY(); ++j) {
            int globalBin = hRecoHalf5->GetBin(i, j);  // Get global bin index
            double error = sqrt(half_cov_matrix5(globalBin - 1, globalBin - 1)); // Convert to zero-based index
            hRecoHalf5->SetBinError(i, j, error);
        }
    }*/

    //RooUnfoldBayes *half_unfold6 = new RooUnfoldBayes(resp_half,hMeasHalf,6);
    //TH2D* hRecoHalf6 = (TH2D*)half_unfold6->Hunfold();
    //hRecoHalf6->SetName("hRecoHalf6");
    /*TMatrixD half_cov_matrix6 = half_unfold6->Ereco(RooUnfold::kCovToy);
    for (int i = 1; i <= hRecoHalf6->GetNbinsX(); ++i) {
        for (int j = 1; j <= hRecoHalf6->GetNbinsY(); ++j) {
            int globalBin = hRecoHalf6->GetBin(i, j);  // Get global bin index
            double error = sqrt(half_cov_matrix6(globalBin - 1, globalBin - 1)); // Convert to zero-based index
            hRecoHalf6->SetBinError(i, j, error);
        }
    }*/
    
    // output histograms 
    std::cout << filename << std::endl;
    TFile *out = new TFile(filename.c_str(),"RECREATE");

    out->cd();
    h_pass_deltaphi->Write();
    h_pass_xj->Write();
    h_pass_spectra->Write();
    h_pass_spectra->Write();
    h_pass_aj_ptavg->Write();
    h_et_towards->Write();
    h_et_transverse->Write();
    h_et_away->Write();
    h_ue_towards->Write();
    h_ue_transverse->Write();
    h_ue_away->Write();
    /*
    for (int i = 0; i < 6; i++) {
        h_ntopo_towards[i]->Write();
        h_ntopo_transverse[i]->Write();
        h_ntopo_away[i]->Write();
        h_sume_topo_towards[i]->Write();
        h_sume_topo_transverse[i]->Write();
        h_sume_topo_away[i]->Write();
    }
    */
    h_ue_xj_towards->Write();
    h_ue_xj_transverse->Write();
    h_ue_xj_away->Write();
    h_ue_pt_towards->Write();
    h_ue_pt_transverse->Write();
    h_ue_pt_away->Write();
    h_pass_truth_deltaphi->Write();
    h_pass_truth_xj->Write();
    h_pass_truth_spectra->Write();
    h_pass_truth_spectra->Write();
    h_pass_truth_aj_ptavg->Write();
    h_et_truth_towards->Write();
    h_et_truth_transverse->Write();
    h_et_truth_away->Write();
    h_ue_truth_towards->Write();
    h_ue_truth_transverse->Write();
    h_ue_truth_away->Write();
    h_ue_xj_truth_towards->Write();
    h_ue_xj_truth_transverse->Write();
    h_ue_xj_truth_away->Write();
    h_ue_pt_truth_towards->Write();
    h_ue_pt_truth_transverse->Write();
    h_ue_pt_truth_away->Write();

    jes_ratio->Write();
    hMeas->Write();
    hTruth->Write();
    hReco->Write();
    hMeasHalf->Write();
    hTruthHalf->Write();
    //hRecoHalf1->Write();
    //hRecoHalf2->Write();
    //hRecoHalf3->Write();
    //hRecoHalf4->Write();
    //hRecoHalf5->Write();
    //hRecoHalf6->Write();
    resp_full->Write();
    resp_half->Write();
    out->Close();

}

bool isInRange(float truthJetPt, float mcWeight)
{
    float ptCutOff1 = -1; 
    float ptCutOff2 = -1;


    if(abs(mcWeight/39.06e-3 - 1) < 1e-7)//Minimum bias
    {
        //std::cout << "MB event found" << std::endl;
        ptCutOff1 = 0;
        ptCutOff2 = 14;
    }
    else if(abs(mcWeight/3.210e-6 - 1) < 1e-7)
    {
        //std::cout << "10GeV event found" << std::endl;
        ptCutOff1 = 14;
        ptCutOff2 = 37;
    }
    else if(abs(mcWeight/2.178e-9 - 1) < 1e-7)
    {
        //std::cout << "30GeV event found" << std::endl;
        ptCutOff1 = 37;
        ptCutOff2 = 3000;
    }

    if(truthJetPt < ptCutOff2 && truthJetPt >= ptCutOff1) return true;
    return false;
    }
 
void getLeadSubleadJet(std::vector<float> *pt, std::vector<float> *eta, int &ind_lead, int &ind_sub)
{
    float temp_lead = -1;
    float temp_sub = -1;
    if (pt->size() < 2 || eta->size() < 2 || pt->size() != eta->size()) { std::cout << "PT and ETA vectors smaller than 2 or not equal, something is wrong!" << std::endl; return; }

    for (int i = 0; i < pt->size(); i++) {
        if (pt->at(i) > temp_lead) {
            if (temp_lead != -1) {
                temp_sub = temp_lead;
                ind_sub = ind_lead;
            }
            temp_lead = pt->at(i);
            ind_lead = i;
        } else if (pt->at(i) > temp_sub) {
            temp_sub = pt->at(i);
            ind_sub = i;
        }
    }
}