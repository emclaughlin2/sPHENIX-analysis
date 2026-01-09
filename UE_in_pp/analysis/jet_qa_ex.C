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
#include "TH1D.h"

using namespace std;

void jet_qa_ex(string output_file = "jet_qa_output.root") {
std::ifstream infile("goodRunList010225.txt");
    if (!infile.is_open()) {
        std::cerr << "Error opening run list file!" << std::endl;
        return;
    }

std::vector<int> run_numbers;
std::string line;
while (std::getline(infile, line)) {
    std::istringstream iss(line);
    int run;
    if (iss >> run) {
        run_numbers.push_back(run);
    }
}

if (run_numbers.empty()) {
    std::cerr << "No valid run numbers found." << std::endl;
    return;
}

TFile* outfile = new TFile(output_file.c_str(), "RECREATE");
if (!outfile->IsOpen()) {
    std::cerr << "Failed to open output file!" << std::endl;
    return;
}

for (int run : run_numbers) {
    std::string path = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/output1.5mrad_ana468_" + std::to_string(run) + "_*.root";
    TChain chain("T");
    chain.Add(path.c_str());

    int nJet;
    vector<float>* eta = nullptr;
    vector<float>* phi = nullptr;
    vector<float>* pt = nullptr;
    vector<float>* e = nullptr;
    float zvtx;
    vector<int>* triggerVector = nullptr;

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("nJet", 1);
    chain.SetBranchStatus("eta", 1);
    chain.SetBranchStatus("phi", 1);
    chain.SetBranchStatus("pt", 1);
    chain.SetBranchStatus("e", 1);
    chain.SetBranchStatus("zvtx", 1);
    chain.SetBranchStatus("triggerVector", 1);

    chain.SetBranchAddress("nJet", &nJet);
    chain.SetBranchAddress("eta", &eta);
    chain.SetBranchAddress("phi", &phi);
    chain.SetBranchAddress("pt", &pt);
    chain.SetBranchAddress("e", &e);
    chain.SetBranchAddress("zvtx", &zvtx);
    chain.SetBranchAddress("triggerVector", &triggerVector);

    TH2F* h_etaphi = new TH2F(Form("etaphi_run_%d", run), "", 100, -1.1, 1.1, 100, -M_PI, M_PI);
    TH1F* h_event_count = new TH1F(Form("event_count_run_%d", run), "", 1, 0, 1);

    Long64_t n_entries = chain.GetEntries();
    std::cout << "run: " << run << " total entries: " << n_entries << std::endl;
    for (Long64_t i = 0; i < n_entries; ++i) {
        chain.GetEntry(i);
        if (std::isnan(zvtx) || zvtx < -30 || zvtx > 30) continue;

        bool jettrig = false;
        if (triggerVector) {
            for (int t : *triggerVector) {
                if (t >= 16 && t <= 23) {
                    jettrig = true;
                    break;
                }
            }
        }
        if (!jettrig) continue;

        for (int j = 0; j < nJet; ++j) {
            if ((*pt)[j] > 5.0) {
                h_etaphi->Fill((*eta)[j], (*phi)[j]);
            }
        }

        h_event_count->Fill(0.5);
    }

    double n_events = h_event_count->GetEntries();
    if (n_events > 0) {
        h_etaphi->Scale(1.0 / n_events);
        h_etaphi->Write();
    }

    h_event_count->Write();
}

outfile->Close();
std::cout << "Done! Output saved in: " << output_file << std::endl;
}
