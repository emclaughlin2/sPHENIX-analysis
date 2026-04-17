#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <TString.h>
#include <iostream>
#include <cmath>

void Hadron_Analysis() {

    TString filenames[3] = {"merged_2sigma.root", "merged_3sigma.root", "merged_4sigma.root"};
    TString tags[3] = {"2s", "3s", "4s"};

    double fit_min[3] = {0.0, 0.0, 0.0};
    double fit_max[3] = {1.4, 2.0, 2.0};

    double init_mean[3]  = {0.442, 0.448, 0.441};
    double init_sigma[3] = {0.305, 0.267, 0.259};

    TFile* outFile = new TFile("Hadron_Isolated_Resolution_3Plots_fit.root", "RECREATE");

    for (int i = 0; i < 3; i++) {

        TFile* f = TFile::Open(filenames[i]);
        if (!f || f->IsZombie()) continue;

        TTree* tree = (TTree*)f->Get("T");
        if (!tree) { f->Close(); continue; }

        Float_t zvtx;
        Int_t clsmult, truthpar_n;
        Float_t cluster_e[450], cluster_eta[450], cluster_phi[450];
        Float_t truthpar_pt[200], truthpar_pz[200], truthpar_eta[200], truthpar_phi[200];
        Int_t truthpar_pid[200];

        tree->SetBranchAddress("zvtx", &zvtx);
        tree->SetBranchAddress("clsmult", &clsmult);
        tree->SetBranchAddress("cluster_e", cluster_e);
        tree->SetBranchAddress("cluster_eta", cluster_eta);
        tree->SetBranchAddress("cluster_phi", cluster_phi);
        tree->SetBranchAddress("truthpar_n", &truthpar_n);
        tree->SetBranchAddress("truthpar_pt", truthpar_pt);
        tree->SetBranchAddress("truthpar_pz", truthpar_pz);
        tree->SetBranchAddress("truthpar_eta", truthpar_eta);
        tree->SetBranchAddress("truthpar_phi", truthpar_phi);
        tree->SetBranchAddress("truthpar_pid", truthpar_pid);

        outFile->cd();

        TH1F* hEP = new TH1F("hEP_hadron_" + tags[i],
                             "Isolated Hadron E/p (" + tags[i] + ");E_{cls}/P_{tru};Entries",
                             100, -2, 4);

        Long64_t nEntries = tree->GetEntries();

        for (Long64_t entry = 0; entry < nEntries; entry++) {

            tree->GetEntry(entry);
            if (fabs(zvtx) > 30) continue;

            for (int itr = 0; itr < truthpar_n; itr++) {

                int pid = truthpar_pid[itr];
                if (abs(pid) != 211) continue;

                if (truthpar_pt[itr] < 3.0 || fabs(truthpar_eta[itr]) >= 0.9) continue;

                bool isIsolated = true;

                for (int jtr = 0; jtr < truthpar_n; jtr++) {

                    if (itr == jtr) continue;

                    if (truthpar_pt[jtr] > 0.2) {

                        double dR_iso = sqrt(pow(truthpar_eta[itr] - truthpar_eta[jtr], 2) +
                                             pow(TVector2::Phi_mpi_pi(truthpar_phi[itr] - truthpar_phi[jtr]), 2));

                        if (dR_iso < 0.6) {
                            isIsolated = false;
                            break;
                        }
                    }
                }

                if (!isIsolated) continue;

                double p_truth = sqrt(pow(truthpar_pt[itr], 2) + pow(truthpar_pz[itr], 2));

                double bestDR = 0.2;
                int bestIdx = -1;

                for (int icl = 0; icl < clsmult; icl++) {

                    double dR_match = sqrt(pow(truthpar_eta[itr] - cluster_eta[icl], 2) +
                                           pow(TVector2::Phi_mpi_pi(truthpar_phi[itr] - cluster_phi[icl]), 2));

                    if (dR_match < bestDR) {
                        bestDR = dR_match;
                        bestIdx = icl;
                    }
                }

                if (bestIdx != -1) {
                    hEP->Fill(cluster_e[bestIdx] / p_truth);
                }
            }
        }

        // ---------- FIT ----------
        TF1* gausFit = new TF1("gausFit_" + tags[i], "gaus", fit_min[i], fit_max[i]);

        gausFit->SetParameters(hEP->GetMaximum(), init_mean[i], init_sigma[i]);

        hEP->Fit(gausFit, "R");

        double mean  = gausFit->GetParameter(1);
        double sigma = gausFit->GetParameter(2);

        // ---------- DRAW ----------
        TCanvas* c = new TCanvas("c_" + tags[i], tags[i], 800, 600);

        hEP->Draw("E");
        gausFit->Draw("same");

        TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
        leg->AddEntry(hEP,"E/p distribution","lep");
        leg->AddEntry(gausFit,
                      Form("Fit: #mu = %.3f, #sigma = %.3f",mean,sigma),
                      "l");

        leg->Draw();

        c->Write();
        hEP->Write();
        gausFit->Write();

        f->Close();
        delete f;
    }

    outFile->Write();
    outFile->Close();

    std::cout << "Successfully saved isolated hadron plots with Gaussian fits." << std::endl;
}
