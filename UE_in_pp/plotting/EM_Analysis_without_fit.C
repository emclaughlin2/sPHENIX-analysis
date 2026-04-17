#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TVector2.h>
#include <TString.h>
#include <iostream>
#include <cmath>

void EM_Analysis_without_fit() {

    TString filenames[3] = {"merged_2sigma.root", "merged_3sigma.root", "merged_4sigma.root"};
    TString tags[3] = {"2s", "3s", "4s"};

    TFile* outFile = new TFile("EM_Isolated_Resolution_3Plots_without_fit.root","RECREATE");

    for(int i=0;i<3;i++){

        TFile* f = TFile::Open(filenames[i]);
        if(!f || f->IsZombie()) continue;

        TTree* tree = (TTree*)f->Get("T");
        if(!tree){f->Close(); continue;}

        Float_t zvtx;
        Int_t clsmult, truthpar_n;

        Float_t cluster_e[450], cluster_eta[450], cluster_phi[450];
        Float_t truthpar_pt[200], truthpar_pz[200], truthpar_eta[200], truthpar_phi[200];
        Int_t truthpar_pid[200];

        tree->SetBranchAddress("zvtx",&zvtx);
        tree->SetBranchAddress("clsmult",&clsmult);
        tree->SetBranchAddress("cluster_e",cluster_e);
        tree->SetBranchAddress("cluster_eta",cluster_eta);
        tree->SetBranchAddress("cluster_phi",cluster_phi);

        tree->SetBranchAddress("truthpar_n",&truthpar_n);
        tree->SetBranchAddress("truthpar_pt",truthpar_pt);
        tree->SetBranchAddress("truthpar_pz",truthpar_pz);
        tree->SetBranchAddress("truthpar_eta",truthpar_eta);
        tree->SetBranchAddress("truthpar_phi",truthpar_phi);
        tree->SetBranchAddress("truthpar_pid",truthpar_pid);

        outFile->cd();

        TH1F* hEP = new TH1F("hEP_em_"+tags[i],
        "Isolated EM E/p ("+tags[i]+");E_{cls}/P_{truth};Entries",100,-2,4);

        Long64_t nEntries = tree->GetEntries();

        for(Long64_t entry=0; entry<nEntries; entry++){

            tree->GetEntry(entry);

            if(fabs(zvtx)>30) continue;

            for(int itr=0; itr<truthpar_n; itr++){

                int pid = truthpar_pid[itr];

                // select electrons, positrons, photons
                if(!(abs(pid)==11 || pid==22)) continue;

                if(truthpar_pt[itr] < 1.0 || fabs(truthpar_eta[itr])>=0.9) continue;

                // ---------- isolation ----------
                bool isIsolated = true;

                for(int jtr=0;jtr<truthpar_n;jtr++){

                    if(itr==jtr) continue;

                    if(truthpar_pt[jtr]>0.2){

                        double dR_iso = sqrt(
                        pow(truthpar_eta[itr]-truthpar_eta[jtr],2) +
                        pow(TVector2::Phi_mpi_pi(truthpar_phi[itr]-truthpar_phi[jtr]),2));

                        if(dR_iso<0.6){
                            isIsolated=false;
                            break;
                        }
                    }
                }

                if(!isIsolated) continue;

                double p_truth = sqrt(
                pow(truthpar_pt[itr],2)+pow(truthpar_pz[itr],2));

                // ---------- matching ----------
                double bestDR=0.1;
                int bestIdx=-1;

                for(int icl=0; icl<clsmult; icl++){

                    double dR_match = sqrt(
                    pow(truthpar_eta[itr]-cluster_eta[icl],2)+
                    pow(TVector2::Phi_mpi_pi(truthpar_phi[itr]-cluster_phi[icl]),2));

                    if(dR_match<bestDR){
                        bestDR=dR_match;
                        bestIdx=icl;
                    }
                }

                if(bestIdx!=-1){

                    double EoverP = cluster_e[bestIdx]/p_truth;
                    hEP->Fill(EoverP);

                }
            }
        }

        // ---------- Draw ----------
        TCanvas* c = new TCanvas("c_em_"+tags[i],tags[i],800,600);

        hEP->Draw("E");

        c->Write();
        hEP->Write();

        f->Close();
        delete f;
    }

    outFile->Write();
    outFile->Close();

    std::cout<<"Successfully saved EM particle plots."<<std::endl;
}
