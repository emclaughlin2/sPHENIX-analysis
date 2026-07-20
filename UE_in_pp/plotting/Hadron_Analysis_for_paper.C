// refit_hadron.C
// Usage: root -l -b -q refit_hadron.C
//
// Reads hEP_hadron_3s from Hadron_Isolated_Resolution_3Plots_0.4.root,
// performs a two-pass Gaussian fit with tight range, and writes the
// histogram + fit to Hadron_Isolated_Resolution_3Plots_fit.root so it
// can be read directly by draw_one() in the plotting macro.

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include <iostream>

void Hadron_Analysis_for_paper()
{
    // ----------------------------------------------------------------
    // 1. Open input file and retrieve histogram
    // ----------------------------------------------------------------
    const char* infile  = "Hadron_Isolated_Resolution_3Plots_0.4.root";
    const char* outfile = "Hadron_Isolated_Resolution_3Plots_paper_fit.root";
    const char* hist_name = "hEP_hadron_3s";
    const char* fit_name  = "gausFit_3s";   // name expected by draw_one()

    TFile* fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "ERROR: cannot open " << infile << std::endl;
        return;
    }

    TH1F* h = nullptr;
    fin->GetObject(hist_name, h);
    if (!h) {
        std::cerr << "ERROR: histogram " << hist_name
                  << " not found in " << infile << std::endl;
        fin->Close();
        return;
    }

    // Detach from input file
    TH1F* hfit = (TH1F*) h->Clone(hist_name);
    hfit->SetDirectory(nullptr);
    fin->Close();

    // ----------------------------------------------------------------
    // 2. First-pass fit: seed range from prior knowledge
    //    hadron peak ~ 0.44, sigma ~ 0.30  -> mean +/- 1.5 sigma
    // ----------------------------------------------------------------
    const double seed_mean  = 0.442;
    const double seed_sigma = 0.305;

    double range_lo = seed_mean - 1.5 * seed_sigma;   // ~ 0.00
    double range_hi = seed_mean + 1.5 * seed_sigma;   // ~ 0.90
    // Clamp low edge to avoid unphysical region
    if (range_lo < 0.05) range_lo = 0.05;

    TF1* gfit = new TF1(fit_name, "gaus", range_lo, range_hi);
    gfit->SetParameters(hfit->GetMaximum(), seed_mean, seed_sigma);

    std::cout << "\n--- Pass 1: range ["
              << range_lo << ", " << range_hi << "] ---" << std::endl;
    hfit->Fit(gfit, "RQ");   // R = use TF1 range, Q = quiet

    // ----------------------------------------------------------------
    // 3. Second-pass fit: re-center range on fitted mean/sigma
    // ----------------------------------------------------------------
    double mean1  = gfit->GetParameter(1);
    double sigma1 = TMath::Abs(gfit->GetParameter(2));

    double range_lo2 = mean1 - 1.5 * sigma1;
    double range_hi2 = mean1 + 1.5 * sigma1;
    if (range_lo2 < 0.05) range_lo2 = 0.05;

    gfit->SetRange(range_lo2, range_hi2);

    std::cout << "--- Pass 2: range ["
              << range_lo2 << ", " << range_hi2 << "] ---" << std::endl;
    hfit->Fit(gfit, "R");   // verbose second pass so fit result prints

    // ----------------------------------------------------------------
    // 4. Print summary
    // ----------------------------------------------------------------
    double mean_f  = gfit->GetParameter(1);
    double sigma_f = TMath::Abs(gfit->GetParameter(2));
    double chi2ndf = (gfit->GetNDF() > 0)
                   ? gfit->GetChisquare() / gfit->GetNDF() : -1.0;

    std::cout << "\n========== Fit Result ==========" << std::endl;
    std::cout << "  Mean  = " << mean_f  << " +/- " << gfit->GetParError(1) << std::endl;
    std::cout << "  Sigma = " << sigma_f << " +/- " << gfit->GetParError(2) << std::endl;
    std::cout << "  Chi2/NDF = " << chi2ndf << std::endl;
    std::cout << "================================" << std::endl;

    // ----------------------------------------------------------------
    // 5. Write histogram + fit to output file
    //    draw_one() calls f->GetObject(hist_name) and f->GetObject(fit_name)
    //    so both must be top-level objects in the file.
    // ----------------------------------------------------------------
    TFile* fout = TFile::Open(outfile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: cannot create " << outfile << std::endl;
        return;
    }

    hfit->Write(hist_name);   // write under the exact name draw_one() expects
    gfit->Write(fit_name);    // write under the exact name draw_one() expects

    fout->Close();
    std::cout << "\nSaved histogram and fit to: " << outfile << std::endl;
}