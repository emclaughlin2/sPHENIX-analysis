#include <TFile.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TF1.h>
#include <TMath.h>
#include <iostream>
#include <vector>

std::vector<double> get_equal_integral_bins(TF1* f, double xmin, double xmax, int nbins, double precision = 1e-4) {
    double total_integral = f->Integral(xmin, xmax);
    double target_bin_integral = total_integral / nbins;
    std::cout << "Total Integral: " << total_integral << ", Target per bin: " << target_bin_integral << std::endl;

    std::vector<double> bin_edges;
    bin_edges.push_back(xmin);
    double current_x = xmin;

    for (int i = 1; i < nbins; ++i) {
        double low = current_x;
        double high = xmax;
        double target_integral = i * target_bin_integral;

        // Binary search for the next bin edge
        while (high - low > precision) {
            double mid = 0.5 * (low + high);
            double integral = f->Integral(xmin, mid);
            if (integral < target_integral) {
                low = mid;
            } else {
                high = mid;
            }
        }

        current_x = 0.5 * (low + high);
        bin_edges.push_back(current_x);
    }

    bin_edges.push_back(xmax);
    return bin_edges;
}

void get_et_binning() {
    TFile* f2 = TFile::Open("analysis_sim_output/output_sim.root");
    if (!f2 || f2->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    TH1D* h_reco = dynamic_cast<TH1D*>(f2->Get("h_reco_et_transverse_record_dijet"));
    TH1D* h_truth = dynamic_cast<TH1D*>(f2->Get("h_truth_et_transverse_record_dijet"));

    if (!h_reco || !h_truth) {
        std::cerr << "Histograms not found!" << std::endl;
        return;
    }

    h_reco->SetDirectory(0);
    h_truth->SetDirectory(0);
    f2->Close();

    // Create TGraphs
    TGraph* g_reco = new TGraph(h_reco);
    TGraph* g_truth = new TGraph(h_truth);

    // Define TF1 using TGraph::Eval via a static wrapper
    TF1* f_reco = new TF1("f_reco", [g_reco](double* x, double*) {
        return g_reco->Eval(x[0], 0, "S");
    }, 0.005, 20, 0);

    TF1* f_truth = new TF1("f_truth", [g_truth](double* x, double*) {
        return g_truth->Eval(x[0], 0, "S");
    }, 0.005, 20, 0);

    double xmin = 0.05;
    double xmax = 10.0;
    int nbins = 10;

    std::vector<double> bin_edges = get_equal_integral_bins(f_reco, xmin, xmax, nbins);

    for (size_t i = 0; i < bin_edges.size(); ++i) {
        std::cout << "Edge " << i << ": " << bin_edges[i] << std::endl;
    }
}

