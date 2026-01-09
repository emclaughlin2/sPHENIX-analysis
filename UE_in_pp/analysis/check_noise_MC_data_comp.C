// check_noise_MC_data_comp.C
// Macro to compare histograms between data and run21_jet10 outputs.
// Produces overlay + ratio plots saved to a folder.

#include <TFile.h>
#include <TKey.h>
#include <TList.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <set>
#include <vector>

// Helper: compare two 1D histograms by overlaying them and drawing a ratio pad.
static void DrawOverlayAndRatio(TH1* hdata, TH1* hmc, const TString &outdir) {
    if (!hdata || !hmc) return;
    // Create clones so axis/attributes don't interfere
    TH1 *d = (TH1*)hdata->Clone(TString(hdata->GetName())+"_clone_for_plot");
    TH1 *m = (TH1*)hmc->Clone(TString(hmc->GetName())+"_clone_for_plot");

    d->SetDirectory(0);
    m->SetDirectory(0);

    d->SetLineColor(kBlack);
    d->SetMarkerColor(kBlack);
    d->SetMarkerStyle(20);
    m->SetLineColor(kRed);
    m->SetMarkerColor(kRed);
    m->SetMarkerStyle(24);

    // Determine histogram name (use data histogram name)
    TString hname = d->GetName();

    TString cname = TString("c_") + d->GetName();
    TCanvas *c = new TCanvas(cname, cname, 800, 800);
    c->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
    pad1->SetBottomMargin(0.02);
    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.35);
    pad1->Draw();
    pad2->Draw();

    pad1->cd();
    // Configure X range and log scale based on histogram type/name
    bool setLogY = false;
    if (hname.BeginsWith("h_ntopo")) {
        // ntopo: x axis 0 -- 30
        d->GetXaxis()->SetRangeUser(0, 30);
        m->GetXaxis()->SetRangeUser(0, 30);
    } else if (hname.Contains("h_sume_topo")) {
        // sume_topo: x axis -3 -- 15 and log y scale
        d->Rebin(8);
        m->Rebin(8);
        d->GetXaxis()->SetRangeUser(-3, 15);
        m->GetXaxis()->SetRangeUser(-3, 15);
        setLogY = true;
    } else if (hname.Contains("h_topo")) {
        // topo: x axis -3 -- 10 and log y scale
        d->Rebin(4);
        m->Rebin(4);
        d->GetXaxis()->SetRangeUser(-3, 10);
        m->GetXaxis()->SetRangeUser(-3, 10);
        setLogY = true;
    } else if (hname.Contains("h_et_spectra")) {
        // et_spectra: x axis -10 -- 80 and log y scale
        d->Rebin(2);
        m->Rebin(2);
        d->GetXaxis()->SetRangeUser(-3, 10);
        m->GetXaxis()->SetRangeUser(-3, 10);
        setLogY = true;
    }
    if (setLogY) pad1->SetLogy(1); else pad1->SetLogy(0);

    // Normalize both histograms to unit integral (self-normalize)
    double id = d->Integral();
    if (id > 0) d->Scale(1.0 / id);
    double im = m->Integral();
    if (im > 0) m->Scale(1.0 / im);

    d->SetStats(0);
    d->Draw("HIST E");
    m->Draw("HIST E SAME");

    TLegend *leg = new TLegend(0.62,0.75,0.88,0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(d, "Data: ana468 0mrad + 1.5mrad", "lep");
    leg->AddEntry(m, "MC: run21 jet10", "lep");
    leg->Draw();

    pad2->cd();
    // ratio = data / mc
    TH1 *r = (TH1*)d->Clone(TString(d->GetName())+"_ratio");
    r->SetDirectory(0);
    r->Divide(m);
    r->SetTitle("");
    r->GetYaxis()->SetTitle("Data / MC");
    r->GetXaxis()->SetTitle(d->GetXaxis()->GetTitle());
    r->GetYaxis()->SetNdivisions(505);
    r->GetYaxis()->SetRangeUser(0, 2);
    r->SetMarkerStyle(20);
    r->Draw("E");

    // Ensure outdir exists
    gSystem->mkdir(outdir, true);
    TString outname = outdir + "/" + d->GetName() + ".png";
    c->SaveAs(outname);

    // cleanup
    delete c;
    delete d;
    delete m;
    delete r;
}

// Find common histogram names in two files that match any of the prefixes in patterns.
static std::vector<TString> FindCommonHistNames(TFile *f1, TFile *f2, const std::vector<TString> &patterns) {
    std::set<TString> names1, names2;

    // Iterate through the list of keys without using TIter (compatibility with older ROOT)
    TList *list1 = f1->GetListOfKeys();
    if (list1) {
        for (int i = 0; i < list1->GetEntries(); ++i) {
            TKey *key = dynamic_cast<TKey*>(list1->At(i));
            if (!key) continue;
            TString name = key->GetName();
            for (auto &p: patterns) {
                if (name.BeginsWith(p)) { names1.insert(name); break; }
            }
        }
    }

    TList *list2 = f2->GetListOfKeys();
    if (list2) {
        for (int i = 0; i < list2->GetEntries(); ++i) {
            TKey *key = dynamic_cast<TKey*>(list2->At(i));
            if (!key) continue;
            TString name = key->GetName();
            for (auto &p: patterns) {
                if (name.BeginsWith(p)) { names2.insert(name); break; }
            }
        }
    }

    std::vector<TString> common;
    for (auto &n: names1) if (names2.find(n) != names2.end()) common.push_back(n);
    return common;
}

void check_noise_MC_data_comp()
{
    const char* sphenixStylePath = "/sphenix/u/egm2153/spring_2023/sPhenixStyle.C";
    // Attempt to load and apply; if file not present ROOT will print an error but code continues
    gROOT->ProcessLine(Form(".L %s+", sphenixStylePath));
    // If the function exists, call it
    // Note: SetsPhenixStyle() is defined in the style file
    gROOT->ProcessLine("SetsPhenixStyle();");

    // Default filenames (edit if your files are located elsewhere)
    const char *dataFile = "dijet_analysis_output/noise_MC_data_comp_data_output.root";
    const char *mcFile   = "dijet_analysis_output/noise_MC_data_comp_run21_jet10_output.root";
    const char *dataTowersFile = "dijet_analysis_output/noise_MC_data_comp_data_towers_output.root";
    const char *mcTowersFile   = "dijet_analysis_output/noise_MC_data_comp_run21_jet10_towers_output.root";

    TFile *fd = TFile::Open(dataFile, "READ");
    TFile *fm = TFile::Open(mcFile, "READ");
    if (!fd || fd->IsZombie()) { std::cerr << "Error opening data file: " << dataFile << std::endl; return; }
    if (!fm || fm->IsZombie()) { std::cerr << "Error opening mc file: " << mcFile << std::endl; return; }

    TString outdir = "check_noise_MC_data_comp_plots";
    gSystem->mkdir(outdir, true);

    // Patterns to search for in the first pair of files
    std::vector<TString> patterns = {"h_ntopo", "h_topo", "h_sume_topo"};

    std::vector<TString> common = FindCommonHistNames(fd, fm, patterns);
    if (common.empty()) {
        std::cout << "No matching histograms found between " << dataFile << " and " << mcFile << " with patterns\n";
    } else {
        std::cout << "Found " << common.size() << " common histograms to compare.\n";
        for (auto &name: common) {
            TObject *od = fd->Get(name);
            TObject *om = fm->Get(name);
            TH1 *hd = dynamic_cast<TH1*>(od);
            TH1 *hm = dynamic_cast<TH1*>(om);
            if (!hd || !hm) {
                std::cout << "Skipping non-TH1 or missing for: " << name << "\n";
                continue;
            }
            DrawOverlayAndRatio(hd, hm, outdir);
        }
    }

    // Now compare the towers ET spectra (single histogram name expected)
    TFile *fdt = TFile::Open(dataTowersFile, "READ");
    TFile *fmt = TFile::Open(mcTowersFile, "READ");
    if (!fdt || fdt->IsZombie()) { std::cerr << "Error opening data towers file: " << dataTowersFile << std::endl; }
    if (!fmt || fmt->IsZombie()) { std::cerr << "Error opening mc towers file: " << mcTowersFile << std::endl; }

    if (fdt && fmt && !fdt->IsZombie() && !fmt->IsZombie()) {
        TH1 *hd_et = dynamic_cast<TH1*>(fdt->Get("h_et_spectra"));
        TH1 *hm_et = dynamic_cast<TH1*>(fmt->Get("h_et_spectra"));
        if (!hd_et || !hm_et) {
            std::cerr << "Could not find h_et_spectra in one of the tower files.\n";
        } else {
            DrawOverlayAndRatio(hd_et, hm_et, outdir);
        }
    }

    // Close files
    fd->Close(); fm->Close();
    if (fdt) fdt->Close();
    if (fmt) fmt->Close();

    std::cout << "Plots written to: " << outdir << std::endl;
}
