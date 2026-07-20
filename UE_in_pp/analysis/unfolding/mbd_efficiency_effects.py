#!/usr/bin/env python3

import ROOT
import sys
import os
from array import array

ROOT.gROOT.LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
ROOT.gROOT.ProcessLine("SetsPhenixStyle()")

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

pxbins = [17, 21, 26, 32.5, 40.5, 63.5, 82]
pybins = [0.0, 0.5, 1.08, 1.97, 3.05, 4.68, 6.2, 15.0, 35.0]

def get_hist(file_handle, name):
    h = file_handle.Get(name)
    if not h:
        raise RuntimeError(f"Could not find histogram: {name}")
    return h


def rebin_th2_variable(hist, x_bins, y_bins, new_name):
    if not hist.InheritsFrom("TH2"):
        raise RuntimeError(f"Histogram {hist.GetName()} is not a TH2")
    if hist.GetNbinsX() != len(x_bins) - 1:
        raise RuntimeError(
            f"X-bin mismatch for {hist.GetName()}: source={hist.GetNbinsX()}, requested={len(x_bins) - 1}"
        )
    if hist.GetNbinsY() != len(y_bins) - 1:
        raise RuntimeError(
            f"Y-bin mismatch for {hist.GetName()}: source={hist.GetNbinsY()}, requested={len(y_bins) - 1}"
        )

    rebinned = ROOT.TH2D(
        new_name,
        hist.GetTitle(),
        len(x_bins) - 1,
        array('d', x_bins),
        len(y_bins) - 1,
        array('d', y_bins),
    )
    rebinned.Sumw2()

    # Copy all bins by index (including under/overflow) so only axis binning changes.
    for ix in range(0, hist.GetNbinsX() + 2):
        for iy in range(0, hist.GetNbinsY() + 2):
            rebinned.SetBinContent(ix, iy, hist.GetBinContent(ix, iy))
            rebinned.SetBinError(ix, iy, hist.GetBinError(ix, iy))

    return rebinned


def main(original_file, corrected_file, output_file):

    hist_name = "h_unfold_calib_dijet_reweight_trim_5_3"
    corrected_name = hist_name + "_etEffCorrected"
    output_root = os.path.basename(output_file).split('.')[0];
    print(output_root)

    f_orig = ROOT.TFile.Open(original_file, "READ")
    f_corr = ROOT.TFile.Open(corrected_file, "READ")

    if not f_orig or not f_corr:
        print("Error opening files.")
        sys.exit(1)

    h_orig_raw = get_hist(f_orig, hist_name)
    h_corr_raw = get_hist(f_corr, corrected_name)

    h_orig = rebin_th2_variable(h_orig_raw, pxbins, pybins, "h_orig_rebinned")
    h_corr = rebin_th2_variable(h_corr_raw, pxbins, pybins, "h_corr_rebinned")

    # ---------------------------------------
    # Create Projections and Profiles
    # ---------------------------------------
    projY_orig = h_orig.ProjectionY("projY_orig")
    projY_corr = h_corr.ProjectionY("projY_corr")

    profX_orig = h_orig.ProfileX("profX_orig")
    profX_corr = h_corr.ProfileX("profX_corr")

    # Styling
    projY_orig.SetLineColor(ROOT.kBlack)
    projY_corr.SetLineColor(ROOT.kRed)

    profX_orig.SetLineColor(ROOT.kBlack)
    profX_corr.SetLineColor(ROOT.kRed)

    # ---------------------------------------
    # Create output file
    # ---------------------------------------
    f_out = ROOT.TFile.Open(output_file, "RECREATE")

    # ---------------------------------------
    # Plot ProjectionY comparison
    # ---------------------------------------
    c1 = ROOT.TCanvas("c_projectionY", "", 800, 600)
    #projY_corr.SetTitle("ProjectionY Before/After Efficiency Correction")
    projY_corr.GetXaxis().SetTitle("#Sigma E_{T} [GeV]")
    projY_corr.GetYaxis().SetTitle("Counts")
    projY_corr.Draw("HIST E")
    projY_orig.Draw("HIST E SAME")

    leg1 = ROOT.TLegend(0.45, 0.7, 0.88, 0.92)
    leg1.AddEntry("","#bf{#it{sPHENIX}} Internal","")
    if "efrac" in original_file:
        leg1.AddEntry("","200 GeV p+p Inclusive Jet","")
    else:
        leg1.AddEntry("","200 GeV p+p Dijet","")
    leg1.AddEntry(projY_orig, "Before MBD Eff Correction", "l")
    leg1.AddEntry(projY_corr, "After MBD Eff Correction", "l")
    leg1.Draw()

    c1.Write()
    c1.SaveAs(f"{output_root}_projectionY.png")

    # ---------------------------------------
    # Plot ProfileX comparison
    # ---------------------------------------
    c2 = ROOT.TCanvas("c_profileX", "", 800, 600)
    #profX_corr.SetTitle("ProfileX Before/After Efficiency Correction")
    profX_corr.GetXaxis().SetTitle("Jet p_{T} [GeV]")
    profX_corr.GetYaxis().SetTitle("<#Sigma E_{T}> [GeV]")
    profX_corr.Draw("E")
    profX_orig.Draw("E SAME")

    leg2 = ROOT.TLegend(0.45, 0.7, 0.88, 0.92)
    leg2.AddEntry("","#bf{#it{sPHENIX}} Internal","")
    if "efrac" in original_file:
        leg2.AddEntry("","200 GeV p+p Inclusive Jet","")
    else:
        leg2.AddEntry("","200 GeV p+p Dijet","")
    leg2.AddEntry(profX_orig, "Before MBD Eff Correction", "l")
    leg2.AddEntry(profX_corr, "After MBD Eff Correction", "l")
    leg2.Draw()

    c2.Write()
    c2.SaveAs(f"{output_root}_profileX.png")

    # ---------------------------------------
    # Save histograms too
    # ---------------------------------------
    projY_orig.Write()
    projY_corr.Write()
    profX_orig.Write()
    profX_corr.Write()

    f_out.Close()
    f_orig.Close()
    f_corr.Close()

    print(f"Comparison plots saved to {output_file}")
    print("PDFs saved as projectionY_comparison.pdf and profileX_comparison.pdf")


if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage:")
        print("  python mbd_efficiency_effects.py original.root corrected.root output.root")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])