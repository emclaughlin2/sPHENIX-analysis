#!/usr/bin/env python3

import ROOT
import sys

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


def get_hist(file_handle, name):
    h = file_handle.Get(name)
    if not h:
        raise RuntimeError(f"Could not find histogram: {name}")
    return h


def main(original_file, corrected_file, output_file):

    hist_name = "h_unfold_calib_dijet_reweight_trim_5_3"
    corrected_name = hist_name + "_etEffCorrected"

    f_orig = ROOT.TFile.Open(original_file, "READ")
    f_corr = ROOT.TFile.Open(corrected_file, "READ")

    if not f_orig or not f_corr:
        print("Error opening files.")
        sys.exit(1)

    h_orig = get_hist(f_orig, hist_name)
    h_corr = get_hist(f_corr, corrected_name)

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
    projY_orig.SetTitle("ProjectionY Before/After Efficiency Correction;#Sigma E_{T} [GeV];")
    projY_corr.Draw("HIST E")
    projY_orig.Draw("HIST E SAME")

    leg1 = ROOT.TLegend(0.6, 0.75, 0.88, 0.88)
    leg1.AddEntry(projY_orig, "Before correction", "l")
    leg1.AddEntry(projY_corr, "After correction", "l")
    leg1.Draw()

    c1.Write()
    c1.SaveAs("mbd_eff_effects_dijet_projectionY.png")

    # ---------------------------------------
    # Plot ProfileX comparison
    # ---------------------------------------
    c2 = ROOT.TCanvas("c_profileX", "", 800, 600)
    profX_orig.SetTitle("ProfileX Before/After Efficiency Correction;Jet p_{T} [GeV];<#Sigma E_{T}> [GeV]")
    profX_corr.Draw("E")
    profX_orig.Draw("E SAME")

    leg2 = ROOT.TLegend(0.6, 0.75, 0.88, 0.88)
    leg2.AddEntry(profX_orig, "Before correction", "l")
    leg2.AddEntry(profX_corr, "After correction", "l")
    leg2.Draw()

    c2.Write()
    c2.SaveAs("mbd_eff_effects_dijet_profileX.png")

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