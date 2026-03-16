#!/usr/bin/env python3

import ROOT
import sys

ROOT.gROOT.SetBatch(True)

def main(input_file, output_file, ratio_file):

    # ----------------------------------------
    # Open input file
    # ----------------------------------------
    f_in = ROOT.TFile.Open(input_file, "READ")
    if not f_in or f_in.IsZombie():
        print(f"Error: could not open {input_file}")
        sys.exit(1)

    # ----------------------------------------
    # Retrieve 2D histograms
    # ----------------------------------------
    h_pass = f_in.Get("h_pass_cut_mbdvertex_truthvertex")
    h_total = f_in.Get("h_total_truthvertex")

    if not h_pass or not h_total:
        print("Error: could not retrieve required histograms.")
        sys.exit(1)

    if not isinstance(h_pass, ROOT.TH2) or not isinstance(h_total, ROOT.TH2):
        print("Error: input histograms are not TH2.")
        sys.exit(1)

    # ----------------------------------------
    # Open scale-factor ratio file
    # ----------------------------------------
    f_ratio = ROOT.TFile.Open(ratio_file, "READ")
    if not f_ratio or f_ratio.IsZombie():
        print(f"Error: could not open {ratio_file}")
        sys.exit(1)

    h_ratio = f_ratio.Get("h_eff_ratio_calo_et_data_over_mc")
    if not h_ratio:
        print("Error: could not retrieve h_eff_ratio_calo_et_data_over_mc.")
        sys.exit(1)

    if not isinstance(h_ratio, ROOT.TH1):
        print("Error: h_eff_ratio_calo_et_data_over_mc is not TH1.")
        sys.exit(1)

    # ----------------------------------------
    # Create Y projections (calo ET axis)
    # ----------------------------------------
    h_pass_et = h_pass.ProjectionY("h_pass_et")
    h_total_et = h_total.ProjectionY("h_total_et")

    # Ensure proper error storage
    h_pass_et.Sumw2()
    h_total_et.Sumw2()

    # ----------------------------------------
    # Create efficiency histogram
    # ----------------------------------------
    h_efficiency_et = h_pass_et.Clone("h_efficiency_et")
    h_efficiency_et.SetTitle(";<#Sigma E_{T}> [GeV];MBD Vertex Efficiency")
    h_efficiency_et.Divide(h_total_et)

    # ----------------------------------------
    # Build efficiency variation histograms
    # scale_up   = ratio_bin + ratio_bin_error
    # scale_down = ratio_bin - ratio_bin_error
    # ----------------------------------------
    h_efficiency_et_var_up = h_efficiency_et.Clone("h_efficiency_et_var_up")
    h_efficiency_et_var_down = h_efficiency_et.Clone("h_efficiency_et_var_down")

    n_bins_eff = h_efficiency_et.GetNbinsX()
    n_bins_ratio = h_ratio.GetNbinsX()
    print(f"n_bins_eff: {n_bins_eff}, n_bins_ratio: {n_bins_ratio}")
    if n_bins_eff != n_bins_ratio:
        print(
            "Warning: bin mismatch between h_efficiency_et "
            f"({n_bins_eff}) and h_eff_ratio_calo_et_data_over_mc ({n_bins_ratio}). "
            "Applying up to min(nbins)."
        )
    n_bins = min(n_bins_eff, n_bins_ratio)

    for ibin in range(1, n_bins + 1):
        eff_val = h_efficiency_et.GetBinContent(ibin)
        eff_err = h_efficiency_et.GetBinError(ibin)
        if ibin < n_bins:
            ratio_val = h_ratio.GetBinContent(ibin+1)
            ratio_err = h_ratio.GetBinError(ibin+1)
        else:
            ratio_val = h_ratio.GetBinContent(ibin)
            ratio_err = h_ratio.GetBinError(ibin)
        scale_up = ratio_val + ratio_err
        scale_down = ratio_val - ratio_err

        h_efficiency_et_var_up.SetBinContent(ibin, eff_val * scale_up)
        h_efficiency_et_var_up.SetBinError(ibin, eff_err * abs(scale_up))
        h_efficiency_et_var_down.SetBinContent(ibin, eff_val * scale_down)
        h_efficiency_et_var_down.SetBinError(ibin, eff_err * abs(scale_down))

    # ----------------------------------------
    # Save output
    # ----------------------------------------
    f_out = ROOT.TFile.Open(output_file, "RECREATE")

    h_efficiency_et.Write()
    h_efficiency_et_var_up.Write()
    h_efficiency_et_var_down.Write()

    f_out.Close()
    f_in.Close()
    f_ratio.Close()

    print(f"Efficiency histogram saved to {output_file}")
    print(f"Added variations from {ratio_file}")


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print("Usage:")
        print("  python make_et_efficiency.py input.root output.root [ratio_scale_file.root]")
        sys.exit(1)

    ratio_file = "efficiency_ratio_scale_factors.root"
    if len(sys.argv) == 4:
        ratio_file = sys.argv[3]

    main(sys.argv[1], sys.argv[2], ratio_file)