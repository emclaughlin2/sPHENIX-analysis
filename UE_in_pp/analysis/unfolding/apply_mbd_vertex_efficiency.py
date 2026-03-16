#!/usr/bin/env python3

import ROOT
import sys
ROOT.gROOT.SetBatch(True)

def apply_efficiency_to_hist2d(hist2d, eff1d):

    hist2d.Sumw2()
    nxbins = hist2d.GetNbinsX()
    nybins = hist2d.GetNbinsY()
    print(f"nxbins: {nxbins}, nybins: {nybins}, eff1d.GetNbinsX(): {eff1d.GetNbinsX()}")

    # Sanity check: matching Y binning
    if nybins != eff1d.GetNbinsX():
        raise RuntimeError(
            f"Y-axis bins ({nybins}) do not match efficiency bins "
            f"({eff1d.GetNbinsX()})"
        )
    for j in range(1, nybins + 1):
        eff = eff1d.GetBinContent(j)
        eff_err = eff1d.GetBinError(j)
        for i in range(1, nxbins + 1):
            content = hist2d.GetBinContent(i, j)
            error = hist2d.GetBinError(i, j)
            if eff > 0:
                corrected = content / eff
                # Proper error propagation
                if content > 0:
                    rel_err2 = (error / content) ** 2 + (eff_err / eff) ** 2
                    corrected_err = corrected * (rel_err2 ** 0.5)
                else:
                    corrected_err = 0.0
                hist2d.SetBinContent(i, j, corrected)
                hist2d.SetBinError(i, j, corrected_err)

            else:
                # Protect against zero efficiency
                hist2d.SetBinContent(i, j, 0.0)
                hist2d.SetBinError(i, j, 0.0)

    return hist2d


def main(eff_file, result_file, output_file):

    f_eff = ROOT.TFile.Open(eff_file, "READ")
    f_res = ROOT.TFile.Open(result_file, "READ")

    if not f_eff or not f_res:
        print("Error opening input files.")
        sys.exit(1)

    h_eff = f_eff.Get("h_efficiency_et")
    h_eff_up = f_eff.Get("h_efficiency_et_var_up")
    h_eff_down = f_eff.Get("h_efficiency_et_var_down")

    if not h_eff:
        print("Error: could not find h_efficiency_et")
        sys.exit(1)

    f_out = ROOT.TFile.Open(output_file, "RECREATE")

    for key in f_res.GetListOfKeys():

        obj = key.ReadObj()

        # --------------------------------------
        # Apply correction ONLY to TH2D
        # --------------------------------------
        if isinstance(obj, ROOT.TH2) and "h_unfold" in obj.GetName():

            print(f"Correcting unfold TH2D: {obj.GetName()}")

            hist_corrected = obj.Clone(obj.GetName() + "_etEffCorrected")
            hist_corrected.SetDirectory(0)

            hist_corrected = apply_efficiency_to_hist2d(
                hist_corrected, h_eff
            )

            hist_corrected_var_up = None
            hist_corrected_var_down = None
            if "h_unfold_calib_dijet_reweight_trim_" in obj.GetName():
                hist_corrected_var_up = hist_corrected.Clone(obj.GetName() + "_var_up_etEffCorrected")
                hist_corrected_var_up = apply_efficiency_to_hist2d(
                    hist_corrected_var_up, h_eff_up
                )
                hist_corrected_var_down = hist_corrected.Clone(obj.GetName() + "_var_down_etEffCorrected")
                hist_corrected_var_down = apply_efficiency_to_hist2d(
                    hist_corrected_var_down, h_eff_down
                )   

            f_out.cd()
            hist_corrected.Write()   
            if hist_corrected_var_up:
                hist_corrected_var_up.Write()
            if hist_corrected_var_down:
                hist_corrected_var_down.Write()

        # --------------------------------------
        # Copy TH1D unchanged
        # --------------------------------------
        else:

            print(f"Copying TH1D/TH2D unchanged: {obj.GetName()}")

            f_out.cd()
            obj.Write()

    f_out.Close()
    f_res.Close()
    f_eff.Close()

    print(f"\nDone. Output written to {output_file}")


if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage:")
        print("  python apply_et_efficiency.py efficiency.root results.root output.root")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])