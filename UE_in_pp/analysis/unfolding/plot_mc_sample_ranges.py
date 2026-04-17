"""
Overlay selected TH1D histograms from ROOT files for each cut setting.

Expected input filenames:
    output_{cut}_bkg_cut_sim_iter_1_jet{trig}.root

with:
    cut  in {"dijet", "none"}
    trig in {10, 15, 20, 30, 50, 70}
"""

# 1.4903e-6 - Jet 12 cross section -> scale 1.4903
# 3.997e-6 - Jet 10 cross section -> scale 3.997

import argparse
import os
import sys

import ROOT


CUTS = ["dijet", "none"]
TRIGS = [12, 20, 30, 40, 50, 60]
#TRIGS = [10,20,30,50]
SCALE = 1.4903 * 0.81608598
#SCALE = 3.997
HIST_NAMES = [
    "h_lead_spectra_record",
    "h_lead_truth_spectra_record",
    #"h_jes_qa",
    "h_et_transverse_record",
    "h_et_truth_transverse_record",
]


def load_sphenix_style() -> None:
    """Load and apply the sPHENIX style macro."""
    ROOT.gROOT.LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C")
    ROOT.gROOT.ProcessLine("SetsPhenixStyle()")
    ROOT.gStyle.SetOptStat(0)


def color_for_index(idx: int) -> int:
    """Return a readable line color for each trigger index."""
    palette = [
        ROOT.kCyan + 2,
        ROOT.kRed + 1,
        ROOT.kBlue + 1,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kOrange + 7,
    ]
    return palette[idx % len(palette)]


def open_root_file(path: str):
    """Open a ROOT file and return the TFile object or None."""
    if not os.path.exists(path):
        print(f"[warn] Missing file: {path}")
        return None

    root_file = ROOT.TFile.Open(path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"[warn] Could not open ROOT file: {path}")
        return None
    return root_file


def fetch_histogram(root_file, hist_name: str):
    """Get histogram by name and detach it from file ownership."""
    hist = root_file.Get(hist_name)
    if not hist:
        return None
    hist.SetDirectory(0)
    return hist


def fit_jes_histogram(hist):
    """
    Fit JES projection histogram with a Gaussian around 1.0.
    Returns (mean, sigma) or (None, None) if fit fails.
    """
    fit_min = 0.6
    fit_max = 1.4
    fit = ROOT.TF1(f"f_gaus_{hist.GetName()}", "gaus", fit_min, fit_max)
    fit.SetParameters(hist.GetMaximum(), 1.0, 0.1)
    fit.SetParLimits(1, 0.8, 1.2)
    fit.SetParLimits(2, 0.01, 0.5)

    fit_status = int(hist.Fit(fit, "RQ0"))
    if fit_status != 0:
        return None, None

    mean = fit.GetParameter(1)
    sigma = fit.GetParameter(2)
    return mean, sigma


def make_overlay_for_cut_and_hist(
    input_dir: str,
    output_dir: str,
    cut: str,
    hist_name: str,
) -> None:
    """Create one overlay canvas for one cut and one histogram name."""
    histograms = []
    total_hist = None

    for i, trig in enumerate(TRIGS):
        filename = f"output_{cut}_bkg_cut_sim_iter_1_jet{trig}.root"
        full_path = os.path.join(input_dir, filename)

        root_file = open_root_file(full_path)
        if root_file is None:
            continue

        hist = fetch_histogram(root_file, hist_name)
        root_file.Close()

        if not hist:
            print(f"[warn] Histogram '{hist_name}' missing in {filename}")
            continue
        hist.Rebin(8)
        hist.Scale(1.0/8.0)
        hist.Scale(SCALE)
        hist.SetLineColor(color_for_index(i))
        hist.SetLineWidth(2)
        hist.SetMarkerColor(color_for_index(i))
        hist.SetMarkerStyle(20 + (i % 10))
        hist.SetMarkerSize(0.9)
        hist.GetYaxis().SetRangeUser(0.0000001,20000)

        fit_mean = None
        fit_sigma = None
        histograms.append((trig, hist, fit_mean, fit_sigma))

        if total_hist is None:
            total_hist = hist.Clone(f"total_{cut}_{hist_name}")
            total_hist.SetDirectory(0)
        else:
            total_hist.Add(hist)

    if not histograms:
        print(f"[warn] No histograms available for cut='{cut}', hist='{hist_name}'")
        return

    canvas_name = f"c_{cut}_{hist_name}"
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 900, 700)
    canvas.SetLogy(1)

    legend = ROOT.TLegend(0.5, 0.60, 0.92, 0.92)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.025)

    first = True
    for trig, hist, fit_mean, fit_sigma in histograms:
        draw_opt = "hist" if first else "hist same"
        hist.Draw(draw_opt)
        legend_label = f"cut={cut}, trig={trig}"

        legend.AddEntry(hist, legend_label, "l")
        first = False

    if total_hist:
        total_hist.SetLineColor(ROOT.kBlack)
        total_hist.SetLineWidth(3)
        total_hist.SetMarkerStyle(0)
        total_hist.Draw("hist same")
        legend.AddEntry(total_hist, f"cut={cut}, trig=all (sum)", "l")

    legend.Draw()
    canvas.Modified()
    canvas.Update()

    out_base = os.path.join(output_dir, f"test_{cut}_{hist_name}_overlay")
    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".pdf")
    print(f"[ok] Wrote: {out_base}.png/.pdf")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Overlay selected TH1D histograms across trigger thresholds."
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        default=".",
        help="Directory containing ROOT files (default: current directory)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="pyroot_overlays",
        help="Directory to store output plots (default: pyroot_overlays)",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    ROOT.gROOT.SetBatch(True)
    load_sphenix_style()

    for cut in CUTS:
        for hist_name in HIST_NAMES:
            make_overlay_for_cut_and_hist(
                input_dir=args.input_dir,
                output_dir=args.output_dir,
                cut=cut,
                hist_name=hist_name,
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
