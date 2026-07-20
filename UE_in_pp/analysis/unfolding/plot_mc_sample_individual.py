"""
Overlay selected TH1D histograms from ROOT files for each cut setting.

Expected input filenames:
    output_{cut}_bkg_cut_sim_iter_1_jet{trig}.root

with:
    cut  in {"dijet", "none"}
    trig in {10, 15, 20, 30, 50, 70}
"""

import argparse
import os
import sys

import ROOT


CUTS = ["dijet", "efrac"]
TRIGS = [12, 20, 30, 40, 50, 60]
HIST_NAMES = [
    "h_lead_spectra_record",
    "h_lead_truth_spectra_record",
    #"h_et_transverse_record",
    #"h_et_truth_transverse_record",
]


def load_sphenix_style() -> None:
    """Load and apply the sPHENIX style macro."""
    ROOT.gROOT.LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C")
    ROOT.gROOT.ProcessLine("SetsPhenixStyle()")
    ROOT.gStyle.SetOptStat(0)


def color_for_index(idx: int) -> int:
    """Return a readable line color for each trigger index."""
    palette = [
        ROOT.kRed + 1,
        ROOT.kBlue + 1,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kOrange + 7,
        ROOT.kCyan + 2,
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


def make_overlay_for_cut_and_hist(
    input_dir: str,
    output_dir: str,
    trig: str,
    hist_name: str,
) -> None:
    """Create one overlay canvas for one cut and one histogram name."""
    histograms = []

    for i, cut in enumerate(CUTS):
        if cut == "efrac":
            filename = f"output_none_bkg_cut_sim_iter_1_jet{trig}.root"
        else:
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
        hist.SetLineColor(color_for_index(i))
        hist.SetLineWidth(2)
        hist.SetMarkerColor(color_for_index(i))
        hist.SetMarkerStyle(20 + (i % 10))
        hist.SetMarkerSize(0.9)
        hist.GetXaxis().SetRangeUser(14,100)
        hist.GetXaxis().SetTitle("Calibrated Lead p_{T} [GeV]")
        #hist.GetYaxis().SetRangeUser(0.0000000001,20000)

        histograms.append((cut, hist))

    if not histograms:
        print(f"[warn] No histograms available for trig='{trig}', hist='{hist_name}'")
        return

    canvas_name = f"c_{trig}_{hist_name}"
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 900, 700)
    canvas.SetLogy(1)

    legend = ROOT.TLegend(0.56, 0.60, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    first = True
    for cut, hist in histograms:
        draw_opt = "hist" if first else "hist same"
        hist.Draw(draw_opt)
        legend.AddEntry(hist, f"cut={cut}, trig={trig}", "l")
        first = False

    legend.Draw()
    canvas.Modified()
    canvas.Update()

    out_base = os.path.join(output_dir, f"test_{trig}_{hist_name}_overlay")
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

    for trig in TRIGS:
        for hist_name in HIST_NAMES:
            make_overlay_for_cut_and_hist(
                input_dir=args.input_dir,
                output_dir=args.output_dir,
                trig=trig,
                hist_name=hist_name,
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
