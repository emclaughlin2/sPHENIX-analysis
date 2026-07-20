"""
Overlay selected TH1 histograms from Herwig MB/Jet trigger samples and
plot trigger/MB ratios in a lower panel.
"""

import argparse
import os
import sys

import ROOT


TRIGS = [0, 5, 12, 20, 30, 40, 50]
XSECS = {
    0: 3.19e10,
    5: 1.8437e08,
    12: 1.0451e06,
    20: 5.2798e03,
    30: 2.0695e03,
    40: 4.4688e02,
    50: 5.044e02,
}
HIST_PATHS = [
    "Jets/R04Jets/Lead_Jets/h_lead_jet_r04_pt",
    "Jets/R04Jets/Lead_Jets/h_lead_jet_r04_eta",
    "Jets/R04Jets/Lead_Jets/h_lead_jet_r04_phi",
    "Jets/R04Jets/Lead_Jets/h_lead_jet_r04_e",
]
MB_FILE = (
    "/sphenix/user/sgross/sphenix_analysis/HerwigToHepMCProduction/"
    "QAHerwigProduction/MB_output/herwig_MB.root"
)
JET_FILE_TEMPLATE = (
    "/sphenix/user/sgross/sphenix_analysis/HerwigToHepMCProduction/"
    "QAHerwigProduction/Jet{trig}_output/herwig_Jet{trig}.root"
)


def load_sphenix_style() -> None:
    """Load and apply the sPHENIX style macro."""
    ROOT.gROOT.LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C")
    ROOT.gROOT.ProcessLine("SetsPhenixStyle()")
    ROOT.gStyle.SetOptStat(0)


def color_for_index(idx: int) -> int:
    """Return a readable line color for each trigger index."""
    palette = [
        ROOT.kBlack,
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


def fetch_histogram(root_file, hist_path: str):
    """Get histogram by path and detach it from file ownership."""
    hist = root_file.Get(hist_path)
    if not hist:
        return None
    hist.SetDirectory(0)
    return hist


def file_path_for_trigger(trig: int) -> str:
    """Return full ROOT filename for a given trigger threshold."""
    if trig == 0:
        return MB_FILE
    return JET_FILE_TEMPLATE.format(trig=trig)


def output_stub_for_hist(hist_path: str) -> str:
    """Convert a ROOT object path into a safe filename stub."""
    return hist_path.replace("/", "_")


def should_use_logy(hist_path: str) -> bool:
    """Use log-y for steeply falling positive distributions."""
    return hist_path.endswith("_pt") or hist_path.endswith("_e") or hist_path.endswith("_eta") or hist_path.endswith("_phi")


def trigger_weight(trig: int) -> float:
    """Return per-trigger weight normalized to trig=0 cross section."""
    if 0 not in XSECS:
        return 1.0
    if trig not in XSECS:
        print(f"[warn] Missing cross section for trig={trig}, using weight=1")
        return 1.0
    return XSECS[trig] / XSECS[0]


def make_overlay_and_ratio_for_hist(
    output_dir: str,
    hist_path: str,
) -> None:
    """Create one two-panel canvas: overlays on top, trigger/MB ratios below."""
    histograms = []
    mb_hist = None

    for i, trig in enumerate(TRIGS):
        full_path = file_path_for_trigger(trig)

        root_file = open_root_file(full_path)
        if root_file is None:
            continue

        hist = fetch_histogram(root_file, hist_path)
        root_file.Close()

        if not hist:
            print(f"[warn] Histogram '{hist_path}' missing in {full_path}")
            continue

        hist.Scale(trigger_weight(trig))

        hist.SetLineColor(color_for_index(i))
        hist.SetLineWidth(2)
        hist.SetMarkerColor(color_for_index(i))
        hist.SetMarkerStyle(20 + (i % 10))
        hist.SetMarkerSize(0.9)
        hist.SetTitle("")

        histograms.append((trig, hist))
        if trig == 0:
            mb_hist = hist

    if not histograms:
        print(f"[warn] No histograms available for hist='{hist_path}'")
        return
    if mb_hist is None:
        print(f"[warn] MB histogram (trig=0) missing for hist='{hist_path}'")
        return

    stub = output_stub_for_hist(hist_path)
    canvas_name = f"c_{stub}"
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 900, 900)
    canvas.Divide(1, 2)

    top_pad = canvas.cd(1)
    top_pad.SetPad(0.0, 0.40, 1.0, 1.0)
    top_pad.SetBottomMargin(0.02)
    top_pad.SetLeftMargin(0.12)
    top_pad.SetRightMargin(0.03)
    if should_use_logy(hist_path):
        top_pad.SetLogy(1)

    ymax = max(hist.GetMaximum() for _, hist in histograms)
    ymin_positive = min(
        (hist.GetBinContent(bin_idx)
         for _, hist in histograms
         for bin_idx in range(1, hist.GetNbinsX() + 1)
         if hist.GetBinContent(bin_idx) > 0.0),
        default=1e-6,
    )

    legend = ROOT.TLegend(0.58, 0.58, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)

    first = True
    for trig, hist in histograms:
        hist.GetYaxis().SetTitleOffset(1.4)
        if should_use_logy(hist_path):
            hist.SetMinimum(max(ymin_positive * 0.5, 1e-10))
            hist.SetMaximum(ymax * 20.0 if ymax > 0 else 1.0)
        else:
            hist.SetMinimum(0.0)
            hist.SetMaximum(ymax * 1.25 if ymax > 0 else 1.0)

        draw_opt = "hist e" if first else "hist e same"
        hist.Draw(draw_opt)
        label = "MB (trig=0)" if trig == 0 else f"Jet{trig}"
        legend.AddEntry(hist, label, "lep")
        first = False

    legend.Draw()

    bottom_pad = canvas.cd(2)
    bottom_pad.SetPad(0.0, 0.0, 1.0, 0.40)
    bottom_pad.SetTopMargin(0.02)
    bottom_pad.SetBottomMargin(0.35)
    bottom_pad.SetLeftMargin(0.12)
    bottom_pad.SetRightMargin(0.03)
    bottom_pad.SetGridy()

    ratio_frame = mb_hist.Clone(f"frame_ratio_{stub}")
    ratio_frame.Reset("ICES")
    ratio_frame.SetTitle("")
    ratio_frame.GetYaxis().SetTitle("Trig N / Trig N-1")
    ratio_frame.GetYaxis().SetNdivisions(505)
    ratio_frame.GetYaxis().SetTitleSize(0.08)
    ratio_frame.GetYaxis().SetTitleOffset(0.6)
    ratio_frame.GetYaxis().SetLabelSize(0.07)
    ratio_frame.GetXaxis().SetTitle(mb_hist.GetXaxis().GetTitle())
    ratio_frame.GetXaxis().SetTitleSize(0.11)
    ratio_frame.GetXaxis().SetTitleOffset(1.1)
    ratio_frame.GetXaxis().SetLabelSize(0.09)
    ratio_frame.SetMinimum(0.0)
    ratio_frame.SetMaximum(10.0)
    ratio_frame.Draw("axis")

    ratio_refs = []
    hist_by_trig = {trig: hist for trig, hist in histograms}
    for idx in range(1, len(TRIGS)):
        trig_curr = TRIGS[idx]
        trig_prev = TRIGS[idx - 1]
        hist_curr = hist_by_trig.get(trig_curr)
        hist_prev = hist_by_trig.get(trig_prev)
        if hist_curr is None or hist_prev is None:
            continue
        ratio = hist_curr.Clone(f"ratio_{stub}_trig{trig_curr}_over_{trig_prev}")
        ratio.SetDirectory(0)
        ratio.Divide(hist_prev)
        ratio.Draw("hist e same")
        ratio_refs.append(ratio)

    one_line = ROOT.TLine(
        mb_hist.GetXaxis().GetXmin(),
        1.0,
        mb_hist.GetXaxis().GetXmax(),
        1.0,
    )
    one_line.SetLineStyle(2)
    one_line.SetLineColor(ROOT.kBlack)
    one_line.Draw()

    canvas.Modified()
    canvas.Update()

    out_base = os.path.join(output_dir, f"{stub}_overlay_ratio")
    canvas.SaveAs(out_base + ".png")
    canvas.SaveAs(out_base + ".pdf")
    print(f"[ok] Wrote: {out_base}.png/.pdf")

    # Keep ratio objects alive until after canvas writeout.
    _ = ratio_refs


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Overlay selected TH1 histograms across trigger thresholds."
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

    for hist_path in HIST_PATHS:
        make_overlay_and_ratio_for_hist(
            output_dir=args.output_dir,
            hist_path=hist_path,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())