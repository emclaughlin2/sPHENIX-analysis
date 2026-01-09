import sys
import re
import ROOT

#!/usr/bin/env python3
# plot_iterations.py
#
# Usage: python plot_iterations.py input.root
#
# Reads TH2D histograms named like
#   h_unfold_calib_dijet_reweight_trim_10_1 ... _20
# projects them to X and Y, and overlays the projections:
# - first iteration: blue
# - middle iterations: black
# - last iteration: red
#
# Requires PyROOT.

# Define the bin edges for the projections
truthptbins = [17, 21, 26, 32.5, 40.5, 63.5, 82]  # 72
truthetbins = [0.0, 0.1, 1.08, 1.97, 3.05, 4.68, 6.2, 15.0, 35.0]
# X-axis range limits for projections and profiles
projx_xmax = 0.82
profx_xmax = 0.82
def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_iterations.py input.root")
        sys.exit(1)

    fname = sys.argv[1]
    f = ROOT.TFile.Open(fname)
    if not f or f.IsZombie():
        print("Failed to open", fname)
        sys.exit(1)

    # collect TH2s matching unfolding pattern
    pattern = re.compile(r'h_unfold_calib_dijet_reweight_trim_10_(\d+)$')
    items = []
    for key in f.GetListOfKeys():
        obj = key.ReadObj()
        if not obj or not obj.InheritsFrom("TH2"):
            continue
        name = obj.GetName()
        m = pattern.search(name)
        if m:
            idx = int(m.group(1))
            items.append((idx, obj))

    if not items:
        print("No unfold histograms found matching pattern.")
        sys.exit(0)

    # sort by iteration index
    items.sort(key=lambda x: x[0])
    n = len(items)
    # wrap ROOT.TCanvas so canvases for projX get log y-scale
    _orig_TCanvas = ROOT.TCanvas
    def _TCanvas_with_logy(*args, **kwargs):
        c = _orig_TCanvas(*args, **kwargs)
        if args and isinstance(args[0], str) and args[0].startswith("c_projX"):
            c.SetLogy(1)
        return c
    ROOT.TCanvas = _TCanvas_with_logy
    projx_list = []
    projy_list = []
    profx_list = [] 
    names = []
    for idx, h2 in items:
        name = h2.GetName()
        names.append(name)
        # create projections and detach from file (SetDirectory(0))
        px = h2.ProjectionX(name + "_projX")
        py = h2.ProjectionY(name + "_projY")
        pfx = h2.ProfileX(name + "_profX", 1, h2.GetNbinsY())
        px.SetDirectory(0)
        py.SetDirectory(0)
        projx_list.append(px)
        projy_list.append(py)
        profx_list.append(pfx)

    # drawing helper
    def draw_overlaid(hist_list, outname_suffix, x_title):
        c = ROOT.TCanvas("c_"+outname_suffix, "", 900, 700)
        ROOT.gStyle.SetOptStat(0)
        # determine global max to set same y-range
        ymax = max(h.GetMaximum() for h in hist_list) * 1.2
        first_idx = 0
        last_idx = len(hist_list) - 1

        drawn = False
        for i, h in enumerate(hist_list):
            if i == first_idx:
                color = ROOT.kBlue
            elif i == last_idx:
                color = ROOT.kRed
            else:
                color = ROOT.kBlack
            
            if x_title == "X":
                h.GetXaxis().SetRangeUser(0, projx_xmax)
            h.SetLineColor(color)
            h.SetLineWidth(2)
            if not drawn:
                h.SetMaximum(ymax)
                h.GetXaxis().SetTitle(x_title)
                h.Draw("HIST")
                drawn = True
            else:
                h.Draw("HIST SAME")

        # legend: representative entries for first, middle (if any), last
        leg = ROOT.TLegend(0.70, 0.70, 0.88, 0.88)
        leg.SetTextSize(0.04)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.AddEntry(hist_list[first_idx], f"Iter 1", "l")
        if len(hist_list) > 2:
            # pick an interior histogram as representative
            mid_idx = 1
            leg.AddEntry(hist_list[mid_idx], "Iter 2-19", "l")
        if len(hist_list) > 1:
            leg.AddEntry(hist_list[last_idx], f"Iter 20", "l")
        leg.Draw()

        outpng = f"{fname}_{outname_suffix}.png"
        outpdf = f"{fname}_{outname_suffix}.pdf"
        c.SaveAs(outpng)
        c.SaveAs(outpdf)
        print("Wrote:", outpng, outpdf)

    # draw projections
    draw_overlaid(projx_list, "projX_overlaid", "X")
    draw_overlaid(projy_list, "projY_overlaid", "Y")
    draw_overlaid(profx_list, "profX_overlaid", "X")

    f.Close()

if __name__ == "__main__":
    main()