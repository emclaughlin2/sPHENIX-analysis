#!/usr/bin/env python3
"""
plot_star_comparison.py
-----------------------
Overlay sPHENIX unfolded dNch/detadphi with STAR transverse-region measurement.

Left  y-axis  (red)   : STAR Transverse Region dNch/detadphi (pT > 0.5 GeV/c)
Right y-axis  (blue)  : sPHENIX result (unfold_hist + syst band)

Usage:
  python3 plot_star_comparison.py
  root -l -b -q plot_star_comparison.py   (if run as a ROOT macro via PyROOT)
"""

import ROOT
import array

ROOT.gROOT.SetBatch(False)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPadTickX(1)
#ROOT.gStyle.SetPadTickY(1)

# ============================================================
# 1.  Load sPHENIX objects from file
# ============================================================
FNAME = ("sphenix_primary_run28_output_files/"
         "plot_result_hist_output_efrac_bkg_cut_run28_iter_3_1000toys.root")

f = ROOT.TFile.Open(FNAME)
if not f or f.IsZombie():
    raise RuntimeError(f"Cannot open file: {FNAME}")

unfold_hist  = f.Get("unfold_hist_calib_dijet_reweight_trim_10_2_etEffCorrected_2")
graph        = f.Get("Graph")          # TGraphAsymmErrors (total syst band)
truth_pythia = f.Get("truth_hist_pythia")
truth_herwig = f.Get("truth_hist_herwig")

# Detach from file so they survive after Close()
for obj in [unfold_hist, truth_pythia, truth_herwig]:
    obj.SetDirectory(0)
f.Close()

# ============================================================
# 2.  STAR data  (Transverse Region, dNch/detadphi, pT > 0.5)
# ============================================================
# fmt: jet pT | value | stat+ | stat- | sys+ | sys-
_jetpt   = [6.0,    8.0,    10.0,   13.0,   17.5,   22.5,   30.0,   40.0  ]
_val     = [0.3355, 0.3368, 0.3261, 0.3073, 0.2547, 0.237,  0.2446, 0.241 ]
_stat_up = [0.0001, 0.0001, 0.0001, 0.0001, 0.0002, 0.0005, 0.001,  0.004 ]
_stat_dn = [0.0001, 0.0001, 0.0001, 0.0001, 0.0002, 0.0005, 0.001,  0.004 ]
_sys_up  = [0.06,   0.04,   0.05,   0.05,   0.05,   0.04,   0.032,  0.032 ]
_sys_dn  = [0.028,  0.04,   0.03,   0.02,   0.018,  0.014,  0.013,  0.015 ]
# Horizontal half-widths (representative bin widths)
_exl     = [1.0,    1.0,    1.5,    2.0,    2.5,    2.5,    5.0,    5.0   ]
_exh     = [1.0,    1.0,    1.5,    2.0,    2.5,    2.5,    5.0,    5.0   ]

n = len(_jetpt)
d = lambda lst: array.array('d', lst)
x, y = d(_jetpt), d(_val)
exl, exh = d(_exl), d(_exh)

star_stat = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, d(_stat_dn), d(_stat_up))
star_sys  = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, d(_sys_dn),  d(_sys_up))

# ============================================================
# 3.  Styling
# ============================================================
SPH_COLOR  = ROOT.kAzure - 2     # blue family for sPHENIX
STAR_COLOR = ROOT.kRed + 1       # red family for STAR

# --- sPHENIX unfolded central values (stat. error bars) ---
unfold_hist.SetLineColor(SPH_COLOR)
unfold_hist.SetLineWidth(2)
unfold_hist.SetMarkerColor(SPH_COLOR)

# --- sPHENIX total systematic band (Graph) ---
graph.SetFillColorAlpha(SPH_COLOR, 0.40)
graph.SetFillStyle(1001)
graph.SetLineWidth(2)
graph.SetMarkerStyle(20)
graph.SetMarkerSize(1)
graph.SetMarkerColor(SPH_COLOR)
graph.SetLineColor(SPH_COLOR)

# --- PYTHIA8 truth ---
truth_pythia.SetLineColor(ROOT.kBlack)
truth_pythia.SetLineWidth(2)
truth_pythia.SetLineStyle(1)
truth_pythia.SetFillStyle(0)

# --- Herwig7 truth ---
truth_herwig.SetLineColor(ROOT.kGray + 1)
truth_herwig.SetLineWidth(2)
truth_herwig.SetLineStyle(7)
truth_herwig.SetFillStyle(0)

# --- STAR statistical error bars ---
star_stat.SetMarkerColor(STAR_COLOR)
star_stat.SetLineColor(STAR_COLOR)
star_stat.SetLineWidth(2)
star_stat.SetMarkerStyle(20)
star_stat.SetMarkerSize(1)

# --- STAR systematic band ---
star_sys.SetFillColorAlpha(STAR_COLOR, 0.30)
star_sys.SetFillStyle(1001)
star_sys.SetLineWidth(2)
star_sys.SetMarkerStyle(20)
star_sys.SetMarkerSize(1)
star_sys.SetMarkerColor(STAR_COLOR)
star_sys.SetLineColor(STAR_COLOR)

# ============================================================
# 4.  Canvas geometry
# ============================================================
LM, RM, BM, TM = 0.14, 0.16, 0.13, 0.08   # margins (left/right/bottom/top)

# Shared x-axis range (covers both experiments)
X_MIN, X_MAX = 5.0, 63.5

# Left  axis range  (sPHENIX)
SPH_YMIN, SPH_YMAX = 0.0,  0.9

# Right axis range  (STAR)
STAR_YMIN, STAR_YMAX = 0.0, 0.5

c = ROOT.TCanvas("c_sphenix_star", "sPHENIX vs STAR UE", 960, 740)

# ============================================================
# 5.  Pad 1 — STAR data (defines left y-axis)
# ============================================================
pad1 = ROOT.TPad("pad1", "STAR", 0, 0, 1, 1)
pad1.SetLeftMargin(LM)
pad1.SetRightMargin(RM)
pad1.SetBottomMargin(BM)
pad1.SetTopMargin(TM)
pad1.Draw()
pad1.cd()

# DrawFrame sets up the STAR coordinate system and renders the left y-axis
hframe = pad1.DrawFrame(X_MIN, STAR_YMIN, X_MAX, STAR_YMAX)
hframe.GetXaxis().SetTitle("p_{T,lead} [GeV]")
hframe.GetYaxis().SetTitle("<#deltaN_{ch}/#delta#eta#delta#phi> [GeV]")

hframe.GetXaxis().SetTitleSize(0.046)
hframe.GetXaxis().SetLabelSize(0.040)
hframe.GetXaxis().SetTitleOffset(1.05)

hframe.GetYaxis().SetTitleSize(0.038)
hframe.GetYaxis().SetLabelSize(0.040)
hframe.GetYaxis().SetTitleOffset(1.55)

# STAR: syst band first, then stat error bars on top
star_sys.Draw("E2 SAME")
star_stat.Draw("EPZ SAME")

pad1.Update()

# ============================================================
# 6.  Pad 2 — sPHENIX data (transparent overlay, right y-axis)
# ============================================================
c.cd()
pad2 = ROOT.TPad("pad2", "sPHENIX", 0, 0, 1, 1)
pad2.SetLeftMargin(LM)
pad2.SetRightMargin(RM)
pad2.SetBottomMargin(BM)
pad2.SetTopMargin(TM)
pad2.SetFillStyle(4000)         # fully transparent background
pad2.SetFrameFillStyle(4000)    # transparent frame
pad2.Draw()
pad2.cd()

# DrawFrame establishes the sPHENIX y-axis coordinate system.
# Suppress all axes — x is already drawn by pad1, right axis drawn via TGaxis.
hframe2 = pad2.DrawFrame(X_MIN, SPH_YMIN, X_MAX, SPH_YMAX)
hframe2.GetXaxis().SetLabelSize(0)
hframe2.GetXaxis().SetTickLength(0)
hframe2.GetYaxis().SetLabelSize(0)
hframe2.GetYaxis().SetTickLength(0)

graph.GetY()[0] = -10;
unfold_hist.SetBinContent(1,-10);

# Drawing order: syst band → data points on top
graph.Draw("E2 SAME")           # total syst band
unfold_hist.Draw("EL X0 SAME")  # stat error bars on top

# --- Right y-axis via TGaxis ---
# "+L" flag draws ticks and labels on the right side
pad2.Update()
right_ax = ROOT.TGaxis(X_MAX, SPH_YMIN, X_MAX, SPH_YMAX,
                       SPH_YMIN, SPH_YMAX, 505, "+L")
right_ax.SetTitle("<#SigmaE_{T}/#delta#eta#delta#phi> [GeV]")
right_ax.SetTitleOffset(1.55)
right_ax.SetTitleSize(0.038)
right_ax.SetLabelFont(42)
right_ax.SetTitleFont(42)
right_ax.SetLabelSize(0.040)
right_ax.Draw()

# ============================================================
# 7.  Legend  (drawn on pad2 so it sits on top of everything)
# ============================================================
pad2.cd()

leg = ROOT.TLegend(0.15, 0.15, 0.6, 0.28)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.033)
leg.SetTextFont(42)

#leg.AddEntry(unfold_hist,  "sPHENIX (stat.)",        "lep")
leg.AddEntry(graph,        "sPHENIX <#kern[0.01]{#SigmaE_{T}/#delta#eta#delta#phi> #kern[0.01]{p_{T}^{ch/neut} > 0.5/0.2 GeV}}",  "fpl")
#leg.AddEntry(truth_pythia, "PYTHIA8 (truth)",         "l")
#leg.AddEntry(truth_herwig, "Herwig7 (truth)",         "l")
#leg.AddEntry(star_stat,    "STAR (stat.)",            "lep")
leg.AddEntry(star_sys,     "STAR <#kern[0.01]{#deltaN_{ch}/#delta#eta#delta#phi> #kern[0.01]{p_{T}^{ch} > 0.5 GeV}}",            "fpl")
leg.Draw()

# ============================================================
# 8.  Labels
# ============================================================
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextFont(42)
tex.SetTextSize(0.038)
tex.DrawLatex(LM + 0.34, 1.0 - TM - 0.07, "#bf{#it{sPHENIX}} Internal")
tex.SetTextFont(42)
tex.SetTextSize(0.033)
tex.DrawLatex(LM + 0.34, 1.0 - TM - 0.12,
              "#sqrt{s} = 200 GeV p+p Inclusive Jet")

# ============================================================
# 9.  Finalize and save
# ============================================================
pad2.Update()
c.Update()
c.SaveAs("sphenix_primary_run28_output_files/sphenix_star_ue_comparison.pdf")
c.SaveAs("sphenix_primary_run28_output_files/sphenix_star_ue_comparison.png")
print("Saved: sphenix_primary_run28_output_files/sphenix_star_ue_comparison.pdf  and  .png")
