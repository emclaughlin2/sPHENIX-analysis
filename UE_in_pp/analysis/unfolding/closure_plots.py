import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency, TH1D, TH2D
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb
from array import array
gROOT.LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
gROOT.ProcessLine("SetsPhenixStyle()")

rgb = [[230, 25, 75], [60, 180, 75], [255, 225, 25], [0, 130, 200], [245, 130, 48], [145, 30, 180], [70, 240, 240], [240, 50, 230], [210, 245, 60], [250, 190, 212], [0, 128, 128], [220, 190, 255], [170, 110, 40], [128, 128, 128], [128, 0, 0], [0, 0, 0], [128, 128, 0], [255, 215, 180], [0, 0, 128], [34, 139, 34]]
colors = [TColor.GetColor(rgb[i][0],rgb[i][1],rgb[i][2]) for i in range(len(rgb))]
full_leg_tags = ["Unfolded spectrum, iter 1", "Unfolded spectrum, iter 2", "Unfolded spectrum, iter 1", "Unfolded spectrum, iter 2"]
half_leg_tags = [f"Unfolded spectrum, iter {i+1}" for i in range(20)]

def draw_unfolded_spectra(truth, measure, unfold, jet, leg_tags, max_iter, x_range = (14,82), output_name=None):
    canvas = ROOT.TCanvas("canvas", "", 600, 800)

    # Top pad for spectra
    pad1 = ROOT.TPad("pad1", "", 0, 0.5, 1, 1.0)
    pad1.SetBottomMargin(0.02)
    pad1.Draw()
    pad1.cd()
    if jet:
        pad1.SetLogy(1)
    else:
        pad1.SetLogx(1)

    # Style histograms
    truth.SetStats(0)
    measure.SetStats(0)
    truth.SetLineColor(2)
    truth.SetMarkerColor(2)
    measure.SetLineColor(1)
    measure.SetMarkerColor(1)
    truth.GetXaxis().SetLabelSize(0)
    measure.GetXaxis().SetLabelSize(0)

    for i in range(max_iter):
        unfold[i].SetStats(0)
        unfold[i].SetLineColor(colors[i + 1])
        unfold[i].SetMarkerColor(colors[i + 1])
        unfold[i].GetXaxis().SetLabelSize(0)

    truth.Draw()
    measure.Draw("same")
    for i in range(max_iter):
        unfold[i].Draw("same")

    leg = ROOT.TLegend(.57, .65, .9, .9)
    leg.AddEntry(truth, "Truth spectrum", "lp")
    leg.AddEntry(measure, "Measured spectrum", "lp")
    for i in range(max_iter):
        leg.AddEntry(unfold[i], leg_tags[i], "lp")
    leg.Draw()
    leg.SetTextSize(0.04)

    # Bottom pad for ratio plots
    canvas.cd()
    pad2 = ROOT.TPad("pad2", "", 0, 0.0, 1, 0.5)
    pad2.SetTopMargin(0.02)
    pad2.SetBottomMargin(0.2)
    pad2.Draw()
    pad2.cd()
    if not jet:
        pad2.SetLogx(1)

    ratios = []
    for i in range(max_iter):
        ratio = unfold[i].Clone(f"ratio{i}")
        ratio.Divide(truth)
        ratios.append(ratio)

    # Format first ratio plot
    r0 = ratios[0]
    r0.GetYaxis().SetTitle("Unfolded/Truth Ratio")
    r0.GetYaxis().SetNdivisions(510)
    r0.GetYaxis().SetRangeUser(0.2, 1.8)
    r0.GetYaxis().SetTitleSize(25)
    r0.GetYaxis().SetTitleFont(43)
    r0.GetYaxis().SetTitleOffset(1.5)
    r0.GetYaxis().SetLabelFont(43)
    r0.GetYaxis().SetLabelSize(25)
    if jet:
        r0.GetXaxis().SetTitle("p_{T} [GeV]")
    else:
        r0.GetXaxis().SetTitle("#SigmaE_{T} [GeV]")
    r0.GetXaxis().SetTitleSize(25)
    r0.GetXaxis().SetTitleFont(43)
    r0.GetXaxis().SetTitleOffset(0)
    r0.GetXaxis().SetLabelFont(43)
    r0.GetXaxis().SetLabelSize(25)

    r0.Draw("ep")
    for i in range(1, 6):
        ratios[i].Draw("ep same")

    # Reference lines
    line0 = ROOT.TLine(x_range[0], 1, x_range[1], 1)
    line1 = ROOT.TLine(x_range[0], 0.95, x_range[1], 0.95)
    line2 = ROOT.TLine(x_range[0], 1.05, x_range[1], 1.05)
    for line in [line0, line1, line2]:
        line.SetLineStyle(1 if line == line0 else 2)
        line.Draw("same")

    canvas.Update()
    canvas.Draw()
    if output_name:
        canvas.SaveAs(output_name)


h_full_unfold = []
h_full_unfold_trim_5 = []
h_full_unfold_trim_10 = []
h_unfold = []
h_unfold_trim_5 = []
h_unfold_trim_10 = []
f2 = ROOT.TFile.Open("/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/analysis/unfolding/output_closure_sim_4bin_newetbin_dijet_1000toys.root")
h_full_truth = TH2D(f2.Get("h_truth_calib_dijet"))
h_full_measure = TH2D(f2.Get("h_measure_calib_dijet"))
h_full_truth.SetDirectory(0)
h_full_measure.SetDirectory(0)
for i in range(2):
    h_full_unfold.append(TH2D(f2.Get(f"h_unfold_calib_full_{i+1}")))
    h_full_unfold[-1].SetDirectory(0)
    h_full_unfold.append(TH2D(f2.Get(f"h_unfold_calib_full_extra_{i+1}")))
    h_full_unfold[-1].SetDirectory(0)
    h_full_unfold_trim_5.append(TH2D(f2.Get(f"h_unfold_calib_trim_5_full_{i+1}")))
    h_full_unfold_trim_5[-1].SetDirectory(0)
    h_full_unfold_trim_5.append(TH2D(f2.Get(f"h_unfold_calib_trim_5_full_extra_{i+1}")))
    h_full_unfold_trim_5[-1].SetDirectory(0)
    h_full_unfold_trim_10.append(TH2D(f2.Get(f"h_unfold_calib_trim_10_full_{i+1}")))
    h_full_unfold_trim_10[-1].SetDirectory(0)
    h_full_unfold_trim_10.append(TH2D(f2.Get(f"h_unfold_calib_trim_10_full_extra_{i+1}")))
    h_full_unfold_trim_10[-1].SetDirectory(0)
h_truth = TH2D(f2.Get("h_truth_calib_dijet_half1"))
h_measure = TH2D(f2.Get("h_measure_calib_dijet_half2"))
h_truth.SetDirectory(0)
h_measure.SetDirectory(0)
for i in range(20):
    h_unfold.append(TH2D(f2.Get(f"h_unfold_calib_half_{i+1}")))
    h_unfold[i].SetDirectory(0)
    h_unfold_trim_5.append(TH2D(f2.Get(f"h_unfold_calib_trim_5_half_{i+1}")))
    h_unfold_trim_5[i].SetDirectory(0)
    h_unfold_trim_10.append(TH2D(f2.Get(f"h_unfold_calib_trim_10_half_{i+1}")))
    h_unfold_trim_10[i].SetDirectory(0)
f2.Close()

hj_full_unfold = []
hj_full_unfold_trim_5 = []
hj_full_unfold_trim_10 = []
hj_unfold = []
hj_unfold_trim_5 = []
hj_unfold_trim_10 = []
hc_full_unfold = []
hc_full_unfold_trim_5 = []
hc_full_unfold_trim_10 = []
hc_unfold = []
hc_unfold_trim_5 = []
hc_unfold_trim_10 = []
hj_full_truth = h_full_truth.ProjectionX("hj_full_truth")
hj_full_measure = h_full_measure.ProjectionX("hj_full_measure")
hc_full_truth = h_full_truth.ProjectionY("hc_full_truth")
hc_full_measure = h_full_measure.ProjectionY("hc_full_measure")
hj_truth = h_truth.ProjectionX("hj_truth")
hj_measure = h_measure.ProjectionX("hj_measure")
hc_truth = h_truth.ProjectionY("hc_truth")
hc_measure = h_measure.ProjectionY("hc_measure")

for i in range(len(h_full_unfold)):
    hj_full_unfold.append(h_full_unfold.ProjectionX(f"hj_full_unfold_{i}"))
    hc_full_unfold.append(h_full_unfold.ProjectionY(f"hc_full_unfold_{i}"))
    hj_full_unfold_trim_5.append(h_full_unfold_trim_5.ProjectionX(f"hj_full_unfold_trim_5_{i}"))
    hc_full_unfold_trim_5.append(h_full_unfold_trim_5.ProjectionY(f"hc_full_unfold_trim_5_{i}"))
    hj_full_unfold_trim_10.append(h_full_unfold_trim_10.ProjectionX(f"hj_full_unfold_trim_10_{i}"))
    hc_full_unfold_trim_10.append(h_full_unfold_trim_10.ProjectionY(f"hc_full_unfold_trim_10_{i}"))

for i in range(len(h_unfold)):
    hj_unfold.append(h_unfold[i].ProjectionX(f"hj_unfold_{i}"))
    hc_unfold.append(h_unfold[i].ProjectionY(f"hc_unfold_{i}"))
    hj_unfold_trim_5.append(h_unfold_trim_5[i].ProjectionX(f"hj_unfold_trim_5_{i}"))
    hc_unfold_trim_5.append(h_unfold_trim_5[i].ProjectionY(f"hc_unfold_trim_5_{i}"))
    hj_unfold_trim_10.append(h_unfold_trim_10[i].ProjectionX(f"hj_unfold_trim_10_{i}"))
    hc_unfold_trim_10.append(h_unfold_trim_10[i].ProjectionY(f"hc_unfold_trim_10_{i}"))

draw_unfolded_spectra(hj_full_truth, hj_full_measure, hj_full_unfold, True, (14,82), full_leg_tags, 4, "figure/h_jet_spectrum_full_closure_1000toys.png")

'''
canvas = ROOT.TCanvas("canvas228", "", 600, 800)
pad1 = ROOT.TPad("pad1", "", 0, 0.5, 1, 1.0)
pad1.SetBottomMargin(0.02)  # Adjust the margin for better separation
pad1.Draw()
pad1.cd()
pad1.SetLogy(1)  # Set logarithmic scale for the spectra plot

# Customize the histograms (color, marker, etc.)
hj_full_truth.SetStats(0)
hj_full_measure.SetStats(0)
hj_full_unfold.SetStats(0)
hj_full_truth.SetLineColor(2)
hj_full_truth.SetMarkerColor(2)
hj_full_measure.SetLineColor(1)
hj_full_measure.SetMarkerColor(1)
hj_full_unfold.SetLineColor(colors[1])
hj_full_unfold.SetMarkerColor(colors[1])
hj_full_truth.GetXaxis().SetLabelSize(0)
hj_full_measure.GetXaxis().SetLabelSize(0)
hj_full_unfold.GetXaxis().SetLabelSize(0)

hj_full_truth.Draw()
hj_full_measure.Draw("same")
hj_full_unfold.Draw("same")

# Add legend
leg = ROOT.TLegend(.57, .65, .9, .9)
leg.AddEntry(hj_full_truth,"Truth spectrum","lp")
leg.AddEntry(hj_full_measure,"Measured spectrum","lp")
leg.AddEntry(hj_full_unfold,f"Unfolded spectrum, iter {1}","lp")
leg.Draw()
leg.SetTextSize(0.04)

canvas.cd()
pad2 = ROOT.TPad("pad2", "", 0, 0.0, 1, 0.5)
pad2.SetTopMargin(0.02)
pad2.SetBottomMargin(0.2)
pad2.Draw()
pad2.cd()
ratios = hj_full_unfold.Clone("ratio")
ratios.Divide(hj_full_truth)

ratios.GetYaxis().SetTitle("Unfolded/Truth Ratio")
ratios.GetYaxis().SetNdivisions(510)
ratios.GetYaxis().SetRangeUser(0.2,1.8)
ratios.GetYaxis().SetTitleSize(25)
ratios.GetYaxis().SetTitleFont(43)
ratios.GetYaxis().SetTitleOffset(1.5)
ratios.GetYaxis().SetLabelFont(43)
ratios.GetYaxis().SetLabelSize(25)
ratios.GetXaxis().SetTitle("p_{T} [GeV]")
ratios.GetXaxis().SetTitleSize(25)
ratios.GetXaxis().SetTitleFont(43)
ratios.GetXaxis().SetTitleOffset(0)
ratios.GetXaxis().SetLabelFont(43)
ratios.GetXaxis().SetLabelSize(25)

# Draw ratio plots
ratios.Draw("ep")

line0 = ROOT.TLine(14,1,82,1)
line0.SetLineStyle(1)
line0.Draw("same")
line1 = ROOT.TLine(14,.95,82,.95)
line1.SetLineStyle(2)
line1.Draw("same")
line2 = ROOT.TLine(14,1.05,82,1.05)
line2.SetLineStyle(2)
line2.Draw("same")

# Update canvas
canvas.Update()
canvas.Draw()
canvas.SaveAs("figure/h_jet_spectrum_full_closure_1000toys.png")

canvas = ROOT.TCanvas("canvas244", "", 600, 800)
pad1 = ROOT.TPad("pad1", "", 0, 0.5, 1, 1.0)
pad1.SetBottomMargin(0.02)  # Adjust the margin for better separation
pad1.Draw()
pad1.cd()
pad1.SetLogy(0)  # Set logarithmic scale for the spectra plot
pad1.SetLogx(1)

# Customize the histograms (color, marker, etc.)
hc_full_truth.SetStats(0)
hc_full_measure.SetStats(0)
hc_full_unfold.SetStats(0)
hc_full_truth.SetLineColor(2)
hc_full_truth.SetMarkerColor(2)
hc_full_measure.SetLineColor(1)
hc_full_measure.SetMarkerColor(1)
hc_full_unfold.SetLineColor(colors[1])
hc_full_unfold.SetMarkerColor(colors[1])

hc_full_truth.GetXaxis().SetLabelSize(0)
hc_full_measure.GetXaxis().SetLabelSize(0)
hc_full_unfold.GetXaxis().SetLabelSize(0)

hc_full_truth.GetXaxis().SetRangeUser(0.1,40)
hc_full_measure.GetXaxis().SetRangeUser(0.1,40)
hc_full_truth.Draw()
hc_full_measure.Draw("same")
hc_full_unfold.GetXaxis().SetRangeUser(0.1,40)
hc_full_unfold.Draw("same")

# Add legend
leg = ROOT.TLegend(.17, .1, .6, .35)
leg.AddEntry(hc_full_truth,"Truth spectrum","lp")
leg.AddEntry(hc_full_measure,"Measured spectrum","lp")
leg.AddEntry(hc_full_unfold,f"Unfolded spectrum, iter {1}","lp")
leg.Draw()
leg.SetTextSize(0.04)

canvas.cd()
pad2 = ROOT.TPad("pad2", "", 0, 0.0, 1, 0.5)
pad2.SetTopMargin(0.02)
pad2.SetBottomMargin(0.2)
pad2.Draw()
pad2.cd()
pad2.SetLogx(1)
ratios = hc_full_unfold.Clone("ratio")
ratios.Divide(hc_full_truth)

ratios.GetYaxis().SetTitle("Unfolded/Truth Ratio")
ratios.GetYaxis().SetNdivisions(510)
ratios.GetYaxis().SetRangeUser(0.2,1.8)
ratios.GetXaxis().SetRangeUser(0.1,40)
ratios.GetYaxis().SetTitleSize(25)
ratios.GetYaxis().SetTitleFont(43)
ratios.GetYaxis().SetTitleOffset(1.5)
ratios.GetYaxis().SetLabelFont(43)
ratios.GetYaxis().SetLabelSize(25)
ratios.GetXaxis().SetTitle("#SigmaE_{T} [GeV]")
ratios.GetXaxis().SetTitleSize(25)
ratios.GetXaxis().SetTitleFont(43)
ratios.GetXaxis().SetTitleOffset(0)
ratios.GetXaxis().SetLabelFont(43)
ratios.GetXaxis().SetLabelSize(25)

# Draw ratio plots
ratios.Draw("ep")

line0 = ROOT.TLine(0.1,1,40,1)
line0.SetLineStyle(1)
line0.Draw("same")
line1 = ROOT.TLine(0.1,.95,40,.95)
line1.SetLineStyle(2)
line1.Draw("same")
line2 = ROOT.TLine(0.1,1.05,40,1.05)
line2.SetLineStyle(2)
line2.Draw("same")

# Update canvas
canvas.Update()
canvas.Draw()
canvas.SaveAs("figure/h_et_spectrum_full_closure_1000toys.png")

canvas = ROOT.TCanvas("canvas", "", 600, 800)
pad1 = ROOT.TPad("pad1", "", 0, 0.5, 1, 1.0)
pad1.SetBottomMargin(0.02)  # Adjust the margin for better separation
pad1.Draw()
pad1.cd()
pad1.SetLogy(1)  # Set logarithmic scale for the spectra plot

# Customize the histograms (color, marker, etc.)
hj_truth.SetStats(0)
hj_measure.SetStats(0)
for i in range(6):
    hj_unfold[i].SetStats(0)
hj_truth.SetLineColor(2)
hj_truth.SetMarkerColor(2)
hj_measure.SetLineColor(1)
hj_measure.SetMarkerColor(1)
for i in range(6):
    hj_unfold[i].SetLineColor(colors[i+1])
    hj_unfold[i].SetMarkerColor(colors[i+1])

hj_truth.GetXaxis().SetLabelSize(0)
hj_measure.GetXaxis().SetLabelSize(0)
for i in range(6):
    hj_unfold[i].GetXaxis().SetLabelSize(0)

hj_truth.Draw()
hj_measure.Draw("same")
for i in range(6):
    hj_unfold[i].Draw("same")

# Add legend
leg = ROOT.TLegend(.57, .65, .9, .9)
leg.AddEntry(hj_truth,"Truth spectrum","lp")
leg.AddEntry(hj_measure,"Measured spectrum","lp")
for i in range(6):
    leg.AddEntry(hj_unfold[i],f"Unfolded spectrum, iter {i+1}","lp")
leg.Draw()
leg.SetTextSize(0.04)

canvas.cd()
pad2 = ROOT.TPad("pad2", "", 0, 0.0, 1, 0.5)
pad2.SetTopMargin(0.02)
pad2.SetBottomMargin(0.2)
pad2.Draw()
pad2.cd()
ratios = []
for i in range(6):
    ratios.append(hj_unfold[i].Clone(f"ratio{i}"))
    ratios[i].Divide(hj_truth)

ratios[0].GetYaxis().SetTitle("Unfolded/Truth Ratio")
ratios[0].GetYaxis().SetNdivisions(510)
ratios[0].GetYaxis().SetRangeUser(0.2,1.8)
ratios[0].GetYaxis().SetTitleSize(25)
ratios[0].GetYaxis().SetTitleFont(43)
ratios[0].GetYaxis().SetTitleOffset(1.5)
ratios[0].GetYaxis().SetLabelFont(43)
ratios[0].GetYaxis().SetLabelSize(25)
ratios[0].GetXaxis().SetTitle("p_{T} [GeV]")
ratios[0].GetXaxis().SetTitleSize(25)
ratios[0].GetXaxis().SetTitleFont(43)
ratios[0].GetXaxis().SetTitleOffset(0)
ratios[0].GetXaxis().SetLabelFont(43)
ratios[0].GetXaxis().SetLabelSize(25)

# Draw ratio plots
ratios[0].Draw("ep")
for i in range(1, 6):
    ratios[i].Draw("ep,same")

line0 = ROOT.TLine(14,1,82,1)
line0.SetLineStyle(1)
line0.Draw("same")
line1 = ROOT.TLine(14,.95,82,.95)
line1.SetLineStyle(2)
line1.Draw("same")
line2 = ROOT.TLine(14,1.05,82,1.05)
line2.SetLineStyle(2)
line2.Draw("same")

# Update canvas
canvas.Update()
canvas.Draw()
#canvas.SaveAs("/sphenix/u/egm2153/"+direct+"/h_jet_spectrum_half_closure_toys.png")

canvas = ROOT.TCanvas("canvas40", "", 600, 800)
pad1 = ROOT.TPad("pad1", "", 0, 0.5, 1, 1.0)
pad1.SetBottomMargin(0.02)  # Adjust the margin for better separation
pad1.Draw()
pad1.cd()
pad1.SetLogy(0)  # Set logarithmic scale for the spectra plot
pad1.SetLogx(1)

# Customize the histograms (color, marker, etc.)
hc_truth.SetStats(0)
hc_measure.SetStats(0)
for i in range(6):
    hc_unfold[i].SetStats(0)
hc_truth.SetLineColor(2)
hc_truth.SetMarkerColor(2)
hc_measure.SetLineColor(1)
hc_measure.SetMarkerColor(1)
for i in range(6):
    hc_unfold[i].SetLineColor(colors[i+1])
    hc_unfold[i].SetMarkerColor(colors[i+1])

hc_truth.GetXaxis().SetLabelSize(0)
hc_measure.GetXaxis().SetLabelSize(0)
for i in range(6):
    hc_unfold[i].GetXaxis().SetLabelSize(0)

hc_truth.GetXaxis().SetRangeUser(0.1,40)
hc_measure.GetXaxis().SetRangeUser(0.1,40)
hc_truth.Draw()
hc_measure.Draw("same")
for i in range(6):
    hc_unfold[i].GetXaxis().SetRangeUser(0.1,40)
    hc_unfold[i].Draw("same")

# Add legend
leg = ROOT.TLegend(.17, .1, .6, .35)
leg.AddEntry(hc_truth,"Truth spectrum","lp")
leg.AddEntry(hc_measure,"Measured spectrum","lp")
for i in range(6):
    leg.AddEntry(hc_unfold[i],f"Unfolded spectrum, iter {i+1}","lp")
leg.Draw()
leg.SetTextSize(0.04)

canvas.cd()
pad2 = ROOT.TPad("pad2", "", 0, 0.0, 1, 0.5)
pad2.SetTopMargin(0.02)
pad2.SetBottomMargin(0.2)
pad2.Draw()
pad2.cd()
pad2.SetLogx(1)
ratios = []
for i in range(6):
    ratios.append(hc_unfold[i].Clone(f"ratio{i}"))
    ratios[i].Divide(hc_truth)

ratios[0].GetYaxis().SetTitle("Unfolded/Truth Ratio")
ratios[0].GetYaxis().SetNdivisions(510)
ratios[0].GetYaxis().SetRangeUser(0.2,1.8)
ratios[0].GetXaxis().SetRangeUser(0.1,40)
ratios[0].GetYaxis().SetTitleSize(25)
ratios[0].GetYaxis().SetTitleFont(43)
ratios[0].GetYaxis().SetTitleOffset(1.5)
ratios[0].GetYaxis().SetLabelFont(43)
ratios[0].GetYaxis().SetLabelSize(25)
ratios[0].GetXaxis().SetTitle("#SigmaE_{T} [GeV]")
ratios[0].GetXaxis().SetTitleSize(25)
ratios[0].GetXaxis().SetTitleFont(43)
ratios[0].GetXaxis().SetTitleOffset(0)
ratios[0].GetXaxis().SetLabelFont(43)
ratios[0].GetXaxis().SetLabelSize(25)

# Draw ratio plots
ratios[0].Draw("ep")
for i in range(1, 6):
    ratios[i].Draw("ep,same")

line0 = ROOT.TLine(0.1,1,35,1)
line0.SetLineStyle(1)
line0.Draw("same")
line1 = ROOT.TLine(0.1,.95,35,.95)
line1.SetLineStyle(2)
line1.Draw("same")
line2 = ROOT.TLine(0.1,1.05,35,1.05)
line2.SetLineStyle(2)
line2.Draw("same")

# Update canvas
canvas.Update()
canvas.Draw()
#canvas.SaveAs("/sphenix/u/egm2153/"+direct+"/h_et_spectrum_half_closure_toys.png")
'''
