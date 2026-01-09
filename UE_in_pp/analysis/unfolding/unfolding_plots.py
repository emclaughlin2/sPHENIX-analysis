import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency, TH1D, TH2D
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb
from array import array
gROOT.LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C");
gROOT.ProcessLine("SetsPhenixStyle()")

direct = 'spring_2025/results_6_2'

rgb = [[230, 25, 75], [60, 180, 75], [255, 225, 25], [0, 130, 200], [245, 130, 48], [145, 30, 180], [70, 240, 240], [240, 50, 230], [210, 245, 60], [250, 190, 212], [0, 128, 128], [220, 190, 255], [170, 110, 40], [128, 128, 128], [128, 0, 0], [0, 0, 0], [128, 128, 0], [255, 215, 180], [0, 0, 128], [34, 139, 34]]
colors = [TColor.GetColor(rgb[i][0],rgb[i][1],rgb[i][2]) for i in range(len(rgb))]

h_unfold = []
f2 = ROOT.TFile.Open("/sphenix/user/egm2153/calib_study/analysis/UE_in_pp/analysis/unfolding/output_unfolded_data_4bin_dijet_trim.root")
h_truth = TH2D(f2.Get("h_truth_calib_dijet"))
h_measure = TH2D(f2.Get("h_calibjet_pt_dijet_eff"))
h_truth.SetDirectory(0)
h_measure.SetDirectory(0)
for i in range(6):
    h_unfold.append(TH2D(f2.Get(f"h_unfold_calib_all_{i+1}")))
    h_unfold[i].SetDirectory(0)
f2.Close()
hj_truth = h_truth.ProjectionX("hj_truth")
hj_measure = h_measure.ProjectionX("hj_measure")
hc_truth = h_truth.ProjectionY("hc_truth")
hc_measure = h_measure.ProjectionY("hc_measure")
hj_unfold = []
hc_unfold = []
for i in range(6):
    hj_unfold.append(h_unfold[i].ProjectionX(f"hj_unfold_{i}"))
    hc_unfold.append(h_unfold[i].ProjectionY(f"hc_unfold_{i}"))
hj_truth.Scale(hj_measure.GetBinContent(1)/hj_truth.GetBinContent(1))
hc_truth.Scale(hc_measure.GetBinContent(1)/hc_truth.GetBinContent(1))
canvas = ROOT.TCanvas("canvas26", "", 600, 800)
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

hj_truth.GetYaxis().SetRangeUser(0.01,500000)
hj_truth.Draw()
hj_measure.Draw("same")
for i in range(5,-1,-1):
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

ratios[5].GetYaxis().SetTitle("Unfolded/Truth Ratio")
ratios[5].GetYaxis().SetNdivisions(510)
ratios[5].GetYaxis().SetRangeUser(0.2,1.8)
ratios[5].GetYaxis().SetTitleSize(25)
ratios[5].GetYaxis().SetTitleFont(43)
ratios[5].GetYaxis().SetTitleOffset(1.5)
ratios[5].GetYaxis().SetLabelFont(43)
ratios[5].GetYaxis().SetLabelSize(25)
ratios[5].GetXaxis().SetTitle("p_{T} [GeV]")
ratios[5].GetXaxis().SetTitleSize(25)
ratios[5].GetXaxis().SetTitleFont(43)
ratios[5].GetXaxis().SetTitleOffset(0)
ratios[5].GetXaxis().SetLabelFont(43)
ratios[5].GetXaxis().SetLabelSize(25)

# Draw ratio plots
ratios[5].Draw("ep")
for i in range(4,-1,-1):
    ratios[i].Draw("ep,same")

line0 = ROOT.TLine(17,1,82,1)
line0.SetLineStyle(1)
line0.Draw("same")
line1 = ROOT.TLine(17,.95,82,.95)
line1.SetLineStyle(2)
line1.Draw("same")
line2 = ROOT.TLine(17,1.05,82,1.05)
line2.SetLineStyle(2)
line2.Draw("same")

# Update canvas
canvas.Update()
canvas.Draw()
canvas.SaveAs("/sphenix/u/egm2153/"+direct+"/h_jet_spectrum_trim_unfolded_toys.png")

canvas = ROOT.TCanvas("canvas51", "", 600, 800)
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
hc_truth.GetYaxis().SetRangeUser(0,15000)
hc_truth.Draw()
hc_measure.Draw('same')
for i in range(6):
    hc_unfold[i].GetXaxis().SetRangeUser(0.1,40)
    hc_unfold[i].Draw("same")

# Add legend
leg = ROOT.TLegend(.17, .65, .6, .9)
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
canvas.SaveAs("/sphenix/u/egm2153/"+direct+"/h_et_spectrum_trim_unfolded_toys.png")
h_truth.GetNbinsY()
truth_prof = h_truth.ProfileX("truth_prof",1,h_truth.GetNbinsY()-1)
measure_prof = h_measure.ProfileX("measure_prof",1,h_measure.GetNbinsY()-1)
unfold_prof = []
for i in range(6):
    unfold_prof.append(h_unfold[i].ProfileX(f"unfold_prof_{i}",1,h_unfold[i].GetNbinsY()-1))
canvas = ROOT.TCanvas("canvas113", "", 600, 500)

# Customize the histograms (color, marker, etc.)
truth_prof.SetStats(0)
measure_prof.SetStats(0)
for i in range(6):
    unfold_prof[i].SetStats(0)
truth_prof.SetLineColor(2)
truth_prof.SetMarkerColor(2)
measure_prof.SetLineColor(1)
measure_prof.SetMarkerColor(1)
for i in range(6):
    unfold_prof[i].SetLineColor(colors[i+1])
    unfold_prof[i].SetMarkerColor(colors[i+1])

truth_prof.GetXaxis().SetLabelSize(0)
measure_prof.GetXaxis().SetLabelSize(0)
for i in range(6):
    unfold_prof[i].GetXaxis().SetLabelSize(0)

truth_prof.GetYaxis().SetRangeUser(0,8)
truth_prof.SetYTitle("<#SigmaE_{T}> [GeV]")
truth_prof.SetLabelSize(0.05)
truth_prof.Draw()
measure_prof.Draw("same")
for i in range(6):
    unfold_prof[i].Draw("same")

# Add legend
leg = ROOT.TLegend(.5, .65, .9, .9)
leg.AddEntry(truth_prof,"Truth spectrum","lp")
leg.AddEntry(measure_prof,"Measured spectrum","lp")
for i in range(6):
    leg.AddEntry(unfold_prof[i],f"Unfolded spectrum, iter {i+1}","lp")
leg.Draw()
leg.SetTextSize(0.04)

# Update canvas
canvas.Update()
canvas.Draw()
canvas.SaveAs("/sphenix/u/egm2153/"+direct+"/h_et_pt_trim_unfolded_toys.png")
for it in range(6):
    for i in range(1, unfold_prof[it].GetNbinsX() + 1):
        print(i, unfold_prof[it].GetBinError(i)/unfold_prof[it].GetBinContent(i))
    print()
for it in range(6):
    for i in range(1, unfold_prof[it].GetNbinsX() + 1):
        print(i, unfold_prof[it].GetBinError(i)/unfold_prof[it].GetBinContent(i))
    print()