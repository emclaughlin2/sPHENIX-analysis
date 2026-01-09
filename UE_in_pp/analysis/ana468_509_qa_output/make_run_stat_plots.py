import pandas as pd
import numpy as np
import sys
import ROOT

ana509 = pd.read_csv("ana509_2024p022_v001_event_production_number.csv")
ana468 = pd.read_csv("ana468_2024p012_v001_event_production_number.csv")
ana509['calofitting_ratio'] = ana509['dst_calofitting'] / ana468['dst_calofitting']
ana509['calo_ratio'] = ana509['dst_calo'] / ana468['dst_calo']
ana509['jet_ratio'] = ana509['dst_jetcalo'] / ana468['dst_jetcalo']
ana509['calofitting_difference'] = ana509['dst_calofitting'] - ana468['dst_calofitting']
ana509['calo_difference'] = ana509['dst_calo'] - ana468['dst_calo']
ana509['jet_difference'] = ana509['dst_jetcalo'] - ana468['dst_jetcalo']

ana509.set_index('runnumber', inplace=True)

print(ana509['dst_calofitting'].sum(),ana468['dst_calofitting'].sum(),ana509['dst_calofitting'].sum()/ana468['dst_calofitting'].sum())
print(ana509['dst_calo'].sum(),ana468['dst_calo'].sum(),ana509['dst_calo'].sum()/ana468['dst_calo'].sum())
print(ana509['DST_Jet'].sum(),ana468['dst_jet'].sum(),ana509['DST_Jet'].sum()/ana468['dst_jet'].sum())
print(ana509['dst_calofitting'].sum(),ana468['dst_calofitting'].sum(),ana509['dst_calofitting'].sum() - ana468['dst_calofitting'].sum())
print(ana509['dst_calo'].sum(),ana468['dst_calo'].sum(),ana509['dst_calo'].sum() - ana468['dst_calo'].sum())
print(ana509['DST_Jet'].sum(),ana468['dst_jet'].sum(),ana509['DST_Jet'].sum() - ana468['dst_jet'].sum())

h_calo_fitting_ratio = ROOT.TH1D("h_calo_fitting_ratio","",4000,47250,51250)
h_calo_ratio = ROOT.TH1D("h_calo_ratio","",4000,47250,51250)
h_jet_ratio = ROOT.TH1D("h_jet_ratio","",4000,47250,51250)



for index, row in ana509.iterrows():
	h_calo_fitting_ratio.Fill(index, row['calofitting_ratio'])
	h_calo_ratio.Fill(index, row['calo_ratio'])
	if not np.isnan(row['jet_ratio']):
		h_jet_ratio.Fill(index, row['jet_ratio'])

file = ROOT.TFile.Open("ana468_509_jet10GeV_qa.root")
h_jet10_ratio = file.Get("h_run_ratio")
h_jet10_ratio.SetDirectory(0)

hists = [h_calo_fitting_ratio, h_jet_ratio, h_jet10_ratio]
legs = ["Ratio_{events,calofitting}", "Ratio_{events,jetskimmed}", "Ratio_{events,jet pT > 10 GeV}"]
colors = [1,2,4]

for hist in hists:
	for i in range(1, hist.GetNbinsX() + 1):
		hist.SetBinError(i, 0)

canvas = ROOT.TCanvas("canvas","",1000,500)
leg = ROOT.TLegend(.6,.7,.9,.9)
leg.SetTextSize(0.035)
for i, hist in enumerate(hists):
	hist.SetMarkerColor(colors[i])
	hist.SetLineColor(colors[i])
	hist.GetYaxis().SetTitle("ana509/ana468 Ratio")
	hist.GetXaxis().SetTitle("runnumber")
	hist.SetMarkerSize(0.5)
	hist.SetMarkerStyle(20)
	hist.SetStats(0)
	leg.AddEntry(hist,legs[i],"p")
	if i == 0:
		hist.Draw("p")
	else:
		hist.Draw("p,same")
leg.Draw()
canvas.SaveAs("ana468_509_run_stats.png")

h_jet10_fit = h_jet10_ratio.Clone("h_jet10_fit")
h_jet10_fit.Divide(h_calo_fitting_ratio)
h_jet10_jet = h_jet10_ratio.Clone("h_jet10_jet")
h_jet10_jet.Divide(h_jet_ratio)
hist_ratios = [h_jet10_fit, h_jet10_jet]
ratio_legs = ["Ratio_{events,jet pT > 10 GeV}/Ratio_{events,calofitting}", "Ratio_{events,jet pT > 10 GeV}/Ratio_{events,jetskimmed}"]

canvas = ROOT.TCanvas("canvas","",1000,500)
leg = ROOT.TLegend(.6,.7,.9,.9)
leg.SetTextSize(0.035)
for i, hist in enumerate(hist_ratios):
	hist.SetMarkerColor(colors[i])
	hist.GetYaxis().SetTitle("(ana509/ana468 Ratio)/(ana509/ana468 Ratio)")
	hist.GetXaxis().SetTitle("runnumber")
	hist.SetMarkerSize(0.5)
	hist.SetMarkerStyle(20)
	hist.SetStats(0)
	leg.AddEntry(hist,ratio_legs[i],"p")
	if i == 0:
		hist.Draw("p")
	else:
		hist.Draw("p,same")
leg.Draw()
canvas.SaveAs("ana468_509_run_ratio_stats.png")

canvas = ROOT.TCanvas("canvas","",1000,500)
leg = ROOT.TLegend(.6,.7,.9,.9)
leg.SetTextSize(0.035)
for i, hist in enumerate(hists):
	hist.SetMarkerColor(colors[i])
	hist.SetLineColor(colors[i])
	hist.GetYaxis().SetRangeUser(0,2)
	hist.GetYaxis().SetTitle("ana509/ana468 Ratio")
	hist.GetXaxis().SetTitle("runnumber")
	hist.SetMarkerSize(0.5)
	hist.SetMarkerStyle(20)
	hist.SetStats(0)
	leg.AddEntry(hist,legs[i],"p")
	if i == 0:
		hist.Draw("p")
	else:
		hist.Draw("p,same")
leg.Draw()
canvas.SaveAs("ana468_509_run_stats_zoom.png")

canvas = ROOT.TCanvas("canvas","",1000,500)
leg = ROOT.TLegend(.6,.7,.9,.9)
leg.SetTextSize(0.035)
for i, hist in enumerate(hist_ratios):
	hist.SetMarkerColor(colors[i])
	hist.GetYaxis().SetTitle("(ana509/ana468 Ratio)/(ana509/ana468 Ratio)")
	hist.GetXaxis().SetTitle("runnumber")
	hist.SetMarkerSize(0.5)
	hist.SetMarkerStyle(20)
	hist.GetYaxis().SetRangeUser(0,2)
	hist.SetStats(0)
	leg.AddEntry(hist,ratio_legs[i],"p")
	if i == 0:
		hist.Draw("p")
	else:
		hist.Draw("p,same")
leg.Draw()
canvas.SaveAs("ana468_509_run_ratio_stats_zoom.png")

'''
canvas = ROOT.TCanvas("canvas","",1000,500)
leg = ROOT.TLegend(.7,.7,.9,.9)
for i, hist in enumerate(hists):
	if i == 2:
		continue
	hist.SetLineColor(i+1)
	leg.AddEntry(hist,legs[i],"l")
	if i == 0:
		hist.Draw("hist")
	else:
		hist.Draw("hist,same")
leg.Draw()
canvas.SaveAs("ana468_509_calofitting_run_stats.png")

canvas1 = ROOT.TCanvas("canvas","",1000,500)
leg1 = ROOT.TLegend(.7,.7,.9,.9)
for i, hist in enumerate(hists):
	if i == 1:
		continue
	hist.SetLineColor(i+1)
	leg1.AddEntry(hist,legs[i],"l")
	if i == 0:
		hist.Draw("hist")
	else:
		hist.Draw("hist,same")
leg1.Draw()
canvas1.SaveAs("ana468_509_jetskim_run_stats.png")
'''
