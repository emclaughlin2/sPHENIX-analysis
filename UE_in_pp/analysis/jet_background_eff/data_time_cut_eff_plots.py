import ROOT
import os

def load_sphenix_style():
    style_path = os.path.expandvars("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C")
    if os.path.exists(style_path):
        ROOT.gROOT.ProcessLine(f'.L {style_path}')
        ROOT.gROOT.ProcessLine('SetsPhenixStyle()')
    else:
        print("Warning: sPHENIX style file not found, using default ROOT style.")

load_sphenix_style()
# Define eta range labels and suffixes
eta_labels = ["All #Delta#eta", "#Delta#eta < 0.4", "#Delta#eta > 0.4"]
eta_suffixes = ["", "_eta_lt_0.4", "_eta_gt_0.4"]

# Open the ROOT file
file = ROOT.TFile.Open("analysis_data_ana509_output/output_dijet_and_time_bkg_cut.root")

# Ensure output directory exists
os.makedirs("time_bkg_cut_plots", exist_ok=True)

# Loop over eta ranges for TH1D plots
for i, (label, suffix) in enumerate(zip(eta_labels, eta_suffixes)):
    h_jettime = file.Get(f"h_jettime{suffix}")
    h_jetunweighttime = file.Get(f"h_jetunweighttime{suffix}")
    h_deltajettime = file.Get(f"h_deltajettime{suffix}")
    h_deltajetunweighttime = file.Get(f"h_deltajetunweighttime{suffix}")

    # Plot Jet Time
    c1 = ROOT.TCanvas(f"c1_{i}", f"Jet Time {label}", 800, 600)
    h_jettime.SetLineColor(ROOT.kBlue+1)
    h_jetunweighttime.SetLineColor(ROOT.kRed+1)
    h_jetunweighttime.GetYaxis().SetRangeUser(0,180)
    h_jettime.SetTitle(f"Jet Time {label}; t [ns]; Entries")
    h_jetunweighttime.SetTitle(f"Jet Time {label}; t [ns]; Entries")
    h_jetunweighttime.Draw("HIST")
    h_jettime.Draw("HIST SAME")

    leg1 = ROOT.TLegend(0.2, 0.8, 0.88, 0.92)
    leg1.AddEntry("", f"{label}, p_{{T,lead reco}} = [30,45] GeV, Dijet requirement, ana509 prod", "")
    leg1.AddEntry(h_jettime, f"Weighted time: Mean={h_jettime.GetMean():.2f}, RMS={h_jettime.GetRMS():.2f}", "l")
    leg1.AddEntry(h_jetunweighttime, f"Unweighted time: Mean={h_jetunweighttime.GetMean():.2f}, RMS={h_jetunweighttime.GetRMS():.2f}", "l")
    leg1.SetTextSize(0.03)
    leg1.Draw()
    c1.SaveAs(f"time_bkg_cut_plots/jet_time{suffix}.png")

    # Plot Delta Jet Time
    c2 = ROOT.TCanvas(f"c2_{i}", f"Delta Jet Time {label}", 800, 600)
    h_deltajettime.SetLineColor(ROOT.kBlue+1)
    h_deltajetunweighttime.SetLineColor(ROOT.kRed+1)
    h_deltajettime.SetTitle(f"Delta Jet Time {label}; #Delta t [ns]; Entries")
    h_deltajetunweighttime.SetTitle(f"Delta Jet Time {label}; #Delta t [ns]; Entries")
    h_deltajettime.GetYaxis().SetRangeUser(0,180)
    h_deltajettime.Draw("HIST")
    h_deltajetunweighttime.Draw("HIST SAME")

    leg2 = ROOT.TLegend(0.2, 0.8, 0.88, 0.92)
    leg2.AddEntry("", f"{label}, p_{{T,lead reco}} = [30,45] GeV, Dijet requirement, ana509 prod", "")
    leg2.AddEntry(h_deltajettime, f"Weighted time: Mean={h_deltajettime.GetMean():.2f}, RMS={h_deltajettime.GetRMS():.2f}", "l")
    leg2.AddEntry(h_deltajetunweighttime, f"Unweighted time: Mean={h_deltajetunweighttime.GetMean():.2f}, RMS={h_deltajetunweighttime.GetRMS():.2f}", "l")
    leg2.SetTextSize(0.03)
    leg2.Draw()
    c2.SaveAs(f"time_bkg_cut_plots/delta_jet_time{suffix}.png")

    # Plot TH2D if available
    h_2D_deltajettime = file.Get(f"h_2D_deltajettime{suffix}")
    h_2D_deltajetunweighttime = file.Get(f"h_2D_deltajetunweighttime{suffix}")
    h_jettimecorr = file.Get(f"h_jettimecorr{suffix}")
    h_deltatcorr = file.Get(f"h_deltatcorr{suffix}")

    if h_2D_deltajettime:
        c3 = ROOT.TCanvas(f"c3_{i}", f"2D Delta Jet Time {label}", 800, 600)
        h_2D_deltajettime.SetTitle(f"2D Delta Jet Time {label}; t_{{lead}} [ns]; #Delta t [ns]")
        h_2D_deltajettime.Draw("COLZ")
        c3.SaveAs(f"time_bkg_cut_plots/2D_delta_jet_time{suffix}.png")

    if h_2D_deltajetunweighttime:
        c4 = ROOT.TCanvas(f"c4_{i}", f"2D Delta Jet Unweighted Time {label}", 800, 600)
        h_2D_deltajetunweighttime.SetTitle(f"2D Delta Jet Unweighted Time {label}; t_{{lead}} [ns]; #Delta t [ns]")
        h_2D_deltajetunweighttime.Draw("COLZ")
        c4.SaveAs(f"time_bkg_cut_plots/2D_delta_jet_unweight_time{suffix}.png")

    if h_jettimecorr:
        c5 = ROOT.TCanvas(f"c5_{i}", f"Jet Time Corr {label}", 800, 600)
        h_jettimecorr.SetTitle(f"Jet Time Corr {label}; Weighted t_{{lead}} [ns]; Unweighted t_{{lead}}")
        h_jettimecorr.Draw("COLZ")
        c5.SaveAs(f"time_bkg_cut_plots/jet_time_corr{suffix}.png")

    if h_deltatcorr:
        c6 = ROOT.TCanvas(f"c6_{i}", f"Delta t Corr {label}", 800, 600)
        h_deltatcorr.SetTitle(f"Delta t Corr {label}; Weighted #Delta t [ns]; Unweighted #Delta t [ns]")
        h_deltatcorr.Draw("COLZ")
        c6.SaveAs(f"time_bkg_cut_plots/delta_t_corr{suffix}.png")
'''
# Retrieve TH1D histograms
h_jettime = file.Get("h_jettime")
h_jetunweighttime = file.Get("h_jetunweighttime")
h_deltajettime = file.Get("h_deltajettime")
h_deltajetunweighttime = file.Get("h_deltajetunweighttime")

# Retrieve TH2D histograms
h_2D_deltajettime = file.Get("h_2D_deltajettime")
h_2D_deltajetunweighttime = file.Get("h_2D_deltajetunweighttime")
h_jettimecorr = file.Get("h_jettimecorr")
h_deltatcorr = file.Get("h_deltatcorr")

# Plot TH1D: t [ns]
c1 = ROOT.TCanvas("c1", "Jet Time", 800, 600)
h_jettime.SetLineColor(ROOT.kBlue+1)
h_jetunweighttime.SetLineColor(ROOT.kRed+1)
h_jetunweighttime.GetYaxis().SetRangeUser(0,180)
h_jettime.SetTitle("Jet Time; t [ns]; Entries")
h_jetunweighttime.SetTitle("Jet Time; t [ns]; Entries")
h_jetunweighttime.Draw("HIST")
h_jettime.Draw("HIST SAME")

leg1 = ROOT.TLegend(0.2, 0.8, 0.88, 0.92)
leg1.AddEntry("","p_{T,lead reco} = [30,45] GeV, Dijet requirement, ana509 prod","")
leg1.AddEntry(h_jettime, f"Weighted time: Mean={h_jettime.GetMean():.2f}, RMS={h_jettime.GetRMS():.2f}", "l")
leg1.AddEntry(h_jetunweighttime, f"Unweighted time: Mean={h_jetunweighttime.GetMean():.2f}, RMS={h_jetunweighttime.GetRMS():.2f}", "l")
leg1.SetTextSize(0.03)
leg1.Draw()
c1.SaveAs("time_bkg_cut_plots/jet_time.png")

# Plot TH1D: delta t [ns]
c2 = ROOT.TCanvas("c2", "Delta Jet Time", 800, 600)
h_deltajettime.SetLineColor(ROOT.kBlue+1)
h_deltajetunweighttime.SetLineColor(ROOT.kRed+1)
h_deltajettime.SetTitle("Delta Jet Time; #Delta t [ns]; Entries")
h_deltajetunweighttime.SetTitle("Delta Jet Time; #Delta t [ns]; Entries")
h_deltajettime.GetYaxis().SetRangeUser(0,180)
h_deltajettime.Draw("HIST")
h_deltajetunweighttime.Draw("HIST SAME")

leg2 = ROOT.TLegend(0.2, 0.8, 0.88, 0.92)
leg2.AddEntry("","p_{T,lead reco} = [30,45] GeV, Dijet requirement, ana509 prod","")
leg2.AddEntry(h_deltajettime, f"Weighted time: Mean={h_deltajettime.GetMean():.2f}, RMS={h_deltajettime.GetRMS():.2f}", "l")
leg2.AddEntry(h_deltajetunweighttime, f"Unweighted time: Mean={h_deltajetunweighttime.GetMean():.2f}, RMS={h_deltajetunweighttime.GetRMS():.2f}", "l")
leg2.SetTextSize(0.03)
leg2.Draw()
c2.SaveAs("time_bkg_cut_plots/delta_jet_time.png")

# Plot TH2D individually
c3 = ROOT.TCanvas("c3", "2D Delta Jet Time", 800, 600)
h_2D_deltajettime.SetTitle("2D Delta Jet Time; t_{lead} [ns]; #Delta t [ns]")
h_2D_deltajettime.Draw("COLZ")
c3.SaveAs("time_bkg_cut_plots/2D_delta_jet_time.png")

c4 = ROOT.TCanvas("c4", "2D Delta Jet Unweighted Time", 800, 600)
h_2D_deltajetunweighttime.SetTitle("2D Delta Jet Unweighted Time; t_{lead} [ns]; #Delta t [ns]")
h_2D_deltajetunweighttime.Draw("COLZ")
c4.SaveAs("time_bkg_cut_plots/2D_delta_jet_unweight_time.png")

c5 = ROOT.TCanvas("c5", "Jet Time Corr", 800, 600)
h_jettimecorr.SetTitle("Jet Time Corr; Weighted t_{lead} [ns]; Unweighted t_{lead}")
h_jettimecorr.Draw("COLZ")
c5.SaveAs("time_bkg_cut_plots/jet_time_corr.png")

c6 = ROOT.TCanvas("c6", "Delta t Corr", 800, 600)
h_deltatcorr.SetTitle("Delta t Corr; Weighted #Delta t [ns]; Unweighted #Delta t [ns]")
h_deltatcorr.Draw("COLZ")
c6.SaveAs("time_bkg_cut_plots/delta_t_corr.png")
'''