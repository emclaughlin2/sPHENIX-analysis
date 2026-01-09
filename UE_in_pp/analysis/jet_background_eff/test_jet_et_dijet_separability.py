import ROOT
import os
from array import array

def load_sphenix_style():
    style_path = os.path.expandvars("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C")
    if os.path.exists(style_path):
        ROOT.gROOT.ProcessLine(f'.L {style_path}')
        ROOT.gROOT.ProcessLine('SetsPhenixStyle()')
    else:
        print("Warning: sPHENIX style file not found, using default ROOT style.")

load_sphenix_style()

colors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+1, ROOT.kBlack, ROOT.kPink+1, ROOT.kRed+3, ROOT.kBlue+3, ROOT.kGreen]
cut_name = ["dijet","none"]
legend_entries = ['dijet','dijet','none','none']

new_files = [
    "analysis_sim_run28_output/output_sim_dijet_bkg_cut.root"
]

# Input ROOT files
old_files = [
    "../jet_calo_analysis_run28_jet10_3sigma_topo_dijetcut_jescorr_notruthmatch_recojet_caloacceptance_vz_lt_60.root",
    "../jet_calo_analysis_run28_jet10_3sigma_topo_nobkgcut_jescorr_notruthmatch_recojet_caloacceptance_vz_lt_60.root"
]

new_hist_names = ["h_pass_cut_measure_tight", "h_total_measure_tight"]
old_hist_names = ["h_ue_pt_transverse","h_et_transverse"]

old_profx = []
'''
files = [ROOT.TFile.Open(path) for path in old_files]
hists = [[f.Get(hname) for hname in old_hist_names] for f in files]   
for file_idx, (f, hlist) in enumerate(zip(files, hists)):
    for hist_idx, h2d in enumerate(hlist):  
        old_profx.append(h2d)
        old_profx[-1].SetDirectory(0)
'''
files = [ROOT.TFile.Open(path) for path in new_files]
hists = [[f.Get(hname) for hname in new_hist_names] for f in files]
for file_idx, (f, hlist) in enumerate(zip(files, hists)):
    for hist_idx, h2d in enumerate(hlist):
        if not h2d:
            print(f"Histogram {hist_names[hist_idx]} not found in file {file_paths[file_idx]}")
            continue
        projy = h2d.ProjectionY()
        projy.SetName(f"{h2d.GetName()}_projy_file{file_idx}")
        profx = h2d.ProfileX()
        profx.SetName(f"{h2d.GetName()}_profx_file{file_idx}")
        old_profx.append(profx)
        old_profx[-1].SetDirectory(0)
        old_profx.append(projy)
        old_profx[-1].SetDirectory(0)

old_profx_hist = []
for i, prof in enumerate(old_profx):
    if i % 2 != 0: 
        continue
    nbins = prof.GetNbinsX()
    xaxis = prof.GetXaxis()
    edges = array('d', [xaxis.GetBinLowEdge(i) for i in range(1, nbins + 2)])
    h1 = ROOT.TH1D(f"{prof.GetName()}_hist", prof.GetTitle(), nbins, edges)
    for ibin in range(1, nbins + 1):
        h1.SetBinContent(ibin, prof.GetBinContent(ibin))
        h1.SetBinError(ibin, prof.GetBinError(ibin))
    old_profx_hist.append(h1)
old_profx_hist[0].Divide(old_profx_hist[0], old_profx_hist[1], 1, 1, "B")

ratio = old_profx[1].Clone("ratio")
ratio.Divide(ratio, old_profx[3], 1.0, 1.0, "B")

canvas = ROOT.TCanvas("c_profx_ratio", "Profiles + Ratio", 800, 900)
pad1 = ROOT.TPad("pad1", "Top pad", 0, 0.35, 1, 1.0)
pad2 = ROOT.TPad("pad2", "Bottom pad", 0, 0.05, 1, 0.35)
pad1.SetBottomMargin(0)  # no X label overlap
pad2.SetTopMargin(0)
pad2.SetBottomMargin(0.3)
pad1.Draw()
pad2.Draw()

# --- Top pad: draw profiles ---
pad1.cd()
legend = ROOT.TLegend(0.6, 0.2, 0.9, 0.35)
legend.SetTextSize(0.03)

for i, prof in enumerate(old_profx):
    if i % 2 != 0:
        continue
    prof.SetLineColor(colors[i])
    prof.SetMarkerColor(colors[i])
    prof.SetMarkerStyle(20 + i)
    drawopt = "E SAME" if i > 0 else "E"
    prof.GetYaxis().SetRangeUser(0, 2.5)
    prof.GetYaxis().SetTitle("Transverse <#SigmaE_{T}> [GeV]")
    prof.Draw(drawopt)
    legend.AddEntry(prof, legend_entries[i], "lep")

legend.Draw()
pad1.RedrawAxis()

# --- Bottom pad: ratio ---
pad2.cd()
old_profx_hist[0].SetLineColor(ROOT.kBlack)
old_profx_hist[0].SetMarkerColor(ROOT.kBlack)
old_profx_hist[0].SetMarkerStyle(20)
old_profx_hist[0].GetYaxis().SetTitle("Ratio")
old_profx_hist[0].GetYaxis().SetNdivisions(505)
old_profx_hist[0].GetYaxis().SetTitleSize(0.1)
old_profx_hist[0].GetYaxis().SetLabelSize(0.09)
old_profx_hist[0].GetXaxis().SetTitleSize(0.1)
old_profx_hist[0].GetXaxis().SetLabelSize(0.09)
old_profx_hist[0].GetYaxis().SetTitleOffset(0.5)
old_profx_hist[0].GetXaxis().SetTitle("p_{T}^{calib lead} [GeV]")
old_profx_hist[0].GetYaxis().SetRangeUser(0, 1.2)
old_profx_hist[0].Draw("E")

# --- Reference line at ratio = 1 ---
line = ROOT.TLine(old_profx_hist[0].GetXaxis().GetXmin(), 1.0,
                  old_profx_hist[0].GetXaxis().GetXmax(), 1.0)
line.SetLineStyle(2)
line.SetLineColor(ROOT.kGray + 2)
line.Draw("SAME")
canvas.SaveAs(f"jet_bkgeff_comparison_binned_dijet.png")

canvas = ROOT.TCanvas("c_etprof_ratio", "Profiles + Ratio", 800, 900)
pad1 = ROOT.TPad("pad1", "Top pad", 0, 0.35, 1, 1.0)
pad2 = ROOT.TPad("pad2", "Bottom pad", 0, 0.05, 1, 0.35)
pad1.SetBottomMargin(0)  # no X label overlap
pad2.SetTopMargin(0)
pad2.SetBottomMargin(0.3)
pad1.Draw()
pad2.Draw()

# --- Top pad: draw profiles ---
pad1.cd()
legend = ROOT.TLegend(0.6, 0.2, 0.9, 0.35)
legend.SetTextSize(0.03)

for i, prof in enumerate(old_profx):
    if i % 2 == 0: 
        continue
    prof.SetLineColor(colors[i])
    prof.SetMarkerColor(colors[i])
    prof.SetMarkerStyle(20 + i)
    drawopt = "E SAME" if i > 1 else "E"
    prof.GetXaxis().SetRangeUser(-1.0,14)
    prof.Draw(drawopt)
    legend.AddEntry(prof, legend_entries[i], "lep")

legend.Draw()
pad1.RedrawAxis()

# --- Bottom pad: ratio ---
pad2.cd()
ratio.SetLineColor(ROOT.kBlack)
ratio.SetMarkerColor(ROOT.kBlack)
ratio.SetMarkerStyle(20)
ratio.GetYaxis().SetTitle("Ratio")
ratio.GetYaxis().SetNdivisions(505)
ratio.GetYaxis().SetTitleSize(0.1)
ratio.GetYaxis().SetLabelSize(0.09)
ratio.GetXaxis().SetTitleSize(0.1)
ratio.GetXaxis().SetLabelSize(0.09)
ratio.GetYaxis().SetTitleOffset(0.5)
ratio.GetXaxis().SetTitle("Transverse #SigmaE_{T} [GeV]")
ratio.GetYaxis().SetRangeUser(0, 1)
ratio.GetXaxis().SetRangeUser(-1.0,14.0)
ratio.Draw("E")

# --- Reference line at ratio = 1 ---
line = ROOT.TLine(ratio.GetXaxis().GetXmin(), 1.0,
                  ratio.GetXaxis().GetXmax(), 1.0)
line.SetLineStyle(2)
line.SetLineColor(ROOT.kGray + 2)
line.Draw("SAME")
canvas.SaveAs(f"etdist_bkgeff_comparison_binned_dijet.png")

# Separate lists for pass_cut and total histograms
projx_pass_cut = []
projx_total = []
profx_pass_cut = []
profx_total = []
projy_pass_cut = []
projy_total = []
# Open files and get histograms
files = [ROOT.TFile.Open(path) for path in new_files]
hists = [[f.Get(hname) for hname in new_hist_names] for f in files]
# For each file and histogram

for file_idx, (f, hlist) in enumerate(zip(files, hists)):
    for hist_idx, h2d in enumerate(hlist):
        if not h2d:
            print(f"Histogram {hist_names[hist_idx]} not found in file {file_paths[file_idx]}")
            continue
        projx = h2d.ProjectionX()
        projx.SetName(f"{h2d.GetName()}_projx_file{file_idx}")
        profx = h2d.ProfileX()
        profx.SetName(f"{h2d.GetName()}_profx_file{file_idx}")
        nxbins = h2d.GetNbinsX()
        temp_projy = []
        for xbin in range(1, nxbins + 1, 2):
            temp_projy.append(h2d.ProjectionY(f"{h2d.GetName()}_projy_xbin{xbin}_file{file_idx}", xbin, xbin))
            temp_projy[-1].SetDirectory(0)
        if hist_idx == 0:
            projx_pass_cut.append(projx)
            projx_pass_cut[-1].SetDirectory(0)
            profx_pass_cut.append(profx)
            profx_pass_cut[-1].SetDirectory(0)
            projy_pass_cut.append(temp_projy)
        else:
            projx_total.append(projx)
            projx_total[-1].SetDirectory(0)
            profx_total.append(profx)
            profx_total[-1].SetDirectory(0)
            projy_total.append(temp_projy)

colors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+1, ROOT.kBlack, ROOT.kPink+1, ROOT.kRed+3, ROOT.kBlue+3, ROOT.kGreen]
cut_name = ["dijet","none"]
legend_entries = ['dijet_new','none_new','dijet_calo_accept','none_calo_accept']

profiles = []
#profiles.append(profx_pass_cut[0])
profiles.append(old_profx[0])

profiles_total = []
#profiles_total.append(profx_total[0])
profiles_total.append(old_profx[1])

profx_hist = []
profx_hist_total = []

for prof in profiles:
    nbins = prof.GetNbinsX()
    xaxis = prof.GetXaxis()
    h1 = ROOT.TH1D(f"{prof.GetName()}_hist", prof.GetTitle(), nbins, xaxis.GetXmin(), xaxis.GetXmax())
    for ibin in range(1, nbins + 1):
        h1.SetBinContent(ibin, prof.GetBinContent(ibin))
        h1.SetBinError(ibin, prof.GetBinError(ibin))
    profx_hist.append(h1)

for prof in profiles_total:
    nbins = prof.GetNbinsX()
    xaxis = prof.GetXaxis()
    h1 = ROOT.TH1D(f"{prof.GetName()}_hist", prof.GetTitle(), nbins, xaxis.GetXmin(), xaxis.GetXmax())
    for ibin in range(1, nbins + 1):
        h1.SetBinContent(ibin, prof.GetBinContent(ibin))
        h1.SetBinError(ibin, prof.GetBinError(ibin))
    profx_hist_total.append(h1)

for i in range(len(profx_hist)):
    profx_hist[i].Divide(profx_hist[i], profx_hist_total[i], 1, 1, "B")