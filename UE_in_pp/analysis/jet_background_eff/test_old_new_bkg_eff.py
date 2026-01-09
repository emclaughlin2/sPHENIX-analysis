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

new_files = [
    "analysis_sim_run28_output/output_sim_jet10_dijet.root",
    "analysis_sim_run28_output/output_sim_jet10_efrac.root"
]

# Input ROOT files
old_files = [
    #"../jet_calo_analysis_run28_jet10_3sigma_topo_dijetcut_truthmatch.root",
    #"../jet_calo_analysis_run28_jet10_3sigma_topo_efraccut_truthmatch.root",
    #"../jet_calo_analysis_run28_jet10_3sigma_topo_nobkgcut_truthmatch.root",
    "../jet_calo_analysis_run28_jet10_3sigma_topo_dijetcut_jescorr_notruthmatch.root",
    "../jet_calo_analysis_run28_jet10_3sigma_topo_efraccut_jescorr_notruthmatch.root",
    "../jet_calo_analysis_run28_jet10_3sigma_topo_nobkgcut_jescorr_notruthmatch.root",
    "../jet_calo_analysis_run28_jet10_3sigma_topo_dijetcut_jescorr_notruthmatch_vz_lt_60.root",
    "../jet_calo_analysis_run28_jet10_3sigma_topo_efraccut_jescorr_notruthmatch_vz_lt_60.root",
    "../jet_calo_analysis_run28_jet10_3sigma_topo_nobkgcut_jescorr_notruthmatch_vz_lt_60.root",
    "../jet_calo_analysis_run28_jet10_3sigma_topo_dijetcut_jescorr_notruthmatch_recojet_caloacceptance_vz_lt_60.root",
    "../jet_calo_analysis_run28_jet10_3sigma_topo_efraccut_jescorr_notruthmatch_recojet_caloacceptance_vz_lt_60.root",
    "../jet_calo_analysis_run28_jet10_3sigma_topo_nobkgcut_jescorr_notruthmatch_recojet_caloacceptance_vz_lt_60.root"
]

new_hist_names = ["h_pass_cut_measure_tight", "h_total_measure_tight"]
old_hist_names = ["h_ue_pt_transverse"]
# Separate lists for pass_cut and total histograms
profx_pass_cut = []
profx_total = []
# Open files and get histograms
files = [ROOT.TFile.Open(path) for path in new_files]
hists = [[f.Get(hname) for hname in new_hist_names] for f in files]
# For each file and histogram
for file_idx, (f, hlist) in enumerate(zip(files, hists)):
    for hist_idx, h2d in enumerate(hlist):
        if not h2d:
            print(f"Histogram {hist_names[hist_idx]} not found in file {file_paths[file_idx]}")
            continue
        profx = h2d.ProfileX()
        profx.SetName(f"{h2d.GetName()}_profx_file{file_idx}")
        if hist_idx == 0:
            profx_pass_cut.append(profx)
            profx_pass_cut[-1].SetDirectory(0)
        else:
            profx_total.append(profx)
            profx_total[-1].SetDirectory(0)

old_profx = []
files = [ROOT.TFile.Open(path) for path in old_files]
hists = [[f.Get(hname) for hname in old_hist_names] for f in files]   
for file_idx, (f, hlist) in enumerate(zip(files, hists)):
    for hist_idx, h2d in enumerate(hlist):  
        old_profx.append(h2d)
        old_profx[-1].SetDirectory(0)

colors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+1, ROOT.kBlack, ROOT.kPink+1, ROOT.kRed+3, ROOT.kBlue+3, ROOT.kGreen]
cut_name = ["dijet","efrac","none"]
legend_entries = ['dijet_new','efrac_new','none_new','dijet_vz_30','efrac_vz_30','none_vz_30','dijet_vz_60','efrac_vz_60','none_vz_60','djiet_calo_accept','efrac_calo_accept','none_calo_accept']

profiles = []
profiles.append(profx_pass_cut[0])
profiles.append(profx_pass_cut[1])
profiles.append(profx_total[0])
for h in old_profx:
    profiles.append(h)
'''
profiles_total = []
profiles_total.append(profx_total[0])
profiles_total.append(old_profx[2])
profiles_total.append(old_profx[5])
profiles_total.append(old_profx[8])
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
    profx_hist[i].Divide(profx_hist[i], profx_hist_total[i//3], 1, 1, "B")
'''
for j in range(4):
    canvas = ROOT.TCanvas("c_profx", "ProfileX", 1000, 800)
    legend = ROOT.TLegend(0.4, 0.2, 0.9, 0.4)
    legend.SetTextSize(0.03)
    for i, hist in enumerate(profiles):
        if i // 3 == j:
            hist.SetLineColor(colors[i])
            hist.SetMarkerColor(colors[i])
            drawopt = "E SAME"
            if i % 3 == 0:
                drawopt = "E"
            legend.AddEntry(hist, legend_entries[i], "l")
            hist.GetYaxis().SetRangeUser(0, 2.5)
            hist.GetXaxis().SetRangeUser(20,70)
            print(drawopt)
            hist.Draw(drawopt)
    legend.Draw("SAME")
    canvas.SaveAs(f"profx_old_new_comparison_{j}.png")




