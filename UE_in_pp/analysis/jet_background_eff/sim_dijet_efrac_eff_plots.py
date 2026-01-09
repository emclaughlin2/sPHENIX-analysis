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

# Input ROOT files
file_paths = [
    "analysis_sim_run28_output/output_sim_efrac_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_dijet_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_dijet_bkg_cut_in_acceptance.root"
]

hist_names = ["h_pass_cut_measure_tight", "h_total_measure_tight"]
# Separate lists for pass_cut and total histograms
projx_pass_cut = []
projx_total = []
profx_pass_cut = []
profx_total = []
projy_pass_cut = []
projy_total = []
# Open files and get histograms
files = [ROOT.TFile.Open(path) for path in file_paths]
hists = [[f.Get(hname) for hname in hist_names] for f in files]
# For each file and histogram
for file_idx, (f, hlist) in enumerate(zip(files, hists)):
    for hist_idx, h2d in enumerate(hlist):
        if not h2d:
            print(f"Histogram {hist_names[hist_idx]} not found in file {file_paths[file_idx]}")
            continue

        #print(f"\nProcessing {hist_names[hist_idx]} from {file_paths[file_idx]}")
        projx = h2d.ProjectionX()
        projx.SetName(f"{h2d.GetName()}_projx_file{file_idx}")
        profx = h2d.ProfileX()
        profx.SetName(f"{h2d.GetName()}_profx_file{file_idx}")
        nxbins = h2d.GetNbinsX()
        temp_projy = []
        for xbin in range(1, nxbins + 1, 2):
            temp_projy.append(h2d.ProjectionY(f"{h2d.GetName()}_projy_xbin{xbin}_file{file_idx}", xbin, xbin))
        if hist_idx == 0:
            projx_pass_cut.append(projx)
            profx_pass_cut.append(profx)
            projy_pass_cut.append(temp_projy)
        else:
            projx_total.append(projx)
            profx_total.append(profx)
            projy_total.append(temp_projy)
print(profx_pass_cut, profx_total, projx_pass_cut, projx_total, projy_pass_cut, projy_total)          

colors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+1, ROOT.kBlack, ROOT.kPink+1]
legend_entries = ["Efrac Cut", "Dijet", "Dijet In Acceptance"]  
cut_name = ["efrac","dijet","dijet_in_acceptance"]

profx_hist_pass_cut = []
profx_hist_total = []

for prof in profx_pass_cut:
    nbins = prof.GetNbinsX()
    xaxis = prof.GetXaxis()
    h1 = ROOT.TH1D(f"{prof.GetName()}_hist", prof.GetTitle(), nbins, xaxis.GetXmin(), xaxis.GetXmax())
    for ibin in range(1, nbins + 1):
        h1.SetBinContent(ibin, prof.GetBinContent(ibin))
        h1.SetBinError(ibin, prof.GetBinError(ibin))
    profx_hist_pass_cut.append(h1)

for prof in profx_total:
    nbins = prof.GetNbinsX()
    xaxis = prof.GetXaxis()
    h1 = ROOT.TH1D(f"{prof.GetName()}_hist", prof.GetTitle(), nbins, xaxis.GetXmin(), xaxis.GetXmax())
    for ibin in range(1, nbins + 1):
        h1.SetBinContent(ibin, prof.GetBinContent(ibin))
        h1.SetBinError(ibin, prof.GetBinError(ibin))
    profx_hist_total.append(h1)

for i in range(len(profx_hist_pass_cut)):
    profx_hist_pass_cut[i].Divide(profx_hist_pass_cut[i], profx_hist_total[i], 1, 1, "B")

canvas = ROOT.TCanvas("c_profx", "ProfileX", 800, 600)
legend = ROOT.TLegend(0.4, 0.2, 0.9, 0.4)
for i, hist in enumerate(profx_hist_pass_cut):
    hist.SetLineColor(colors[i])
    hist.SetMarkerColor(colors[i])
    drawopt = "E" if i == 0 else "E SAME"
    legend.AddEntry(hist, legend_entries[i], "l")
    hist.GetYaxis().SetRangeUser(0, 1.2)
    hist.Draw(drawopt)
legend.Draw("SAME")
canvas.SaveAs("profx.png")

for i in range(len(projx_pass_cut)):
    projx_pass_cut[i].Divide(projx_pass_cut[i], projx_total[i], 1, 1, "B")

canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.2, 0.9, 0.4)
for i, hist in enumerate(projx_pass_cut):
    hist.SetLineColor(colors[i])
    hist.SetMarkerColor(colors[i])
    drawopt = "E" if i == 0 else "E SAME"
    legend.AddEntry(hist, legend_entries[i], "l")
    hist.GetYaxis().SetRangeUser(0, 1.2)
    hist.Draw(drawopt)
legend.Draw("SAME")
canvas.SaveAs("projx.png")

for i in range(len(projy_pass_cut)):
    for j in range(len(projy_pass_cut[i])):
        projy_pass_cut[i][j].Divide(projy_pass_cut[i][j], projy_total[i][j], 1, 1, "B")
        # Draw all ProjectionY histograms for each x-bin on the same canvas, one canvas per x-bin

for i in range(len(projy_pass_cut)):
    canvas = ROOT.TCanvas(f"c_projy_xbin{i+1}", f"ProjectionY xbin {i+1}", 800, 600)
    legend = ROOT.TLegend(0.4, 0.2, 0.9, 0.4)
    for j in range(len(projy_pass_cut[i])):
        hist = projy_pass_cut[i][j]
        hist.SetLineColor(colors[j])
        hist.SetMarkerColor(colors[j])
        hist.GetYaxis().SetRangeUser(0, 1.2)
        drawopt = "E" if j == 0 else "E SAME"
        hist.Draw(drawopt)
        legend_label = f"xbin {j}"
        legend.AddEntry(hist, legend_label, "l")
        first_drawn = True
    legend.Draw("SAME")
    canvas.SaveAs(f"projy_{cut_name[i]}.png")