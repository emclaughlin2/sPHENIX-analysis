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

file_paths = [     
    "analysis_sim_run28_output/output_sim_zvtx_lt_10_none_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_zvtx_lt_10_efrac_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_zvtx_lt_10_dijet_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_zvtx_lt_20_none_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_zvtx_lt_20_efrac_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_zvtx_lt_20_dijet_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_zvtx_lt_30_none_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_zvtx_lt_30_efrac_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_zvtx_lt_30_dijet_bkg_cut.root"
]

# Input ROOT files
file_paths2 = [
    "analysis_sim_run28_output/output_sim_zvtx_lt_60_none_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_zvtx_lt_60_efrac_bkg_cut.root",
    "analysis_sim_run28_output/output_sim_zvtx_lt_60_dijet_bkg_cut.root"
]

hist_names = ["h_pass_cut_mbdvertex_truthvertex","h_total_truthvertex"]
hist_names2 = ["h_pass_cut_truth", "h_total_truth"]
# Separate lists for pass_cut and total histograms
projx_pass_cut = []
projx_total = []
profx_pass_cut = []
profx_total = []
projy_pass_cut = []
projy_total = []
all_projy_pass_cut = []
all_projy_total = []

for file_idx, path in enumerate(file_paths):
    file = ROOT.TFile.Open(path)
    file.ls()
    for hist_idx, hname in enumerate(hist_names):
        h2d = ROOT.TH2D(file.Get(hname))
        print(f"Found histogram: {h2d.GetName()} in file: {path}")
        projx = ROOT.TH1D(h2d.ProjectionX())
        projx.SetName(f"{h2d.GetName()}_projx_file{file_idx}")
        profx = ROOT.TProfile(h2d.ProfileX())
        profx.SetName(f"{h2d.GetName()}_profx_file{file_idx}")
        all_projy = ROOT.TH1D(h2d.ProjectionY("projy",2,6))
        all_projy.SetName(f"{h2d.GetName()}_all_projy_file{file_idx}")
        nxbins = h2d.GetNbinsX()
        temp_projy = []
        for xbin in range(2, nxbins):
            temp_projy.append(ROOT.TH1D(h2d.ProjectionY(f"{h2d.GetName()}_projy_xbin{xbin}_file{file_idx}", xbin, xbin)))
            temp_projy[-1].SetDirectory(0)
        if hist_idx == 0:
            projx_pass_cut.append(projx)
            profx_pass_cut.append(profx)
            projy_pass_cut.append(temp_projy)
            all_projy_pass_cut.append(all_projy)
            projx_pass_cut[-1].SetDirectory(0)
            profx_pass_cut[-1].SetDirectory(0)
            all_projy_pass_cut[-1].SetDirectory(0)
        else:
            projx_total.append(projx)
            profx_total.append(profx)
            projy_total.append(temp_projy)
            all_projy_total.append(all_projy)
            projx_total[-1].SetDirectory(0)
            profx_total[-1].SetDirectory(0)
            all_projy_total[-1].SetDirectory(0)

for file_idx, path in enumerate(file_paths2):
    file = ROOT.TFile.Open(path)
    file.ls()
    for hist_idx, hname in enumerate(hist_names2):
        h2d = ROOT.TH2D(file.Get(hname))
        print(f"Found histogram: {h2d.GetName()} in file: {path}")
        projx = ROOT.TH1D(h2d.ProjectionX())
        projx.SetName(f"{h2d.GetName()}_projx_file{file_idx}")
        profx = ROOT.TProfile(h2d.ProfileX())
        profx.SetName(f"{h2d.GetName()}_profx_file{file_idx}")
        all_projy = ROOT.TH1D(h2d.ProjectionY("projy",2,6))
        all_projy.SetName(f"{h2d.GetName()}_all_projy_file{file_idx}")
        nxbins = h2d.GetNbinsX()
        temp_projy = []
        for xbin in range(2, nxbins):
            temp_projy.append(ROOT.TH1D(h2d.ProjectionY(f"{h2d.GetName()}_projy_xbin{xbin}_file{file_idx+9}", xbin, xbin)))
            temp_projy[-1].SetDirectory(0)
        if hist_idx == 0:
            projx_pass_cut.append(projx)
            profx_pass_cut.append(profx)
            projy_pass_cut.append(temp_projy)
            all_projy_pass_cut.append(all_projy)
            projx_pass_cut[-1].SetDirectory(0)
            profx_pass_cut[-1].SetDirectory(0)
            all_projy_pass_cut[-1].SetDirectory(0)
        else:
            projx_total.append(projx)
            profx_total.append(profx)
            projy_total.append(temp_projy)
            all_projy_total.append(all_projy)
            projx_total[-1].SetDirectory(0)
            profx_total[-1].SetDirectory(0)
            all_projy_total[-1].SetDirectory(0)

print(profx_pass_cut, profx_total, projx_pass_cut, projx_total, projy_pass_cut, projy_total)

colors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+1, ROOT.kBlack, ROOT.kPink+1, ROOT.kRed-4, ROOT.kBlue-4, ROOT.kGreen-4]
legend_entries = ["|v_{z,truth}| < 10cm No Cut","|v_{z,truth}| < 10cm Efrac Cut", "|v_{z,truth}| < 10cm Dijet",
                  "|v_{z,truth}| < 20cm No Cut","|v_{z,truth}| < 20cm Efrac Cut", "|v_{z,truth}| < 20cm Dijet",
                  "|v_{z,truth}| < 30cm No Cut","|v_{z,truth}| < 30cm Efrac Cut", "|v_{z,truth}| < 30cm Dijet",
                  "|v_{z,truth}| < 60cm No Cut","|v_{z,truth}| < 60cm Efrac Cut", "|v_{z,truth}| < 60cm Dijet"]  
cut_name = ["vz_lt_10_none","vz_lt_10_efrac","vz_lt_10_dijet",
            "vz_lt_20_none","vz_lt_20_efrac","vz_lt_20_dijet",
            "vz_lt_30_none","vz_lt_30_efrac","vz_lt_30_dijet",
            "vz_lt_60_none","vz_lt_60_efrac","vz_lt_60_dijet"]
ptbin_vals = ['22-27','27-34','34-40','40-52','52-63']

for i in range(len(projx_pass_cut)):
    projx_pass_cut[i].Divide(projx_pass_cut[i], projx_total[i], 1, 1, "B")

first = True
canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
for i, hist in enumerate(projx_pass_cut):
    if i % 3 == 0:
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        hist.SetXTitle("Jet p_{T}^{truth,lead} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        drawopt = "E" if first else "E SAME"
        legend.AddEntry(hist, legend_entries[i], "l")
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.Draw(drawopt)
        first = False
legend.Draw("SAME")
canvas.SaveAs("mc_plots/projx_all_vertex.png")

first = True
canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
for i, hist in enumerate(projx_pass_cut):
    if i < 3:
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        hist.SetXTitle("Jet p_{T}^{truth,lead} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        drawopt = "E" if first else "E SAME"
        legend.AddEntry(hist, legend_entries[i], "l")
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.Draw(drawopt)
        first = False
legend.Draw("SAME")
canvas.SaveAs("mc_plots/projx_vertex_lt_10.png")

first = True
canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
for i, hist in enumerate(projx_pass_cut):
    if i >= 3 and i < 6:
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        hist.SetXTitle("Jet p_{T}^{truth,lead} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        drawopt = "E" if first else "E SAME"
        legend.AddEntry(hist, legend_entries[i], "l")
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.Draw(drawopt)
        first = False
legend.Draw("SAME")
canvas.SaveAs("mc_plots/projx_vertex_lt_20.png")

first = True
canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
for i, hist in enumerate(projx_pass_cut):
    if i >= 6 and i < 9:
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        hist.SetXTitle("Jet p_{T}^{truth,lead} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        drawopt = "E" if first else "E SAME"
        legend.AddEntry(hist, legend_entries[i], "l")
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.Draw(drawopt)
        first = False
legend.Draw("SAME")
canvas.SaveAs("mc_plots/projx_vertex_lt_30.png")

first = True
canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
for i, hist in enumerate(projx_pass_cut):
    if i >= 9:
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        hist.SetXTitle("Jet p_{T}^{truth,lead} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        drawopt = "E" if first else "E SAME"
        legend.AddEntry(hist, legend_entries[i], "l")
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.Draw(drawopt)
        first = False
legend.Draw("SAME")
canvas.SaveAs("mc_plots/projx_vertex_lt_60.png")

for i in range(len(projy_pass_cut)):
    for j in range(len(projy_pass_cut[i])):
        projy_pass_cut[i][j].Divide(projy_pass_cut[i][j], projy_total[i][j], 1, 1, "B")
        # Draw all ProjectionY histograms for each x-bin on the same canvas, one canvas per x-bin

#for i in range(len(projy_pass_cut)):
#    for j in range(len(projy_pass_cut[i])):
#        projy_pass_cut[i][j].Rebin(2)
#        projy_pass_cut[i][j].Scale(1.0/2)

for i in range(len(projy_pass_cut)):
    canvas = ROOT.TCanvas(f"c_projy_xbin{i+1}", f"ProjectionY xbin {i+1}", 800, 600)
    legend = ROOT.TLegend(0.4, 0.6, 0.9, 0.9)
    legend.SetTextSize(0.035)
    legend.AddEntry("",legend_entries[i],"")
    for j in range(len(projy_pass_cut[i])):
        print(i,j)
        hist = projy_pass_cut[i][j]
        hist.SetXTitle("Transverse E_{T}^{truth} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        hist.SetLineColor(colors[j])
        hist.SetMarkerColor(colors[j])
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.GetXaxis().SetRangeUser(0,15)
        drawopt = "E" if j == 0 else "E SAME"
        hist.Draw(drawopt)
        legend.AddEntry(hist, "p_{T}^{jet} = "+ptbin_vals[j]+"GeV", "l")
    legend.Draw("SAME")
    canvas.SaveAs(f"mc_plots/projy_{cut_name[i]}.png")


for i in range(len(projy_pass_cut)):
    for j in range(len(projy_pass_cut[i])):
        #print(i,j,ptbin_vals[j],projx_pass_cut[i].GetBinContent(j+2),projy_pass_cut[i][j].GetBinContent(2))
        projy_pass_cut[i][j].Scale(1.0/projx_pass_cut[i].GetBinContent(j+2))

for i in range(len(projy_pass_cut)):
    canvas = ROOT.TCanvas(f"c_norm_projy_xbin{i+1}", f"ProjectionY xbin {i+1}", 800, 600)
    legend = ROOT.TLegend(0.4, 0.6, 0.9, 0.9)
    legend.SetTextSize(0.035)
    legend.AddEntry("",legend_entries[i],"")
    for j in range(len(projy_pass_cut[i])):
        print(i,j)
        hist = projy_pass_cut[i][j]
        hist.SetXTitle("Transverse E_{T}^{truth} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency (Jet p_{T} Norm)")
        hist.SetLineColor(colors[j])
        hist.SetMarkerColor(colors[j])
        hist.GetYaxis().SetRangeUser(0.6, 1.5)
        hist.GetXaxis().SetRangeUser(0,15)
        drawopt = "E" if j == 0 else "E SAME"
        hist.Draw(drawopt)
        legend.AddEntry(hist, "p_{T}^{jet} = "+ptbin_vals[j]+"GeV", "l")
    legend.Draw("SAME")
    canvas.SaveAs(f"mc_plots/norm_projy_{cut_name[i]}.png")

for i in range(len(all_projy_pass_cut)):
    all_projy_pass_cut[i].Divide(all_projy_pass_cut[i], all_projy_total[i], 1, 1, "B")
    #all_projy_pass_cut[i].Rebin(2)
    #all_projy_pass_cut[i].Scale(1.0/2)

first = True
canvas = ROOT.TCanvas("c_all_projy", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
for i, hist in enumerate(all_projy_pass_cut):
    if i % 3 == 0:
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        hist.SetXTitle("Transverse E_{T}^{truth} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        drawopt = "E" if first else "E SAME"
        legend.AddEntry(hist, legend_entries[i], "l")
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.GetXaxis().SetRangeUser(0,15)
        hist.Draw(drawopt)
        first = False
legend.Draw("SAME")
canvas.SaveAs("mc_plots/all_projy_all_vertex.png")

first = True
canvas = ROOT.TCanvas("c_all_projy", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
for i, hist in enumerate(all_projy_pass_cut):
    if i < 3:
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        hist.SetXTitle("Transverse E_{T}^{truth} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        drawopt = "E" if first else "E SAME"
        legend.AddEntry(hist, legend_entries[i], "l")
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.GetXaxis().SetRangeUser(0,15)
        hist.Draw(drawopt)
        first = False
legend.Draw("SAME")
canvas.SaveAs("mc_plots/all_projy_vertex_lt_10.png")

first = True
canvas = ROOT.TCanvas("c_all_projy", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
for i, hist in enumerate(all_projy_pass_cut):
    if i >= 3 and i < 6:
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        hist.SetXTitle("Transverse E_{T}^{truth} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        drawopt = "E" if first else "E SAME"
        legend.AddEntry(hist, legend_entries[i], "l")
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.GetXaxis().SetRangeUser(0,15)
        hist.Draw(drawopt)
        first = False
legend.Draw("SAME")
canvas.SaveAs("mc_plots/all_projy_vertex_lt_20.png")

first = True
canvas = ROOT.TCanvas("c_all_projy", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
for i, hist in enumerate(all_projy_pass_cut):
    if i >= 6 and i < 9:
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        hist.SetXTitle("Transverse E_{T}^{truth} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        drawopt = "E" if first else "E SAME"
        legend.AddEntry(hist, legend_entries[i], "l")
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.GetXaxis().SetRangeUser(0,15)
        hist.Draw(drawopt)
        first = False
legend.Draw("SAME")
canvas.SaveAs("mc_plots/all_projy_vertex_lt_30.png")

first = True
canvas = ROOT.TCanvas("c_all_projy", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
for i, hist in enumerate(all_projy_pass_cut):
    if i >= 9:
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        hist.SetXTitle("Transverse E_{T}^{truth} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        drawopt = "E" if first else "E SAME"
        legend.AddEntry(hist, legend_entries[i], "l")
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.GetXaxis().SetRangeUser(0,15)
        hist.Draw(drawopt)
        first = False
legend.Draw("SAME")
canvas.SaveAs("mc_plots/all_projy_vertex_lt_60.png")

