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
    "analysis_sim_run28_output/output_sim_zvtx_lt_60_60_none_bkg_cut.root"
]

hist_names = ["h_pass_cut_mbdvertex_truthvertex","h_total_truthvertex"]
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
        all_projy = ROOT.TH1D(h2d.ProjectionY("projy",2,5))
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

print(profx_pass_cut, profx_total, projx_pass_cut, projx_total, projy_pass_cut, projy_total)

colors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+1, ROOT.kBlack, ROOT.kPink+1, ROOT.kRed-4, ROOT.kBlue-4, ROOT.kGreen-4]
legend_entries = ["|v_{z,truth}| < 60cm"]  
cut_name = ["vz_lt_60_none"]
ptbin_vals = ['21-26','26-32.5','32.5-40.5','40.5-63.5']

for i in range(len(projx_pass_cut)):
    projx_pass_cut[i].Divide(projx_pass_cut[i], projx_total[i], 1, 1, "B")

first = True
canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
legend.AddEntry("","#bf{#it{sPHENIX}} Internal","")
legend.AddEntry("","Pythia8 200 GeV p+p","")
for i, hist in enumerate(projx_pass_cut):
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
canvas.SaveAs("new_mc_plots/projx_all_vertex.png")

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
    legend = ROOT.TLegend(0.2, 0.7, 0.9, 0.9)
    legend.SetTextSize(0.035)
    legend.SetNColumns(2)
    legend.AddEntry("","#bf{#it{sPHENIX}} Internal","")
    legend.AddEntry("","","")
    legend.AddEntry("","Pythia8 200 GeV p+p","")
    legend.AddEntry("",legend_entries[i],"")
    for j in range(len(projy_pass_cut[i])):
        print(i,j)
        hist = projy_pass_cut[i][j]
        hist.SetXTitle("Transverse #Sigma E_{T}^{truth} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency")
        hist.SetLineColor(colors[j])
        hist.SetMarkerColor(colors[j])
        hist.GetYaxis().SetRangeUser(0, 1.2)
        hist.GetXaxis().SetRangeUser(0,15)
        drawopt = "E" if j == 0 else "E SAME"
        hist.Draw(drawopt)
        legend.AddEntry(hist, "p_{T}^{jet} = "+ptbin_vals[j]+"GeV", "l")
    legend.Draw("SAME")
    canvas.SaveAs(f"new_mc_plots/projy_{cut_name[i]}.png")


for i in range(len(projy_pass_cut)):
    for j in range(len(projy_pass_cut[i])):
        #print(i,j,ptbin_vals[j],projx_pass_cut[i].GetBinContent(j+2),projy_pass_cut[i][j].GetBinContent(2))
        projy_pass_cut[i][j].Scale(1.0/projx_pass_cut[i].GetBinContent(j+2))

for i in range(len(projy_pass_cut)):
    canvas = ROOT.TCanvas(f"c_norm_projy_xbin{i+1}", f"ProjectionY xbin {i+1}", 800, 600)
    legend = ROOT.TLegend(0.2, 0.7, 0.9, 0.9)
    legend.SetTextSize(0.035)
    legend.SetNColumns(2)
    legend.AddEntry("","#bf{#it{sPHENIX}} Internal","")
    legend.AddEntry("","","")
    legend.AddEntry("","Pythia8 200 GeV p+p","")
    legend.AddEntry("",legend_entries[i],"")
    for j in range(len(projy_pass_cut[i])):
        print(i,j)
        hist = projy_pass_cut[i][j]
        hist.SetXTitle("Transverse #Sigma E_{T}^{truth} [GeV]")
        hist.SetYTitle("MBD Vertex Efficiency (Jet p_{T} Norm)")
        hist.SetLineColor(colors[j])
        hist.SetMarkerColor(colors[j])
        hist.GetYaxis().SetRangeUser(0.6, 1.5)
        hist.GetXaxis().SetRangeUser(0,15)
        drawopt = "E" if j == 0 else "E SAME"
        hist.Draw(drawopt)
        legend.AddEntry(hist, "p_{T}^{jet} = "+ptbin_vals[j]+"GeV", "l")
    legend.Draw("SAME")
    canvas.SaveAs(f"new_mc_plots/norm_projy_{cut_name[i]}.png")

for i in range(len(all_projy_pass_cut)):
    all_projy_pass_cut[i].Divide(all_projy_pass_cut[i], all_projy_total[i], 1, 1, "B")
    #all_projy_pass_cut[i].Rebin(2)
    #all_projy_pass_cut[i].Scale(1.0/2)

first = True
canvas = ROOT.TCanvas("c_all_projy", "ProjectionX", 800, 600)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
legend.AddEntry("","#bf{#it{sPHENIX}} Internal","")
legend.AddEntry("","Pythia8 200 GeV p+p","")
for i, hist in enumerate(all_projy_pass_cut):
    hist.SetLineColor(colors[i])
    hist.SetMarkerColor(colors[i])
    hist.SetXTitle("Transverse #Sigma E_{T}^{truth} [GeV]")
    hist.SetYTitle("MBD Vertex Efficiency")
    drawopt = "E" if first else "E SAME"
    legend.AddEntry(hist, legend_entries[i], "l")
    hist.GetYaxis().SetRangeUser(0, 1.2)
    hist.GetXaxis().SetRangeUser(0,15)
    hist.Draw(drawopt)
    first = False
legend.Draw("SAME")
canvas.SaveAs("new_mc_plots/all_projy_all_vertex.png")
