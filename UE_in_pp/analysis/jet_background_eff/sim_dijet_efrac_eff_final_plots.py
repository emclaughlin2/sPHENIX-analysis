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
    #"analysis_sim_run28_output/output_sim_dijet_bkg_cut.root",
    #"analysis_sim_run28_output/output_sim_dijet_bkg_cut_in_acceptance.root"
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
legend_entries = ["Efrac Cut"]#, "Dijet", "Dijet In Acceptance"]  
cut_name = ["efrac"]#,"dijet","dijet_in_acceptance"]

#for i in range(len(projx_pass_cut)):
#    projx_pass_cut[i].Divide(projx_pass_cut[i], projx_total[i], 1, 1, "B")

canvas = ROOT.TCanvas("c_projx", "ProjectionX", 600, 800)
pad1 = ROOT.TPad("pad1", "Top Pad", 0.05, 0.30, 0.95, 1.0)
pad2 = ROOT.TPad("pad2", "Bottom Pad", 0.05, 0.0, 0.95, 0.30)
pad1.Draw()
pad2.Draw()

pad1.cd()
pad1.SetBottomMargin(0.01)
pad1.SetTopMargin(0.05)
pad1.SetLogy(1)
legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
legend.SetTextSize(0.04)
legend.AddEntry("","#bf{#it{sPHENIX}} Internal","")
legend.AddEntry("","Pythia8 200 GeV p+p","")
projx_total[0].GetYaxis().SetTitleFont(43)
projx_total[0].GetYaxis().SetTitleSize(20)
projx_total[0].GetYaxis().SetLabelFont(43)
projx_total[0].GetYaxis().SetLabelSize(20)
projx_total[0].SetLineColor(colors[1])
projx_total[0].SetMarkerColor(colors[1])
projx_pass_cut[0].SetLineColor(colors[0])
projx_pass_cut[0].SetMarkerColor(colors[0])
legend.AddEntry(projx_total[0], "No bkg cut", "lp")
legend.AddEntry(projx_pass_cut[0],"Energy frac bkg cut","lp")
projx_total[0].GetYaxis().SetTitle("Counts")
projx_total[0].GetXaxis().SetLabelSize(0)
projx_total[0].Draw()
projx_pass_cut[0].Draw("same")
legend.Draw("SAME")
pad2.cd()
pad2.SetBottomMargin(0.33)
ratio = projx_pass_cut[0].Clone("ratio")
ratio.Divide(projx_pass_cut[0],projx_total[0],1,1,"B")
ratio.SetLineColor(1)
ratio.SetMarkerColor(1)
ratio.GetYaxis().SetTitleFont(43)
ratio.GetYaxis().SetTitleSize(20)
ratio.GetYaxis().SetLabelFont(43)
ratio.GetYaxis().SetLabelSize(20)
ratio.GetXaxis().SetTitleFont(43)
ratio.GetXaxis().SetTitleSize(20)
ratio.GetXaxis().SetLabelFont(43)
ratio.GetXaxis().SetLabelSize(20)
ratio.GetYaxis().SetRangeUser(0, 1.2)
ratio.GetXaxis().SetTitle("p_{T}^{calib lead} [GeV]")
ratio.GetYaxis().SetTitle("Efficiency")
ratio.Draw()
canvas.SaveAs("projx_efrac_final.png")


