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

zvertex_list = [10,20]

for zvertex in zvertex_list:

        # Input ROOT files
        file_paths = [ f"analysis_sim_run28_output/output_sim_zvtx_lt_{zvertex}_100_none_bkg_cut.root",
                       f"analysis_data_ana509_output/output_zvtx_lt_{zvertex}_100_efrac.root" ]

        hist_names = [
            "h_pass_cut_measure_mbdvertex_siliconvertex","h_total_measure_siliconvertex",
            "h_pass_cut_measure_siliconvertex_mbdvertex","h_total_measure_mbdvertex"]

        # Separate lists for pass_cut and total histograms
        projx_pass_cut = []
        projx_total = []
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
                projy = h2d.ProjectionY()
                projy.SetName(f"{h2d.GetName()}_projy_file{file_idx}")
                if hist_idx % 2 == 0:
                    projx_pass_cut.append(projx)
                    projy_pass_cut.append(projy)
                else:
                    projx_total.append(projx)
                    projy_total.append(projy)

        colors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kRed+2, ROOT.kBlue+2]
        legend_entries = ["MC MBD/silicon","MC silicon/MBD","Data MBD/silicon","Data silicon/MBD"]  

        for i in range(len(projx_pass_cut)):
            projx_pass_cut[i].Divide(projx_pass_cut[i], projx_total[i], 1, 1, "B")

        for i in range(len(projy_pass_cut)):
            projy_pass_cut[i].Divide(projy_pass_cut[i], projy_total[i], 1, 1, "B")

        canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
        legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
        legend.SetTextSize(0.03)
        legend.AddEntry("","|z_{silicon}| < "+f"{zvertex}cm,"+" |z_{MBD}|"+f" < 100cm","")
        for i, hist in enumerate(projx_pass_cut):
            if i%2 == 0:
                hist.SetLineColor(colors[i])
                hist.SetMarkerColor(colors[i])
                if i < 2:
                    hist.SetMarkerStyle(24)
                if i > 2:
                    hist.SetMarkerStyle(20)
                hist.SetXTitle("Jet p_{T}^{lead} [GeV]")
                hist.SetYTitle("Vertex Efficiency")
                drawopt = "E" if i == 0 else "E SAME"
                legend.AddEntry(hist, legend_entries[i], "lp")
                hist.GetYaxis().SetRangeUser(0, 1.5)
                hist.Draw(drawopt)
            legend.Draw("SAME")
        canvas.SaveAs(f"data_sim_mbd_silicon_comp_projx_lt_{zvertex}_100.png")

        canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
        legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
        legend.SetTextSize(0.03)
        legend.AddEntry("","|z_{silicon}| < "+f"100cm,"+" |z_{MBD}|"+f" < {zvertex}cm","")
        for i, hist in enumerate(projx_pass_cut):
            if i%2 == 1:
                hist.SetLineColor(colors[i])
                hist.SetMarkerColor(colors[i])
                if i < 2:
                    hist.SetMarkerStyle(24)
                if i > 2:
                    hist.SetMarkerStyle(20)
                hist.SetXTitle("Jet p_{T}^{lead} [GeV]")
                hist.SetYTitle("Vertex Efficiency")
                drawopt = "E" if i == 1 else "E SAME"
                legend.AddEntry(hist, legend_entries[i], "lp")
                hist.GetYaxis().SetRangeUser(0, 1.5)
                hist.Draw(drawopt)
            legend.Draw("SAME")
        canvas.SaveAs(f"data_sim_silicon_mbd_comp_projx_lt_{zvertex}_100.png")

        canvas = ROOT.TCanvas("c_projy", "ProjectionY", 800, 600)
        legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
        legend.SetTextSize(0.03)
        legend.AddEntry("","|z_{silicon}| < "+f"{zvertex}cm,"+" |z_{MBD}|"+f" < 100cm","")
        for i, hist in enumerate(projy_pass_cut):
            if i%2 == 0:
                hist.SetLineColor(colors[i])
                hist.SetMarkerColor(colors[i])
                if i < 2:
                    hist.SetMarkerStyle(24)
                if i > 2:
                    hist.SetMarkerStyle(20)
                hist.SetXTitle("Transverse Region #SigmaE_{T} [GeV]")
                hist.SetYTitle("Vertex Efficiency")
                drawopt = "E" if i == 0 else "E SAME"
                legend.AddEntry(hist, legend_entries[i], "lp")
                hist.GetYaxis().SetRangeUser(0, 1.5)
                hist.Draw(drawopt)
            legend.Draw("SAME")
        canvas.SaveAs(f"data_sim_mbd_silicon_comp_projy_lt_{zvertex}_100.png")

        canvas = ROOT.TCanvas("c_projy", "ProjectionY", 800, 600)
        legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
        legend.SetTextSize(0.03)
        legend.AddEntry("","|z_{silicon}| < "+f"100cm,"+" |z_{MBD}|"+f" < {zvertex}cm","")
        for i, hist in enumerate(projy_pass_cut):
            if i%2 == 1:
                hist.SetLineColor(colors[i])
                hist.SetMarkerColor(colors[i])
                if i < 2:
                    hist.SetMarkerStyle(24)
                if i > 2:
                    hist.SetMarkerStyle(20)
                hist.SetXTitle("Transverse Region #SigmaE_{T} [GeV]")
                hist.SetYTitle("Vertex Efficiency")
                drawopt = "E" if i == 1 else "E SAME"
                legend.AddEntry(hist, legend_entries[i], "lp")
                hist.GetYaxis().SetRangeUser(0, 1.5)
                hist.Draw(drawopt)
            legend.Draw("SAME")
        canvas.SaveAs(f"data_sim_silicon_mbd_comp_projy_lt_{zvertex}_100.png")
