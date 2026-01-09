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

sim = True
zvertex_list = [10,20]
bkg_cut_list = ["efrac","dijet","none"]

for zvertex in zvertex_list:
    for bkg_cut in bkg_cut_list:

        # Input ROOT files
        if sim:
            file_paths = [ f"analysis_sim_run28_output/output_sim_zvtx_lt_{zvertex}_100_{bkg_cut}_bkg_cut.root" ]
        else: 
            file_paths = [ f"analysis_data_ana509_output/output_zvtx_lt_{zvertex}_100_{bkg_cut}.root" ]

        if sim:
            hist_names = [
            "h_pass_cut_measure_mbdvertex_trackvertex","h_total_measure_trackvertex",
            "h_pass_cut_measure_trackvertex_mbdvertex","h_total_measure_mbdvertex",
            "h_pass_cut_mbdvertex_truthvertex","h_total_truthvertex",
            "h_pass_cut_trackvertex_truthvertex","h_total_truthvertex",
            "h_pass_cut_mbdvertex_trackvertex","h_total_trackvertex",
            "h_pass_cut_trackvertex_mbdvertex","h_total_mbdvertex"]
        else:
            hist_names = [
            "h_pass_cut_measure_mbdvertex_trackvertex","h_total_measure_trackvertex",
            "h_pass_cut_measure_trackvertex_mbdvertex","h_total_measure_mbdvertex"]

        # Separate lists for pass_cut and total histograms
        projx_pass_cut = []
        projx_total = []
        profx_pass_cut = []
        profx_total = []
        projy_pass_cut = []
        projy_total = []
        all_projy_pass_cut = []
        all_projy_total = []
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
                projy = h2d.ProjectionY()
                projy.SetName(f"{h2d.GetName()}_projy_file{file_idx}")
                nxbins = h2d.GetNbinsX()
                temp_projy = []
                for xbin in range(1, nxbins + 1, 2):
                    temp_projy.append(h2d.ProjectionY(f"{h2d.GetName()}_projy_xbin{xbin}_file{file_idx}", xbin, xbin))
                if hist_idx % 2 == 0:
                    projx_pass_cut.append(projx)
                    profx_pass_cut.append(profx)
                    projy_pass_cut.append(temp_projy)
                    all_projy_pass_cut.append(projy)
                else:
                    projx_total.append(projx)
                    profx_total.append(profx)
                    projy_total.append(temp_projy)
                    all_projy_total.append(projy)

        print([prof.GetName() for prof in profx_pass_cut])
        print([prof.GetName() for prof in profx_total])

        if sim:
            colors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1, ROOT.kOrange-2, ROOT.kCyan+1, ROOT.kGray+1, ROOT.kBlack, ROOT.kPink+1]
        else: 
            colors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1, ROOT.kOrange+7, ROOT.kCyan+1, ROOT.kGray+1, ROOT.kBlack, ROOT.kPink+1]
        if sim:
            legend_entries = ["Measure MBD/Track","Measure Track/MBD","MBD/Truth","Track/Truth", "MBD/Track","Track/MBD"]  
            cut_name = ["mbd_truth","track_truth","mbd_track","track_mbd","measure_mbd_track","measure_track_mbd"]
        else:
            legend_entries = ["Measure MBD/Track","Measure Track/MBD"]  
            cut_name = ["measure_mbd_track","measure_track_mbd"]

        for i in range(len(projx_pass_cut)):
            projx_pass_cut[i].Divide(projx_pass_cut[i], projx_total[i], 1, 1, "B")

        for i in range(len(all_projy_pass_cut)):
            all_projy_pass_cut[i].Divide(all_projy_pass_cut[i], all_projy_total[i], 1, 1, "B")
        '''
        canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
        legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
        legend.SetTextSize(0.03)
        for i, hist in enumerate(projx_pass_cut):
            hist.SetLineColor(colors[i])
            hist.SetMarkerColor(colors[i])
            hist.SetXTitle("Jet p_{T}^{lead} [GeV]")
            hist.SetYTitle("Vertex Efficiency")
            drawopt = "E" if i == 0 else "E SAME"
            legend.AddEntry(hist, legend_entries[i], "l")
            hist.GetYaxis().SetRangeUser(0, 1.5)
            hist.Draw(drawopt)
        legend.Draw("SAME")
        if sim:
            canvas.SaveAs(f"sim_projx_lt_{zvertex}_100_{bkg_cut}.png")
        else: 
            canvas.SaveAs(f"data_projx_lt_{zvertex}_100_{bkg_cut}.png")
        '''
        canvas = ROOT.TCanvas("c_projy", "ProjectionX", 800, 600)
        legend = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
        legend.SetTextSize(0.03)
        for i, hist in enumerate(all_projy_pass_cut):
            hist.SetLineColor(colors[i])
            hist.SetMarkerColor(colors[i])
            hist.SetXTitle("Transverse Region #SigmaE_{T} [GeV]")
            hist.SetYTitle("Vertex Efficiency")
            drawopt = "E" if i == 0 else "E SAME"
            legend.AddEntry(hist, legend_entries[i], "l")
            hist.GetYaxis().SetRangeUser(0, 1.5)
            hist.Draw(drawopt)
        legend.Draw("SAME")
        if sim:
            canvas.SaveAs(f"sim_projy_lt_{zvertex}_100_{bkg_cut}.png")
        else: 
            canvas.SaveAs(f"data_projy_lt_{zvertex}_100_{bkg_cut}.png")
