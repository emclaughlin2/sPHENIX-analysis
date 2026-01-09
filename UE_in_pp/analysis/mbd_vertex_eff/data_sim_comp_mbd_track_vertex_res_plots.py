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

zvertex_list = [10,20,10,20]
wide_zvertex_list = [10,20,100,100]

for zvertex, wide_zvertex in zip(zvertex_list, wide_zvertex_list):

    # Input ROOT files
    file_paths = [ f"analysis_sim_run28_output/output_sim_zvtx_lt_{zvertex}_{wide_zvertex}_none_bkg_cut.root",
                   f"analysis_data_ana509_output/output_zvtx_lt_{zvertex}_{wide_zvertex}_dijet.root"]

    hist_names = [
        "h_measure_mbd_silicon_diff","h_measure_silicon_mbd_diff",
        "h_track_mult_pass_cut_measure_siliconvertex_mbdvertex","h_track_mult_no_pass_cut_measure_siliconvertex_mbdvertex",
        "h_track_pt_spectra_pass_cut_measure_siliconvertex_mbdvertex","h_track_pt_spectra_no_pass_cut_measure_siliconvertex_mbdvertex"]

    # Separate lists for pass_cut and total histograms
    mbd_track_diff = []
    track_mult = []
    track_spectra = []
    # Open files and get histograms
    files = [ROOT.TFile.Open(path) for path in file_paths]
    hists = [[f.Get(hname) for hname in hist_names] for f in files]
    # For each file and histogram
    for file_idx, (f, hlist) in enumerate(zip(files, hists)):
        for hist_idx, h2d in enumerate(hlist):
            if not h2d:
                print(f"Histogram {hist_names[hist_idx]} not found in file {file_paths[file_idx]}")
                continue
            h2d.Scale(1.0/h2d.Integral())
            if hist_idx < 2:
                mbd_track_diff.append(h2d)
                mbd_track_diff[-1].SetDirectory(0)
            elif hist_idx < 4:
                track_mult.append(h2d)
                track_mult[-1].SetDirectory(0)
            else:
                track_spectra.append(h2d)
                track_spectra[-1].SetDirectory(0)

    colors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kRed+2, ROOT.kBlue+2]
    legend_entries = ["MC MBD/silicon","MC silicon/MBD","Data MBD/silicon","Data silicon/MBD"]  

    for i, hist in enumerate(mbd_track_diff):
        hist.Fit("gaus","","",-5,5)
        fit = hist.GetFunction("gaus")
        fit.SetLineColor(colors[i])
        legend_entries[i] = legend_entries[i] + ": #mu = "+f"{fit.GetParameter(1):.2f} "+"#sigma = "+f"{fit.GetParameter(2):.2f}"

    canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
    legend = ROOT.TLegend(0.17, 0.7, 0.5, 0.9)
    legend.SetTextSize(0.03)
    legend.AddEntry("","|z_num| < "+f"{wide_zvertex}cm,"+" |z_dem|"+f" < {zvertex}cm","")
    for i, hist in enumerate(mbd_track_diff):
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        if i < 2:
            hist.SetMarkerStyle(24)
        if i > 2:
            hist.SetMarkerStyle(20)
        drawopt = "E" if i == 0 else "E SAME"
        legend.AddEntry(hist, legend_entries[i], "lp")
        hist.GetYaxis().SetRangeUser(0,0.06)
        hist.GetXaxis().SetRangeUser(-10,10)
        hist.Draw(drawopt)
    legend.Draw("SAME")
    canvas.SaveAs(f"data_sim_mbd_silicon_vz_res_lt_{zvertex}_{wide_zvertex}.png")

    track_legends = ['MC silicon Vertex Reco','MC No silicon Vertex Reco','Data silicon Vertex Reco','Data No silicon Vertex Reco']

    canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
    #canvas.SetLogy(1)
    legend = ROOT.TLegend(0.55, 0.7, 0.9, 0.9)
    legend.SetTextSize(0.03)
    legend.AddEntry("","|z_{silicon}| < "+f"{wide_zvertex}cm,"+" |z_{MBD}|"+f" < {zvertex}cm","")
    for i, hist in enumerate(track_mult):
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        if i < 2:
            hist.SetMarkerStyle(24)
        if i > 2:
            hist.SetMarkerStyle(20)
        drawopt = "HIST" if i == 0 else "HIST SAME"
        legend.AddEntry(hist, track_legends[i], "l")
        hist.GetXaxis().SetRangeUser(0,30)
        hist.GetYaxis().SetRangeUser(0.0,1.0)
        hist.GetXaxis().SetTitle("N_{trk}")
        hist.Draw(drawopt)
    legend.Draw("SAME")
    canvas.SaveAs(f"data_silicon_mult_lt_{zvertex}_{wide_zvertex}.png")

    canvas = ROOT.TCanvas("c_projx", "ProjectionX", 800, 600)
    canvas.SetLogy(1)
    legend = ROOT.TLegend(0.55, 0.7, 0.9, 0.9)
    legend.SetTextSize(0.03)
    legend.AddEntry("","|z_{silicon}| < "+f"{wide_zvertex}cm,"+" |z_{MBD}|"+f" < {zvertex}cm","")
    for i, hist in enumerate(track_spectra):
        hist.SetLineColor(colors[i])
        hist.SetMarkerColor(colors[i])
        if i < 2:
            hist.SetMarkerStyle(24)
        if i > 2:
            hist.SetMarkerStyle(20)
        drawopt = "HIST" if i == 0 else "HIST SAME"
        legend.AddEntry(hist, track_legends[i], "l")
        hist.GetXaxis().SetRangeUser(0,30)
        hist.GetXaxis().SetTitle("p_{T,trk} [GeV]")
        hist.Draw(drawopt)
    legend.Draw("SAME")
    canvas.SaveAs(f"data_silicon_spectra_lt_{zvertex}_{wide_zvertex}.png")
