import ROOT
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1I, TH1F, TH2F, TH3F, TColor, TEfficiency
from ROOT import gROOT, gBenchmark, gRandom, gSystem
import numpy as np
import pdb
from array import array

gROOT.LoadMacro("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C")
gROOT.ProcessLine("SetsPhenixStyle()")

topo_thres = ['-9999','0','100','200','300','500']
thres_string = ['All E_{topo}','E_{topo} > 0 MeV','E_{topo} > 100 MeV','E_{topo} > 200 MeV','E_{topo} > 300 MeV','E_{topo} > 500 MeV']

leg_tags = ['0-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9-10']
rgb = [[230, 25, 75], [60, 180, 75], [255, 225, 25], [0, 130, 200], [245, 130, 48], [145, 30, 180], [70, 240, 240], [240, 50, 230], [210, 245, 60], [250, 190, 212], [0, 128, 128], [220, 190, 255], [170, 110, 40], [128, 128, 128], [128, 0, 0], [0, 0, 0], [128, 128, 0], [255, 215, 180], [0, 0, 128], [34, 139, 34]]
colors = [TColor.GetColor(rgb[i][0],rgb[i][1],rgb[i][2]) for i in range(len(rgb))]

mb_files = ['minbias_analysis_output/minbias_analysis_run22.root',
         'minbias_analysis_output/minbias_analysis_1.5mrad_ana468_3sigma_topo.root']

files = ['pileup_0mrad_and_1.5mrad_UE_analysis_leadjet_15_20_GeV_zvtx_lt_10cm_no_low_ET_towers_rate_range_0.000000_0.020000.root',
                'pileup_0mrad_and_1.5mrad_UE_analysis_leadjet_15_20_GeV_zvtx_lt_10cm_no_low_ET_towers_rate_range_0.020000_0.030000.root',
                'pileup_0mrad_and_1.5mrad_UE_analysis_leadjet_15_20_GeV_zvtx_lt_10cm_no_low_ET_towers_rate_range_0.030000_0.040000.root',
                'pileup_0mrad_and_1.5mrad_UE_analysis_leadjet_15_20_GeV_zvtx_lt_10cm_no_low_ET_towers_rate_range_0.040000_0.050000.root',
                'pileup_0mrad_and_1.5mrad_UE_analysis_leadjet_15_20_GeV_zvtx_lt_10cm_no_low_ET_towers_rate_range_0.050000_0.060000.root',
                'pileup_0mrad_and_1.5mrad_UE_analysis_leadjet_15_20_GeV_zvtx_lt_10cm_no_low_ET_towers_rate_range_0.060000_0.070000.root',
                'pileup_0mrad_and_1.5mrad_UE_analysis_leadjet_15_20_GeV_zvtx_lt_10cm_no_low_ET_towers_rate_range_0.070000_0.080000.root',
                'pileup_0mrad_and_1.5mrad_UE_analysis_leadjet_15_20_GeV_zvtx_lt_10cm_no_low_ET_towers_rate_range_0.080000_0.090000.root',
                'pileup_0mrad_and_1.5mrad_UE_analysis_leadjet_15_20_GeV_zvtx_lt_10cm_no_low_ET_towers_rate_range_0.090000_0.100000.root']

h_ue_total = []
h_ntopo_total = []
h_topo_total = []

mb_leg_tags = ['Pythia8 MB sim','1.5mrad MB data']

for i, file in enumerate(mb_files):
    f2 = ROOT.TFile.Open(file)
    events = f2.Get("h_vz").GetEntries()
    h_ue_total.append(TH1F(f2.Get("h_et_total")))
    h_ue_total[i].Rebin(2)
    h_ue_total[i].Scale(1.0/h_ue_total[i].Integral())
    h_ue_total[i].SetDirectory(0)
    ntopo_total = []
    topo_total = []
    for j, t in enumerate(topo_thres):
        ntopo_total.append(f2.Get('h_ntopo'+t+'_total'))
        ntopo_total[j].Scale(1.0/ntopo_total[j].Integral())
        ntopo_total[j].SetDirectory(0)
        topo_total.append(f2.Get('h_topo'+t+'_total'))
        topo_total[j].Rebin(2)
        topo_total[j].Scale(1.0/events)
        topo_total[j].SetDirectory(0)
    h_ntopo_total.append(ntopo_total)
    h_topo_total.append(topo_total)
    f2.Close()
    

c = ROOT.TCanvas("c", "", 800, 600)
leg = ROOT.TLegend(.3,.78,.92,.92)
for i, h in enumerate(h_ue_total):
    h.SetLineColor(i+1)
    h.SetMarkerColor(i+1)
    leg.AddEntry(h,mb_leg_tags[i]+" mean = "+str(h.GetMean()),"pl")
    if i == 0:
        h.GetXaxis().SetRangeUser(-1,15)
        h.SetXTitle("#SigmaE_{T} [GeV]")
        h.Draw()
    else:
        h.Draw("same")
leg.Draw("same")
c.Draw()
c.SaveAs("minbias_ET.png")

'''
tm, te = [], [], [], [], [], []
for i in range(len(h_ue_total)):
    tm.append(h_ue_total[i].GetMean())
    te.append(h_ue_total[i].GetMeanError())

tpd = []
tpe = []
trpd = []
trpe = []
apd = []
ape = []

for i in range(1, len(tm)):
    tpd.append((tm[i] - tm[0]) /tm[0] * 100)
    tpe.append(100 * np.sqrt((te[i] / tm[0]) ** 2 + ((tm[i] * te[0]) / (tm[0] ** 2)) ** 2))
    trpd.append((trm[i] - trm[0]) /trm[0] * 100)
    trpe.append(100 * np.sqrt((tre[i] / trm[0]) ** 2 + ((trm[i] * tre[0]) / (trm[0] ** 2)) ** 2))
    apd.append((am[i] - am[0]) /am[0] * 100)
    ape.append(100 * np.sqrt((ae[i] / am[0]) ** 2 + ((am[i] * ae[0]) / (am[0] ** 2)) ** 2))

for i in range(len(tm)):
    print(f'{tm[i]:.3f} +/- {te[i]:.3f}')
print()
for i in range(len(trm)):
    print(f'{trm[i]:.3f} +/- {tre[i]:.3f}')
print()
for i in range(len(am)):
    print(f'{am[i]:.3f} +/- {ae[i]:.3f}')
print()
for i in range(len(tpd)):
    print(f'{tpd[i]:.3f} +/- {tpe[i]:.3f}')
print()
for i in range(len(trpd)):
    print(f'{trpd[i]:.3f} +/- {trpe[i]:.3f}')
print()
for i in range(len(apd)):
    print(f'{apd[i]:.3f} +/- {ape[i]:.3f}')

x_vals = array('f', [1.0, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
n_points = len(x_vals)

tm_arr = array('f', tm)
trm_arr = array('f', trm)
am_arr = array('f', am)
te_arr = array('f', te)
tre_arr = array('f', tre)
ae_arr = array('f', ae)

x_err = array('f', [0]*n_points)  # No x errors

# Draw the graphs
c = ROOT.TCanvas("c8", "TGraphErrors Example", 800, 600)

g_tm = ROOT.TGraphErrors(n_points, x_vals, tm_arr, x_err, te_arr)
g_trm = ROOT.TGraphErrors(n_points, x_vals, trm_arr, x_err, tre_arr)
g_am = ROOT.TGraphErrors(n_points, x_vals, am_arr, x_err, ae_arr)

# Style the graphs
g_tm.SetMarkerStyle(20)
g_tm.SetMarkerColor(ROOT.kBlue)
g_tm.SetLineColor(ROOT.kBlue)

g_trm.SetMarkerStyle(20)
g_trm.SetMarkerColor(ROOT.kRed)
g_trm.SetLineColor(ROOT.kRed)

g_am.SetMarkerStyle(20)
g_am.SetMarkerColor(ROOT.kGreen+2)
g_am.SetLineColor(ROOT.kGreen+2)

#g_tm.GetYaxis().SetRangeUser(3.74,3.8)
g_tm.GetXaxis().SetTitle("<PU> [%]")
g_tm.GetYaxis().SetTitle("#SigmaE_{T}/#delta#eta#delta#phi [GeV]") 
g_tm.GetYaxis().SetTitleOffset(1.6)
g_tm.Draw("AP")

# Add legend
legend = ROOT.TLegend(0.17, 0.8, 0.47, 0.92)
legend.AddEntry(g_tm, "Towards Region 0mrad + 1.5mrad", "lp")
legend.Draw()

c.Update()
c.Draw()
c.SaveAs("0mrad_plots/h_ue_pileup_towards_0mrad_and_1.5mrad_15_20GeV_Topoclusters_no_low_ET_towers.png")
# Draw the graphs
c = ROOT.TCanvas("c6", "TGraphErrors Example", 800, 600)

f_linear = ROOT.TF1("f_linear", "[0] + [1]*x", 0, 10)
f_linear.SetParameters(0, 0.004)  # Optional: initial guess (intercept=0, slope=1)
g_trm.Fit(f_linear, "R")  # "R" restricts fit to function range

g_trm.GetXaxis().SetTitle("<PU> [%]")
g_trm.GetYaxis().SetTitle("#SigmaE_{T}/#delta#eta#delta#phi [GeV]") 
g_trm.GetYaxis().SetTitleOffset(1.6)
g_trm.Draw("AP")
f_linear.Draw("same")

# Add legend
legend = ROOT.TLegend(0.17, 0.8, 0.47, 0.92)
legend.AddEntry(g_trm, "Transverse Region 0mrad + 1.5mrad", "lp")
legend.Draw()

c.Update()
c.Draw()
c.SaveAs("0mrad_plots/h_ue_pileup_transverse_0mrad_and_1.5mrad_15_20GeV_Topoclusters_no_low_ET_towers.png")

# Draw the graphs
c = ROOT.TCanvas("c11", "TGraphErrors Example", 800, 600)

g_am.GetXaxis().SetTitle("<PU> [%]")
g_am.GetYaxis().SetTitle("#SigmaE_{T}/#delta#eta#delta#phi [GeV]") 
g_am.GetYaxis().SetTitleOffset(1.6)
g_am.GetYaxis().SetRangeUser(2.7,2.9)
g_am.Draw("AP")

# Add legend
legend = ROOT.TLegend(0.17, 0.8, 0.47, 0.92)
legend.AddEntry(g_am, "Away Region 0mrad + 1.5mrad", "lp")
legend.Draw()

c.Update()
c.Draw()
c.SaveAs("0mrad_plots/h_ue_pileup_away_0mrad_and_1.5mrad_15_20GeV_Topoclusters_no_low_ET_towers.png")
'''


'''

xlist0 = [-125,-25,75,175,275,475]
xlist1 = [-115,-15,85,185,285,485]
xlist2 = [-105,-5,95,195,295,495]
xlist3 = [-95,5,105,205,305,505]
xlist4 = [-85,15,115,215,315,515]
xlist5 = [-75,25,125,225,325,525]
xerr = [0,0,0,0,0,0]
x0 = np.array(xlist0, dtype='float64')
x1 = np.array(xlist1, dtype='float64')
x2 = np.array(xlist2, dtype='float64')
x3 = np.array(xlist3, dtype='float64')
x4 = np.array(xlist4, dtype='float64')
x5 = np.array(xlist5, dtype='float64')
# Create a canvas
canvas = ROOT.TCanvas("canvas2", "Graph", 500, 400)

# Define marker styles and colors
marker_styles = {
    "towards": ROOT.kFullCircle,
    "transverse": ROOT.kFullCircle,
    "away": ROOT.kFullCircle
}

colors = {
    "ntopo": ROOT.kBlack,
    "mc_ntopo": ROOT.kRed,
    "clus_ntopo": ROOT.kBlue,
    "nz_ntopo": ROOT.kMagenta,
    #"dt_ntopo": ROOT.kCyan
}

# Prepare a list of all data groups
data_groups = [
    ("mc_ntopo", x0, mean_mc_ntopo_towards, mean_mc_ntopo_transverse, mean_mc_ntopo_away, std_mc_ntopo_towards, std_mc_ntopo_transverse, std_mc_ntopo_away),
    ("clus_ntopo", x1, mean_clus_ntopo_towards, mean_clus_ntopo_transverse, mean_clus_ntopo_away, std_clus_ntopo_towards, std_clus_ntopo_transverse, std_clus_ntopo_away),
    ("nz_ntopo", x2, mean_nz_ntopo_towards, mean_nz_ntopo_transverse, mean_nz_ntopo_away, std_nz_ntopo_towards, std_nz_ntopo_transverse, std_nz_ntopo_away),
    #("dt_ntopo", x3, mean_dt_ntopo_towards, mean_dt_ntopo_transverse, mean_dt_ntopo_away, std_dt_ntopo_towards, std_dt_ntopo_transverse, std_dt_ntopo_away),
    ("ntopo", x4, mean_ntopo_towards, mean_ntopo_transverse, mean_ntopo_away, std_ntopo_towards, std_ntopo_transverse, std_ntopo_away)
]

graphs = []
    
# Create TGraphs for each data group
for group, x, towards, transverse, away, towards_width, transverse_width, away_width in data_groups:
    n_points = len(towards)

    # Create graphs for towards, transverse, and away
    graph_towards = ROOT.TGraphErrors(n_points, x, np.array(towards, dtype='float64'), np.array(xerr, dtype='float64'), np.array(towards_width, dtype='float64'))
    graph_transverse = ROOT.TGraphErrors(n_points, x, np.array(transverse, dtype='float64'), np.array(xerr, dtype='float64'), np.array(transverse_width, dtype='float64'))
    graph_away = ROOT.TGraphErrors(n_points, x, np.array(away, dtype='float64'), np.array(xerr, dtype='float64'), np.array(away_width, dtype='float64'))

    # Set marker styles and colors
    graph_towards.SetMarkerStyle(marker_styles["towards"])
    graph_transverse.SetMarkerStyle(marker_styles["transverse"])
    graph_away.SetMarkerStyle(marker_styles["away"])

    graph_towards.SetMarkerColor(colors[group])
    graph_transverse.SetMarkerColor(colors[group])
    graph_away.SetMarkerColor(colors[group])

    graph_towards.SetLineColor(colors[group])
    graph_transverse.SetLineColor(colors[group])
    graph_away.SetLineColor(colors[group])

    graphs.extend([graph_towards, graph_transverse, graph_away])

# Draw all graphs on the same canvas
for i, graph in enumerate(graphs):
    if i % 3 == 0:
        if i == 0:
            graph.Draw("AP")
            graph.GetYaxis().SetRangeUser(5,12)
            graph.GetXaxis().SetTitle("E_{topo} thres [MeV]")
            graph.GetYaxis().SetTitle("<N_{topo}>")
        else:
             graph.Draw("P SAME")

# Adding a legend
legend = ROOT.TLegend(0.53, 0.65, 0.92, 0.92)
legend.SetNColumns(1)
legend.SetTextSize(0.035)
legend.AddEntry("","#bf{Towards Region}","")
legend.AddEntry("","15 < p_{T,lead} < 20 GeV","")
legend.AddEntry(graphs[10], "<PU> < 1.5%", "pe")
legend.AddEntry(graphs[1], "1.5 < <PU> < 2.0%", "pe")
legend.AddEntry(graphs[4], "2.0 < <PU> < 2.8%", "pe")
legend.AddEntry(graphs[7], "<PU> > 2.8%", "pe")
#legend.AddEntry(graphs[10], "<PU> > 2.8%", "ep")
legend.Draw()

# Update and display the canvas
canvas.Update()
canvas.Draw()
canvas.SaveAs("/sphenix/u/egm2153/spring_2025/"+direct+"/h_mean_ntopo_15_20GeV_towards.png")
# Create a canvas
canvas = ROOT.TCanvas("canvas6", "Graph", 500, 400)

# Define marker styles and colors
marker_styles = {
    "towards": ROOT.kFullCircle,
    "transverse": ROOT.kFullCircle,
    "away": ROOT.kFullCircle
}

colors = {
    "ntopo": ROOT.kBlack,
    "mc_ntopo": ROOT.kRed,
    "clus_ntopo": ROOT.kBlue,
    "nz_ntopo": ROOT.kMagenta,
    #"dt_ntopo": ROOT.kCyan
}

# Prepare a list of all data groups
data_groups = [
    ("mc_ntopo", x0, mean_mc_ntopo_towards, mean_mc_ntopo_transverse, mean_mc_ntopo_away, std_mc_ntopo_towards, std_mc_ntopo_transverse, std_mc_ntopo_away),
    ("clus_ntopo", x1, mean_clus_ntopo_towards, mean_clus_ntopo_transverse, mean_clus_ntopo_away, std_clus_ntopo_towards, std_clus_ntopo_transverse, std_clus_ntopo_away),
    ("nz_ntopo", x2, mean_nz_ntopo_towards, mean_nz_ntopo_transverse, mean_nz_ntopo_away, std_nz_ntopo_towards, std_nz_ntopo_transverse, std_nz_ntopo_away),
    #("dt_ntopo", x3, mean_dt_ntopo_towards, mean_dt_ntopo_transverse, mean_dt_ntopo_away, std_dt_ntopo_towards, std_dt_ntopo_transverse, std_dt_ntopo_away),
    ("ntopo", x4, mean_ntopo_towards, mean_ntopo_transverse, mean_ntopo_away, std_ntopo_towards, std_ntopo_transverse, std_ntopo_away)
]

graphs = []
    
# Create TGraphs for each data group
for group, x, towards, transverse, away, towards_width, transverse_width, away_width in data_groups:
    n_points = len(towards)

    # Create graphs for towards, transverse, and away
    graph_towards = ROOT.TGraphErrors(n_points, x, np.array(towards, dtype='float64'), np.array(xerr, dtype='float64'), np.array(towards_width, dtype='float64'))
    graph_transverse = ROOT.TGraphErrors(n_points, x, np.array(transverse, dtype='float64'), np.array(xerr, dtype='float64'), np.array(transverse_width, dtype='float64'))
    graph_away = ROOT.TGraphErrors(n_points, x, np.array(away, dtype='float64'), np.array(xerr, dtype='float64'), np.array(away_width, dtype='float64'))

    # Set marker styles and colors
    graph_towards.SetMarkerStyle(marker_styles["towards"])
    graph_transverse.SetMarkerStyle(marker_styles["transverse"])
    graph_away.SetMarkerStyle(marker_styles["away"])

    graph_towards.SetMarkerColor(colors[group])
    graph_transverse.SetMarkerColor(colors[group])
    graph_away.SetMarkerColor(colors[group])

    graph_towards.SetLineColor(colors[group])
    graph_transverse.SetLineColor(colors[group])
    graph_away.SetLineColor(colors[group])

    graphs.extend([graph_towards, graph_transverse, graph_away])

# Draw all graphs on the same canvas
for i, graph in enumerate(graphs):
    if (i-1) % 3 == 0:
        if i == 1:
            graph.Draw("AP")
            graph.GetYaxis().SetRangeUser(0.5,5)
            graph.GetXaxis().SetTitle("E_{topo} thres [MeV]")
            graph.GetYaxis().SetTitle("<N_{topo}>")
        else:
             graph.Draw("P SAME")

# Adding a legend
legend = ROOT.TLegend(0.53, 0.65, 0.92, 0.92)
legend.SetNColumns(1)
legend.SetTextSize(0.035)
legend.AddEntry("","#bf{Transverse Region}","")
legend.AddEntry("","15 < p_{T,lead} < 20 GeV","")
legend.AddEntry(graphs[10], "<PU> < 1.5%", "pe")
legend.AddEntry(graphs[1], "1.5 < <PU> < 2.0%", "pe")
legend.AddEntry(graphs[4], "2.0 < <PU> < 2.8%", "pe")
legend.AddEntry(graphs[7], "<PU> > 2.8%", "pe")
#legend.AddEntry(graphs[10], "Run 150 (Detroit)", "ep")
legend.Draw()

# Update and display the canvas
canvas.Update()
canvas.Draw()
canvas.SaveAs("/sphenix/u/egm2153/spring_2025/"+direct+"/h_mean_ntopo_15_20GeV_transverse.png")
# Create a canvas
canvas = ROOT.TCanvas("canvas7", "Graph", 500, 400)

# Define marker styles and colors
marker_styles = {
    "towards": ROOT.kFullCircle,
    "transverse": ROOT.kFullCircle,
    "away": ROOT.kFullCircle
}

colors = {
    "ntopo": ROOT.kBlack,
    "mc_ntopo": ROOT.kRed,
    "clus_ntopo": ROOT.kBlue,
    "nz_ntopo": ROOT.kMagenta,
    #"dt_ntopo": ROOT.kCyan
}

# Prepare a list of all data groups
data_groups = [
    ("mc_ntopo", x0, mean_mc_ntopo_towards, mean_mc_ntopo_transverse, mean_mc_ntopo_away, std_mc_ntopo_towards, std_mc_ntopo_transverse, std_mc_ntopo_away),
    ("clus_ntopo", x1, mean_clus_ntopo_towards, mean_clus_ntopo_transverse, mean_clus_ntopo_away, std_clus_ntopo_towards, std_clus_ntopo_transverse, std_clus_ntopo_away),
    ("nz_ntopo", x2, mean_nz_ntopo_towards, mean_nz_ntopo_transverse, mean_nz_ntopo_away, std_nz_ntopo_towards, std_nz_ntopo_transverse, std_nz_ntopo_away),
    #("dt_ntopo", x3, mean_dt_ntopo_towards, mean_dt_ntopo_transverse, mean_dt_ntopo_away, std_dt_ntopo_towards, std_dt_ntopo_transverse, std_dt_ntopo_away),
    ("ntopo", x4, mean_ntopo_towards, mean_ntopo_transverse, mean_ntopo_away, std_ntopo_towards, std_ntopo_transverse, std_ntopo_away)
]

graphs = []
    
# Create TGraphs for each data group
for group, x, towards, transverse, away, towards_width, transverse_width, away_width in data_groups:
    n_points = len(towards)

    # Create graphs for towards, transverse, and away
    graph_towards = ROOT.TGraphErrors(n_points, x, np.array(towards, dtype='float64'), np.array(xerr, dtype='float64'), np.array(towards_width, dtype='float64'))
    graph_transverse = ROOT.TGraphErrors(n_points, x, np.array(transverse, dtype='float64'), np.array(xerr, dtype='float64'), np.array(transverse_width, dtype='float64'))
    graph_away = ROOT.TGraphErrors(n_points, x, np.array(away, dtype='float64'), np.array(xerr, dtype='float64'), np.array(away_width, dtype='float64'))

    # Set marker styles and colors
    graph_towards.SetMarkerStyle(marker_styles["towards"])
    graph_transverse.SetMarkerStyle(marker_styles["transverse"])
    graph_away.SetMarkerStyle(marker_styles["away"])

    graph_towards.SetMarkerColor(colors[group])
    graph_transverse.SetMarkerColor(colors[group])
    graph_away.SetMarkerColor(colors[group])

    graph_towards.SetLineColor(colors[group])
    graph_transverse.SetLineColor(colors[group])
    graph_away.SetLineColor(colors[group])

    graphs.extend([graph_towards, graph_transverse, graph_away])

# Draw all graphs on the same canvas
for i, graph in enumerate(graphs):
    if (i-2) % 3 == 0:
        if i == 2:
            graph.Draw("AP")
            graph.GetYaxis().SetRangeUser(4,14)
            graph.GetXaxis().SetTitle("E_{topo} thres [MeV]")
            graph.GetYaxis().SetTitle("<N_{topo}>")
        else:
             graph.Draw("P SAME")

# Adding a legend
legend = ROOT.TLegend(0.53, 0.65, 0.92, 0.92)
legend.SetNColumns(1)
legend.SetTextSize(0.035)
legend.AddEntry("","#bf{Away Region}","")
legend.AddEntry("","15 < p_{T,lead} < 20 GeV","")
legend.AddEntry(graphs[1], "1.5 < <PU> < 2.0%", "pe")
legend.AddEntry(graphs[4], "2.0 < <PU> < 2.8%", "pe")
legend.AddEntry(graphs[7], "<PU> > 2.8%", "pe")
#legend.AddEntry(graphs[10], "Run 150 (Detroit)", "ep")
legend.AddEntry(graphs[10], "<PU> < 1.5%", "pe")
legend.Draw()

# Update and display the canvas
canvas.Update()
canvas.Draw()
canvas.SaveAs("/sphenix/u/egm2153/spring_2025/"+direct+"/h_mean_ntopo_15_20GeV_away.png")


# Create a canvas
canvas = ROOT.TCanvas("canvas", "Graph", 500, 400)

# Define marker styles and colors
marker_styles = {
    "towards": ROOT.kFullCircle,
    "transverse": ROOT.kFullCircle,
    "away": ROOT.kFullCircle
}

colors = {
    "ntopo": ROOT.kBlack,
    "mc_ntopo": ROOT.kRed,
    "clus_ntopo": ROOT.kBlue,
    "nz_ntopo": ROOT.kMagenta,
    "dt_ntopo": ROOT.kCyan
}

# Prepare a list of all data groups
data_groups = [
    ("mc_ntopo", x0, mean_mc_topo_towards, mean_mc_topo_transverse, mean_mc_topo_away, std_mc_topo_towards, std_mc_topo_transverse, std_mc_topo_away),
    ("clus_ntopo", x1, mean_clus_topo_towards, mean_clus_topo_transverse, mean_clus_topo_away, std_clus_topo_towards, std_clus_topo_transverse, std_clus_topo_away),
    ("nz_ntopo", x2, mean_nz_topo_towards, mean_nz_topo_transverse, mean_nz_topo_away, std_nz_topo_towards, std_nz_topo_transverse, std_nz_topo_away),
    #("dt_ntopo", x3, mean_dt_topo_towards, mean_dt_topo_transverse, mean_dt_topo_away, std_dt_topo_towards, std_dt_topo_transverse, std_dt_topo_away),
    ("ntopo", x4, mean_topo_towards, mean_topo_transverse, mean_topo_away, std_topo_towards, std_topo_transverse, std_topo_away)
]

graphs = []
    
# Create TGraphs for each data group
for group, x, towards, transverse, away, towards_width, transverse_width, away_width in data_groups:
    n_points = len(towards)

    # Create graphs for towards, transverse, and away
    graph_towards = ROOT.TGraphErrors(n_points, x, np.array(towards, dtype='float64'), np.array(xerr, dtype='float64'), np.array(towards_width, dtype='float64'))
    graph_transverse = ROOT.TGraphErrors(n_points, x, np.array(transverse, dtype='float64'), np.array(xerr, dtype='float64'), np.array(transverse_width, dtype='float64'))
    graph_away = ROOT.TGraphErrors(n_points, x, np.array(away, dtype='float64'), np.array(xerr, dtype='float64'), np.array(away_width, dtype='float64'))

    # Set marker styles and colors
    graph_towards.SetMarkerStyle(marker_styles["towards"])
    graph_transverse.SetMarkerStyle(marker_styles["transverse"])
    graph_away.SetMarkerStyle(marker_styles["away"])

    graph_towards.SetMarkerColor(colors[group])
    graph_transverse.SetMarkerColor(colors[group])
    graph_away.SetMarkerColor(colors[group])

    graph_towards.SetMarkerSize(1)
    graph_transverse.SetMarkerSize(1)
    graph_away.SetMarkerSize(1)

    graph_towards.SetLineColor(colors[group])
    graph_transverse.SetLineColor(colors[group])
    graph_away.SetLineColor(colors[group])

    graphs.extend([graph_towards, graph_transverse, graph_away])

# Draw all graphs on the same canvas
for i, graph in enumerate(graphs):
    if i % 3 == 0:
        if i == 0:
            graph.Draw("AP")
            graph.GetYaxis().SetRangeUser(1,4)
            graph.GetXaxis().SetTitle("E_{topo} thres [MeV]")
            graph.GetYaxis().SetTitle("<E_{T,topo}> [GeV]")
        else:
             graph.Draw("P SAME")

# Adding a legend
legend = ROOT.TLegend(0.17, 0.65, 0.5, 0.92)
legend.SetNColumns(1)
legend.SetTextSize(0.035)
legend.AddEntry("","#bf{Towards Region}","")
legend.AddEntry("","15 < p_{T,lead} < 20 GeV","")
legend.AddEntry(graphs[10], "<PU> < 1.5%", "pe")
legend.AddEntry(graphs[1], "1.5 < <PU> < 2.0%", "pe")
legend.AddEntry(graphs[4], "2.0 < <PU> < 2.8%", "pe")
legend.AddEntry(graphs[7], "<PU> > 2.8%", "pe")
#legend.AddEntry(graphs[10], "Run 150 (Detroit)", "ep")
legend.Draw()

# Update and display the canvas
canvas.Update()
canvas.Draw()
canvas.SaveAs("/sphenix/u/egm2153/spring_2025/"+direct+"/h_mean_etopo_15_20GeV_towards.png")
# Create a canvas
canvas = ROOT.TCanvas("canvas", "Graph", 500, 400)

# Define marker styles and colors
marker_styles = {
    "towards": ROOT.kFullCircle,
    "transverse": ROOT.kFullCircle,
    "away": ROOT.kFullCircle
}

colors = {
    "ntopo": ROOT.kBlack,
    "mc_ntopo": ROOT.kRed,
    "clus_ntopo": ROOT.kBlue,
    "nz_ntopo": ROOT.kMagenta,
    "dt_ntopo": ROOT.kCyan
}

# Prepare a list of all data groups
data_groups = [
    ("mc_ntopo", x0, mean_mc_topo_towards, mean_mc_topo_transverse, mean_mc_topo_away, std_mc_topo_towards, std_mc_topo_transverse, std_mc_topo_away),
    ("clus_ntopo", x1, mean_clus_topo_towards, mean_clus_topo_transverse, mean_clus_topo_away, std_clus_topo_towards, std_clus_topo_transverse, std_clus_topo_away),
    ("nz_ntopo", x2, mean_nz_topo_towards, mean_nz_topo_transverse, mean_nz_topo_away, std_nz_topo_towards, std_nz_topo_transverse, std_nz_topo_away),
    #("dt_ntopo", x3, mean_dt_topo_towards, mean_dt_topo_transverse, mean_dt_topo_away, std_dt_topo_towards, std_dt_topo_transverse, std_dt_topo_away),
    ("ntopo", x4, mean_topo_towards, mean_topo_transverse, mean_topo_away, std_topo_towards, std_topo_transverse, std_topo_away)
]

graphs = []
    
# Create TGraphs for each data group
for group, x, towards, transverse, away, towards_width, transverse_width, away_width in data_groups:
    n_points = len(towards)

    # Create graphs for towards, transverse, and away
    graph_towards = ROOT.TGraphErrors(n_points, x, np.array(towards, dtype='float64'), np.array(xerr, dtype='float64'), np.array(towards_width, dtype='float64'))
    graph_transverse = ROOT.TGraphErrors(n_points, x, np.array(transverse, dtype='float64'), np.array(xerr, dtype='float64'), np.array(transverse_width, dtype='float64'))
    graph_away = ROOT.TGraphErrors(n_points, x, np.array(away, dtype='float64'), np.array(xerr, dtype='float64'), np.array(away_width, dtype='float64'))

    # Set marker styles and colors
    graph_towards.SetMarkerStyle(marker_styles["towards"])
    graph_transverse.SetMarkerStyle(marker_styles["transverse"])
    graph_away.SetMarkerStyle(marker_styles["away"])

    graph_towards.SetMarkerColor(colors[group])
    graph_transverse.SetMarkerColor(colors[group])
    graph_away.SetMarkerColor(colors[group])
    
    graph_towards.SetMarkerSize(1)
    graph_transverse.SetMarkerSize(1)
    graph_away.SetMarkerSize(1)

    graph_towards.SetLineColor(colors[group])
    graph_transverse.SetLineColor(colors[group])
    graph_away.SetLineColor(colors[group])

    graphs.extend([graph_towards, graph_transverse, graph_away])

# Draw all graphs on the same canvas
for i, graph in enumerate(graphs):
    if (i-1) % 3 == 0:
        if i == 1:
            graph.Draw("AP")
            graph.GetYaxis().SetRangeUser(0,1)
            graph.GetXaxis().SetTitle("E_{topo} thres [MeV]")
            graph.GetYaxis().SetTitle("<E_{T,topo}> [GeV]")
        else:
             graph.Draw("P SAME")

# Adding a legend
legend = ROOT.TLegend(0.17, 0.65, 0.5, 0.92)
legend.SetNColumns(1)
legend.SetTextSize(0.035)
legend.AddEntry("","#bf{Transverse Region}","")
legend.AddEntry("","15 < p_{T,lead} < 20 GeV","")
legend.AddEntry(graphs[10], "<PU> < 1.5%", "pe")
legend.AddEntry(graphs[1], "1.5 < <PU> < 2.0%", "pe")
legend.AddEntry(graphs[4], "2.0 < <PU> < 2.8%", "pe")
legend.AddEntry(graphs[7], "<PU> > 2.8%", "pe")
#legend.AddEntry(graphs[10], "Run 150 (Detroit)", "ep")
legend.Draw()

# Update and display the canvas
canvas.Update()
canvas.Draw()
canvas.SaveAs("/sphenix/u/egm2153/spring_2025/"+direct+"/h_mean_etopo_15_20GeV_transverse.png")
# Create a canvas
canvas = ROOT.TCanvas("canvas", "Graph", 500, 400)

# Define marker styles and colors
marker_styles = {
    "towards": ROOT.kFullCircle,
    "transverse": ROOT.kFullCircle,
    "away": ROOT.kFullCircle
}

colors = {
    "ntopo": ROOT.kBlack,
    "mc_ntopo": ROOT.kRed,
    "clus_ntopo": ROOT.kBlue,
    "nz_ntopo": ROOT.kMagenta,
    "dt_ntopo": ROOT.kCyan
}

# Prepare a list of all data groups
data_groups = [
    ("mc_ntopo", x0, mean_mc_topo_towards, mean_mc_topo_transverse, mean_mc_topo_away, std_mc_topo_towards, std_mc_topo_transverse, std_mc_topo_away),
    ("clus_ntopo", x1, mean_clus_topo_towards, mean_clus_topo_transverse, mean_clus_topo_away, std_clus_topo_towards, std_clus_topo_transverse, std_clus_topo_away),
    ("nz_ntopo", x2, mean_nz_topo_towards, mean_nz_topo_transverse, mean_nz_topo_away, std_nz_topo_towards, std_nz_topo_transverse, std_nz_topo_away),
    #("dt_ntopo", x3, mean_dt_topo_towards, mean_dt_topo_transverse, mean_dt_topo_away, std_dt_topo_towards, std_dt_topo_transverse, std_dt_topo_away),
    ("ntopo", x4, mean_topo_towards, mean_topo_transverse, mean_topo_away, std_topo_towards, std_topo_transverse, std_topo_away)
]

graphs = []
    
# Create TGraphs for each data group
for group, x, towards, transverse, away, towards_width, transverse_width, away_width in data_groups:
    n_points = len(towards)

    # Create graphs for towards, transverse, and away
    graph_towards = ROOT.TGraphErrors(n_points, x, np.array(towards, dtype='float64'), np.array(xerr, dtype='float64'), np.array(towards_width, dtype='float64'))
    graph_transverse = ROOT.TGraphErrors(n_points, x, np.array(transverse, dtype='float64'), np.array(xerr, dtype='float64'), np.array(transverse_width, dtype='float64'))
    graph_away = ROOT.TGraphErrors(n_points, x, np.array(away, dtype='float64'), np.array(xerr, dtype='float64'), np.array(away_width, dtype='float64'))

    # Set marker styles and colors
    graph_towards.SetMarkerStyle(marker_styles["towards"])
    graph_transverse.SetMarkerStyle(marker_styles["transverse"])
    graph_away.SetMarkerStyle(marker_styles["away"])

    graph_towards.SetMarkerColor(colors[group])
    graph_transverse.SetMarkerColor(colors[group])
    graph_away.SetMarkerColor(colors[group])

    graph_towards.SetLineColor(colors[group])
    graph_transverse.SetLineColor(colors[group])
    graph_away.SetLineColor(colors[group])

    graphs.extend([graph_towards, graph_transverse, graph_away])

# Draw all graphs on the same canvas
for i, graph in enumerate(graphs):
    if (i-2) % 3 == 0:
        if i == 2:
            graph.Draw("AP")
            graph.GetYaxis().SetRangeUser(0.5,2.5)
            graph.GetXaxis().SetTitle("E_{topo} thres [MeV]")
            graph.GetYaxis().SetTitle("<E_{T,topo}> [GeV]")
        else:
             graph.Draw("P SAME")

# Adding a legend
legend = ROOT.TLegend(0.17, 0.65, 0.5, 0.92)
legend.SetNColumns(1)d
legend.SetTextSize(0.035)
legend.AddEntry("","#bf{Away Region}","")
legend.AddEntry("","15 < p_{T,lead} < 20 GeV","")
legend.AddEntry(graphs[10], "<PU> < 1.5%", "pe")
legend.AddEntry(graphs[1], "1.5 < <PU> < 2.0%", "pe")
legend.AddEntry(graphs[4], "2.0 < <PU> < 2.8%", "pe")
legend.AddEntry(graphs[7], "<PU> > 2.8%", "pe")
#legend.AddEntry(graphs[10], "Run 150 (Detroit)", "ep")
legend.Draw()

# Update and display the canvas
canvas.Update()
canvas.Draw()
canvas.SaveAs("/sphenix/u/egm2153/spring_2025/"+direct+"/h_mean_etopo_15_20GeV_away.png")
'''