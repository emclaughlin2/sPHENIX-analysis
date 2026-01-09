import ROOT
import os

outdir = "jet_time_plots"
os.makedirs(outdir, exist_ok=True)

# ----------------------------------------------------------------------
# Input files and labels
# ----------------------------------------------------------------------
files = [
    #"analysis_data_ana509_output/output_low_pu_runs_efrac_and_time_bkg_cut.root",
    "analysis_data_ana509_output/output_none_and_time_bkg_cut.root",
    "analysis_data_ana509_output/output_efrac_and_time_bkg_cut.root",
    "analysis_data_ana509_output/output_efrac_bkg_only_and_time_bkg_cut.root"
]

#labels = ["Efrac cut 1.5mrad", 
labels = ["No bkg cut all", "Efrac cut all", "Only efrac cut bkg all"]
labels_string = ["efrac_cut_1.5mrad", "efrac_cut_all", "no_bkg_cut_all","only_bkg_all"]
em_frac_labels = ["0-20%","20-40%","40-60%","60-80%","80-100%"]
jet_pt_labels = ['20-30 GeV','30-40 GeV','40-50 GeV','50-60 GeV']

# Histograms you want to extract:
hist_names = [
    "h_zvertex",
    "h_lead_spectra_record",
    "h_et_transverse_record",
    "h_total_measure_tight",
    "h_pass_time_cut_measure_tight",
    "h_mbd_time",
    "h_jet_mbd_time",
    "h_jet_mbd_timediff",
    "h_jettime",
    "h_em_jet_frac_jettime",
    "h_em_jet_frac_jet_mbd_timediff",
    "h_em_jet_frac_mbd_time",
    "h_jet_pt_jettime",
    "h_jet_pt_jet_mbd_timediff",
    "h_jet_pt_mbd_time"
]

full_hist_names = [
    "h_zvertex",
    "h_lead_spectra_record",
    "h_et_transverse_record",
    "h_total_measure_tight",
    "h_pass_time_cut_measure_tight",
    "h_mbd_time",
    "h_jet_mbd_time",
    "h_jet_mbd_timediff",
    "h_jettime",
    "h_em_jet_frac_jettime",
    "h_em_jet_frac_jet_mbd_timediff",
    "h_em_jet_frac_mbd_time",
    "h_em_jet_frac_jettime_1",
    "h_em_jet_frac_jet_mbd_timediff_1",
    "h_em_jet_frac_mbd_time_1",
    "h_em_jet_frac_jettime_2",
    "h_em_jet_frac_jet_mbd_timediff_2",
    "h_em_jet_frac_mbd_time_2",
    "h_em_jet_frac_jettime_3",
    "h_em_jet_frac_jet_mbd_timediff_3",
    "h_em_jet_frac_mbd_time_3",
    "h_em_jet_frac_jettime_4",
    "h_em_jet_frac_jet_mbd_timediff_4",
    "h_em_jet_frac_mbd_time_4",
    "h_em_jet_frac_jettime_5",
    "h_em_jet_frac_jet_mbd_timediff_5",
    "h_em_jet_frac_mbd_time_5",
    "h_em_jet_frac_jettime_efrac_cut_1.5mrad",
    "h_em_jet_frac_jettime_efrac_cut_all",
    "h_em_jet_frac_jettime_no_bkg_cut_all",
    "h_em_jet_frac_jettime_only_bkg_all",
    "h_em_jet_frac_jet_mbd_timediff_efrac_cut_1.5mrad",
    "h_em_jet_frac_jet_mbd_timediff_efrac_cut_all",
    "h_em_jet_frac_jet_mbd_timediff_no_bkg_cut_all",
    "h_em_jet_frac_jet_mbd_timediff_only_bkg_all",
    "h_em_jet_frac_mbd_time_efrac_cut_1.5mrad",
    "h_em_jet_frac_mbd_time_efrac_cut_all",
    "h_em_jet_frac_mbd_time_no_bkg_all",
    "h_em_jet_frac_mbd_time_only_bkg_all",
    "h_jet_pt_jettime",
    "h_jet_pt_jet_mbd_timediff",
    "h_jet_pt_mbd_time",
    "h_jet_pt_jettime_1",
    "h_jet_pt_jet_mbd_timediff_1",
    "h_jet_pt_mbd_time_1",
    "h_jet_pt_jettime_2",
    "h_jet_pt_jet_mbd_timediff_2",
    "h_jet_pt_mbd_time_2",
    "h_jet_pt_jettime_3",
    "h_jet_pt_jet_mbd_timediff_3",
    "h_jet_pt_mbd_time_3",
    "h_jet_pt_jettime_4",
    "h_jet_pt_jet_mbd_timediff_4",
    "h_jet_pt_mbd_time_4",
    "h_jet_pt_jettime_efrac_cut_1.5mrad",
    "h_jet_pt_jettime_efrac_cut_all",
    "h_jet_pt_jettime_no_bkg_cut_all",
    "h_jet_pt_jettime_only_bkg_all",
    "h_jet_pt_jet_mbd_timediff_efrac_cut_1.5mrad",
    "h_jet_pt_jet_mbd_timediff_efrac_cut_all",
    "h_jet_pt_jet_mbd_timediff_no_bkg_cut_all",
    "h_jet_pt_jet_mbd_timediff_only_bkg_all",
    "h_jet_pt_mbd_time_efrac_cut_1.5mrad",
    "h_jet_pt_mbd_time_efrac_cut_all",
    "h_jet_pt_mbd_time_no_bkg_cut_all",
    "h_jet_pt_mbd_time_only_bkg_all"
]

# ----------------------------------------------------------------------
# Load ROOT files and histograms
# ----------------------------------------------------------------------
root_files = []
all_histos = {name: [] for name in full_hist_names}

for f in files:
    rf = ROOT.TFile.Open(f)
    if not rf or rf.IsZombie():
        print(f"Error opening {f}")
        continue
    root_files.append(rf)

    for name in hist_names:
        h = rf.Get(name)
        if not h:
            print(f"Histogram {name} not found in file {f}")
            all_histos[name].append(None)
        else:
            h.SetDirectory(0)  # detach from file so file can close
            all_histos[name].append(h.Clone())

def save_canvas(canvas, name):
    """Save a canvas to PNG inside jet_time_plots/."""
    filename = os.path.join(outdir, f"{name}.png")
    canvas.SaveAs(filename)
    print("Saved:", filename)

# ----------------------------------------------------------------------
# Helper: overlay 1D histograms
# ----------------------------------------------------------------------
def plot_overlay(hists, labels, title, ymax = 0):
    c = ROOT.TCanvas(title, title, 900, 700)
    c.SetGrid()

    first = True
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2, ROOT.kMagenta, ROOT.kBlack]

    for i, h in enumerate(hists):
        if not h: 
            continue
        h.SetLineColor(colors[i])
        h.SetLineWidth(2)
        if h.GetEntries() != 0:
            h.Scale(1.0/h.Integral())
        if ymax > 0:
            h.GetYaxis().SetRangeUser(0,ymax)
        if first:
            h.Draw("hist")
            h.SetTitle(title)
            h.SetXTitle(title)
            h.GetYaxis().SetTitle("Counts")
            first = False
        else:
            h.Draw("hist same")

    # Legend
    leg = ROOT.TLegend(0.65, 0.7, 0.88, 0.88)
    for i, h in enumerate(hists):
        if h:
            leg.AddEntry(h, labels[i], "l")
    leg.Draw()

    c.Update()
    save_canvas(c, title.replace(" ", "_"))

# ----------------------------------------------------------------------
# Helper: show 2D histograms side-by-side
# ----------------------------------------------------------------------
def plot_side_by_side(hists, labels, title):
    c = ROOT.TCanvas(title, title, 1800, 600)
    c.Divide(3, 1)

    for i, h in enumerate(hists):
        c.cd(i+1)
        if h:
            h.Draw("COLZ")
            h.SetXTitle("Leading Jet Time [ns]")
            h.SetYTitle("MBD Time [ns]")
            h.SetTitle(f"{title} {labels[i]}")
    c.Update()
    save_canvas(c, title.replace(" ", "_"))

# ----------------------------------------------------------------------
# Example plots:
# ----------------------------------------------------------------------
'''
# === Overlaid 1D histograms ===
plot_overlay(all_histos["h_zvertex"], labels, "Z Vertex Comparison")
plot_overlay(all_histos["h_lead_spectra_record"], labels, "Leading Jet Spectra")
plot_overlay(all_histos["h_mbd_time"], labels, "MBD Time [ns]")
plot_overlay(all_histos["h_jettime"], labels, "Leading Jet Time [ns]")
plot_overlay(all_histos["h_jet_mbd_timediff"], labels, "Leading Jet - MBD Time [ns]")

# === 2D histograms (side-by-side) ===
plot_side_by_side(all_histos["h_total_measure_tight"], labels, "Total Measure Tight")
plot_side_by_side(all_histos["h_jet_mbd_time"], labels, "Jet MBD Time")

for h, label in zip(all_histos["h_em_jet_frac_jettime"], labels):
    for i in range(1, h.GetNbinsX() + 1):
        hist = h.ProjectionY(f"h_em_jet_frac_jettime_{label}_{i}",i,i)
        all_histos[f"h_em_jet_frac_jettime_{i}"].append(hist.Clone())

for h, label in zip(all_histos["h_em_jet_frac_jet_mbd_timediff"], labels):
    for i in range(1, h.GetNbinsX() + 1):
        hist = h.ProjectionY(f"h_em_jet_frac_jet_mbd_timediff_{label}_{i}",i,i)
        all_histos[f"h_em_jet_frac_jet_mbd_timediff_{i}"].append(hist.Clone())

for h, label in zip(all_histos["h_em_jet_frac_mbd_time"], labels):
    for i in range(1, h.GetNbinsX() + 1):
        hist = h.ProjectionY(f"h_em_jet_frac_mbd_time_{label}_{i}",i,i)
        all_histos[f"h_em_jet_frac_mbd_time_{i}"].append(hist.Clone())

for i in range(5,0,-1):
    labels_temp = [label+" "+em_frac_labels[i-1] for label in labels]
    plot_overlay(all_histos[f"h_em_jet_frac_jettime_{i}"], labels_temp, f"Leading Jet Time [ns] - {em_frac_labels[i-1]} EM frac")

for i in range(5,0,-1):
    labels_temp = [label+" "+em_frac_labels[i-1] for label in labels]
    plot_overlay(all_histos[f"h_em_jet_frac_jet_mbd_timediff_{i}"], labels_temp, f"Leading Jet - MBD Time [ns] - {em_frac_labels[i-1]} EM frac")

for i in range(5,0,-1):
    labels_temp = [label+" "+em_frac_labels[i-1] for label in labels]
    plot_overlay(all_histos[f"h_em_jet_frac_mbd_time_{i}"], labels_temp, f"MBD Time [ns] - {em_frac_labels[i-1]} EM frac")

for l, label in enumerate(labels):
    for i in range(5,0,-1):
        all_histos[f'h_em_jet_frac_jettime_{labels_string[l]}'].append(all_histos[f"h_em_jet_frac_jettime_{i}"][l])
        all_histos[f'h_em_jet_frac_jet_mbd_timediff_{labels_string[l]}'].append(all_histos[f"h_em_jet_frac_jet_mbd_timediff_{i}"][l])
        all_histos[f'h_em_jet_frac_mbd_time_{labels_string[l]}'].append(all_histos[f"h_em_jet_frac_mbd_time_{i}"][l])
    
    plot_overlay(all_histos[f'h_em_jet_frac_jettime_{labels_string[l]}'], em_frac_labels, f"Leading Jet Time [ns] - {label}",0.05)
    plot_overlay(all_histos[f'h_em_jet_frac_jet_mbd_timediff_{labels_string[l]}'], em_frac_labels, f"Leading Jet - MBD Time [ns] - {label}",0.15)
    plot_overlay(all_histos[f'h_em_jet_frac_mbd_time_{labels_string[l]}'], em_frac_labels, f"MBD Time [ns] - {label}",0.05)
'''
for h, label in zip(all_histos["h_jet_pt_jettime"], labels):
    for i in range(1, h.GetNbinsX() + 1):
        hist = h.ProjectionY(f"h_jet_pt_jettime_{label}_{i}",i,i)
        all_histos[f"h_jet_pt_jettime_{i}"].append(hist.Clone())

for h, label in zip(all_histos["h_jet_pt_jet_mbd_timediff"], labels):
    for i in range(1, h.GetNbinsX() + 1):
        hist = h.ProjectionY(f"h_jet_pt_jet_mbd_timediff_{label}_{i}",i,i)
        all_histos[f"h_jet_pt_jet_mbd_timediff_{i}"].append(hist.Clone())

for h, label in zip(all_histos["h_jet_pt_mbd_time"], labels):
    for i in range(1, h.GetNbinsX() + 1):
        hist = h.ProjectionY(f"h_jet_pt_mbd_time_{label}_{i}",i,i)
        all_histos[f"h_jet_pt_mbd_time_{i}"].append(hist.Clone())

for i in range(1,5):
    labels_temp = [label+" "+jet_pt_labels[i-1] for label in labels]
    plot_overlay(all_histos[f"h_jet_pt_jettime_{i}"], labels_temp, f"Leading Jet Time [ns] - Uncalib Jet pT {jet_pt_labels[i-1]}")

for i in range(1,5):
    labels_temp = [label+" "+jet_pt_labels[i-1] for label in labels]
    plot_overlay(all_histos[f"h_jet_pt_jet_mbd_timediff_{i}"], labels_temp, f"Leading Jet - MBD Time [ns] - Uncalib Jet pT {jet_pt_labels[i-1]}")

for i in range(1,5):
    labels_temp = [label+" "+jet_pt_labels[i-1] for label in labels]
    plot_overlay(all_histos[f"h_jet_pt_mbd_time_{i}"], labels_temp, f"MBD Time [ns] - Uncalib Jet pT {jet_pt_labels[i-1]}")

for l, label in enumerate(labels):
    for i in range(1,5):
        all_histos[f'h_jet_pt_jettime_{labels_string[l]}'].append(all_histos[f"h_jet_pt_jettime_{i}"][l])
        all_histos[f'h_jet_pt_jet_mbd_timediff_{labels_string[l]}'].append(all_histos[f"h_jet_pt_jet_mbd_timediff_{i}"][l])
        all_histos[f'h_jet_pt_mbd_time_{labels_string[l]}'].append(all_histos[f"h_jet_pt_mbd_time_{i}"][l])

    plot_overlay(all_histos[f'h_jet_pt_jettime_{labels_string[l]}'], jet_pt_labels, f"Leading Jet Time [ns] - {label}")
    plot_overlay(all_histos[f'h_jet_pt_jet_mbd_timediff_{labels_string[l]}'], jet_pt_labels, f"Leading Jet - MBD Time [ns] - {label}")
    plot_overlay(all_histos[f'h_jet_pt_mbd_time_{labels_string[l]}'], jet_pt_labels, f"MBD Time [ns] - {label}")
    
