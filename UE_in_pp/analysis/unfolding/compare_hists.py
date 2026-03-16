import ROOT
ROOT.gROOT.SetBatch(True)

types = ["dijet", "efrac"]
labels = ["lowpuruns", "midpuruns", "lowpuruns_notiming_njet_cut"]

hist_names = [
    "h_zvertex",
    "h_lead_spectra_record",
    "h_et_transverse_record",
    "h_ue_pt_transverse_record",
    "h_calibjet_pt_dijet_eff",
    "h_calibjet_pt_dijet_pu_correct_et"
]

# ------------------------------------------------------
def open_file(jet_type, label):
    fname = f"analysis_data_run28_output/output_pu_correct_data_{jet_type}_bkg_cut_{label}.root"
    f = ROOT.TFile.Open(fname)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open {fname}")
    return f

# ------------------------------------------------------
def get_hist(file, name, jet_type, label):
    h = file.Get(name)
    if not h:
        raise RuntimeError(
            f"Histogram '{name}' not found in {file.GetName()}"
        )

    if name == "h_ue_pt_transverse_record" or name == "h_calibjet_pt_dijet_eff" or name == "h_calibjet_pt_dijet_pu_correct_et":
        unique_name = f"{name}_jet{jet_type}_{label}"
        h2_clone = h.Clone(f"{unique_name}_2")
        h2_clone.SetDirectory(0)
        h_clone = h2_clone.ProfileX(unique_name)
        h_clone.SetDirectory(0)
    else: 
        unique_name = f"{name}_jet{jet_type}_{label}"
        h_clone = h.Clone(unique_name)
        h_clone.SetDirectory(0)

    return h_clone

# ------------------------------------------------------
def plot_overlay(h1, h2, label1, label2, title, outname):
    c = ROOT.TCanvas(outname, "", 800, 600)

    h1.SetLineColor(ROOT.kRed)
    h2.SetLineColor(ROOT.kBlue)

    h1.SetLineWidth(2)
    h2.SetLineWidth(2)

    h1.SetTitle(title)
    h1.Draw("hist")
    h2.Draw("hist same")

    leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    leg.AddEntry(h1, label1, "l")
    leg.AddEntry(h2, label2, "l")
    leg.Draw()

    c.SaveAs("compare_reco_plots/" + outname + ".pdf")

def plot_ratio(h1, h2, label1, label2, title, outname):
    # Clone to avoid modifying originals
    ratio = h2.Clone(outname + "_ratio")
    ratio.SetDirectory(0)

    # Safe division
    ratio.Divide(h1)

    c = ROOT.TCanvas(outname + "_ratio", "", 800, 600)

    ratio.SetTitle(title)
    ratio.SetLineWidth(2)
    ratio.SetLineColor(ROOT.kBlack)
    ratio_label = label2 + " / " + label1
    ratio.GetYaxis().SetTitle(ratio_label)
    ratio.GetYaxis().SetRangeUser(0, 2)  # adjust if needed

    ratio.Draw("hist")

    # Reference line at 1
    line = ROOT.TLine(
        ratio.GetXaxis().GetXmin(),
        1.0,
        ratio.GetXaxis().GetXmax(),
        1.0,
    )
    line.SetLineStyle(2)
    line.Draw()

    c.SaveAs("compare_reco_plots/" + outname + "_ratio.pdf")

# ======================================================
# 1) Plot 21 vs 28 for each type
# ======================================================

for jet_type in types:
    f1 = open_file(jet_type, labels[0])
    f2 = open_file(jet_type, labels[1])
    f3 = open_file(jet_type, labels[2])

    for hname in hist_names:
        h1 = get_hist(f1, hname, jet_type, labels[0])
        h2 = get_hist(f2, hname, jet_type, labels[1])
        h3 = get_hist(f3, hname, jet_type, labels[2])

        #h1.Scale(1.0/h1.Integral())
        #h2.Scale(1.0/h2.Integral())
        #h3.Scale(1.0/h3.Integral())

        #h1.Rebin(10)
        #h2.Rebin(10)
        #h3.Rebin(10)
        #h1.Scale(1.0/10)
        #h2.Scale(1.0/10)
        #h3.Scale(1.0/10)

        plot_overlay(h1,h2,labels[0],labels[1],f"{hname} ({jet_type} bkg cut)",f"{hname}_{jet_type}_bkg_cut_{labels[0]}_vs_{labels[1]}_binned")
        plot_ratio(h1,h2,labels[0],labels[1],f"{hname} Ratio ({jet_type} bkg cut)",f"{hname}_{jet_type}_bkg_cut_{labels[0]}_vs_{labels[1]}_binned")
        plot_overlay(h1,h3,labels[0],labels[2],f"{hname} ({jet_type} bkg cut)",f"{hname}_{jet_type}_bkg_cut_{labels[0]}_vs_{labels[2]}_binned")
        plot_ratio(h3,h1,labels[2],labels[0],f"{hname} Ratio ({jet_type} bkg cut)",f"{hname}_{jet_type}_bkg_cut_{labels[0]}_vs_{labels[2]}_binned")
        
    f1.Close()
    f2.Close()
    f3.Close()
