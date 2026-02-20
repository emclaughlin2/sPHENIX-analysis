import ROOT
ROOT.gROOT.SetBatch(True)

#types = ["10", "15", "20", "30", "50"]
types = ["jet20", "jet30", "jet50"]
labels = ["21", "28"]

hist_names = [
    "h_leadingunsubjet_pt",
    "h_leadingtruthjet_pt",
    "h_zvertex",
    #"h_leadingunsubjet_eta",
    #"h_leadingtruthjet_eta",
    #"h_leadingtruthjet_zvertex",
    #"h_leadingunsubjet_phi",
    #"h_leadingtruthjet_phi",
]

cross_sections = {
    "10": 3.997,
    "15": 4.073e-1,
    "jet20": 6.218e-2,
    "jet30": 2.502e-3,
    "jet50": 7.2695e-6,
}

# ------------------------------------------------------
def open_file(jet_type, label):
    fname = f"output_after_bug_fix/emma_{jet_type}_output_r04_{label}.root"
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

    unique_name = f"{name}_jet{jet_type}_r04_{label}"
    h_clone = h.Clone(unique_name)
    h_clone.SetDirectory(0)

    return h_clone

# ------------------------------------------------------
def plot_overlay(h21, h28, title, outname):
    c = ROOT.TCanvas(outname, "", 800, 600)

    h21.SetLineColor(ROOT.kRed)
    h28.SetLineColor(ROOT.kBlue)

    h21.SetLineWidth(2)
    h28.SetLineWidth(2)

    h21.SetTitle(title)
    h21.Draw("hist")
    h28.Draw("hist same")

    leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    leg.AddEntry(h21, "r04_21", "l")
    leg.AddEntry(h28, "r04_28", "l")
    leg.Draw()

    c.SaveAs("compare_run21_run28_plots/bug_fix_test_" + outname + ".pdf")

def plot_ratio(h21, h28, title, outname):
    # Clone to avoid modifying originals
    ratio = h28.Clone(outname + "_ratio")
    ratio.SetDirectory(0)

    # Safe division
    ratio.Divide(h21)

    c = ROOT.TCanvas(outname + "_ratio", "", 800, 600)

    ratio.SetTitle(title)
    ratio.SetLineWidth(2)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.GetYaxis().SetTitle("28 / 21")
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

    c.SaveAs("compare_run21_run28_plots/bug_fix_test_" + outname + "_ratio.pdf")

# ======================================================
# 1) Plot 21 vs 28 for each type
# ======================================================

for jet_type in types:
    f21 = open_file(jet_type, "21")
    f28 = open_file(jet_type, "28")

    for hname in hist_names:
        h21 = get_hist(f21, hname, jet_type, "21")
        h28 = get_hist(f28, hname, jet_type, "28")

        plot_overlay(
            h21,
            h28,
            f"{hname} (Jet {jet_type} GeV)",
            f"{hname}_Jet{jet_type}_21_vs_28",
        )

        plot_ratio(
            h21,
            h28,
            f"{hname} Ratio (Jet {jet_type} GeV)",
            f"{hname}_Jet{jet_type}_28_over_21",
        )


    f21.Close()
    f28.Close()

# ======================================================
# 2) Scale and sum histograms
# ======================================================

combined = {
    "21": {},
    "28": {},
}

# initialize containers
combined = {label: {hname: None for hname in hist_names} for label in labels}

for jet_type in types:
    scale = cross_sections[jet_type]

    for label in labels:
        f = open_file(jet_type, label)

        for hname in hist_names:
            h = get_hist(f, hname, jet_type, label)  # already cloned + detached
            h.Scale(scale)

            if combined[label][hname] is None:
                combined[label][hname] = h  # TAKE OWNERSHIP
            else:
                combined[label][hname].Add(h)

        f.Close()

# ======================================================
# 3) Plot combined 21 vs 28
# ======================================================

for hname in hist_names:
    h21 = combined["21"][hname]
    h28 = combined["28"][hname]

    plot_overlay(
        h21,
        h28,
        f"{hname} (Combined, scaled)",
        f"{hname}_combined_21_vs_28",
    )

    plot_ratio(
        h21,
        h28,
        f"{hname} Ratio (Combined, scaled)",
        f"{hname}_combined_28_over_21",
    )

