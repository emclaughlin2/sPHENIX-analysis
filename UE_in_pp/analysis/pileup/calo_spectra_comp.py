import ROOT
import os
import array

# Load sPHENIX style if available
def load_sphenix_style():
    style_path = os.path.expandvars("/sphenix/u/egm2153/spring_2023/sPhenixStyle.C")
    if os.path.exists(style_path):
        ROOT.gROOT.ProcessLine(f'.L {style_path}')
        ROOT.gROOT.ProcessLine('SetsPhenixStyle()')
    else:
        print("Warning: sPHENIX style file not found, using default ROOT style.")

load_sphenix_style()

# File info: (filename, label)
files = [
    ("pileup_inclusive_jet_UE_analysis_test_Oct2025_calo_tower_sum_dijet_bkg_cut_maxrate_0.020000.root", "dijet, <PU> < 0.02"),
    ("pileup_inclusive_jet_UE_analysis_test_Oct2025_calo_tower_sum_dijet_bkg_cut_maxrate_0.050000.root", "dijet, 0.03 < <PU> < 0.05"),
    ("pileup_inclusive_jet_UE_analysis_test_Oct2025_calo_tower_sum_efrac_bkg_cut_maxrate_0.020000.root", "efrac, <PU> = 0.02"),
    ("pileup_inclusive_jet_UE_analysis_test_Oct2025_calo_tower_sum_efrac_bkg_cut_maxrate_0.050000.root", "efrac, 0.03 < <PU> < 0.05"),
]

calos = ["emcal", "ihcal", "ohcal"]
regions = ["towards", "transverse", "away"]

# Colors and styles for distinction
colors = [ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kOrange+7]
styles = [1, 2, 3, 4]

# Read histograms
histos = {}  # (calo, region): [h1, h2, h3, h4]
for calo in calos:
    for region in regions:
        key = (calo, region)
        histos[key] = []
        hname = f"h_{calo}_et_{region}_spectra"
        for fname, label in files:
            f = ROOT.TFile.Open(fname)
            h = f.Get(hname)
            if not h:
                print(f"Warning: {hname} not found in {fname}")
                histos[key].append(None)
            else:
                h.SetDirectory(0)  # Detach from file
                histos[key].append(h)
            f.Close()

# Normalize histograms by number of entries
for key in histos:
    for h in histos[key]:
        if h and h.GetEntries() > 0:
            h.Scale(1.0 / h.GetEntries())

ROOT.gStyle.SetOptStat(0)
for calo in calos:
    for region in regions:
        c = ROOT.TCanvas(f"c_{calo}_{region}", f"{calo} {region} spectra", 800, 600)
        c.SetLogy()
        leg = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)
        leg.SetTextSize(0.03)
        leg.AddEntry(0, f"{calo.upper()} {region.capitalize()} spectra", "")
        first = True
        for i, (h, (fname, label)) in enumerate(zip(histos[(calo, region)], files)):
            if not h:
                continue
            h.SetLineColor(colors[i])
            #h.SetLineStyle(styles[i])
            h.SetLineWidth(2)
            h.SetTitle(f"{calo.upper()} {region.capitalize()} spectra")
            h.GetXaxis().SetTitle("E_{T} [GeV]")
            h.GetYaxis().SetTitle("Entries")
            if first:
                h.Draw("hist")
                first = False
            else:
                h.Draw("hist same")
            leg.AddEntry(h, label, "l")
        leg.Draw()
        c.Update()
        c.SaveAs(f"1.5mrad_jet_bkg_plots/1.5mrad_pu_bkg_testing_{calo}_{region}_spectra_comp.png")

        # Plot ratio of ohcal_towards for the two efrac distributions

# Indices for efrac files in 'files' list
efrac_indices = [2, 3]
key = ("ohcal", "towards")
h1 = histos[key][efrac_indices[0]]
h2 = histos[key][efrac_indices[1]]

if h1 and h2:
    ratio = h1.Clone("ratio_ohcal_towards_efrac")
    ratio.SetTitle("OHCal Towards: efrac, PU=0.02 / PU=0.05")
    ratio.Divide(h2)
    c_ratio = ROOT.TCanvas("c_ratio_ohcal_towards_efrac", "OHCal Towards Ratio (efrac)", 800, 600)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(2)
    ratio.GetXaxis().SetTitle("E_{T} [GeV]")
    ratio.GetYaxis().SetTitle("Ratio (PU=0.02 / PU=0.05)")
    ratio.Draw("hist")
    c_ratio.Update()
    c_ratio.SaveAs("1.5mrad_jet_bkg_plots/1.5mrad_pu_bkg_testing_ohcal_towards_efrac_ratio.png")
else:
    print("Warning: Could not create ratio plot for ohcal_towards efrac distributions.")