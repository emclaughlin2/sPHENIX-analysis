import ROOT
import numpy as np

ROOT.gStyle.SetOptStat(0)

# --------------------------------------------------
# Dijet cut efficiency
# --------------------------------------------------

cuts = [0.3, 0.4, 0.5]

c_eff = ROOT.TCanvas("c_eff", "", 800, 600)

leg = ROOT.TLegend(0.60, 0.20, 0.85, 0.40)
leg.SetBorderSize(0)
leg.SetFillStyle(0)

colors = [
    ROOT.kBlack,
    ROOT.kRed+1,
    ROOT.kBlue+1,
    ROOT.kGreen+2
]

first = True
h_pass = []
h_total = []
h_eff = []

for i, cut in enumerate(cuts):

    fname = (
        f"sphenix_primary_analysis_sim_run28_output/"
        f"output_pt3_lt_{cut}_dijet_bkg_cut_sim_iter_1.root"
    )

    f = ROOT.TFile.Open(fname)

    h_pass.append(f.Get("h_exclusive_dijet_passcut"))
    h_total.append(f.Get("h_exclusive_dijet_total"))
    h_pass[i].SetDirectory(0)
    h_total[i].SetDirectory(0)

    h_eff.append(h_pass[i].Clone(f"h_eff_{cut}"))
    h_eff[-1].SetDirectory(0)
    h_eff[-1].Divide(h_pass[i], h_total[i], 1.0, 1.0, "B")

    h_eff[i].SetLineColor(colors[i])
    h_eff[i].SetMarkerColor(colors[i])
    h_eff[i].SetMarkerStyle(20+i)

    h_eff[i].GetYaxis().SetTitle("p_{T,3} Cut Efficiency")
    h_eff[i].GetXaxis().SetTitle("p_{T,lead} [GeV]")
    h_eff[i].GetYaxis().SetRangeUser(0.5, 1.5)

    if first:
        h_eff[i].Draw("E1")
        first = False
    else:
        h_eff[i].Draw("E1 SAME")
    legvalue = "p_{T,3} < "+f"{cut}"+"p_{T,1}"
    leg.AddEntry(h_eff[i], legvalue, "lep")

leg.Draw()

c_eff.SaveAs("dijet_cut_efficiency_variations.pdf")

c_eff_ue = ROOT.TCanvas("c_eff_ue", "", 800, 600)

leg1 = ROOT.TLegend(0.60, 0.65, 0.85, 0.85)
leg1.SetBorderSize(0)
leg1.SetFillStyle(0)

colors = [
    ROOT.kBlack,
    ROOT.kRed+1,
    ROOT.kBlue+1,
    ROOT.kGreen+2
]

first = True
h2_ue_pt = []
h_ue_pt = []

for i, cut in enumerate(cuts):

    fname = (
        f"sphenix_primary_analysis_sim_run28_output/"
        f"output_pt3_lt_{cut}_dijet_bkg_cut_sim_iter_1.root"
    )

    f = ROOT.TFile.Open(fname)

    h2_ue_pt.append(f.Get("h_ue_pt_truth_transverse_record"))
    h2_ue_pt[i].SetDirectory(0)

    h_ue_pt.append(h2_ue_pt[i].ProfileX(f"h_ue_pt_cut_{cut}"))
    h_ue_pt[i].SetDirectory(0)

    h_ue_pt[i].Scale(3.0/(2.2*2*np.pi))

    h_ue_pt[i].SetLineColor(colors[i])
    h_ue_pt[i].SetMarkerColor(colors[i])
    h_ue_pt[i].SetMarkerStyle(20+i)

    h_ue_pt[i].GetYaxis().SetTitle("<#SigmaE_{T}^{truth}> [GeV]")
    h_ue_pt[i].GetXaxis().SetTitle("p_{T,lead}^{truth} [GeV]")
    h_ue_pt[i].GetYaxis().SetRangeUser(0, 0.6)

    if first:
        h_ue_pt[i].Draw("E1")
        first = False
    else:
        h_ue_pt[i].Draw("E1 SAME")
    legvalue = "p_{T,3} < "+f"{cut}"+"p_{T,1}"
    leg1.AddEntry(h_ue_pt[i], legvalue, "lep")

leg1.Draw()

c_eff_ue.SaveAs("dijet_cut_efficiency_ue_distribution.pdf")

# --------------------------------------------------
# UE ProfileX comparison
# --------------------------------------------------

f_pythia = ROOT.TFile.Open(
    "sphenix_primary_analysis_sim_run28_output/"
    "output_none_bkg_cut_sim_iter_1.root"
)

f_herwig = ROOT.TFile.Open(
    "sphenix_primary_analysis_sim_run28_output/"
    "output_none_bkg_cut_sim_iter_1_herwig.root"
)

h2_pythia = f_pythia.Get("h_ue_pt_truth_transverse_record_full")
h2_herwig = f_herwig.Get("h_ue_pt_truth_transverse_record_full")

p_pythia = h2_pythia.ProfileX("p_pythia")
p_herwig = h2_herwig.ProfileX("p_herwig")

p_pythia.Scale(3.0)
p_herwig.Scale(3.0)

p_pythia.SetLineColor(ROOT.kBlue+1)
p_pythia.SetMarkerColor(ROOT.kBlue+1)
p_pythia.SetMarkerStyle(20)

p_herwig.SetLineColor(ROOT.kRed+1)
p_herwig.SetMarkerColor(ROOT.kRed+1)
p_herwig.SetMarkerStyle(24)

c_prof = ROOT.TCanvas("c_prof", "", 800, 600)

p_pythia.GetYaxis().SetRangeUser(0,13)
p_pythia.GetYaxis().SetTitle("<#SigmaE_{T}> Total Calo Accept. [GeV]")
p_pythia.GetXaxis().SetTitle("p_{T,lead}^{truth} [GeV]")
p_pythia.Draw("E1")
p_herwig.Draw("E1 SAME")

leg2 = ROOT.TLegend(0.60, 0.75, 0.85, 0.88)
leg2.SetBorderSize(0)
leg2.SetFillStyle(0)
leg2.AddEntry(p_pythia, "Pythia", "lep")
leg2.AddEntry(p_herwig, "Herwig", "lep")
leg2.Draw()

c_prof.SaveAs("pythia_herwig_ue_pt_profile_comparison.pdf")