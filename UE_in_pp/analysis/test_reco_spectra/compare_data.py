import ROOT
import glob
import os
import re
import math

ROOT.gROOT.SetBatch(True)

# tolerance for floating point comparison of means
MEAN_TOL = 1e-4

# Directory containing the ROOT files
indir = "output"

# Find all hanpu files
hanpu_files = glob.glob(os.path.join(indir, "hanpu_jet10_output_r04_*_0_0.root"))

# Regex to extract runnumber
pattern = re.compile(r"hanpu_jet10_output_r04_(\d+)_0_0\.root")

print(f"{'Run':>8} | {'Hist':<22} | "
      f"{'Hanpu Entries':>14} {'Hanpu Mean':>12} | "
      f"{'Emma Entries':>14} {'Emma Mean':>12}")
print("-" * 90)

for hanpu_file in hanpu_files:
    m = pattern.search(os.path.basename(hanpu_file))
    if not m:
        continue

    runnumber = m.group(1)
    emma_file = os.path.join(
        indir, f"emma_jet10_output_r04_{runnumber}_0_0.root"
    )

    # Require both files to exist
    if not os.path.exists(emma_file):
        continue

    f_hanpu = ROOT.TFile.Open(hanpu_file)
    f_emma  = ROOT.TFile.Open(emma_file)

    if not f_hanpu or f_hanpu.IsZombie():
        continue
    if not f_emma or f_emma.IsZombie():
        continue

    for hname in ["h_leadingunsubjet_pt", "h_zvertex"]:
        h_hanpu = f_hanpu.Get(hname)
        h_emma  = f_emma.Get(hname)

        if not h_hanpu or not h_emma:
            continue

        hanpu_entries = h_hanpu.GetEntries()
        hanpu_mean    = h_hanpu.GetMean()

        emma_entries  = h_emma.GetEntries()
        emma_mean     = h_emma.GetMean()

        entries_diff = hanpu_entries != emma_entries
        mean_diff    = not math.isclose(hanpu_mean, emma_mean, rel_tol=0, abs_tol=MEAN_TOL)

        if entries_diff or mean_diff:
            print(f"{runnumber:>8} | {hname:<22} | "
              f"{hanpu_entries:14.0f} {hanpu_mean:12.4f} | "
              f"{emma_entries:14.0f} {emma_mean:12.4f}")

    f_hanpu.Close()
    f_emma.Close()
