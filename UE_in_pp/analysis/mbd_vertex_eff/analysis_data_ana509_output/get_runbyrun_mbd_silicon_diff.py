import ROOT
import re
import glob

# List of your files (you can also use glob if stored in a folder)
file_list = [
    "output_zvtx_lt_10.000000_100.000000_efrac_53018.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53046.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53079.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53081.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53195.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53196.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53513.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53517.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53530.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53531.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53532.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53571.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53586.root",
    "output_zvtx_lt_10.000000_100.000000_efrac_53590.root"
]

for fname in file_list:
    # Extract the last 5 digits (run number)
    match = re.search(r'(\d{5})\.root$', fname)
    runnumber = match.group(1) if match else "Unknown"

    f = ROOT.TFile.Open(fname)
    if not f or f.IsZombie():
        print(f"Could not open {fname}")
        continue

    hist = f.Get("h_measure_mbd_silicon_diff")
    if not hist:
        print(f"Histogram not found in {fname}")
        f.Close()
        continue

    mean = hist.GetMean()
    rms = hist.GetRMS()
    entries = hist.GetEntries()
    print(f"Run {runnumber}: Mean = {mean:.4f}, RMS = {rms:.4f}, Entries = {entries}")

    f.Close()

