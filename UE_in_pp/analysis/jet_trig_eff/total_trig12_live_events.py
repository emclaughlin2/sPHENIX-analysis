import ROOT
import glob

# Pattern matching your files
#file_pattern = "/sphenix/tg/tg01/jets/egm2153/UEinppOutput/sim_run28_jet10_total_events_zvtx_lt_60_output_*.root"
file_pattern = "ana509_analysis_output/dijet_and_mbdtimecut_output_4*.root"
file_list = glob.glob(file_pattern)

total_entries = 0

for fname in file_list:
    f = ROOT.TFile.Open(fname)
    if not f or f.IsZombie():
        print(f"Could not open {fname}")
        continue
    h = f.Get("h_vz_mb_wzcut_online")
    if not h:
        print(f"h_vz_mb_wzcut_online not found in {fname}")
        continue
    entries = h.GetEntries()
    print(f"{fname}: {entries} entries")
    total_entries += entries

    f.Close()

print(f"\nTotal number of entries across all files: {total_entries}")