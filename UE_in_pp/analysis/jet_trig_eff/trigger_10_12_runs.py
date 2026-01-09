import glob
import re
import csv

# Pattern matching your files
file_pattern = "output_logs/output_*.out"
file_list = glob.glob(file_pattern)

pattern = "trigger 10 entries ([0-9]+) trigger 12 entries ([0-9]+)"
run_pattern = r"Processing jet_trig_eff_analysis.C\(([0-9]+)\)" 
#file_list = ['output_50000.out']

run_trigger = [];

for fname in file_list:
    with open(fname) as my_file:
        run = '-9999'
        for line in my_file:
            run_match = re.findall(run_pattern, line)
            matches = re.findall(pattern, line)
            if run_match:
                run = run_match[0]
            if matches: 
                run_trig = {'run': run, 'trig10': matches[0][0], 'trig12': matches[0][1]};
                run_trigger.append(run_trig)
            
with open("run_trigger_info.csv", 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=run_trigger[0].keys())
    writer.writeheader()
    writer.writerows(run_trigger)
