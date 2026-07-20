runs = []
lumis = []

with open('runlist_for_lumi.list', 'r') as f:
    for line in f:
        runs.append(line.strip())

with open('60cmlumi_allxing.list', 'r') as f:
    for line in f:
        if len(line.strip().split()) > 0:
            if line.strip().split()[0] in runs:
                print(line.strip())
                lumis.append(float(line.strip().split()[5]))

for run, lumi in zip(runs, lumis):
    print(run, lumi)

print(sum(lumis))