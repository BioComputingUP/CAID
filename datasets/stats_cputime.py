

import json
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt


# sudo mount 172.21.2.103:/volume1/Pleiadi/projects/CAID/qsub /home/damiano/Projects/caid/data/qsub/
# sudo mount 172.21.2.103:/volume1/Pleiadi/projects/CAID/run /home/damiano/Projects/caid/data/run/


targets = {}
with open("../data/qsub_hhblits_out") as f:
    for line in f:
        if line.startswith("Start") or line.startswith("End"):
            if len(line.split()) == 9:
                # Start chimera Tue Oct 2 19:14:07 CEST 2018 DP00085.fasta
                # End ortro Tue Oct 2 19:22:05 CEST 2018 DP00080.fasta
                _, _, _, m, d, tm, _, _, disprot_id = line.strip().split()
                disprot_id = disprot_id.split('.')[0]

                targets.setdefault(disprot_id, [None, None])
                date = datetime.datetime.strptime("{},{},{}".format(m, d, tm), "%b,%d,%H:%M:%S")

                if line.startswith("Start") and targets[disprot_id][0] is None:
                    targets[disprot_id][0] = date
                elif line.startswith("End") and targets[disprot_id][1] is None:
                    targets[disprot_id][1] = date

time_hhblits = [(targets[target][1] - targets[target][0]).total_seconds() for target in targets if targets[target][1] is not None and targets[target][0] is not None]

# qacct -j qsub_psiblast.sh | grep '_time\|jobnumber'|paste - - - - > qsub_psiblast_time.txt
time_psiblast = []
with open("../data/qsub_psiblast_time.txt") as f:
    for line in f:
        # jobnumber    31154               	qsub_time    Mon Oct  1 12:17:06 2018	start_time   Thu Oct  4 22:44:38 2018	end_time     Thu Oct  4 22:51:30 2018
        _, _, _, _, _, _, _, _, _, _, s_m, s_d, s_time, _, _, _, e_m, e_d, e_time, _ = line.strip().split()
        start = datetime.datetime.strptime("{},{},{}".format(s_m, s_d, s_time), "%b,%d,%H:%M:%S")
        end = datetime.datetime.strptime("{},{},{}".format(e_m, e_d, e_time), "%b,%d,%H:%M:%S")
        time_psiblast.append((end - start).total_seconds())


with open("../data/method_names.json") as f:
    methods_list = json.load(f)

time_methods = []
time_methods_labels = []
for method in sorted(methods_list, key=lambda k: k['name']):
    outfile = None
    if os.path.isfile("../data/qsub/out/{}".format(method['dir'])):
         outfile = "../data/qsub/out/{}".format(method['dir'])
    elif os.path.isfile("../data/run/out/{}".format(method['dir'])):
        outfile = "../data/run/out/{}".format(method['dir'])
    else:
        print("Missing file {}".format(method['dir']))

    if outfile is not None:
        targets = {}
        with open(outfile) as f:
            for line in f:
                if line.startswith("job-"):
                    # print(line)
                    _, disprot_id, _, start_date, start_time = line.strip().split()
                    # 2019:01:15 15:10:55.139
                    date = datetime.datetime.strptime("{},{}".format(start_date, start_time), "%Y:%m:%d,%H:%M:%S.%f")
                    targets.setdefault(disprot_id, [None, None])

                    if line.startswith("job-start") and targets[disprot_id][0] is None:
                        targets[disprot_id][0] = date
                    elif line.startswith("job-end") and targets[disprot_id][1] is None:
                        targets[disprot_id][1] = date

        time_methods.append(np.array([(targets[target][1] - targets[target][0]).total_seconds() for target in targets if targets[target][1] is not None and targets[target][0] is not None]))
        time_methods_labels.append(method)


time_colors = ["white" for method in time_methods_labels]
time_psiblast_mean = np.mean(time_psiblast)
time_hhblits_mean = np.mean(time_hhblits)

for i, method in enumerate(time_methods_labels):
    # if method['dir'] == 'mobidblite':
    #     time_colors[i] = 'blue'

    if 'evolution' in method:
        time_colors[i] = 'red'

    # Add Psiblast time
    if "input" in method and "pssm" in method["input"]:
        time_methods[i] += time_psiblast_mean

    # Add HHblits time
    if "input" in method and "hhblits" in method["input"]:
        time_methods[i] += time_hhblits_mean


fig, ax = plt.subplots()
bplot = ax.boxplot(time_methods, labels=["{}{}".format(method['name'], '*' if method["dir"] == "mobidblite" else '') for method in time_methods_labels], sym="", patch_artist=True)
# ax.set_ylim(0, 4000)
ax.set_ylabel("Seconds")
ax.set_yscale("log")
plt.xticks(rotation=90)

for patch, color in zip(bplot['boxes'], time_colors):
    patch.set_facecolor(color)

plt.tight_layout()
plt.savefig('cputime.png'.format(), dpi=300)
