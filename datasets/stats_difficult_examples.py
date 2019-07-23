import numpy as np
import matplotlib.pyplot as plt


reference = {}
c = 0
with open("../data/disorder/new-disprot-all_simple.txt") as f:
    for line in f:
        if line[0] != '#':
            if line[0] == '>':
                name = line.strip().split()[0][1:]
                c = 0
            elif c == 1:
                pass
            elif c == 2:
                reference[name] = np.array(list(line.strip()), dtype=np.int)
            c += 1
# print(reference)

consensus = []
with open("../data/new-disprot-all_simple_predstack.csv") as f:
    header = next(f)
    for line in f:
        if line.startswith("D"):
            pass
        else:
            consensus.append(list(line.strip().split()[1]))
consensus = np.array(consensus)
consensus = consensus.astype(np.int)  # convert string into int
mean = np.mean(consensus, axis=0)  # calculate consensus score
# print(consensus.shape)
# print(consensus)
# print(mean)


examples = ["DP01196", "DP01128", "DP01339", "DP01883", "DP01971", "DP01248"]

for head in header.split("#"):
    if head:

        disprot_id, start, end = head.split(",")
        start = int(start)
        end = int(end)
        x = np.arange(end - start)

        # print(head)
        # print(reference[disprot_id])
        # print(mean[start:end])
        # print(x)

        if disprot_id in examples:

            fig, axes = plt.subplots(figsize=(7, 4))

            axes.bar(x, reference[disprot_id], width=1.0, label="DisProt", color='tomato')
            axes.bar(x, mean[start:end], width=1.0, label="Predictions", alpha=0.8, color='navy')
            axes.set_xlabel("Sequence")
            axes.set_ylabel("Consensus")
            axes.set_title(disprot_id)

            plt.legend()

            plt.tight_layout()
            plt.savefig('../data/example_{}.png'.format(disprot_id), dpi=300)


