

# Combine output files to speed up the parsing
# find . -maxdepth 1 -type f -name '*' -print0 | xargs -0 cat > ../pairwise.txt

new_targets = []
with open("../data/new_id_list.txt") as f:
    for line in f:
        new_targets.append(line.strip())

old_targets = []
with open("../data/old_id_list.txt") as f:
    for line in f:
        old_targets.append(line.strip())

all_targets = new_targets + old_targets

new_targets = frozenset(new_targets)
old_targets = frozenset(old_targets)
all_targets = frozenset(all_targets)


seq_len = {}
with open("../data/disprot8_all.fasta") as f:
    for line in f:
        if line.startswith(">"):
            name = line.split()[0][1:]
        else:
            seq_len[name] = len(line.strip())


"""
# 1: DP00003
# 2: DP00028
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 596
# Identity:      21/596 ( 3.5%)
# Similarity:    25/596 ( 4.2%)
# Gaps:         545/596 (91.4%)
# Score: 49.5
"""
data = {}
with open("../data/pairwise.txt") as f:
    for i, line in enumerate(f):
        if line.startswith("# 1: "):
            p1 = line.strip().split()[2]
        elif line.startswith("# 2: "):
            p2 = line.strip().split()[2]
        elif line.startswith("# Length: "):
            l = int(line.strip().split()[2])
        elif line.startswith("# Identity: "):
            identity = int(line.strip().split()[2].split("/")[0])
        elif line.startswith("# Similarity: "):
            similarity = int(line.strip().split()[2].split("/")[0])
        elif line.startswith("# Gaps: "):
            gaps = int(line.strip().split()[2].split("/")[0])
        elif line.startswith("# Score: "):
            score = float(line.strip().split()[2])

            if p1 != p2 and p1 in all_targets and p2 in all_targets:  # Just a check, should not happen
                identity_1 = float(identity) / seq_len[p1]
                identity_2 = float(identity) / seq_len[p2]
                labels = ["all"]  # dataset to consider

                # Inside new
                if p1 in new_targets and p2 in new_targets:
                    labels.append("new")
                # Inside old
                if p1 not in new_targets and p2 not in new_targets:
                    labels.append("old")
                # New vs old
                if p1 in new_targets and p2 not in new_targets:
                    labels.append("new_old")
                # New vs old
                if p2 in new_targets and p1 not in new_targets:
                    labels.append("new_old")

                for label in labels:
                    data.setdefault(label, {})

                    data[label].setdefault(p1, 0.0)
                    if identity_1 > data[label][p1]:
                        data[label][p1] = identity_1

                    data[label].setdefault(p2, 0.0)
                    if identity_2 > data[label][p2]:
                        data[label][p2] = identity_2

        if i % 1000000 == 0:
            print(i)

with open("../data/pairwise_identity.txt", "w") as fout:
    for label in data:
        for p in data[label]:
            if label != 'new_old' or p in new_targets:  # Exclude old from the new_old count
                fout.write("{} {} {}\n".format(label, p, data[label][p]))
