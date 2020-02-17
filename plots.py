from pathlib import Path

# obtain from cli arg
resultdir = Path("/home/marnec/Projects/CAID/caid/results")


# iterate over file in dir (foreach reference)
reference = "/home/marnec/Projects/CAID/caid/data/disorder/new-disprot-all_simple.txt"
# for body start
reference = Path(reference)
refname = reference.stem

dataset_metrics = resultdir / "{}.analysis.all.dataset.default.metrics.csv".format(refname)

print(dataset_metrics)



