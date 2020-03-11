
import json
import gzip
import logging.config

# Set logging
logging.basicConfig(format='%(asctime)s - %(process)d - %(name)s - %(message)s',
                    level=logging.getLevelName("INFO"))

logging.info("parse_mobidb.py started")


caid_obj = {}
with open("data/caid.json") as f:
    for doc in json.load(f)["data"]:
        if doc["name"] != "Genome polyprotein":  # Filter virus "Genome polyproteins"
            caid_obj[doc["acc"]] = doc

with gzip.open("data/derived_disorder.mjson.gz") as f:
    for line in f:
        obj = json.loads(line)
        if obj["acc"] in caid_obj:
            if obj["sequence"] == caid_obj[obj["acc"]]["sequence"]:
                regions = []
                for ele in obj["mobidb_consensus"]["disorder"]["derived"]:
                    # Consider all types of experiments
                    if ele["method"] == "missing_residues":
                        for start, end, label in ele["regions"]:
                            # if label in ["S"]:  # high confidence 90% agreement as defined in MobiDB
                            if label in ["S", "C"]:  # at least
                                regions.append((start, end, label))
                if regions:
                    seq = ['-'] * len(obj["sequence"])
                    for start, end, type in regions:
                        for i in range(start - 1, end):
                            seq[i] = "S"

                    print(">{} {}\n{}".format(caid_obj[obj["acc"]]["disprot_id"], obj["acc"], "".join(seq)))
            else:
                logging.warning("Different sequences {} {}: Disprot {}, Mobidb {}".format(caid_obj[obj["acc"]]["disprot_id"], obj["acc"], len(caid_obj[obj["acc"]]["sequence"]), len(obj["sequence"])))

# python3 parse_mobidb.py > data/pdb-atleast.fasta