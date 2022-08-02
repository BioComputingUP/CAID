'''python parse_mobidb.py ../data/annotations/disprot-2018-11.json ../data/annotations/derived_disorder.mjson.gz > ../data/annotations/pdb-atleast.fasta'''

import argparse
import json
import gzip
import logging.config


def parse_mobidb(disprot, derived):
    logging.info("parse_mobidb.py started")
    caid_obj = {}
    with open(disprot) as f:
        for doc in json.load(f)["data"]:
            if doc["name"] != "Genome polyprotein":  # Filter virus "Genome polyproteins"
                caid_obj[doc["acc"]] = doc

    with gzip.open(derived) as f:
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


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('disprot', help='directory where CAID predictors results are saved')
    parser.add_argument('derived', help='directory where CAID predictors results are saved')

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # Set logging
    logging.basicConfig(format='%(asctime)s - %(process)d - %(name)s - %(message)s',
                        level=logging.getLevelName("INFO"))
    args = parse_args()
    parse_mobidb(args.disprot, args.derived)