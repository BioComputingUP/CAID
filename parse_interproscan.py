"""
InterProScan (5.38-76.0) output generated with the following command:
interproscan.sh -f tsv -dp -iprlookup -T /tmp/ -i caid.fasta -o caid_interproscan -T /tmp/
python parse_interproscan.py ../data/annotations/disprot-2018-11-disorder.fasta ../data/annotations/disprot-2018-11.json caid_interproscan  > ../data/annotations/data/gene3d.fasta
"""

import argparse
import logging.config
import json


def gen_block(f):
    """
    Parse and split the input.
    The input must be sorted by target name, first column.
    """
    name, old_name = None, None
    chunk = []

    for line in f:
        if line and line[0] != '#':
            # name, term, ontology, score, identities (optional), uniprot_ids (optional)
            name = line.split()[0]
            if name != old_name and old_name:
                logging.debug("block {}\n{}".format(old_name, chunk))
                yield (old_name, chunk)
                chunk = []
            old_name = name
            chunk.append(line.strip())
    if old_name:
        logging.debug("block {}\n{}".format(old_name, chunk))
        yield (old_name, chunk)


def parse_interproscan(disprot_fasta, disprot_json, interproscan):
    logging.info("parse_interproscan.py started")

    caid_fasta = {}
    with open(disprot_fasta) as f:
        for line in f:
            if line[0] == ">":
                name = line.split()[0][1:]
            else:
                caid_fasta[name] = line.strip()

    # The JSON as downloaded from DisProt
    caid_obj = {}
    with open(disprot_json) as f:
        for doc in json.load(f)["data"]:
            if doc["name"] != "Genome polyprotein":  # Filter virus "Genome polyproteins"
                caid_obj[doc["disprot_id"]] = doc

    # Check the sequnces used as input for predictors are the same as those available from the DisProt website
    for disprot_id in caid_fasta:
        if caid_fasta[disprot_id] != caid_obj[disprot_id]["sequence"]:
            logging.warning("{} different sequence".format(disprot_id))

    # Print the gene3d merged regions
    with open(interproscan) as f:
        for n, lines in gen_block(f):
            if n in caid_obj:  # Filter virus "Genome polyproteins"
                regions = []
                for line in lines:
                    target, _, length, method, label, name, start, end, score, _, _ = line.split("\t")[:11]
                    if method == "Gene3D":
                        regions.append((int(start), int(end), label))
                if regions:
                    seq = ['-'] * len(caid_obj[n]["sequence"])
                    for start, end, type in regions:
                        for i in range(start - 1, end):
                            seq[i] = "S"

                    print(">{} {}\n{}".format(n, caid_obj[n]["acc"], "".join(seq)))


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('disprotfasta', help='directory where CAID predictors results are saved')
    parser.add_argument('disprotjson', help='directory where CAID predictors results are saved')
    parser.add_argument('interproscan', help='directory where CAID predictors results are saved')

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    # Set logging
    logging.basicConfig(format='%(asctime)s - %(process)d - %(name)s - %(message)s',
                        level=logging.getLevelName("INFO"))

    parse_interproscan(args.disprotfasta, args.disprotjson, args.interproscan)
