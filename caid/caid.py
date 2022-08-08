# module imports
import argparse
import logging
import os
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from vectorized_metrics.vectorized_metrics import bvaluation


def parse_args():
    parser = argparse.ArgumentParser(
            prog='caid-assess', description="CAID: Critical Assessment of Intrinsic Disorder",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('references_path', help='Path to the references file')

    parser.add_argument('predictions', help="directory containing prediction file(s)")

    parser.add_argument('-r', '--refList', help='reference csv file with associations between predictor and reference', required=True)

    parser.add_argument('-o', '--outputDir', default='.',
                        help='directory where the output will be written (default: cwd)')

    parser.add_argument('-b', '--labels', default=None, help='filename with labels')

    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="WARNING",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels >= choice will be displayed')

    args = parser.parse_args()
    return args


# def parse_config(config_file):
#     config_parser = configparser.ConfigParser()
#     config_parser.optionxform = str
#     config_parser.read(config_file)
#
#     return config_parser


def set_logger(logfile, level):
    logging.basicConfig(level=level,
                        format='%(asctime)s | %(module)-13s | %(levelname)-8s | %(message)s',
                        stream=open(logfile) if logfile else sys.stderr)


if __name__ == '__main__':
    args = parse_args()
    set_logger(args.log, args.logLevel)
    all_pred_paths = list(Path(args.predictions).glob('*.caid'))
    all_ref_paths = list(Path(args.references_path).glob('*.fasta'))

    pred_references = pd.read_csv(args.refList, sep=',')

    references = pred_references.columns[1:]

    os.makedirs(args.outputDir, exist_ok=True)

    for reference in tqdm(references, desc='Bvaluation on references'):
        pred_for_reference = pred_references[pred_references[reference] == 1].Method.to_list()
        pred_paths = list(filter(lambda x: x.stem in pred_for_reference, all_pred_paths))

        if len(pred_paths) == 0:
            raise ValueError(f'No predictors found for reference {reference}')

        ref_path = list(filter(lambda x: x.stem == reference, all_ref_paths))

        if len(ref_path) == 0:
            raise ValueError(f'No reference found for {reference}')

        bvaluation(str(ref_path[0]), pred_paths, outpath=args.outputDir, dataset=True, bootstrap=True, target=True)
