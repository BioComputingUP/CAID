import os
import logging
import argparse
import configparser


def load_names(labels_file):
    with open(labels_file) as f:
        labels = dict(zip(*line.strip().split()) for line in f)
    return labels

def parse_config(config_file):
    config_parser = configparser.ConfigParser()
    config_parser.optionxform = str
    config_parser.read(config_file)

    return config_parser


def set_logger(logfile, level):
    handlers = list()
    log_formatter = logging.Formatter('%(asctime)s | %(module)-13s | %(levelname)-8s | %(message)s')

    if logfile:
        file_handler = logging.FileHandler(logfile, 'a')
        file_handler.setFormatter(log_formatter)
        handlers.append(file_handler)
    else:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(log_formatter)
        handlers.append(console_handler)

    logging.basicConfig(level=level,
                        format=log_formatter,
                        handlers=handlers)





def parse_args_plots(wd):
    parser = argparse.ArgumentParser(
        prog='caid-assess', description="CAID: Critical Assessment of Intrinsic Disorder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument('mode', choices=['disorder', 'binding', 'test'],
    #                     help='category for which to produce caid-uniformed prediction files')

    parser.add_argument('csv',
                        help='csv file')
    parser.add_argument('-s', '--secondaryCsv', default=None)
    parser.add_argument('-m', '--metricName', default='acc')

    # parser.add_argument('-s', '--strictNegatives', action='store_true', default=False,
    #                     help="path to an alternative configuration file.")
    parser.add_argument('-o', '--outputDir', help='directory where the output will be written',
                        default=None)
    parser.add_argument('-d', '--codes', action='store_true', default=False)
    parser.add_argument('-a', '--hideAfter', default=0, type=int)
    # config file
    parser.add_argument('-c', '--conf', type=str,
                        default=os.path.join(wd, 'config.ini'),
                        help="path to an alternative configuration file.")
    # log options
    parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    parser.add_argument("-ll", "--logLevel", default="ERROR",
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help='log level filter. All levels <= choice will be displayed')

    args = parser.parse_args()
    return args
