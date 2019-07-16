import os
import abc
import sys
import logging
import argparse
import configparser
from launchers.job import Job


class Wrapper(object):
    def __init__(self, arglist):
        self.name = 'genericWrapper'
        self._construct_argparser()
        self.args = self.argparser.parse_args(arglist)
        __metaclass__ = abc.ABCMeta


    def _construct_argparser(self):
        self.cd = os.path.join(os.path.dirname(os.path.realpath(__file__)))
        self.argparser = argparse.ArgumentParser(
            prog=self.__class__.__name__, description="placeholder",
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            add_help=False)

        # manage config
        self.argparser.add_argument(
            "-c", "--config",
            default=os.path.join(self.cd, '{}.ini'.format(
                os.path.splitext(sys.modules[self.__class__.__module__].__file__)[0])))
        self.argparser.add_argument('-pi', '--packageInfo', default=None)
        self.argparser.add_argument("-ll", "--logLevel", default="INFO",
                                    choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
        self.argparser.add_argument("-lf", "--logfile", default=None,
                                    help="filename of the log file")

    def _set_config(self):
        # parse config file
        if os.path.isfile(self.args.config):
            config = configparser.ConfigParser()
            config.optionxform = str
            config.read(self.args.config)
        else:
            logging.warning('File not found: %s', self.args.config)
            config = None

        return config

    def sigint_handler(self, signal, frame):
        logging.info('SIGINT detected. Attempting to stop \'%s\'', self.name)
        Job('qdel {}'.format(self.name)).submit_local()
        sys.exit(0)

    @abc.abstractmethod
    def build(self, **kwargs):
        """Method that should build JSON from output files"""
