import os
import sys
import logging

# relative imports
from launchers.job import Job
from launchers.wrapper import Wrapper


class FessLauncher(Wrapper):
    def __init__(self, arglist):
        super(FessLauncher, self).__init__(arglist)
        self._construct_argparser()
        self.args = self._build_args(arglist)
        self.config = self._set_config()

        self.basename, _ = os.path.splitext(
            os.path.basename(os.path.normpath(self.args.fasta)))
        self.outfile = os.path.join(self.args.outDir, '{}.txt'.format(self.basename))

    # Private methods
    # - Init methods
    def _construct_argparser(self):

        super()._construct_argparser()
        dirs = self.argparser.add_argument_group('dirs')
        # manage input
        self.argparser.add_argument("fasta", default=None,
                                    help="Input fasta file (multi fasta expected)")
        # manage output
        dirs.add_argument('-o', '--outDir', default='.')
        # manage options
        self.argparser.add_argument("-b", "--bin", default='.')
        # self.argparser.add_argument('-e', '--environModule', default='fess')

    def _build_args(self, arglist=None):
        args = self.argparser.parse_args(arglist if arglist else sys.argv[1:])
        # make paths absolute
        args.fasta = os.path.abspath(args.fasta) if args.fasta else None
        args.outDir = os.path.abspath(args.outDir)
        # make dirs if missing
        os.makedirs(args.outDir, exist_ok=True)

        return args

    def make_command(self, bindir='', infile=None):
        infile = infile if infile else self.args.fasta
        bindir = os.path.abspath(bindir) if bindir else bindir
        if infile:
            cline = '{} {} {}'.format(os.path.join(bindir, 'fess'), infile, self.outfile)
        else:
            logging.error('Input file missing. Required input file')

        logging.debug('cmd -- %s', cline)
        return cline

    def run(self, infile=None):
        cmd = self.make_command(self.args.bin, infile=infile)
        job = Job(cmd)
        proc = job.submit_local(cwd=self.args.bin)

        out , err = proc.communicate()

        if err:
            logging.error('%s', ''.join(err.decode('utf-8').split('\n')[1:-1]))

        return out, err