import os
import sys
import copy
import logging

# relative imports
from launchers.job import Job
from launchers.wrapper import Wrapper


class InterProScanLauncher(Wrapper):
    def __init__(self, arglist):
        super(InterProScanLauncher, self).__init__(arglist)
        self._construct_argparser()
        self.args = self._build_args(arglist)
        self.config = self._set_config()
        self.tools = self._get_tools()

        self.basename, _ = os.path.splitext(
            os.path.basename(os.path.normpath(self.args.inputFasta)))
        self.outfile = os.path.join(self.args.outDir, '{}.tsv'.format(self.basename))

    # Private methods
    # - Init methods
    def _construct_argparser(self):

        super()._construct_argparser()
        dirs = self.argparser.add_argument_group('dirs')
        # manage input
        self.argparser.add_argument("inputFasta", default=None,
                                    help="Input fasta file (multi fasta expected)")
        # manage output
        dirs.add_argument('-o', '--outDir', default='.')
        # manage tools to activate
        self.argparser.add_argument("-t", "--tools", default=None,
                                    help="Tools to run")
        self.argparser.add_argument('-e', '--environModule', default=None)

    def _build_args(self, arglist=None):
        args = self.argparser.parse_args(arglist if arglist else sys.argv[1:])
        # make paths absolute
        args.inputFasta = os.path.abspath(args.inputFasta) if args.inputFasta else None
        args.outDir = os.path.abspath(args.outDir)
        # make dirs if missing
        os.makedirs(args.outDir, exist_ok=True)

        return args

    def _get_tools(self):
        if self.args.tools:
            tools = self.args.tools
        else:
            if self.config:
                tools = list()
                for tool, switch in dict(self.config.items('tools')).items():
                    if switch == '1':
                        tools.append(tool)
                        logging.debug('%s included', tool)
                    else:
                        logging.debug('%s excluded', tool)
            else:
                logging.warning('Config file not found. Tools to launch not set.')
                tools = None

        return tools

    # - Static methods
    @staticmethod
    def _parse_progress(raw_string):
        for iprs_msg in raw_string:
            iprs_msg = iprs_msg.decode('utf-8')
            if '%' in iprs_msg:
                yield iprs_msg.strip('\n').rstrip().lstrip().split()[2]

    @staticmethod
    def _parse_tool_list(raw_list):
        tools = list()
        for tool in raw_list.split('\n'):
            if tool and not tool.startswith('Deactivated'):
                ts = tool.split(':')
                if len(ts) == 2:
                    tool_desc = map(lambda s: s.lstrip().rstrip().strip(')').replace(' (', '-'), ts)

                else:
                    tool_desc = map(lambda s: s.lstrip().rstrip(), ts)

                tools.append(tool_desc)

        return tools

    # Public methods
    def check_tools_available(self, bindir='', envmod=None):
        job = Job(self.make_command(bindir=bindir, just_check=True))
        proc = job.submit_local(envmod)
        output, err = proc.communicate()
        available, unavailable = None, None

        for i, chunk in enumerate(output.decode('utf-8').split('analyses:\n')):

            if i == 1:
                available = list()
                for tname, tdesc in self._parse_tool_list(chunk):
                    logging.debug('%-*s -- Available' % (28, tname))
                    available.append(tname)

            if i == 2:
                unavailable = list()
                for tname, tdesc, tbin in self._parse_tool_list(chunk):
                    logging.debug('%-*s -- Unavailable' % (28, tname))
                    unavailable.append(tname)

        return available, unavailable

    def make_command(self, bindir='', just_check=False, avail=None):
        bindir = os.path.abspath(bindir) if bindir else bindir
        cline = os.path.join(bindir, 'interproscan.sh')

        tools_to_run = copy.copy(self.tools)

        if not just_check:
            if tools_to_run and avail:
                apps = list()
                for app in avail:
                    if app.split('-')[0] in tools_to_run:
                        apps.append(app)
                        tools_to_run.remove(app.split('-')[0])
                for tool in tools_to_run:
                    logging.warning('A selected tool is not available: %s', tool)
                if apps:
                    appl = '-appl {}'.format(','.join(apps))
                else:
                    appl = ''
                    logging.warning('No selected tool was available. '
                                    'Running with default setting')
            else:
                logging.warning('No tool specified. Running with default setting')
                appl = ''

            cline = '{} -i {} -f TSV -goterm -dp {}'.format(
                cline, self.args.inputFasta, appl)

        logging.debug('cmd -- %s', cline)
        return cline

    def run(self, available=None, bindir=None):
        if not available:
            logging.warning('Available tools not specified. Attempting to retrieve')
            available, _ = self.check_tools_available(envmod=self.args.environModule, bindir=bindir)

        cmd = self.make_command(avail=available, bindir=bindir)
        job = Job(cmd)
        proc = job.submit_local(envmod=self.args.environModule)

        try:
            for progress in self._parse_progress(iter(proc.stdout.readline, b'')):
                logging.info('progress %s', progress)
        except Exception as e:
            proc.kill()
            raise Exception(e)

        out, err = proc.communicate()
        if err:
            logging.error('%s', ''.join(err.decode('utf-8').split('\n')[1:-1]))

        return out, err
