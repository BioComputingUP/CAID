import os
import stat
import time
import shlex
import logging
import subprocess
# from modulecmd import Modulecmd


class Queue(object):
    def __init__(self, _jobname):
        self.jobname = _jobname
        self.jobs_count = None

    def get_jobs_count(self):
        time.sleep(3)
        qstat = Job('qstat -f').submit_local()
        grep_cmd = 'grep -c "{}"'.format(self.jobname)
        grep = subprocess.Popen(shlex.split(grep_cmd), stdin=qstat.stdout, stdout=subprocess.PIPE)
        qstat.stdout.close()
        jobs = grep.communicate()[0]

        try:
            jobs = int(jobs.decode('utf-8').strip('\n'))
        except ValueError as exp:
            logging.error("SGE  Queue %s", exp)
            jobs = None

        self.jobs_count = jobs

    def wait_completion(self, check_each_seconds=60):
        continue_checking = True

        while continue_checking is True:
            self.get_jobs_count()

            if self.jobs_count and self.jobs_count != 0:
                logging.info('SGE   Jid %s', self.jobname)
                logging.info('SGE   Queued Jobs %s', self.jobs_count)
                time.sleep(check_each_seconds)

            else:
                continue_checking = False


class Job(object):
    def __init__(self, _cmd):
        self.cmd = _cmd

    def submit_sge(self, script, mdl='', stde='', stdo='', reserve_cores=None, queues=None):

        qs_stat = os.stat(script)
        os.chmod(script, qs_stat.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        stde = '-e {}'.format(stde) if stde else ''
        stdo = '-o {}'.format(stdo) if stdo else ''
        pe = '-pe multithread {}'.format(reserve_cores) if reserve_cores else ''
        queue = '-q {}'.format(','.join(queues)) if queues else ''

        cmd = "qsub {} {} {} {} {} {} '{}'".format(pe, queue, stde, stdo, script, mdl, self.cmd)
        logging.debug(cmd)

        proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        logging.debug('Popen Pid %s', proc.pid)

        return proc

    def submit_sge_bin(self, stde, stdo, reserve_cores=None, queues=None):

        stde = '-e {}'.format(stde) if stde else ''
        stdo = '-o {}'.format(stdo) if stdo else ''
        pe = '-pe multithread {}'.format(reserve_cores) if reserve_cores else ''
        queue = '-q {}'.format(','.join(queues)) if queues else ''

        cmd = "qsub {} {} {} {} -b y {}".format(pe, queue, stde, stdo, self.cmd)
        logging.debug(cmd)

        proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        logging.debug('Popen Pid %s', proc.pid)

        return proc

    def submit_local(self, envmod=None, bufsize=-1, cwd=None):
        # m = Modulecmd()
        # if envmod:
        #     m.load(envmod)

        proc = subprocess.Popen(shlex.split(self.cmd),
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=bufsize,
                                cwd=cwd)
        # if envmod:
        #     m.unload(envmod)
        # del m

        return proc
