#!/usr/bin/env python
# _*_ coding:utf-8 _*_

import commands
import getpass
import glob
import os
import stat
import re
import warnings

from Settings import BioConfig


class AtomJob(object):
    """
    用于根据参数构建job文件中shell相关区块的内容
    """

    def __init__(self, job_name, shell_path, p = None,fun=None, memory=None, status='waiting', sched=None,
                 is_run_local=False):
        self.job_name          = job_name
        self._memory           = memory
        self.status            = status
        self.function_name     = fun
        self.p                 = p
        self.default_memory_setting = BioConfig.MEMORY
        self.sched = sched + ' -l p={}'.format(self.p) if sched else self._get_default_sched()
        self.shell_path = shell_path
        self.is_run_local = is_run_local

    _template = '\njob_begin\n' \
                '  name {job_name}\n' \
                '  memory {memory}\n' \
                '  status {status}\n' \
                '  sched_options {sched}\n' \
                '  cmd sh {shell_path}\n' \
                'job_end\n'

    _local_template = '\njob_begin\n' \
                      '  name {job_name}\n' \
                      '  host localhost\n' \
                      '  cmd sh {shell_path}\n' \
                      'job_end\n'

    def _get_default_sched(self):
        """
        construct default sched options useing system qselect command
        """
        user = getpass.getuser()
        try:
            temp = commands.getoutput("qselect -U {user}".format(user=user))
            temp = 'test.q'
            queues = set(i.split('@')[0] for i in temp.strip().split(os.linesep) if 'tjsmp' not in i)
            if 'joyce.q' in queues:
                queues.remove('joyce.q')
            default_queues = " -q " + ' -q '.join(
                queues
            )
            default_sched = "-V -cwd " + ' -l p={}'.format(self.p) + default_queues
            return default_sched
        except IOError as e:
            exit(e)

    def set_sched(self, value, action="r"):
        """
        modify sched for default "-S /bin/bash -cwd -V -q disease1.q -q novo.q -q all.q"
        must called before job writing
        :param value: your sjm sched setting, string type
        :param action: must be "r" or "+" or "-", if r, default will be replaced, if +,
        value will add to default, if -, remove value from default
        """
        if action not in ['r', '+', '-']:
            print "invalid action, ached will not be change!!!"
        if action == 'r':
            self.sched = value
        if action == '+':
            self.sched += " "
            self.sched += value
        if action == "-":
            self.sched = self.sched.replace(value, '')

    def set_status(self, value):
        if value.strip().lower() in ['done', 'waiting', ]:
            self.status = value.strip().lower()
        else:
            print("Invalid status value, so the status will not be changed!!! status should be done or waiting")

    def _get_default_memory(self):
        return self.default_memory_setting.get(self.function_name,'2G')

    @property
    def memory(self):
        return self._memory or self._get_default_memory()

    @memory.setter
    def memory(self, value):

        if re.match(r'^\d+[gGmM]$', value):
            self._memory = value
        else:
            raise ValueError('Invalid memory settings')

    @property
    def context(self):
        if self.is_run_local:
            return self._local_template.format(
                **self.__dict__
            )
        else:
            return self._template.format(
                memory=self.memory, **self.__dict__
            )

    def __str__(self):
        return self.job_name


class Job(object):
    """
    job 包含三个部分内容：
        1. job 运行脚本相关区块, 称之为atom_job,
        2. job log 目录
        3. job orders区块
    """

    def __init__(self, log_dir):
        self.log_dir = log_dir
        self.atom_jobs = list()
        self.orders = list()

        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)

    @property
    def atom_job_dict(self):
        return {atom_job.job_name: atom_job for atom_job in self.atom_jobs}

    @property
    def order_job_names(self):
        return set(order.strip().split()[1] for order in self.orders).update(
            set(order.strip().split()[3] for order in self.orders)
        )

    def check_orders(self):
        pass

    def add_atom_job(self, atom_job):
        if isinstance(atom_job, AtomJob):
            self.atom_jobs.append(atom_job)
        else:
            exit('Only can add AtomJob instance!!!!')

    def add_order(self, job_name, previous_jobs=None, after_jobs=None):
        if isinstance(previous_jobs, str):
            previous_jobs = [previous_jobs]
        if isinstance(after_jobs, str):
            after_jobs = [after_jobs]
        if previous_jobs:
            order = '\n'.join(
                'order {} after {}'.format(job_name, previous_job) for previous_job in previous_jobs
            )
            self.orders.append(order)
        if after_jobs:
            order = '\n'.join(
                'order {} before {}'.format(job_name, after_job) for after_job in after_jobs
            )
            self.orders.append(order)

    def write(self, job_file):
        with open(job_file, 'w') as f:
            for atom_job in self.atom_jobs:
                f.write(atom_job.context)

            for order in self.orders:
                f.write(order + '\n')
            
            f.write("\nlog_dir {log_dir}\n".format(log_dir=self.log_dir))
                

def make_dir(analydir):
    if not os.path.exists(analydir):
        os.makedirs(analydir)


def write_shell(code,shell_path):
    shell_directory = os.path.dirname(shell_path)
    if not os.path.exists(shell_directory):
        os.makedirs(shell_directory)
    with open(shell_path,'w+') as ouf:
        ouf.writelines(code)
    
              
if __name__ == '__main__':
    pass
