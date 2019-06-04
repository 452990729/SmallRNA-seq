#!/ure/bin/env python
# _*_coding:utf-8_*_

import os
import getpass
import argparse
import grp
from datetime import datetime
from mongoengine import *
connect('qcstatcommon',host='172.17.244.172', port=8085,username='usercommon',password='novo2017')

__author__ = "hz"
__concat__ = "huzhigang@novogene.com"

user = getpass.getuser()
group = grp.getgrgid(os.getgid()).gr_name
date = datetime.now()

class qcstat_coll(DynamicDocument):
    LIBID          = StringField(max_length=50)
    RUNID          = StringField()
    USERNAME       = StringField(max_length=50)
    READS_RAW      = FloatField()
    READS_CLEAN    = FloatField()
    DATA_RAW       = FloatField()
    DATA_CLEAN     = FloatField()
    RAW_N_FQ1      = FloatField()
    N_FQ1          = FloatField()
    RAW_LOWQ_FQ1   = FloatField()
    LOWQ_FQ1       = FloatField()
    RAW_Q20_FQ1    = FloatField()
    Q20_FQ1        = FloatField()
    RAW_Q30_FQ1    = FloatField()
    Q30_FQ1        = FloatField()
    RAW_GC_FQ1     = FloatField()
    GC_FQ1         = FloatField()
    RAW_ERROR_FQ1  = FloatField()
    ERROR_FQ1      = FloatField()
    RAW_DUP_FQ1    = FloatField()
    DUP_FQ1        = FloatField()
    DIS_N          = FloatField()
    DIS_LOWQ       = FloatField()
    PROJPATH       = StringField()
    PROJECT_NUM    = StringField()
    SAMPLE_NAME    = StringField()
    GROUPNAME      = StringField()
    QCRPORT_TIME   = DateTimeField()


def store_qc_info_by_sample(projpath,qc_dir, sample_name, lib_id, run_id, project_num=''):
    qc_stat_file = os.path.join(qc_dir, sample_name, 'clean_data', '{}.stat'.format(sample_name))
    assert os.path.isfile(qc_stat_file)
    with open(qc_stat_file) as f:
        qcstats = qcstat_coll.objects(
            LIBID=lib_id, RUNID=run_id, SAMPLENAME=sample_name,
            USERNAME=user, GROUPNAME=group)
        context = f.readlines()
        records = dict()
        records['PROJPATH']                          = projpath
        records['QCREPORT_TIME']                     = date
        records['PROJECT_NUM']                       = project_num
        records['READS_RAW'], records['READS_CLEAN'] = [int(i) for i in context[1].strip().split('\t')[1:]]
        records['DATA_RAW'], records['DATA_CLEAN']   = [int(i.split('(')[0]) for i in context[2].strip().split('\t')[1:]]
        records['RAW_N_FQ1'], records['N_FQ1']         = [
            float(i.strip('%'))/100 for i in context[3].strip().split('\t')[1:]
        ]
        records['RAW_LOWQ_FQ1'], records['LOWQ_FQ1']   = [
            float(i.strip('%')) / 100 for i in context[4].strip().split('\t')[1:]
        ]
        records['RAW_Q20_FQ1'], records['Q20_FQ1']     = [
            float(i.strip('%')) / 100 for i in context[5].strip().split('\t')[1:]
        ]
        records['RAW_Q30_FQ1'], records['Q30_FQ1']     = [
            float(i.strip('%')) / 100 for i in context[6].strip().split('\t')[1:]
        ]
        records['RAW_GC_FQ1'], records['GC_FQ1']       = [
            float(i.strip('%')) / 100 for i in context[7].strip().split('\t')[1:]
        ]
        records['RAW_ERROR_FQ1'], records['ERROR_FQ1'] = [
            float(i.strip('%')) / 100 for i in context[8].strip().split('\t')[1:]
        ]
        records['RAW_DUP_FQ1'], records['DUP_FQ1']     = [
            float(i.strip('%')) / 100 for i in context[9].strip().split('\t')[1:]
        ]
        records['DIS_N'] = float(context[10].strip().split('\t')[1].strip('%')) / 100
        records['DIS_LOWQ'] = float(context[11].strip().split('\t')[1].strip('%')) / 100
        qcstats.update(upsert=True, **records)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='store qc statistic records into mongoen db')
    parser.add_argument('--projpath',help="sRNA projpath Directory")
    parser.add_argument('--qc_dir', help="sRNA QC Directory", default=os.getcwd())
    parser.add_argument('--mapfile', help="abspath of project mapfile", required=True)
    parser.add_argument('--project', help="project information")
    args = parser.parse_args()

    mapfile  = args.mapfile
    assert os.path.isfile(mapfile)
    projpath = args.projpath 
    qc_dir   = args.qc_dir
    project  = args.project

    with open(mapfile) as f:
        for row in f:
            raw_path, lib_id, sample_name = row.strip().split('\t')
            run_id = os.path.split(raw_path)[-1]
            store_qc_info_by_sample(projpath,qc_dir, sample_name, lib_id, run_id, project)
