#!/ure/bin/env python
# _*_coding:utf-8_*_

# import block
import os
import argparse
from string import Template
from collections import defaultdict
import re


WS = r'(?P<ws>\s+)'   # write space
PUNC = r'(?P<pnuc>[,\-])'      # punctuation
DROP_WORD = r'(?P<drop_word>\w+)'
NAME = r'name\s+(?P<job_name>\w+)'
SLOTS = r'p=(?P<slots_1>\d+)|slots\s+(?P<slots_2>\d+)'
MEMORY = r'vf=(?P<memory_1>\d+[GgMm])|memory\s+(?P<memory_2>\d+)'
QUEUE = r'queue\s+(?P<queue_1>\S+)|-q\s+(?P<queue_2>\S+)|-P\s+(?P<queue_3>\S+)'
CMD = r'(?P<cmd>sh\s+\S+\.sh)'


# print("|".join([NAME, SLOTS, MEMORY, QUEUE, CMD, DROP_WORD, WS, NL, COMMA]))
process_patten = re.compile("|".join([NAME, SLOTS, MEMORY, QUEUE, CMD, DROP_WORD, WS, PUNC]))
order_patten = re.compile(
    r'order\s+(?P<prev_1>\w+)\s+before\s+(?P<after_1>\w+)|'
    r'order\s+(?P<after_2>\w+)\s+after\s+(?P<prev_2>\w+)'
)


# nextflow lfs resource manager job file template
NEXTFLOW_LFS_PROCESS_TEMPLATE = Template("""
process $job_name {
    $input_str
    $output
    $executor
    $cpus
    $memory
    $queues
    $cmd
}
""")



#nextflow sge resource manager job file template
NEXTFLOW_SGE_PROCESS_TEMPLATE = Template("""
process $job_name {
    $input_str
    $output
    $cmd
}
""")


def read_sjm_job(sjm_job_file):
    with open(sjm_job_file) as f:
        process_block = ''
        for row in f:
            if row.startswith('order'):
                flag = "order"
                yield row, flag
            elif row.startswith('job_begin'):
                flag = 'process'
                while True:
                    row = f.next()
                    if row.startswith('job_end'):
                        yield process_block, flag
                        process_block = ''
                        break
                    else:
                        process_block += row


def sjm2nextflow(sjm_job_file, nextflow_job_file, executor):
    order_tree = defaultdict(list)
    priority_dict = defaultdict(int)
    process_contexts = dict()
    for ctx, flag in read_sjm_job(sjm_job_file):
        if flag == 'order':
            scanner = order_patten.scanner(ctx)
            for m in iter(scanner.match, None):
                groupdict = m.groupdict()
                prev = groupdict.get('prev_1', None) or groupdict.get('prev_2', None)
                after = groupdict.get('after_1', None) or groupdict.get('after_2', None)
                order_tree[after].append(prev)
                priority_dict[after] = max(priority_dict[after], priority_dict[prev]+1)
        else:
            context = dict()
            scanner = process_patten.scanner(ctx)
            context['executor'] = "" if "localhost" in ctx else "executor\t'{}'".format(executor)
            for m in iter(scanner.match, None):
                if m.lastgroup == 'job_name':
                    context['output'] = "output:\n\tstdout {}".format(m.groupdict().get(m.lastgroup))
                    context['job_name'] = m.groupdict().get(m.lastgroup)
                elif m.lastgroup == 'cmd':
                    context['cmd'] = '"""\n\t{}\n    """'.format(m.groupdict().get(m.lastgroup))
                '''
                elif m.lastgroup.startswith('memory'):
                    if m.group().endswith(("m", 'M')):
                        context['memory'] = 'memory\t"1 GB"'
                    elif m.group().endswith(('g', 'G')):
                        context['memory'] = 'memory\t"{} GB"'.format(m.groupdict().get(m.lastgroup)[:-1])
                    else:
                        context['memory'] = 'memory\t"{} GB"'.format(m.groupdict().get(m.lastgroup))
                elif m.lastgroup.startswith('slots'):
                    context['cpus'] = "cpus\t{}".format(m.groupdict().get(m.lastgroup))
                elif m.lastgroup.startswith('queue'):
                    context.setdefault('queues', '')
                    context['queues'] += 'queue\t"{}"\n\t'.format(m.groupdict().get(m.lastgroup))
                '''
            '''
            if not context.get('cpus', None):
                context['cpus'] = 'cpus\t1'
            if not context.get('memory', None):
                context['memory'] = 'memory\t"1 GB"'
            if not context.get('queues', None):
                context['queues'] = ''
            #context['clusterOptions'] = 'clusterOptions \'-cwd -V -S /bin/bash -q {queues} -l p={cpus},vf={memory}\''.format(
            #        queues = context['queues'].strip('queue\t\n').strip('"'),cpus = context['cpus'].strip().split()[1],
            #        memory = context['memory'].split()[1].strip('"')+'G'
            #)
            '''
            if executor != 'sge':
                process_contexts[context['job_name']] = NEXTFLOW_LFS_PROCESS_TEMPLATE.safe_substitute(
                    **context
                )
            else:
                process_contexts[context['job_name']] = NEXTFLOW_SGE_PROCESS_TEMPLATE.safe_substitute(
                    **context
                )                

    for k in priority_dict:
        if order_tree.get(k):
            priority_dict[k] = max(priority_dict[pk] + 1 for pk in order_tree[k])
    with open(nextflow_job_file, 'w') as out:
        for job_name in sorted(priority_dict, key=lambda k: priority_dict[k]):
            process_ctx = process_contexts[job_name]
            #print process_ctx
            if order_tree.get(job_name, None):
                input_str = "input:"
                for prev in order_tree.get(job_name):
                    input_str += "\n\tstdin {}".format(prev)
            else:
                input_str = ''
            #print job_name
            #print input_str
            #print Template(process_ctx).substitute(input_str=input_str)
            out.write(
                Template(process_ctx).substitute(input_str=input_str)
            )
            #print process_ctx 
       # for job_name, process_ctx in process_contexts.items():
        #     # print(job_name)
        #     # print(order_tree.get(job_name))
        #     if order_tree.get(job_name, None):
        #         input_str = "input:"
        #         for prev in order_tree.get(job_name):
        #             input_str += "\n\tstdin {}".format(prev)
        #     else:
        #         input_str = ''
        #     out.write(
        #         Template(process_ctx).substitute(input_str=input_str)
        #     )


def commandline():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--sjm_job', help="sjm job file to be converted")
    parser.add_argument(
        '-o', '--nextflow_job', help="nextflow job file to be converted to, default is "
                                     "the basename of sjm job file wiht 'nf' appendix "
    )
    parser.add_argument('-e', '--executor', choices=['lsf', 'sge'], default='lsf')
    return parser.parse_args()


def main():
    args = commandline()
    sjm_job_file = args.sjm_job
    nextflow_job_file = args.nextflow_job if args.nextflow_job else '{}.nr'.format(
        os.path.splitext(sjm_job_file)[0]
    )
    executor = args.executor
    sjm2nextflow(sjm_job_file, nextflow_job_file, executor)

if __name__ == '__main__':
    main()
