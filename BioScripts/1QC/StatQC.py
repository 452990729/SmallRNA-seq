#!/usr/bin/env python2


import sys
import re
from glob import glob

def GetInfo(file_in):
    with open(file_in, 'r') as f:
        for line in f:
            list_split = re.split('\t', line.strip())
            if list_split[0] == 'total reads:':
                TOTAL = list_split[1]
            elif list_split[0] == 'clean reads:':
                CLEAN = list_split[1]
                RATIO = list_split[2].lstrip('(').rstrip(' %)')
    return str(int(TOTAL)-int(CLEAN))+'\t'+TOTAL+'\t'+CLEAN+'\t'+str(round(float(RATIO), 2))

def main():
    sample_path = glob(sys.argv[1]+'/*/*_clean_process_overview.txt')
    print 'Sample\tFiltered_Flagment\tTotal_Flagment\tPassed_Flagment\tPassRate(%)'
    for sample in sample_path:
        lb = re.split('\/', sample)[-2]
        info = GetInfo(sample)
        print lb+'\t'+info


if __name__ == '__main__':
    main()
