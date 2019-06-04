#/usr/bin/env python
import re
import sys
import os.path

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    T1='\033[96m'
    F1='\033[91m'
    END1='\033[90m'
    UNDERLINE = '\033[4m'
    FLE  = '\033[1;31;40m'
    TUE  = '\033[1;32;40m'
pattern=re.compile(r'/\S*')
with open(sys.argv[1],"r") as F1:
    for eachline in F1.readlines():
        if re.findall(pattern,eachline):
            path=pattern.findall(eachline)[0].split()[0].split("'")[0]
            if os.path.exists(path):
                print "%s %s%s%s" %(path,bcolors.TUE,os.path.exists(path),bcolors.ENDC)
            else:
                print "%s %s%s%s" %(path,bcolors.FLE,os.path.exists(path),bcolors.ENDC)
F1.close()
