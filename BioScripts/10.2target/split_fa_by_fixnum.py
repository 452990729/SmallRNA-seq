#!/usr/bin/env python
#-*- coding:utf8 -*-
# Power by wangyunkai@novogene.com 2018-05-15 10:22:17

import os
import sys
import linecache

def split_fa_by_fixnum(file,num):
   
    num = int(num)
    Quotients,Remainders = divmod(len(linecache.getlines(file))/2,num)
    for linenum in xrange(1,num+1):
        analydir    = os.path.dirname(file)
        filename    = os.path.basename(file)+'_'+str(linenum)+'.fasta'
        if Remainders > 0 :
            if linenum <= (num-2):
                index_start = (Quotients+1)*(linenum-1)*2
                index_end   = (Quotients+1)*linenum*2
            elif linenum == (num-1):
                average     = (len(linecache.getlines(file)) - (Quotients+1)*(num-2)*2)/4
                index_start = (Quotients+1)*(num-2)*2
                index_end   = (Quotients+1)*(num-2)*2 + average*2
            else:
                index_start = (Quotients+1)*(num-2)*2 + average*2
                index_end   = len(linecache.getlines(file))
        else:
            index_start = Quotients*(linenum-1)*2
            index_end   = Quotients*linenum*2
        open(os.path.join(analydir,filename),'w+').writelines(linecache.getlines(file)[index_start:index_end])

split_fa_by_fixnum(sys.argv[1],sys.argv[2])
