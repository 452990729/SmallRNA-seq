#!/usr/bin/env python
# _*_ coding:utf-8 _*_
#用于删除和压缩已结题并已释放数据项目中冗余数据,该脚本适用于已执行byebye后的项目，是对byebye脚本的补充
#目前sRNA项目路径为/TJPROJ1/RNA/ncRNA/sRNA201707；/TJPROJ1/RNA/ncRNA/sRNA；/TJPROJ1/RNA/ncRNA201704
import sys
import glob
import re
import os 
import os.path
import argparse
parser = argparse.ArgumentParser(description="delete redundant data of project directory after byebye ")
parser.add_argument('--outdir',help='filename',required=True)
argv = vars(parser.parse_args())
pro_workdir=argv['outdir'].strip()
os.chdir(pro_workdir)
log_del=open("log_del.txt","w")
sign1_byebye=glob.glob("shell/byebye.sh.e*")
#sign1_release_data=glob.glob("shell/release_data.sh.e*")
sign2_release_data=glob.glob("*data_give")
if sign1_byebye and sign2_release_data:
    log_del.writelines("%s 该项目已经释放数据，已执行byebye,执行过程文件压缩和删除工作\n" % pro_workdir)
    #known
    #p2=glob.glob("*.known/*/*firstbase.png")
    #p3=glob.glob("*.known/*/*position.png")
    os.system("rm -rf *.known/*/*firstbase.png")
    os.system("rm -rf *.known/*/*position.png")
    os.system("rm -rf *.known/*/*.known/image")
    if glob.glob("*.known/all.map.collapse.fa"):
        p1=glob.glob("*.known/all.map.collapse.fa")[0]
        os.system("gzip %s" % p1)
    if glob.glob("*.known/*/pdfs_*.known"):
        p4=glob.glob("*.known/*/pdfs_*.known")[0]
        os.system("tar -zvcf pdfs_known.tar.gz %s --remove-file " % p4)
    if glob.glob("*.known/*/*.known/known_miRNA.unmap.fas"):
       p5=glob.glob("*.known/*/*.known/known_miRNA.unmap.fas")[0]
       os.system("gzip %s" % p5)
    #if glob.glob("*.known/*/*.known/image"):
    #   p6=glob.glob("*.known/*/*.known/image")[0]
    #   os.system("tar -zvcf known.image.tar.gz %s --remove-file " % p6)
    #ncRNA
    if glob.glob("*.ncRNA/*/output/ncRNA.unmap.fas"):
        p7=glob.glob("*.ncRNA/*/output/ncRNA.unmap.fas")[0]
        os.system("gzip %s" % p7)
    #repeat
    os.system("rm -rf *.repeat/*bar.png")
    os.system("rm -rf *.repeat/*bar.pdf")
    os.system("rm -rf *.repeat/repeat.map.fas")
    os.system("rm -rf *.repeat/repeat.unmap.fas")
    if glob.glob("*.repeat/mdm/repeat.map.fas"):
        p12=glob.glob("*.repeat/mdm/repeat.map.fas")[0]
        os.system("gzip %s" % p12)
    if glob.glob("*.repeat/mdm/repeat.map.fas"):
        p13=glob.glob("*.repeat/mdm/repeat.map.fas")[0]
        os.system("gzip %s" % p13)
    #NAT
    if glob.glob("*.NAT/*/output/*.unmap.fas"):
        p14=glob.glob("*.NAT/*/output/*.unmap.fas")[0]
        os.system("gzip %s" % p14)
    #gene
    if glob.glob("*.gene/*/output/gene.unmap.fas"):
        p15=glob.glob("*.gene/*/output/gene.unmap.fas")[0]
        os.system("gzip %s" % p15)
    #novel
    #p16=glob.glob("*.novel/*/*firstbase.png")
    #p17=glob.glob("*.novel/*/*position.png")
    #p18=glob.glob("*.novel/*/pdfs_*.novel")[0]
    #p19=glob.glob("*.novel/*/*.novel/novel_miRNA.unmap.fas")[0]
    #p20=glob.glob("*.novel/*/*.novel/image")[0]
    os.system("rm -rf *.known/*/*firstbase.png")
    os.system("rm -rf *.known/*/*position.png")
    os.system("rm -rf *.known/*/*.novel/image")
    if glob.glob("*.novel/*/*.novel/novel_miRNA.unmap.fas"):
        p19=glob.glob("*.novel/*/*.novel/novel_miRNA.unmap.fas")[0]
        os.system("gzip %s" % p19) 
    if glob.glob("*.novel/*/pdfs_*.novel"):
        p18=glob.glob("*.novel/*/pdfs_*.novel")[0]
        os.system("tar -zvcf novel.image.tar.gz %s --remove-file " % p18)
    if glob.glob("*.novel/*/*.novel/image"):
        p20=glob.glob("*.novel/*/*.novel/image")[0]
        os.system("tar -zvcf image.tar.gz %s --remove-file " % p20) 
    #TAS
    if glob.glob("*.TAS/*/output/TAS.unmap.fas"):
        p21=glob.glob("*.TAS/*/output/TAS.unmap.fas")[0]
        os.system("gzip %s" % p21)
    #Category 
    #p22=glob.glob("*Category/*pie.png")
    #p23=glob.glob(".Category/*pie.pdf")
    os.system("rm -rf *Category/*pie.png")
    os.system("rm -rf .Category/*pie.pdf")
    #edit_family
    #target
    #diff
    #enrich
else:
#    print " 请check是否执行释放数据，执行byebye"
    log_del.writelines("%s 该项目未执行释放数据，未执行byebye，请核对\n" % pro_workdir) 
log_del.close()
