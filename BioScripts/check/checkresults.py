#coding=utf-8
#sRNA check results  scripts
import sys
import os
import os.path
import re
import argparse
outdir=os.getcwd()
parser = argparse.ArgumentParser(description="sRNA check results  scripts")		
parser.add_argument('--outdir',help="the outputdir of analysis result",default=None)
parser.add_argument('--sample',help="sample names(sample1,sample2,...)",required=True)
parser.add_argument('--project',help='project name, maybe same with the name of root dir, which will be displayed in the final report title',required=True)
parser.add_argument('--venn_cluster',help="venn ploting way; suit for 2~4 compare groups; '_' split compare group in the same one plot;  ',' split different plots, e.g. '2:1_1:3_2:3,1:3_2:3'",required=False,default=None)
parser.add_argument('--compare',help="group compare way, e.g. '2:1,1:3,2:3'  1,2,3 were groupname order, 1:2(treat:control), intra-group splited by ':',  inter-group by ','",required=False,default=None)
parser.add_argument('--groupname',help='group name, splited by "," , e.g. TR1,TR2,TS1,TS2,TR,TS',required=False,default=None)
parser.add_argument('--group',help="sample grouping way, e.g. TR1,TR2,TS1,TS2,TR1:TR2,TS1:TS2, ':' used in intra-group;  ','used in inter-group",default=None)
parser.add_argument('--org',help='organism type',required=True)
parser.add_argument('--abbr', help='abbreviation of species from miRBase', required=True)
argv = vars(parser.parse_args())
sample=argv['sample'].strip()
samples=sample.split(',')

project=argv['project'].strip()
outdir=argv['outdir'].strip()
out_dir=outdir+'/'+project+'_sRNA_result'+'/'+project+'_results'
abbr=argv['abbr']
abbrs=[]
if abbr == "all" or abbr == "plant" or abbr == "animal":
    abbrs.append('hsa')
else:
    for each in abbr.split(','):
        abbrs.append(each)
#====================group ================================
if argv['group'] == None:
    groups=samples
    groups_iter=samples
    group=sample
    flag_repeat=False
else:
    groups_iter=[]
    if ':' in argv['group'].strip(':'):
        flag_repeat=True
        groups=[each.strip().split(':') for each in argv['group'].strip().strip(':').split(',') if each.strip() != '']
        for each in groups:
            groups_iter+=each
            group_iter_n=[]
        for each in groups_iter:
            if each not in group_iter_n:   # uniq
                group_iter_n.append(each)
                groups_iter=group_iter_n   # uniq
                group=','.join([':'.join(each) for each in groups])
    else:   # group name  not equal sample name, but is not repeat!
        flag_repeat=False
        groups=[each.strip() for each in argv['group'].strip().split(',') if each.strip() != '']
        for each in groups:
            groups_iter.append(each)
        group=','.join(groups)
	#assert len(groups_iter) == len(set(groups_iter))
    assert set(groups_iter).issubset(samples)
group_iter=','.join(groups_iter)

if argv['groupname'] == None:

    if flag_repeat == False:
        groupnames=groups
    else:
        groupnames=['group'+str(k+1) for k in range(len(groups))]		
else:
    groupnames=[each.strip() for each in argv['groupname'].split(',') if each.strip() != '']
    assert len(groupnames) == len(groups)
groupname=','.join(groupnames)


#==================add com.txt  compares ======================
if argv['compare'] != None:
    groupname_dict = {}
    compares_temp=[]
    if os.path.isfile(argv['compare']):
        for g_index,g in enumerate(groupnames):
	    groupname_dict[g] = g_index + 1
            with open(argv['compare']) as f:
	        for line in f:
		    temp = line.split()
		    compare_temp = '%d:%d' %( groupname_dict[ temp[1] ], groupname_dict[ temp[2] ] )
		    compares_temp.append(compare_temp)
	    icompare = ','.join(compares_temp)
    else:
	icompare=argv['compare']
else:
    icompare='1:1'
#print 'icompare is %s' %(icompare)
compares=[each.strip().split(':') for each in icompare.strip().split(',') if each.strip() != '']

M=[]
for each1 in compares:
    assert len(each1) == 2
    for each2 in each1:
        assert each2.isdigit()
        M.append(int(each2))
assert max(M) <= len(groupnames)
assert min(M) > 0	
compare=','.join([':'.join(each) for each in compares])
temp2=[]
for each1 in compares:
    temp1=[]
    for each2 in each1:
        temp1.append(groupnames[int(each2)-1])
    temp2.append(':'.join(temp1))
compare_name=','.join(temp2)    # compare_name
compare_name2=compare_name.replace(':','vs').split(",")
print compare_name2
#==========================================================================
#====================add venn.txt =========================================
if argv['venn_cluster'] != None:
    if os.path.isfile(argv['venn_cluster']):
        compare_temp=compare.split(',')
        ven_clust = []
        with open(argv['venn_cluster']) as f:     
            for line in f:
                element = []
                temp = re.findall(r'\d+',line)[1:]
                transfm = [int(i)-1 for i in temp]
                for i in transfm:
                    element.append(compare_temp[i])
                    element = '_'.join(element)
                    ven_clust.append(element)
                    venn_cluster = ','.join(ven_clust)
    else:
        venn_cluster=argv['venn_cluster'].strip()
        com_pairs=compare.split(',')
        venn_clusters=[each.split('_') for each in venn_cluster.split(',')]
        temp1=[]
        for each1 in venn_clusters:
            temp2=[]
            for each2 in each1:
                assert each2 in com_pairs
                temp3=each2.split(':')
                assert len(temp3) == 2
                temp2.append(groupnames[int(temp3[0])-1]+':'+groupnames[int(temp3[1])-1])
                temp1.append('_'.join(temp2))	
                venn_cluster_name=','.join(temp1)	#venn_cluster_name	
else:
    venn_cluster = None
    #venn_cluster=compare.replace(',','_')
    #venn_cluster_name=compare_name.replace(',','_')  # default plot all comparegroup venn
    venn_cluster_name=None

venn_cluster_vs_names=[]
if venn_cluster_name != None:
    for each in venn_cluster_name.split(','):
        if each.count(':') >= 2 and each.count(':') <= 5:
            flag_venn=True
            venn_cluster_vs_names.append(each.replace(':','vs'))
#print venn_cluster_vs_names
#=============================================================================================

def checkfile(file):
    if os.path.exists(file):
        if os.path.getsize(file):
            message=" is OK"
            return message
        else:
            message=" is empty,failed"
            return message
    else:
        message=" not exist,failed"
        return message
checklog=open(outdir+"/checkresultlog","w")  
###########check 2.QC result ####################
path_QC=os.path.join(out_dir,"2.QualityControl")
checklog.writelines("RawData_Stat.xls "+checkfile(path_QC+"/2.2.RawData_Stat/RawData_Stat.xls")+"\n")
assert len(open(path_QC+"/2.2.RawData_Stat/RawData_Stat.xls","r").readlines())-1 == len(samples)
checklog.writelines("clean_process_overview.xls"+checkfile(path_QC+"/2.3.ReadsClassification/clean_process_overview.xls")+"\n")
assert len(open(path_QC+"/2.3.ReadsClassification/clean_process_overview.xls","r").readlines())-1 == len(samples)
for eachsample in samples:
    checklog.writelines(eachsample+"_error_rate_distribution.png"+checkfile(path_QC+"/2.1.RawData_ErrorRate/"+eachsample+"_error_rate_distribution.png")+"\n")
    checklog.writelines(eachsample+"_seq_len_distribution.pdf"+checkfile(path_QC+"/2.4.Length_Filter/"+eachsample+"_seq_len_distribution.pdf")+"\n")     
###########check 3.map result ####################
for eachsample in samples:
    eachpath_map=os.path.join(out_dir,"3.Mapping_Stat")  
    checklog.writelines("3.Mapping_Stat/"+eachsample+".mapping.stat"+checkfile(eachpath_map+"/"+eachsample+".mapping.stat")+"\n")
    checklog.writelines("3.Mapping_Stat/"+eachsample+".circos.svg"+checkfile(eachpath_map+"/"+eachsample+".circos.svg")+"\n")
    checklog.writelines("3.Mapping_Stat/"+eachsample+".circos.png"+checkfile(eachpath_map+"/"+eachsample+".circos.png")+"\n")
    checklog.writelines("3.Mapping_Stat/"+"reference.mapping.stat"+checkfile(eachpath_map+"/"+"reference.mapping.stat")+"\n")    
    assert len(open(eachpath_map+"/reference.mapping.stat","r").readlines())-1 == len(samples)
###########check 4.Known_miRNA result ####################
known_path=os.path.join(out_dir,"4.Known_miRNA")
#checklog.writelines("4.Known_miRNA/"+"all.map.collapse.fa"+checkfile(known_path+"/all.map.collapse.fa")+"\n")
checklog.writelines("4.Known_miRNA/"+"hairpin.fa"+checkfile(known_path+"/hairpin.fa")+"\n")
checklog.writelines("4.Known_miRNA/"+"hairpin_mature.fa"+checkfile(known_path+"/hairpin_mature.fa")+"\n")
checklog.writelines("4.Known_miRNA/"+"hairpin_mature.pairs"+checkfile(known_path+"/hairpin_mature.pairs")+"\n")
checklog.writelines("4.Known_miRNA/"+"known_miRNA.map.fas"+checkfile(known_path+"/known_miRNA.map.fas")+"\n")
checklog.writelines("4.Known_miRNA/"+"known_miRNA.map.stat"+checkfile(known_path+"/known_miRNA.map.stat")+"\n")
checklog.writelines("4.Known_miRNA/"+"mature.fa"+checkfile(known_path+"/mature.fa")+"\n")
checklog.writelines("4.Known_miRNA/"+"mature.readcount"+checkfile(known_path+"/mature.readcount")+"\n")
checklog.writelines("4.Known_miRNA/"+"miRBase.mrd"+checkfile(known_path+"/miRBase.mrd")+"\n")
for eachsample in samples:
    path_3=os.path.join(out_dir,"4.Known_miRNA")
    checklog.writelines("/4.Known_miRNA/"+eachsample+".firstbase"+checkfile(path_3+"/"+eachsample+".firstbase")+"\n")
    checklog.writelines("/4.Known_miRNA/"+eachsample+".firstbase.png"+checkfile(path_3+"/"+eachsample+".firstbase.png")+"\n") 
    checklog.writelines("/4.Known_miRNA/"+eachsample+".position"+checkfile(path_3+"/"+eachsample+".position")+"\n")
    checklog.writelines("/4.Known_miRNA/"+eachsample+".position.png"+checkfile(path_3+"/"+eachsample+".position.png")+"\n")
###########check 5.ncRNA ####################
ncRNA_path=os.path.join(out_dir,"5.ncRNA")
checklog.writelines("/5.ncRNA/ncRNA.map.fas"+checkfile(ncRNA_path+"/ncRNA.map.fas")+"\n")
checklog.writelines("/5.ncRNA/rRNA.map.fas"+checkfile(ncRNA_path+"/rRNA.map.fas")+"\n")
checklog.writelines("/5.ncRNA/snoRNA.map.fas"+checkfile(ncRNA_path+"/snoRNA.map.fas")+"\n")
checklog.writelines("/5.ncRNA/snRNA.map.fas"+checkfile(ncRNA_path+"/snRNA.map.fas")+"\n")
checklog.writelines("/5.ncRNA/tRNA.map.fas"+checkfile(ncRNA_path+"/tRNA.map.fas")+"\n")
checklog.writelines("/5.ncRNA/rc.stat"+checkfile(ncRNA_path+"/rc.stat")+"\n")
checklog.writelines("/5.ncRNA/uc.stat"+checkfile(ncRNA_path+"/uc.stat")+"\n")
###########check 6.Repeat for ref  ####################
if argv['org']=='refplant' or argv['org']=='refanimal':
    repeat_path=os.path.join(out_dir,"6.Repeat")
    checklog.writelines("/6.Repeat/repeat.rc.stat"+checkfile(repeat_path+"/repeat.rc.stat")+"\n")
    checklog.writelines("/6.Repeat/repeat.uc.stat"+checkfile(repeat_path+"/repeat.uc.stat")+"\n")
    checklog.writelines("/6.Repeat/repeat.map.fas"+checkfile(repeat_path+"/repeat.map.fas")+"\n")
    for eachsample in samples:
        checklog.writelines("/6.Repeat/rc.stat_"+eachsample+"_bar.png"+checkfile(repeat_path+"/rc.stat_"+eachsample+"_bar.png")+"\n")
        checklog.writelines("/6.Repeat/rc.stat_"+eachsample+"_bar.pdf"+checkfile(repeat_path+"/rc.stat_"+eachsample+"_bar.pdf")+"\n")
        checklog.writelines("/6.Repeat/uc.stat_"+eachsample+"_bar.png"+checkfile(repeat_path+"/uc.stat_"+eachsample+"_bar.png")+"\n")
        checklog.writelines("/6.Repeat/uc.stat_"+eachsample+"_bar.pdf"+checkfile(repeat_path+"/uc.stat_"+eachsample+"_bar.pdf")+"\n")
###########check 7.NAT for refplant  ####################        
if argv['org']=='refplant':
    nat_path=os.path.join(out_dir,"7.NAT")
    checklog.writelines("/7.NAT/cis-NAT.map.fas"+checkfile(nat_path+"/cis-NAT.map.fas")+"\n")
    checklog.writelines("/7.NAT/NAT.map.fas"+checkfile(nat_path+"/NAT.map.fas")+"\n")
    checklog.writelines("/7.NAT/trans-NAT.map.fas"+checkfile(nat_path+"/trans-NAT.map.fas")+"\n")
    checklog.writelines("/7.NAT/rc.stat"+checkfile(nat_path+"/rc.stat")+"\n")
    checklog.writelines("/7.NAT/uc.stat"+checkfile(nat_path+"/uc.stat")+"\n")
    checklog.writelines("/7.NAT/NAT.rc.stat"+checkfile(nat_path+"/NAT.rc.stat")+"\n")
    checklog.writelines("/7.NAT/NAT.uc.stat"+checkfile(nat_path+"/NAT.uc.stat")+"\n")
###########check Novel_miRNA  ####################  
if argv['org']:
    if argv['org']=='refplant':
        num="9"
    elif argv['org']=='refanimal':
        num="8"
    elif argv['org']=='norefanimal' or argv['org']=='norefplant':
        num="6"
    novel_path=os.path.join(out_dir,num+".Novel_miRNA")
    checklog.writelines("Novel_miRNA/hairpin.fa"+checkfile(novel_path+"/hairpin.fa")+"\n")
    checklog.writelines("Novel_miRNA/hairpin_mature.fa"+checkfile(novel_path+"/hairpin_mature.fa")+"\n")
    checklog.writelines("Novel_miRNA/hairpin_mature.pairs"+checkfile(novel_path+"/hairpin_mature.pairs")+"\n")
    checklog.writelines("Novel_miRNA/novel_miRNA.map.fas"+checkfile(novel_path+"/novel_miRNA.map.fas")+"\n")
    checklog.writelines("Novel_miRNA/novel_miRNA.map.stat"+checkfile(novel_path+"/novel_miRNA.map.stat")+"\n")
    checklog.writelines("Novel_miRNA/mature.readcount"+checkfile(novel_path+"/mature.readcount")+"\n")
    checklog.writelines("Novel_miRNA/miRBase.mrd"+checkfile(novel_path+"/miRBase.mrd")+"\n")
    for eachsample in samples:
        checklog.writelines("Novel_miRNA/"+eachsample+".firstbase"+checkfile(novel_path+"/"+eachsample+".firstbase")+"\n")
        checklog.writelines("Novel_miRNA/"+eachsample+".firstbase.png"+checkfile(novel_path+"/"+eachsample+".firstbase.png")+"\n")
        checklog.writelines("Novel_miRNA/"+eachsample+".position"+checkfile(novel_path+"/"+eachsample+".position")+"\n")
        checklog.writelines("Novel_miRNA/"+eachsample+".position.png"+checkfile(novel_path+"/"+eachsample+".position.png")+"\n")
###########check gene for ref ####################        
if argv['org']=='refplant' or argv['org']=='refanimal':  
    if argv['org']=='refplant':
        num="8"
    elif argv['org']=='refanimal':
        num="7"
    gene_path=os.path.join(out_dir,num+".gene")
    checklog.writelines("gene/"+"intron.map.fas"+checkfile(gene_path+"/intron.map.fas")+"\n")
    checklog.writelines("gene/"+"gene.map.fas"+checkfile(gene_path+"/gene.map.fas")+"\n")
    checklog.writelines("gene/"+"rc.stat"+checkfile(gene_path+"/rc.stat")+"\n")
    checklog.writelines("gene/"+"exon.map.fas"+checkfile(gene_path+"/exon.map.fas")+"\n")
    checklog.writelines("gene/"+"uc.stat"+checkfile(gene_path+"/uc.stat")+"\n")
###########check Category ####################
if argv['org']:
    if argv['org']=='refplant':
        num="11"
    elif argv['org']=='refanimal':
        num="9"
    elif argv['org']=='norefanimal':
        num="7"
    elif argv['org']=='norefplant':
        num="8"
    Category_path=os.path.join(out_dir,num+".Category")
    checklog.writelines("Category/category_rc_full.txt"+checkfile(Category_path+"/category_rc_full.txt")+"\n")
    checklog.writelines("Category/category_uc_full.txt"+checkfile(Category_path+"/category_uc_full.txt")+"\n")
    for eachsample in samples:
        checklog.writelines("Category/"+eachsample+".category_rc_pie.pdf"+checkfile(Category_path+"/"+eachsample+".category_rc_pie.pdf")+"\n")
        checklog.writelines("Category/"+eachsample+".category_rc_pie.png"+checkfile(Category_path+"/"+eachsample+".category_rc_pie.png")+"\n")
        checklog.writelines("Category/"+eachsample+".category_uc_pie.pdf"+checkfile(Category_path+"/"+eachsample+".category_uc_pie.pdf")+"\n")
        checklog.writelines("Category/"+eachsample+".category_uc_pie.png"+checkfile(Category_path+"/"+eachsample+".category_uc_pie.png")+"\n")
###########check miRNA_editing ####################
if argv['org']:
    if argv['org']=='refplant':
        num="12"
    elif argv['org']=='refanimal':
        num="10"
    elif argv['org']=='norefanimal':
        num="8"
    elif argv['org']=='norefplant':
        num="9" 
    edit_path=os.path.join(out_dir,num+".miRNA_editing")
    for eachsample in samples:
        checklog.writelines("miRNA_editing/"+eachsample+"_editing_stats.example.txt"+checkfile(edit_path+"/"+eachsample+"_editing_stats.example.txt")+"\n")
        checklog.writelines("miRNA_editing/"+eachsample+"_editing_stats.example.txt"+checkfile(edit_path+"/"+eachsample+"_editing_stats.txt")+"\n")
###########check miRNA_family ####################
if argv['org']:
    if argv['org']=='refplant':
        num="13"
    elif argv['org']=='refanimal':
        num="11"
    elif argv['org']=='norefanimal':
        num="9"
    elif argv['org']=='norefplant':
        num="10"
    family_path=os.path.join(out_dir,num+".miRNA_family")
    checklog.writelines("miRNA_family/"+abbrs[0]+"_miRNA_family.detail.txt"+checkfile(family_path+"/"+abbrs[0]+"_miRNA_family.detail.txt")+"\n")
    checklog.writelines("miRNA_family/"+abbrs[0]+"_miRNA_family.mir_num.txt"+checkfile(family_path+"/"+abbrs[0]+"_miRNA_family.mir_num.txt")+"\n")
    checklog.writelines("miRNA_family/"+abbrs[0]+"_miRNA_family.mir_sign.txt"+checkfile(family_path+"/"+abbrs[0]+"_miRNA_family.mir_sign.txt")+"\n")
    checklog.writelines("miRNA_family/"+abbrs[0]+"_miRNA_family.mir_sign.example.txt"+checkfile(family_path+"/"+abbrs[0]+"_miRNA_family.mir_sign.example.txt")+"\n")  
###########check TAS for plant ####################
if argv['org']=='refplant' or argv['org']=='norefplant':
    if argv['org']=='refplant':
        num="10"
    elif argv['org']=='norefplant':
        num="7"
    tas_path=os.path.join(out_dir,num+".TAS")
    checklog.writelines("TAS/TAS.uc.stat"+checkfile(tas_path+"/TAS.uc.stat")+"\n")
    checklog.writelines("TAS/TAS.map.fas"+checkfile(tas_path+"/TAS.map.fas")+"\n")
    checklog.writelines("TAS/TAS.rc.stat"+checkfile(tas_path+"/TAS.rc.stat")+"\n")
###########check DiffExprAnalysis ####################
if argv['org']:
    if argv['org']=='refplant':
        num="14"
    elif argv['org']=='refanimal':
        num="12"
    elif argv['org']=='norefanimal':
        num="10"
    elif argv['org']=='norefplant':
        num="11"
    diff_path=os.path.join(out_dir,num+".DiffExprAnalysis")
    checklog.writelines("DiffExprAnalysis/miRNAExp/Readcount_TPM.xls"+checkfile(diff_path+"/"+num+".1.miRNAExp/Readcount_TPM.xls")+"\n")
    checklog.writelines("DiffExprAnalysis/miRNAExp/TPM_interval.xls"+checkfile(diff_path+"/"+num+".1.miRNAExp/TPM_interval.xls")+"\n")
    checklog.writelines("DiffExprAnalysis/miRNAExpdensity/TPM_boxplot.pdf"+checkfile(diff_path+"/"+num+".2.miRNAExpdensity/TPM_boxplot.pdf")+"\n")
    checklog.writelines("DiffExprAnalysis/miRNAExpdensity/TPM_boxplot.png"+checkfile(diff_path+"/"+num+".2.miRNAExpdensity/TPM_boxplot.png")+"\n")
    checklog.writelines("DiffExprAnalysis/miRNAExpdensity/TPM_density_distribution.pdf"+checkfile(diff_path+"/"+num+".2.miRNAExpdensity/TPM_density_distribution.pdf")+"\n")
    checklog.writelines("DiffExprAnalysis/miRNAExpdensity/TPM_density_distribution.png"+checkfile(diff_path+"/"+num+".2.miRNAExpdensity/TPM_density_distribution.png")+"\n")
    checklog.writelines("DiffExprAnalysis/CorAnalysis/cor_kendall.xls"+checkfile(diff_path+"/"+num+".3.CorAnalysis/cor_kendall.xls")+"\n")
    checklog.writelines("DiffExprAnalysis/CorAnalysis/cor_pearson.pdf"+checkfile(diff_path+"/"+num+".3.CorAnalysis/cor_pearson.pdf")+"\n")
    checklog.writelines("DiffExprAnalysis/CorAnalysis/cor_pearson.png"+checkfile(diff_path+"/"+num+".3.CorAnalysis/cor_pearson.png")+"\n")
    checklog.writelines("DiffExprAnalysis/CorAnalysis/cor_pearson.xls"+checkfile(diff_path+"/"+num+".3.CorAnalysis/cor_pearson.xls")+"\n")
    checklog.writelines("DiffExprAnalysis/CorAnalysis/cor_spearman.xls"+checkfile(diff_path+"/"+num+".3.CorAnalysis/cor_spearman.xls")+"\n")
#        for eachcompare in compare:
#            checklog.writelines("DiffExprAnalysis/14.3.CorAnalysis/"+eachcompare+".scatter.pdf"+checkfile(diff_path+"/14.3.CorAnalysis/"+eachcompare+".scatter.pdf"+"\n")
#            checklog.writelines("DiffExprAnalysis/14.3.CorAnalysis/"+eachcompare+".scatter.png"+checkfile(diff_path+"/14.3.CorAnalysis/"+eachcompare+".scatter.png"+"\n")
    for eachcompare in compare_name2:
        checklog.writelines("DiffExprAnalysis/DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DElist_up.txt"+checkfile(diff_path+"/"+num+".4.DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DElist_up.txt")+"\n")
        checklog.writelines("DiffExprAnalysis/DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DE_up.xls"+checkfile(diff_path+"/"+num+".4.DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DE_up.xls")+"\n")
        checklog.writelines("DiffExprAnalysis/DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DE.xls"+checkfile(diff_path+"/"+num+".4.DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DE.xls")+"\n")
        checklog.writelines("DiffExprAnalysis/DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DE_down.xls"+checkfile(diff_path+"/"+num+".4.DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DE_down.xls")+"\n")
        checklog.writelines("DiffExprAnalysis/DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DElist.txt"+checkfile(diff_path+"/"+num+".4.DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DElist.txt")+"\n")
        checklog.writelines("DiffExprAnalysis/DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".Differential_analysis_results.xls"+checkfile(diff_path+"/"+num+".4.DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".Differential_analysis_results.xls")+"\n")
        checklog.writelines("DiffExprAnalysis/DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DElist_down.txt"+checkfile(diff_path+"/"+num+".4.DiffExprAnalysis/"+eachcompare+"/"+eachcompare+".DElist_down.txt")+"\n")
        checklog.writelines("DiffExprAnalysis/DEsFilter/"+eachcompare+"/"+".Volcanoplot.pdf"+checkfile(diff_path+"/"+num+".5.DEsFilter/"+eachcompare+".Volcanoplot.pdf")+"\n")
        checklog.writelines("DiffExprAnalysis/DEsFilter/"+eachcompare+"/"+".Volcanoplot.png"+checkfile(diff_path+"/"+num+".5.DEsFilter/"+eachcompare+".Volcanoplot.png")+"\n")
    checklog.writelines("DiffExprAnalysis/DEcluster/DE_union_for_cluster"+checkfile(diff_path+"/"+num+".6.DEcluster/DE_union_for_cluster")+"\n")
    checklog.writelines("DiffExprAnalysis/DEcluster/Hcluster_heatmap.detail.pdf"+checkfile(diff_path+"/"+num+".6.DEcluster/Hcluster_heatmap.detail.pdf")+"\n")
    checklog.writelines("DiffExprAnalysis/DEcluster/Hcluster_heatmap.pdf"+checkfile(diff_path+"/"+num+".6.DEcluster/Hcluster_heatmap.pdf")+"\n")
    checklog.writelines("DiffExprAnalysis/DEcluster/Hcluster_heatmap.png"+checkfile(diff_path+"/"+num+".6.DEcluster/Hcluster_heatmap.png")+"\n")
    checklog.writelines("DiffExprAnalysis/DEcluster/K_means_cluster/K_means_cluster.png"+checkfile(diff_path+"/"+num+".6.DEcluster/K_means_cluster/K_means_cluster.png")+"\n")
    checklog.writelines("DiffExprAnalysis/DEcluster/K_means_cluster/K_means_cluster.pdf"+checkfile(diff_path+"/"+num+".6.DEcluster/K_means_cluster/K_means_cluster.pdf")+"\n")
    checklog.writelines("DiffExprAnalysis/DEcluster/SOM_cluster/SOM_cluster.png"+checkfile(diff_path+"/"+num+".6.DEcluster/SOM_cluster/SOM_cluster.png")+"\n")    
    if argv['venn_cluster'] != None:
        for eachvenn in venn_cluster_vs_names:
            checklog.writelines("DEvenn/venn/"+eachvenn+"/"+eachvenn+".venn.png"+checkfile(diff_path+"/"+num+".7.DEvenn/venn/"+eachvenn+"/"+eachvenn+".venn.png")+"\n")
            checklog.writelines("DEvenn/venn/"+eachvenn+"/"+eachvenn+".venn.pdf"+checkfile(diff_path+"/"+num+".7.DEvenn/venn/"+eachvenn+"/"+eachvenn+".venn.pdf")+"\n")
            checklog.writelines("DEvenn/venn/"+eachvenn+"/"+eachvenn+".venn.xls"+checkfile(diff_path+"/"+num+".7.DEvenn/venn/"+eachvenn+"/"+eachvenn+".venn.xls")+"\n")

###########check Enrichment ####################
if argv['org']:
    if argv['org']=='refplant':
        num="16"
    elif argv['org']=='refanimal':
        num="14"
    elif argv['org']=='norefanimal':
        num="12"
    elif argv['org']=='norefplant':
        num="13"
    if argv['compare']:
        enrich_path=os.path.join(out_dir,num+".Enrichment")
        for eachcompare in compare_name2:
            checklog.writelines("Enrichment/"+eachcompare+"/"+eachcompare+".diffmiRNA-geneid"+checkfile(enrich_path+"/"+eachcompare+"/"+eachcompare+".diffmiRNA-geneid")+"\n")
            checklog.writelines("Enrichment/"+eachcompare+"/"+eachcompare+".diffmiRNAID"+checkfile(enrich_path+"/"+eachcompare+"/"+eachcompare+".diffmiRNAID")+"\n")
            checklog.writelines("GOenrichment/"+eachcompare+".DEG_Enriched_GO_bp_DAG.png"+checkfile(enrich_path+"/"+eachcompare+"/GOenrichment/"+eachcompare+".DEG_Enriched_GO_bp_DAG.png")+"\n")
            checklog.writelines("GOenrichment/"+eachcompare+".DEG_Enriched_GO_bp_DAG.pdf"+checkfile(enrich_path+"/"+eachcompare+"/GOenrichment/"+eachcompare+".DEG_Enriched_GO_bp_DAG.pdf")+"\n")
            checklog.writelines("GOenrichment/"+eachcompare+".DEG_Enriched_GO_classification_gene_count.txt"+checkfile(enrich_path+"/"+eachcompare+"/GOenrichment/"+eachcompare+".DEG_Enriched_GO_classification_gene_count.txt")+"\n")
            checklog.writelines("GOenrichment/"+eachcompare+".DEG_Enriched_GO_classification.pdf"+checkfile(enrich_path+"/"+eachcompare+"/GOenrichment/"+eachcompare+".DEG_Enriched_GO_classification.pdf")+"\n")
            checklog.writelines("GOenrichment/"+eachcompare+".DEG_Enriched_GO_classification.png"+checkfile(enrich_path+"/"+eachcompare+"/GOenrichment/"+eachcompare+".DEG_Enriched_GO_classification.png")+"\n")
            checklog.writelines("GOenrichment/"+eachcompare+".DEG_Enriched_GO_mf_DAG.pdf"+checkfile(enrich_path+"/"+eachcompare+"/GOenrichment/"+eachcompare+".DEG_Enriched_GO_mf_DAG.pdf")+"\n")
            checklog.writelines("GOenrichment/"+eachcompare+".DEG_Enriched_GO_mf_DAG.png"+checkfile(enrich_path+"/"+eachcompare+"/GOenrichment/"+eachcompare+".DEG_Enriched_GO_mf_DAG.png")+"\n")
            checklog.writelines("GOenrichment/"+eachcompare+".DEG_GO_enrichment_result.xls"+checkfile(enrich_path+"/"+eachcompare+"/GOenrichment/"+eachcompare+".DEG_GO_enrichment_result.xls")+"\n")
            checklog.writelines("KEGGenrichment/"+eachcompare+".DEG_enriched_KEGG_pathway_scatterplot.pdf"+checkfile(enrich_path+"/"+eachcompare+"/KEGGenrichment/"+eachcompare+".DEG_enriched_KEGG_pathway_scatterplot.pdf")+"\n")
            checklog.writelines("KEGGenrichment/"+eachcompare+".DEG_enriched_KEGG_pathway_scatterplot.png"+checkfile(enrich_path+"/"+eachcompare+"/KEGGenrichment/"+eachcompare+".DEG_enriched_KEGG_pathway_scatterplot.png")+"\n")
            checklog.writelines("KEGGenrichment/"+eachcompare+".DEG_enriched_KEGG_pathway_API.html"+checkfile(enrich_path+"/"+eachcompare+"/KEGGenrichment/"+eachcompare+".DEG_enriched_KEGG_pathway_API.html")+"\n")
            checklog.writelines("KEGGenrichment/"+eachcompare+".DEG_KEGG_pathway_enrichment_result.xls"+checkfile(enrich_path+"/"+eachcompare+"/KEGGenrichment/"+eachcompare+".DEG_KEGG_pathway_enrichment_result.xls")+"\n")
            checklog.writelines("KEGGenrichment/"+eachcompare+".DEG_enriched_KEGG_pathway_top20.xls"+checkfile(enrich_path+"/"+eachcompare+"/KEGGenrichment/"+eachcompare+".DEG_enriched_KEGG_pathway_top20.xls")+"\n")
###########check miRNA_target ####################
if argv['org']:
    if argv['org']=='refplant':
        num="15"
        target_path=os.path.join(out_dir,num+".miRNA_target")
        checklog.writelines("miRNA_target/"+abbrs[0]+"_targets.pairs.annotate"+checkfile(target_path+"/"+abbrs[0]+"_targets.pairs.annotate")+"\n")
        checklog.writelines("miRNA_target/"+abbrs[0]+"_targets.txt"+checkfile(target_path+"/"+abbrs[0]+"_targets.txt")+"\n")
        checklog.writelines("miRNA_target/"+abbrs[0]+"_targets.pairs"+checkfile(target_path+"/"+abbrs[0]+"_targets.pairs")+"\n")
        checklog.writelines("miRNA_target/"+abbrs[0]+"_targets.pairs_example.xls"+checkfile(target_path+"/"+abbrs[0]+"_targets.pairs_example.xls")+"\n")
    elif argv['org']=='refanimal':
        num="13"
        target_path=os.path.join(out_dir,num+".miRNA_target")
        checklog.writelines("miRNA_target/commom_target.xls"+checkfile(target_path+"/commom_target.xls")+"\n")
        checklog.writelines("miRNA_target/all_target_gene.xls"+checkfile(target_path+"/all_target_gene.xls")+"\n")
        checklog.writelines("miRNA_target/commom_target_example.xls"+checkfile(target_path+"/commom_target_example.xls")+"\n")
        checklog.writelines("miRNA_target/miranda_targets_out.fmt.gz"+checkfile(target_path+"/miranda_targets_out.fmt.gz")+"\n")
        checklog.writelines("miRNA_target/RNAhybrid_miRNA_target_pairs.gz"+checkfile(target_path+"/RNAhybrid_miRNA_target_pairs.gz")+"\n")
        checklog.writelines("miRNA_target/PITA_pita_results_targets.tab.gz"+checkfile(target_path+"/PITA_pita_results_targets.tab.gz")+"\n")
    elif argv['org']=='norefanimal':
        num="11"
        target_path=os.path.join(out_dir,num+".miRNA_target")
        checklog.writelines("miRNA_target/"+abbrs[0]+"_mature.miranda_targets_gene.pairs"+checkfile(target_path+"/"+abbrs[0]+"_mature.miranda_targets_gene.pairs")+"\n")
        checklog.writelines("miRNA_target/"+abbrs[0]+"_mature.miranda_targets.pairs"+checkfile(target_path+"/"+abbrs[0]+"_mature.miranda_targets.pairs")+"\n")
        checklog.writelines("miRNA_target/"+abbrs[0]+"_mature.miranda_targets_example.xls"+checkfile(target_path+"/"+abbrs[0]+"_mature.miranda_targets_example.xls")+"\n")
    elif argv['org']=='norefplant':
        num="12"
        target_path=os.path.join(out_dir,num+".miRNA_target")
        checklog.writelines("miRNA_target/"+abbrs[0]+"_targets.pairs.annotate"+checkfile(target_path+"/"+abbrs[0]+"_targets.pairs.annotate")+"\n")
        checklog.writelines("miRNA_target/"+abbrs[0]+"_targets.txt"+checkfile(target_path+"/"+abbrs[0]+"_targets.txt")+"\n")
        checklog.writelines("miRNA_target/"+abbrs[0]+"_targets.pairs"+checkfile(target_path+"/"+abbrs[0]+"_targets.pairs")+"\n")
        checklog.writelines("miRNA_target/"+abbrs[0]+"_targets.pairs_example.xls"+checkfile(target_path+"/"+abbrs[0]+"_targets.pairs_example.xls")+"\n") 
checklog.close()
os.chdir(outdir)
os.system('grep failed checkresultlog >failed_resultlog')
