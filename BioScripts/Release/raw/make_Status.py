#-*- coding:utf-8 -*-
#sun


import re
import os
import glob
import time

root_dir=os.getcwd()
projectTotalName=root_dir.split('/')[-1]
contractNumber = projectTotalName.split('_')[0]
data_give_gz=glob.glob(root_dir+'/*_give/*tar.gz')

projectNumber='_ _'
yunyingjingli='_ _'
xinxifenxi='_ _'

fenqi = projectTotalName.split('_')[-1]
fenqinumber = fenqi if fenqi.isdigit() else '1'

xinxifenxi = os.popen("grep '\-ownername' run.sh").read().strip().split()[1]
yunyingjingli = os.popen("grep '\-yunying' run.sh").read().strip().split()[1]

project = os.popen("grep '\-contract' run.sh").read().strip().split()[1]
projectNumber =project.split('_')[0]

###############  Get the create date of tar.gz   #######################

if len(data_give_gz)!=0:
	filetime=os.stat(data_give_gz[0])
        arrytime=time.localtime(filetime.st_ctime)
#	print arrytime[1]
	year = arrytime[0]
	month = arrytime[1]
	date = arrytime[2]
	date = str(date) if len(str(date)) >1 else str(0)+str(date)

	if month + 3 > 12 :
		endtime = str(year)+str(month)+date
		deadline = str(year+1)+str(0)+str(month-9)+date
	else:
		endtime = str(year)+str(0)+str(month)+date
		if len(str(month+3))==1:
			deadline = str(year)+str(0)+str(month+3)+date
		else:
			deadline = str(year)+str(month+3)+date
	if month + 6 > 12 :
		deadline6 = str(year+1)+str(0)+str(month-6)+date
	else:
		if len(str(month+6))==1:
			deadline6 = str(year)+str(0)+str(month+6)+date
		else:
			deadline6 = str(year)+str(month+6)+date

else:
	exit('please cheak tar.gz')


statusFileName = projectTotalName+'_'+str(deadline)+'_'+str(deadline6)+'_Status'
fileStatus=open('/TJPROJ1/RNA/WORK/project_status/'+statusFileName,'w')


###############  Grab information from oms!    #######################

fileStatus.write("交付日期	合同编号	项目编号	操作步骤（3-5）	删除期限	是否为合作项目或普通项目	信息负责人	运营经理	路径	执行脚本step3路径	执行脚本step4路径	执行脚本step5路径\n")
fileStatus.write('\t'.join([endtime,contractNumber,projectNumber,'2',deadline,'N',xinxifenxi,yunyingjingli,root_dir,root_dir+'/shell/backup.sh',root_dir+'/shell/del1.sh',root_dir+'/shell/delraw.sh'])+'\n')


fileStatus.close()
#os.system('cp {statusFile} {status_dir}'.format(statusFile=statusFileName,status_dir=total_status_dir))

