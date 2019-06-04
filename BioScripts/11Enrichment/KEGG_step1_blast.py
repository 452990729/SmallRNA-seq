#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = "Liang Peiquan (liangpeiquan@novogene.com)"
__date__ = "$Date: 2014/12/16 $"

import sys
import os
import os.path
import glob
import argparse
import ConfigParser

BASIC_DIR = os.path.dirname(os.path.abspath(__file__))
config = ConfigParser.ConfigParser()
config.read('{}/../../BioModule/Settings/config.ini'.format(BASIC_DIR))
kobas_v20140801 = config.get('database','kobas_v20140801')
blast = config.get('software','blast2')

sp_list = "{}/abbr_list.txt".format(kobas_v20140801)
sp_data = "{}/seq_pep/".format(kobas_v20140801)

if len( sys.argv ) != 5:
	sys.stderr.write("Script to prepare KEGG Enrichment.\n\n")
	sys.stderr.write("Usage: python %s <in.diffgene.fasta> <species> <out.blastout> <out.script>\n\n" % sys.argv[0] )
	sys.stderr.write( "This script is used to run KEGG blast\n" )
	sys.exit(1)

diffgene = sys.argv[1]
species = sys.argv[2]
blastout = sys.argv[3]
script = sys.argv[4]

number_of_processors = '4'
e_value = '1e-5'

abbr = []
sp_file = open(sp_list)
all_text = sp_file.read()
list = all_text.strip().split('\n')
if species not in list:
	print "Your species is not in sqlite!!!\n"
	print "Please check %s !!!\n" % (sp_list)
else:
	for each in glob.glob('{}/seq_pep/*.fasta'.format(kobas_v20140801)):
		filename = os.path.basename(each).split('.')
		filename_abbr = filename[0]
		abbr.append(filename_abbr)
	if species not in abbr:
		print "Your species is not in seq_pep!!!\n"
		print "Please check %s !!!" % (sp_data)
	else:
		sp_pep_seq = sp_data + species + '.pep.fasta'

		code = '''
%s/blastx -query %s -db %s -evalue %s -outfmt 5 -max_target_seqs 1 -num_threads %s -out %s\n
''' % (blast,diffgene, sp_pep_seq, e_value, number_of_processors, blastout)

		open(script,'w').write(code)
