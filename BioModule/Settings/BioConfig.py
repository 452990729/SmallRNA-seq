#!/usr/bin/env python
#_*_ coding:utf-8 _*_

import os
import sys
import ConfigParser

BASE_DIR = os.path.dirname(__file__)


#get configini
config_ini = os.path.join(BASE_DIR,'config.ini')
config     = ConfigParser.ConfigParser()
config.read(config_ini)


#Analysis code

SYMBOL_CODE = {1.1:'quality_control',
               1.2:'quality_control_and_sequence_sharedif',
               2.1:'mapping_ref',
               2.2:'mapping_noref'
              } 

#default analysis content
ORG_ANALYSIS_CODE = {'refplant': [1.1,2.1,3.1],
                     'norefplant': [1.1,2.2,3.1],
                     'refanimal': [1.1,2.1,3.1],
                     'norefanimal': [1.1,2.2,3.1]
                    }


# default memory settings

MEMORY = {
        'rawdata_zcat': '2G',
        'raw2clean': '2G',
        'adapter2clean' : '1G',
        'category' : '1G',
        'sharespecific' : '2G',
        'MapStat': '2G',
        'qc_report' : '1G',
        'circos' : '1G',
        'get_hairpin_mature_seq': '1G',
        'known_miRNA_quantify': '1G',
        'known_miRNA_csv2pdf': '1G',
        'known_miRNA_stat': '1G',
        'ncRNA_map_and_stat': '1G',
        'repeat_mapping': '2G',
        'NAT_siRNA_plant': '2G',
        'ncRNA_map_and_stat':'2G',
        'exon_and_intron_build_index':'4G',
        'predict_novel_miRNA': '10G',
        'novel_miRNA_quantify':'2G',
        'novel_miRNA_stat':'2G',
        'known_TAS_and_novel_TAS_analysis': '2G',
        'TAS_mapping_and_TAS_stat': '2G',
        'sRNA_Category':'2G',
        'edit_family_analy':'1G',
        'diffAnalysis':'1G',
        'mirna_target_refplant':'1G',
        'mirna_target_norefplant':'1G',
        'kegg_blast_refplant':'3G',
        'kegg_blast_refanimal':'2G',
        'GO_and_kegg_enrich':'2G',
}
