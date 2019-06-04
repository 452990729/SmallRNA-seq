# Embedded file name: pathway_annotation_flow_parallel_simple_tolerant.py
import re
import os
import sys
import os.path
import urllib
import argparse
import ConfigParser
from bs4 import BeautifulSoup as bs
from multiprocessing import Pool, cpu_count
from django.template import Template, Context, loader
from django.conf import settings
settings.configure(DEBUG=True, TEMPLATE_DEBUG=True, TEMPLATE_DIRS=(sys.path[0],))

skipp_patn  = re.compile(r'.*ko(01100|01110).*')
uncomplete_pathway_list = []
line_num = 0

config = ConfigParser.ConfigParser()
config.read('{}/../../Pipeline/config.ini'.format(sys.path[0]))
Pathway = config.get('database','pathway')

def file_wash(filename):
    washed = []
    k = 0
    tmp = []
    for eachLine in open(filename):
        tmp.append(eachLine.strip())
    for line_stat in tmp:
        if(line_stat.startswith('#Term')):
	    global line_num
            line_num = len(line_stat.split("\t"))
            break
    for line in tmp:
        if line != '' and not line.startswith('#Term') and len(line.split('\t')) == line_num :
            washed.append(line.strip() + '\t' + str(k))
            k = k + 1
    return washed

def read_html(abbr,pathway):
    filepath = '{}/'.format(Pathway)+ abbr + "/" + pathway + ".html"
    file_object = open(filepath)
    all_the_text = file_object.read( )
    file_object.close( )
    return all_the_text
def copy_img(abbr,pathway,target_path):
    filepath = '{}/'.format(Pathway) + abbr + "/" +  pathway + ".png"
    assert not os.system('cp ' + filepath + " " + target_path)
def main_flow(abbr,row):
    global line_num
    
    each = row.strip().split('\t')
    pathway = each[2].strip().lower()
    for col in range(7, line_num-1):
        temp = [ each_temp.strip() for each_temp in each[col].split('|') if each_temp.strip() != '' ]
        if list(set(temp)) == [] or list(set(temp)) == ['NA']:
            each[col] = 'NA'
        else:
            each[col] = ' '.join(temp)

    ko = [ one_ko.strip() for one_ko in each[8].strip().split() ]
    gene = [ one_gene.strip() for one_gene in each[7].strip().split() ]
    try:
        #raise len(ko) == len(gene) or AssertionError
        assert len(ko) == len(gene)
    except:
        sys.exit('assert error!' + pathway)

    ko_gene = {}
    for i, each_ko in enumerate(ko):
        if each_ko not in ko_gene:
            ko_gene[each_ko] = []
        ko_gene[each_ko].append(gene[i])
    try:
	
        content = read_html(abbr,pathway)
	
        soup = bs(content)
        copy_img(abbr,pathway,'src/')
        map = soup.map
        html_content = ""
        if pathway.endswith("01200") or pathway.endswith("01210") or pathway.endswith("01212") or pathway.endswith("01230") or pathway.endswith("01220") or pathway.endswith("01502"):
                for each_area in map.find_all('area'):
                    each_area['href'] = 'http://www.kegg.jp' + each_area['href']
        else:
            index = 0
            for each_area in map.find_all('area'):
                ko_set = [ each_ko.strip() for each_ko in re.search('\\?(.*)', each_area['href']).group(1).split('+') ]
                num = 0
                for each_ko in ko_gene:
		    #print each_ko
                    if each_ko in ko_set:
                        num = num + 1
                each_area['href'] = 'http://www.kegg.jp' + each_area['href']
		
                if num != 0:
		    #print pathway + ":" + str(num)
                    coords = each_area['coords'].split(",")
                    top = int(coords[1]) + 8
                    left = int(coords[0]) + 8
                    each_area['coords'] = "0,0,47,18"
                    html_content += "<div><span class='maps'><img class='maps' src='bg_pink.png' style='top:"+ str(top) +"px;left:"+ str(left) +"px;width:47px;height:18px;position:absolute;' usemap='#map_"+ str(index)  +"'>" + "<map name='map_"+ str(index) +"'>" + str(each_area.prettify(formatter=None)) + "</map></span> <span style='display:none;'>"
                    
                    inner_html = '<ul>'
                    for each_ko in ko_gene:
                        if each_ko in ko_set:
                            inner_html += '<li>%s</li>' % each_ko
                            inner_html += '<ul><li>%s</li></ul>' % ' '.join(ko_gene[each_ko])
                    inner_html += '</ul></span></div>'
                    html_content += inner_html
                index = index + 1
        t = loader.get_template('Kegg_map_template.html')
        c = Context({'title': pathway,'add_content':html_content,
         'map_content': str(map.prettify(formatter=None)),
         'image': pathway + '.png'})
        html = t.render(c)
        open('src/' + pathway + '.html', 'w').write(html)
        link = 'src/' + pathway + '.html'
        term = [link, each[0]] + each[3:-2] + [each[-1]]
    except:
        print pathway + ' detail map fail...'
        term = ['#', each[0]] + each[3:-2] + [each[-1]]

    return term

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='KEGG pathway enrichment web visualization, diff results not needed, used for sRNA lncRNA and alike')
    parser.add_argument('--table', required=True, help='standard formated pathway enrichment input file')
    parser.add_argument('--abbr', required=True, help='species abbr')
    argv = vars(parser.parse_args())
    filename = argv['table']
    abbr = argv['abbr']
    name = os.path.basename(filename)
    assert not os.system('cp -r %s .' % (sys.path[0] + '/src'))
    if not os.path.exists('src'): os.system('mkdir src')

    parallel_result = []
    row_pathway = file_wash(filename)
    pool = Pool(processes=cpu_count())
    for eachrow in row_pathway:
        parallel_result.append(pool.apply_async(main_flow, (abbr,eachrow,)))

    pool.close()
    pool.join()
    result = [ '' for i in range(len(row_pathway)) ]
    for eachone in parallel_result:
        try:
            temp = eachone.get(timeout=10)
            result[int(temp[-1])] = temp[0:-1]
        except:
            uncomplete_pathway = row_pathway[parallel_result.index(eachone)]
            sys.stderr.write( '%s not complete W/I 2mins' % uncomplete_pathway )
            if not re.match(skipp_patn, uncomplete_pathway):
                uncomplete_pathway_list.append(uncomplete_pathway)

    t = loader.get_template('Table_template.html')
    c = Context({'terms': result})
    html = t.render(c)
    open(name + '_rendered_html_detail.html', 'w').write(html)

if uncomplete_pathway_list:
    open('unfinished_pathway.txt','w').writelines(uncomplete_pathway_list)
