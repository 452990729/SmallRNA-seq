#!/usr/bin/python

import os, sys
from optparse import OptionParser

from norefkbs import annot, config, dbutils, discover, exception, output

def config_option():
    usage = 'Usage: %prog -f fgfile [-b bgfile] [-d databases] [-m test] [-n fdr] [-o outfile] [-c cutoff]'
    p = OptionParser(usage)
    p.add_option(
        '-f', '--fgfile', dest = 'fgfile', action = 'store',
        type = 'string', help = 'foreground file, the output of annotate')
    p.add_option(
        '-b', '--bgfile', dest = 'bgfile', default = 'same', action = 'store',
        type = 'string', help = 'background file, the output of annotate (3 or 4-letter file name is not allowed), or species abbreviation (for example: hsa for Homo sapiens, mmu for Mus musculus, dme for Drosophila melanogaster, ath for Arabidopsis thaliana, sce for Saccharomyces cerevisiae and eco for Escherichia coli K-12 MG1655), default same species as annotate')
    p.add_option(
        '-d', '--db', dest = 'db', default = 'K/n/b/R/B/p/o/k/f/g/N/G', action = 'store',
        type = 'string', help = 'databases for selection, 1-letter abbreviation separated by "/": K for KEGG PATHWAY, n for PID, b for BioCarta, R for Reactome, B for BioCyc, p for PANTHER, o for OMIM, k for KEGG DISEASE, f for FunDO, g for GAD, N for NHGRI GWAS Catalog and G for Gene Ontology, default K/n/b/R/B/p/o/k/f/g/N/G')
    p.add_option(
        '-m', '--method', dest = 'method', default = 'h', action = 'store',
        type = 'string', help = "choose statistical test method: b for binomial test, c for chi-square test, h for hypergeometric test / Fisher's exact test, and x for frequency list, default hypergeometric test / Fisher's exact test")
    p.add_option(
        '-n', '--fdr', dest = 'fdr', default = 'BH', action = 'store',
        type = 'string', help = 'choose false discovery rate (FDR) correction method: BH for Benjamini and Hochberg, BY for Benjamini and Yekutieli, QVALUE, and None, default BH')
    p.add_option(
        '-o', '--outfile', dest = 'outfile', action = 'store',
        type = 'string', help = 'output file for identification result, default stdout')
    p.add_option(
        '-c', '--cutoff', dest = 'cutoff', default = 5, action = 'store',
        type = 'int', help = 'terms with less than cutoff number of genes are not used for statistical tests, default 5')
    p.add_option(
	'-t', '--type', dest = 'type', default = 'all', action = 'store',
        type = 'string', help = 'choose KEGG pathway type')
    opt, args = p.parse_args()
    return (p, opt, args)

METHODS = {'b': 'BinomTest', 'c': 'ChiTest', 'h': 'FisherTest', 'x': 'FreqList'}
FULLNAMES = {'b': 'binomial test', 'c': 'chi-square test', 'h': "hypergeometric test / Fisher's exact test",
             'BH': 'Benjamini and Hochberg', 'BY': 'Benjamini and Yekutieli', 'QVALUE': 'QVALUE'}

def verify_file(file_name):
    ##verify file existence and accession, and return file handle
    if os.access(file_name, os.F_OK):
        return open(file_name)
    else:
        print 'file %s does not exist' % file_name
        sys.exit(1)

def get_items(handle):
    items = set()
    for item in annot.Iterator(handle):
        if item.has_links():
            items.add((item.query, item.links))
    return items

def get_class(test_method):
    return getattr(discover, METHODS[test_method])

def find_terms(res, test_method, fdistr, bdistr):
    ##find significant terms of input from backgroud distribution
    test_class = get_class(test_method)
    test = test_class(fdistr, bdistr)
    return test(res)

opt_parser, opt, args = config_option()
##define ptype
type = opt.type
if type == 'nr':
	type = 'n'
##process opt.fgfile
if not opt.fgfile:
    opt_parser.error('Option -f must be assigned.\n')
else:
    fg_handle = verify_file(opt.fgfile)
    abbr, idmapping = annot.get_args(fg_handle)
    fg_items = get_items(fg_handle)

##process opt.method
if opt.method not in METHODS.keys():
    opt_parser.error('%s method is not supported yet, only %s are supported now.\n' % (opt.method, ', '.join(METHODS.keys())))

##process opt.bgfile
if len(opt.bgfile) in [3, 4]:
    default_bg = True
    if opt.bgfile == 'same':
        if opt.method != 'x' and abbr == 'ko':
            opt_parser.error('A background file or a species abbreviation must be assigned for statistical test for KO.\n')
        else:
            opt.bgfile = abbr
else:
    default_bg = False
    bg_handle = verify_file(opt.bgfile)
    bg_items = get_items(bg_handle)

##process opt.outfile
if opt.outfile:
    global old_stdout
    old_stdout = sys.stdout
    sys.stdout = open(opt.outfile, 'w')

##KOBAS environment configuration
kobasrc = config.getrc()

##open KOBASDB
speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + abbr + '.db')

##process opt.db
input_dbs = set(opt.db.split('/'))
avail_dbs = {}
for database in speciesdb.databases_from_abbr(abbr):
    avail_dbs[database[0]] = database[1]
dbs = input_dbs.intersection(set(avail_dbs.keys()))
if dbs:
    print '##Databases: %s' % ', '.join([avail_dbs[db] for db in dbs])
else:
    opt_parser.error('No supported databases are selected. Supported databases are %s, but your input databases are: %s.' % \
        ('/'.join(avail_dbs.keys()), opt.db))

odistr = discover.Distr()

if opt.method == 'x':
    res = discover.TestResult(
        ['Term', 'Database', 'ID', 'Number', 'Proportion', 'Input'])
    for db in dbs:
        if db == 'K':
            fdistr, odistr = discover.distr_from_file(fg_items, db, speciesdb, abbr, type)
        else:
            fdistr = discover.distr_from_file(fg_items, db, speciesdb, abbr, type)
        total_num = fdistr.size()
        for (term, queries) in fdistr.items():
            num = len(queries)
            proportion = float(num) / total_num
            res.add_row(list(term) + ['%d' % num, '%.3g' % proportion, '|'.join(list(queries))])
    res.sort(key = -2, order = 1)

else:
    res = discover.TestResult(
        ['Term', 'Database', 'ID', 'Input number', 'Background number', 'P-Value', 'Input'])
    for db in dbs:
        ##foreground distr and background distr
        if db == 'K':
            fdistr, odistr = discover.distr_from_file(fg_items, db, speciesdb, abbr, type)
        else:
            fdistr = discover.distr_from_file(fg_items, db, speciesdb, abbr, type)
        if default_bg:
            bdistr = discover.distr_from_default(opt.bgfile, db, speciesdb, abbr, type, opt.cutoff, idmapping)
        else:
            bdistr = discover.distr_from_file(bg_items, db, speciesdb, abbr, type, opt.cutoff)
        ##do statistical test
        try:
            res = find_terms(res, opt.method, fdistr, bdistr)
        except ArithmeticError, msg:
            exception.error(msg)
            sys.exit(1)
    res.sort(key = -2)
    ##do FDR correction
    if opt.fdr != 'None':
        try:
            res.add_fdr(method = opt.fdr)
        except ArithmeticError, msg:
            exception.error(msg)
            sys.exit(1)
        except TypeError, msg:
            exception.error(msg)
            sys.exit(1)

##report identification result
    print '##Statistical test method: %s' % FULLNAMES[opt.method]
    if opt.fdr != 'None':
        print '##FDR correction method: %s' % FULLNAMES[opt.fdr]

output.identify(res, abbr, odistr)
