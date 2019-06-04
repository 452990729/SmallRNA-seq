'''Identify enriched pathways or human diseases by frequency of terms or statistical significance of terms.
Species abbreviation can be designated to indicate the preexistent background file.
'''

import os, sys
#sys.path.append('{}/../../BioDB/KEGG/kobas2.0-20120208/src'.format(sys.path[0]))

from optparse import OptionParser

from kobas2 import config, dbutils, discover, exception, output

def config_option():
    '''configure options and arguments for identify.py
    '''
    p = OptionParser()
    p.add_option(
        '-f', '--fgfile', dest = 'fgfile', action = 'store',
        type = 'string', help = 'foreground file')
    p.add_option(
        '-k', '--isko', dest = 'isko', action = 'store_true',
        help = 'whether the foreground file is ko annotation or not, default not ko')
    p.add_option(
        '-b', '--bgfile', dest = 'bgfile', action = 'store',
        type = 'string', help = 'background file, can be species abbreviation')
    p.add_option(
        '-d', '--db', dest = 'db', default = '/'.join(all_dbs), action = 'store',
        type = 'string', help = 'selected databases, one letter abbreviation separated by "/", \
K is KEGG PATHWAY, n is PID Curated, b is PID BioCarta, r is PID Reactome, B is BioCyc, R is Reactome, p is PANTHER, \
k is KEGG DISEASE, f is FunDO, o is OMIM, g is GAD, N is NHGRI, \
default %s' % '/'.join(all_dbs))
    p.add_option(
        '-m', '--method', dest = 'method', default = 'h', action = 'store',
        type = 'string', help = 'choose statistic method, \
b is binomial test, c is chi-square test, f is fisher exact test, h is hypergeometric test and x is frequency list, \
default hypergeometric test')
    p.add_option(
        '-n', '--fdr', dest = 'fdr', default = 'QVALUE', action = 'store',
        type = 'string', help = 'choose false discovery rate (FDR) correction method: QVALUE, BH, BY or none, default QVALUE')
    p.add_option(
        '-o', '--outfile', dest = 'outfile', action = 'store',
        type = 'string', help = 'output file for identification result, default stdout')
##  new options
    p.add_option(
        '-s', '--standard', dest = 'standard', action = 'store_true',
        help = 'whether the number of tests is standard or not, default not standard')
    p.add_option(
        '-c', '--cutoff', dest = 'cutoff', default = 5, action = 'store',
        type = 'int', help = 'the cutoff of the gene numbers in a term, default 5')
    opt, args = p.parse_args()
    return (p, opt, args)

METHODS = {'b': 'BinomTest', 'c': 'ChiTest', 'f': 'FisherTest', 'h': 'HyperTest', 'x': 'FreqList'}
FULLNAME = {'b': 'binomial test', 'c': 'chi-square test', 'f': 'fisher exact test', 'h': 'hypergeometric test', \
    'QVALUE': 'QVALUE', 'BH': 'Benjamini and Hochberg (1995)', 'BY': 'Benjamini and Yekutieli (2001)'}
all_dbs = set(dbutils.DATABASES.keys())

def file2distr(file, dbs, kobas2db, bg = None, small_num = 5):
    '''get distribution from annot file or species abbreviation
    '''
    if not os.access(file, os.F_OK):
        print 'file %s does not exist' % file
        sys.exit(1)
    return discover.distr_from_annot_file(open(file), dbs, kobas2db, bg, small_num)

def find_terms(res, test_method, fdistr, bdistr):
    '''find significant terms of sample from backgroud distribution
    test_method: method abbreviation
    fdistr: {term: amount}, sample distribution
    bdistr: {term: amount}, background distribution
    '''
    test_class = get_class(test_method)
    test = test_class(fdistr, bdistr)
    return test(res)

def get_class(test_method):
    return getattr(discover, METHODS[test_method])

if __name__ == '__main__':
    opt_parser, opt, args = config_option()

    if not opt.fgfile:
        opt_parser.print_help()
        sys.exit(1)

    if opt.method not in METHODS.keys():
        sys.exit('%s method is not supported yet, only %s is supported now' %
            (opt.method, ', '.join(METHODS.keys())))
    elif opt.method != 'x':
        if not opt.bgfile:
            opt_parser.print_help()
            sys.exit(1)

    if opt.outfile:
        global old_stdout
        old_stdout = sys.stdout
        sys.stdout = open(opt.outfile, 'w')
    
    # process databases option
    dbs = set(opt.db.split('/'))
    dbs.discard('') # for weblab web server
    other_dbs = dbs - all_dbs
    if len(other_dbs) != 0:
        sys.exit('%s database(s) is(are) not supported yet, only %s are supported now' % (', '.join(other_dbs), ', '.join(all_dbs)))

    # KOBAS environment configuration
    kobas2rc = config.getrc()

    # open kobas2db
    kobas2db = dbutils.kobas2db(kobas2rc['kobas2db'])    

    if opt.method == 'x':
        res = discover.TestResult(
            ['Term', 'Database', 'Id', 'Genes', 'Number', 'Proportion'])
        for db in dbs:
            tqdistr = file2distr(opt.fgfile, [db], kobas2db)[1]
            total_querys = tqdistr.size()
            for (term, querys) in tqdistr.items():
                number = len(querys)
                proportion = float(number) / total_querys
                res.add_row([term[2], term[1], term[0], '|'.join(querys), '%d / %d' % (number, total_querys), proportion])
        res.sort(order = 1)
        print '##The most frequent pathway or disease terms: \n'

    else:
        res = discover.TestResult(
            ['Term', 'Database', 'Id', 'Genes', 'Sample number', 'Background number', 'P-Value'])

        bg_term_num = 0
        if opt.bgfile in kobas2db.get_species_by_databases(only_abbr = True):
            FILE = 'default_fg'
        else:
            FILE = 'annotated_bg'

        for db in dbs:
            # foreground distr and background distr
            if FILE == 'default_fg':
                ftqdistr = file2distr(opt.fgfile, [db], kobas2db, bg = FILE, small_num = opt.cutoff)[1]
                btqdistr = discover.distr_from_species_abbr(opt.bgfile, [db], kobas2db, opt.isko, opt.cutoff)[1]
            else:
                ftqdistr = file2distr(opt.fgfile, [db], kobas2db, bg = None, small_num = opt.cutoff)[1]
                btqdistr = file2distr(opt.bgfile, [db], kobas2db, bg = FILE, small_num = opt.cutoff)[1]

            bg_term_num = bg_term_num + len(btqdistr)

            # Do statistical test
            try:
                res = find_terms(res, opt.method, ftqdistr, btqdistr)
            except ArithmeticError, msg:
                exception.error(msg)
                sys.exit(1)

        # Do FDR correction
        if opt.fdr == 'none':
            res.sort(key = -1)
        else:
            try:
                res.add_fdr(method = opt.fdr, test_num = (opt.standard and bg_term_num))
            except ArithmeticError, msg:
                exception.error(msg)
                res.sort(key = -1)
            except TypeError, msg:
                raise
            else:
                res.sort(key = -2)
        print '##The most enriched pathway or disease terms: '
        print '##Statistic method: %s ' % FULLNAME[opt.method]
        if opt.fdr != 'none':
            print '##FDR correction method: %s ' % FULLNAME[opt.fdr]
        print 

    # get species_abbr
    file = open(opt.fgfile)
    line = file.readline()
    file.close()
    species_abbr = line.split(' ')[1]

    output.identify(res, species_abbr)

