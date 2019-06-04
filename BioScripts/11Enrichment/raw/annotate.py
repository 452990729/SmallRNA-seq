#!/PUBLIC/software/public/System/Python-2.7.6/bin/python
'''Annotate a set of sequences with KEGG GENES/ORTHOLOGY.
'''

import sys
from optparse import OptionParser

from kobas2 import annot, config, dbutils, exception, fasta, output

def config_option():
    '''configure options and arguments for annotate.py
    '''
    p = OptionParser()
    # options for information
    p.add_option(
        '-l', '--list', dest = 'list', action = 'store_true',
        help = 'list available species, or list available databases for specific species')
    # basic options
    p.add_option(
        '-i', '--infile', dest = 'infile', action = 'store',
        type = 'string', help = 'input data file')
    p.add_option(
        '-t', '--intype', dest = 'intype', default = 'fasta:pro', action = 'store',
        type = 'string', help = 'input type (%s, blastout:xml, blastout:tab, %s), default fasta:pro' %
            (', '.join(PROGRAMS.keys()), ', '.join(DBLINKS.keys())))
    p.add_option(
        '-s', '--species', dest = 'species', action = 'store',
        type = 'string', help = 'species abbr for databases')
    p.add_option(
        '-o', '--outfile', dest = 'outfile', action = 'store',
        type = 'string', help = 'output file for annotation result, default stdout')
    # options for blast and parsing blast result
    p.add_option(
        '-e', '--evalue', dest = 'evalue', default = 1e-5, action = 'store',
        type = 'float', help = 'expect threshold for BLAST, default 1e-5')
    p.add_option(
        '-r', '--rank', dest = 'rank', default = 5, action = 'store',
        type = 'int', help = 'rank cutoff for valid hits from BLAST result, default 5')
    p.add_option(
        '-n', '--nCPUs', dest = 'nCPUs', default = 1, action = 'store',
        type = 'int', help = 'number of CPUs to be used by BLAST, default 1')
    # option for id mapping
    p.add_option(
        '-S', '--inspecies', dest = 'inspecies', action = 'store_true',
        help = 'map id to other species')
    # option for reviewer
    p.add_option(
        '-c', '--coverage', dest = 'coverage', default = 0.0, action = 'store',
        type = 'float', help = 'subject coverage cutoff for BLAST, default 0')
    p.add_option(
        '-z', '--ortholog', dest = 'ortholog', default = 'NO', action = 'store',
        type = 'string', help = 'whether only use ortholog for cross-species annotation or not, default NO (If only use ortholog, give species abbr)')
    opt, args = p.parse_args()
    return (p, opt, args)

PROGRAMS = {'fasta:pro': 'blastp', 'fasta:nuc': 'blastx'}
DBLINKS = {'id:ncbigene': 'ng', 'id:ncbigi': 'gi', 'id:uniprot': 'up'}

if __name__ == '__main__':
    opt_parser, opt, args = config_option()

    # KOBAS environment configuration
    kobas2rc = config.getrc()

    # open kobas2db
    kobas2db = dbutils.kobas2db(kobas2rc['kobas2db'])

    if opt.list:
        if opt.species:
            print 'Available databases for %s: ' % opt.species
            databases = kobas2db.get_databases_by_species_abbr(species_abbr = opt.species)
            for line in databases:
                print line
        else:
            if opt.inspecies:
                print 'Available species for gene annotations:'
                species = kobas2db.get_species()
            else:
                print 'Available species for database annotations: '
                species = kobas2db.get_species_by_databases()
            for line in species:
                print '\t'.join(line)
        sys.exit(1)

    if opt.infile:
        args.insert(0,opt.infile)

    if len(args) != 1:
        opt_parser.print_help()
        sys.exit(1)

    if opt.outfile:
        global old_stdout
        old_stdout = sys.stdout
        sys.stdout = open(opt.outfile, 'w')

    # get species name
    species_name = kobas2db.get_species_name_by_species_abbr(species_abbr = opt.species)

    if opt.intype in PROGRAMS.keys():
        # set program for blast
        program = PROGRAMS[opt.intype]
        # verify fasta file
        try:
            f = open(args[0])
            try:
                fasta.verify(f)
            finally:
                f.close()
        except exception.FastaIOError, msg:
            exception.error(msg)
            sys.exit(1)
        # key step
        annotator = annot.Annotator(
            reader = annot.BlastProgReader(
                kobas2rc['blast'], program, kobas2rc['blastdb'] + species_name + '.pep.fasta', args[0], opt.nCPUs),
            selector = annot.BlastoutXMLSelector(kobas2db, opt.species, [opt.rank, opt.evalue, opt.coverage])) # add coverage
    elif opt.intype == 'blastout:xml':
        annotator = annot.Annotator( 
            reader = annot.BlastoutXMLReader(open(args[0])),
            selector = annot.BlastoutXMLSelector(kobas2db, opt.species, [opt.rank, opt.evalue, opt.coverage])) # add coverage
    elif opt.intype == 'blastout:tab': 
        annotator = annot.Annotator(
            reader = annot.BlastoutTabReader(open(args[0])),
            selector = annot.BlastoutTabSelector(kobas2db, opt.species, [opt.rank, opt.evalue]))
    elif opt.intype in DBLINKS.keys():
        db = DBLINKS[opt.intype]
        annotator = annot.Annotator(
            reader = annot.IdMappingReader(kobas2db, args[0], db),
            selector = annot.IdMappingSelector(kobas2db, opt.species, opt.inspecies))
    else:
        sys.exit('%s input is not supported yet, only %s, blastout:xml, blastout:tab, %s' %
            (opt.intype, ', '.join(PROGRAMS.keys()), ', '.join(DBLINKS.keys())))

    items = [item for item in annotator.annotate()]

    # add ortholog
    if (opt.intype in PROGRAMS.keys() and opt.species != 'ko' and opt.ortholog != 'NO'):
        o_species_name = kobas2db.get_species_name_by_species_abbr(species_abbr = opt.ortholog)
        # run BLAST for the species itself
        annotator = annot.Annotator(
            reader = annot.BlastProgReader(
                kobas2rc['blast'], PROGRAMS[opt.intype], kobas2rc['blastdb'] + o_species_name + '.pep.fasta', args[0], opt.nCPUs),
            selector = annot.BlastoutXMLSelector(kobas2db, opt.species, [opt.rank, opt.evalue, opt.coverage])) # add coverage
        # filter with ortholog
        orthologs = {}
        for item in annotator.annotate():
            orthologs[item.query] = item.links
        for item in items:
            if item.links != set() and orthologs[item.query] != set() and kobas2db.is_ortholog(list(orthologs[item.query])[0][0], list(item.links)[0][0], opt.ortholog, opt.species):
                pass
            else:
                items[items.index(item)].links = set()
                
    # Report annotation result
    num_genes_has_annot, num_genes_hasnot_annot = 0, 0
    for item in items:
        if item.has_links():
            num_genes_has_annot += 1
        else:
            num_genes_hasnot_annot += 1

    print '##Species: %s (%s)' % (opt.species, species_name)
    if opt.intype in ('fasta:pro', 'fasta:nuc', 'blastxml', 'blasttab'):
        if opt.species == 'ko':
            print '##Method: BLAST\tOptions: evalue <= %s; rank <= %s' % (str(opt.evalue), str(opt.rank))
        else:
            print '##Method: BLAST\tOptions: evalue <= %s' % str(opt.evalue)
    elif opt.intype in DBLINKS.keys():
        print '##Method: Id mapping'
    print '##Summary:\t%d succeed, %d fail' % (num_genes_has_annot, num_genes_hasnot_annot)
    if opt.species == 'ko':
        print '\n#Query\tKo id|Ko name|hyperlink'
    else:
        print '\n#Query\tGene id|Gene name|hyperlink'

    output.annotate_table(items, opt.species)
    output.annotate_text(items, kobas2db)

