#! /usr/bin/env python
# $URL: https://amberfrog.googlecode.com/svn/trunk/pyensembl/code/gene.py $
# $Rev: 130 $

"""gene.py Genus [species] stable-id

A command line "gene browser".  Prints out DNA sequence.
"""

import ensembl

def geneout(out, species, geneid):
    db = ensembl.Binomial(species)
    gene = db.fetch_gene_stable_id(geneid)
    print >>out, repr(gene)
    print >>out, gene

def main(argv=None):
    import sys
    if argv is None:
        argv = sys.argv

    arg = argv[1:]
    if len(arg) == 3:
        name = "%s %s" % (arg[:2])
        geneid = arg[2]
    elif len(arg) == 2:
        name = arg[0]
        geneid = arg[1]

    geneout(sys.stdout, name, geneid)


if __name__ == '__main__':
    main()
