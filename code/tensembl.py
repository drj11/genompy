#! /usr/bin/env python
# $URL: https://amberfrog.googlecode.com/svn/trunk/pyensembl/code/tensembl.py $
# $Rev: 153 $
# tensembl.py

"""Tests for Ensembl."""

# http://docs.python.org/release/2.5.4/lib/module-random.html
import random
# http://docs.python.org/release/2.5.4/lib/module-re.html
import re
# http://docs.python.org/release/2.5.4/lib/module-urllib.html
import urllib

# Local modules
import ensembl
import util

# == Tests using Danio rerio

# Note that ensembl has a sort of REST interface:
# http://www.ensembl.org/Danio_rerio/Export/Output/Location?db=core;flank3_display=0;flank5_display=0;output=fasta;r=25:61-1000000;strand=1;coding=yes;cdna=yes;utr3=yes;peptide=yes;intron=yes;exon=yes;genomic=unmasked;utr5=yes;_format=Text
def fetchREST(species, chromosome, begin, end):
    """Fetch a DNA sequence using the URL based REST interface to the
    Ensembl server.  *species* names the species using its binomial
    name.  *chromosome* specifies the chromosome (as
    a string); *begin* and *end* specifies the base positions.
    Following the usual Ensembl convention the end points *begin* and
    *end* are inclusive.
    """
    url = urlREST(species, chromosome, begin, end)
    return util.FASTAseq(urllib.urlopen(url))

def urlREST(species, chromosome, begin, end):
    """Just the URL for *fetchREST*."""
    template = """http://www.ensembl.org/%(species)s/Export/Output/Location?db=core;output=fasta;r=%(chromosome)s:%(begin)d-%(end)d;strand=1;genomic=unmasked;_format=Text"""

    species = speciesURL(species)

    url = template % locals()
    return url

def speciesURL(n):
    """Canonicalise the binomial name n by capitalising the Genus and
    separating the two words with an underscore (for the URL).
    """
    g,s = re.compile('[ _]+').split(n)
    g = g.capitalize()
    species = '_'.join([g,s])
    return species

def geneREST(species, geneid):
    """Fetch DNA sequence.  *geneid* can be either a stable id (like
    'ENSG00000186442') or a name (like 'KRT3')."""

    template = """http://www.ensembl.org/%(species)s/Export/Output/Gene?db=core;g=%(geneid)s;output=fasta;strand=feature;genomic=unmasked;_format=Text"""

    species = speciesURL(species)
    url = template % locals()

    return util.FASTAseq(urllib.urlopen(url))

def transcriptREST(species, transcriptid):
    """Fetch DNA sequence, given a transcript id."""

    template = """http://www.ensembl.org/%(species)s/Export/Output/Transcript?db=core;output=fasta;strand=feature;t=%(transcriptid)s;param=coding;_format=Text"""

    species = speciesURL(species)
    url = template % locals()

    return util.FASTAseq(urllib.urlopen(url))

def DrerioRange(chromoname, begin, end):
    print "D.rerio %s %d %d" % (chromoname, begin, end)
    b,e = map(int, [begin, end])
    fragrest = fetchREST('Danio rerio', chromoname, b, e)
    fragensembl = str(ensembl.Binomial('danio').chromosome(chromoname)[b:e+1])
    util.asFASTA(fragrest, open('D.rerio%srest.fas' % chromoname, 'w'),
      title='D.rerio chromosome %s %d to %d from REST'%(chromoname, b, e))
    util.asFASTA(fragensembl, open('D.rerio%sensembl.fas' % chromoname, 'w'),
      title='D.rerio chromosome %s %d to %d from Ensembl'%(chromoname, b, e))
    assert fragrest == fragensembl

def DrerioA():
    """Test a particular segment."""

    return DrerioRange(25, 1e6, 2e6)

def DrerioRandom():
    """Test a random segment."""

    import math
    import random

    danio = ensembl.Binomial('Danio')
    chromosome = random.choice(danio.chromosomes())

    # log of sequence length follows a uniform distribution
    x = math.log10(chromosome.length)
    p,q = sorted(int(10**(x*random.random())) for _ in range(2))
    # Clamp segment length to 1e6 (because the REST interface rejects
    # large requests (and it take a long time to transfer)).
    if chromosome.length > 10**6:
        q = p+int(10**(6*random.random()))
    return DrerioRange(chromosome.name, p, q)

# == Tests using Arabidopsis thaliana

def chloroplastFTPArabidopsis():
    """Get chloroplast sequence straight from FTP site."""

    import urllib

    return util.FASTAseq(urllib.urlopen(
      "ftp://ftp.arabidopsis.org"
      "/home/tair/Sequences/whole_chromosomes/TAIR10_chrC.fas"))

def Athaliana():
    print "A.thaliana"
    s = ensembl.Server(host='mysql.ebi.ac.uk', port=4157)
    db = s.database('Arabidopsis thaliana')
    chloroplast_seq = str(db.chromosome('Pt')[:])
    assert chloroplast_seq == chloroplastFTPArabidopsis()

# == Tests using Mus musculus

def musMT():
    """Mus mitochondrial DNA."""

    print "musMT"
    mus = ensembl.Binomial('mus')
    str(mus.chromosome('MT')[:])

def musY():
    """Mus Y chromosome."""

    print "musY"
    mus = ensembl.Binomial('Mus')
    str(mus.chromosome('Y')[1:2000])
    # Includes part of a gap.
    str(mus.chromosome('Y')[194000:195000])
    # Overlaps two contigs.
    str(mus.chromosome('Y')[2665000:2666000])

def homo20():
    """Homo chromosome 20."""

    # This used to select the wrong database, so worth keeping as a
    # test:
    print "Homo 20"
    homo = ensembl.Binomial('homo')
    c20 = homo.chromosome('20')
    # See Issue 2
    str(c20[1596200:1596400])

# == Miscellaneous tests

def musChid():
    """(Somewhat internal test) Check chromosome coord system id."""

    import sys
    sys.stdout.write("Mus Chromosome Coord System")

    mus = ensembl.database('mus_musculus_core_56_37i')
    coordsys = mus.chromosome_coord_system_id()
    sys.stdout.write(" %d\n" % coordsys)
    assert 1 == coordsys

def Hseq():
    """(used to) provoke a pymysql bug (that was originally detected
    by chromosome 20.
    """

    dbid = 18559
    print "Homo assembly 37a seq_region_id %d" % dbid
    # Note: We require a particular release and assembly.
    homo = ensembl.database('homo_sapiens_core_56_37a')
    cursor=homo.conn.cursor()
    cursor.execute('''select sequence from dna where
      seq_region_id = %s;''', dbid)
    result, = list(cursor)
    result, = result
    assert 37972 == len(result)

def revN():
    """Includes a reversed contig that has an 'N' nucleotide.  Which
    caused problems when we didn't complement this properly.
    """

    print "revN"

    drerio = ensembl.Binomial('Danio').chromosome('25')[1000000:2000000]

def gene99889():
    """Gene from the Ensembl tutorial."""

    print "gene 99889"

    homo = ensembl.Binomial('Homo')
    gid = 'ENSG00000099889'
    gene = homo.fetch_gene_stable_id(gid)
    dbseq = str(gene[:])
    restseq = geneREST('Homo sapiens', gid)
    assert dbseq == restseq

def tut1():
    """Examples from the Ensembl Perl API tutorial."""
    # http://www.ensembl.org/info/docs/api/core/core_tutorial.html
    ensembl.Binomial('Homo').fetch_region('chromosome', 'X')
    ensembl.Binomial('Homo').fetch_region('clone', 'AL359765.6')
    # Not stable, doesn't work with this release.
    # ensembl.Binomial('Homo').fetch_region('supercontig', 'NT_011333')

def catalogue():
    """Basic test of catalogue feature."""
    import sys
    sys.stdout.write("catalogue ")
    sys.stdout.flush()
    s = ensembl.Server()
    p = len(s.catalogue())
    sys.stdout.write("%d " % p)
    sys.stdout.flush()
    # See http://plants.ensembl.org/index.html
    g = ensembl.Server(host='mysql.ebi.ac.uk', port=4157)
    q = len(g.catalogue())
    print "+", q, "species"

def KRT3():
    """A particular gene."""

    print "Homo KRT3"

    homo = ensembl.Binomial('Homo')
    gname = 'KRT3'
    gene = homo.fetch_gene_name(gname)
    dbseq = str(gene[:])
    restseq = geneREST('Homo sapiens', gname)
    assert dbseq == restseq
    return gene

def KRT3transCCDS():
    """One transcript of the KRT3 gene has a CCDS xref."""

    gname = 'KRT3'
    tname='ENST00000417996'
    print "Homo %s %s" % (gname, tname)
    homo = ensembl.Binomial('homo sapiens')
    gene = homo.fetch_gene_name(gname)
    ts = gene.Transcripts()[tname]
    xrefs = ensembl.externals(ts)
    assert 'CCDS' in xrefs

def BRCA2trans():
    """BRCA2 Transcripts."""

    print "BRCA2 Transcripts"""

    homo = ensembl.database('homo_sapiens_core_62_37g')
    ts = homo.fetch_gene_name('BRCA2').Transcripts()
    assert 6 == len(ts)
    return ts

def BRCA2_380152():
    """One particular BRCA2 transcript."""

    tname = 'ENST00000380152'
    print "BRCA2", tname

    homo = ensembl.Binomial('Homo')
    ts = homo.fetch_gene_name('BRCA2').Transcripts()
    dbseq = ts[tname].translated_seq()
    restseq = transcriptREST('Homo sapiens', tname)
    assert dbseq == restseq

def LRRN2_367175():
    """One particular LRRN2 transcript.  Of note because it starts and
    ends in the same exon."""

    tname = 'ENST00000367175'
    print "LRRN2", tname

    homo = ensembl.Binomial('Homo')
    ts = homo.fetch_gene_name('LRRN2').Transcripts()
    dbseq = ts[tname].translated_seq()
    restseq = transcriptREST('Homo sapiens', tname)
    assert dbseq == restseq


def main():
    KRT3()
    KRT3transCCDS()
    DrerioRandom()
    LRRN2_367175()
    catalogue()
    Hseq()
    gene99889()
    musChid()
    musMT()
    musY()
    homo20()
    revN()
    BRCA2trans()
    BRCA2_380152()
    DrerioA()
    Athaliana()
    tut1()

if __name__ == '__main__':
    main()
