#! /usr/bin/env python
# $URL: https://amberfrog.googlecode.com/svn/trunk/pyensembl/code/util.py $
# $Rev: 118 $

def FASTAseq(inp):
    """Convert FASTA file to sequence data as a string."""

    return ''.join(
      l.strip() for l in inp if not l.startswith('>'))

def asFASTA(s, out, title='some FASTA sequence', linelength=79):
    """Write string as FASTA file."""
    out.write('> %s\n' % title)
    for i in range(0,len(s),linelength):
        out.write(s[i:i+linelength])
        out.write('\n')
