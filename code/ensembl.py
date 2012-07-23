#!/usr/bin/env python
# $URL: https://amberfrog.googlecode.com/svn/trunk/pyensembl/code/ensembl.py $
# $Rev: 156 $
# An interface to the Ensembl database
#
# [EBI 2006-11-13] "The Ensembl Database Schema";
# http://medblog.stanford.edu/lane-faq/archives/EnsSchema.pdf
#
# :detuple:magic: This trick is useful enough to be referred to multiple
# times.  If a (SQL) query result has only one column (for example,
# "select o.ensembl_id from ..." then each row in the result is a
# 1-tuple.  Often instead of a list of 1-tuples, you want a list of
# whatever the 1-tuples contain.  If l is a list of 1-tuples then the
# following expression detuplizes it: "l, = zip(*l)".
#
# Facts:
#
# (By inspection)
# Each organism gets a different mysql database.  In fact, each assembly
# (I think that's what they're called) gets a different database.  EG
# "mus_musculus_core_56_37i".  Example mysql command line:
# mysql -h ensembldb.ensembl.org -P 5306 -u anonymous mus_musculus_core_56_37i
# (later I found: http://metazoa.ensembl.org/info/docs/api/versions.html )
#
# (from [EBI 2006-11-13]) Sequence data ([GACT]*) is stored in the "dna"
# table; the "seq_region_id" will specify a "coord_system_id" that has
# "sequence_level" in its "attrib" column.
#
# (By inspection, and [EBI 2006-11-13]) the coord_system_id with name
# "chromosome" and attrib "default_version" identifies the chromosomes.
# Assemblies with that coord_system_id correspond to chromosomes.  This
# will print them out:
#   select s.seq_region_id,s.name,s.length from seq_region s join
#     coord_system c on s.coord_system_id = c.coord_system_id where c.attrib
#     like '%default_version%' and c.name = 'chromosome';
# or, if you know the coord_system_id for chromosome:
#   select * from seq_region where coord_system_id = 1;
#
# (By inspection)
# The sequences that together comprise an assembly can be got with:
#   select * from assembly a join seq_region s on a.cmp_seq_region_id =
#   s.seq_region_id join coord_system c on s.coord_system_id =
#   c.coord_system_id where a.asm_seq_region_id = 3 and c.attrib like
#   '%sequence_level%';
# (This gets seq_region_id 3 which in the mus_musculus_core_56_37i database
# is the "Y" chromosome, see above for list of chromosomes)
# The following is similar, but only extracting a critical set of fields
# and presenting in assembly order:
#   select
#   a.cmp_seq_region_id,a.asm_start,a.asm_end,a.cmp_start,a.cmp_end,a.ori
#   from assembly a join seq_region s on a.cmp_seq_region_id =
#   s.seq_region_id join coord_system c on s.coord_system_id =
#   c.coord_system_id where a.asm_seq_region_id = 3 and c.attrib like
#   '%sequence_level%' order by a.asm_start;
#
# (By inspection, and [EBI 2006-11-13]) the following almost gets you
# the gene known as KRT1:
#   select * from gene g join object_xref o on o.ensembl_id = g.gene_id join
#   xref x on o.xref_id = x.xref_id where x.display_label='krt1' and
#   o.ensembl_object_type = 'Gene';
# (see also the source to method "_type_by_external_id" of perl module
# "Bio::EnsEMBL::DBSQL::DBEntryAdaptor" in the Ensembl perl software
# release.)
# The above query returns multiple rows with the same o.ensembl_id (=
# g.gene_id when joining on the "gene" tale); bung them in a set().
# The perl code also performs a second separate query, joining "gene",
# "object_xref", and "external_synonym".
#
# CCDS.
#   select * from object_xref join xref on object_xref.xref_id =
#   xref.xref_id where ensembl_id = 189216;
# xref.display_label starts "CCDS".
# Bug: Note [EBI 2006-11-13] has "dbname" instead of "db_name" as a
# column in table "external_db".
# external_db.db_name = "CCDS"


# requires drj's locally hacked version of pymysql for which quoting
# works.
import pymysql
# Not (yet) documented
if 0:
  import pymysql.connections
  pymysql.connections.DEBUG = True

# http://docs.python.org/release/2.5.4/lib/module-itertools.html
import itertools
# http://docs.python.org/release/2.5.4/lib/module-re.html
import re

# Local module
import util

class Error(Exception):
    pass

class IntegrityError(Error):
    """An error in the expected format / content of database tables or
    similar."""

class Server:
    """Models an internet server using the ensembl protocol.  The
    archetypal server is the one at enembldb.ensembl.org:5306.
    """

    host = "ensembldb.ensembl.org"
    port = 5306

    def __init__(self, **keys):
        if 'host' in keys:
            self.host = keys['host']
        if 'port' in keys:
            self.port = keys['port']

    def catalogue(self, versions=False):
        """Retrieve a list of species stored on the server.  If
        *versions* is true then the list will include the Ensembl
        software and genome versions (and therefore have many entries
        per species, and be much longer).  Normally version is false in
        which case each species appears only once."""

        cursor = self.cursor()
        cursor.execute(
           """select schema_name from information_schema.schemata
           where schema_name like %s;""", "%_core%")
        result = list(cursor)
        cursor.close()
        # :detuple:magic
        result, = zip(*result)
        if not versions:
            result = set(
              re.match(r'^((?:.*?)_(?:.*?))_', x).group(1) for x in result)
        return sorted(result)

    def databaseList(self, name):
        """List of database names."""

        return BinomialToEnsemblDBNameList(name, self)
    
    def database(self, name):
        """Return a database instance corresponding to the species with
        binomial name *name*.
        """

        return Binomial(name, self)

    def cursor(self):
        """Return a database cursor connected to the server.  Mostly an
        internal convenience."""

        conn = pymysql.connect(user='anonymous', passwd='',
          host=self.host, port=self.port)
        return conn.cursor()


def BinomialToEnsemblDBNameList(name, server=None):
    """Return a list of database names for the various assemblies of the
    species named by *name*.  This function is an internal helper for the
    :meth:`Binomial` function, exposed in case you want to use it.  See
    :meth:`Binomial` for documentation.

    When *name* specifies a genus for which there are more than one
    species, then None is returned.

    *server* should be an instance of the Server class that specifies
    a host and port; if it's None, the default, then a default
    server is used.
    """

    if server is None:
        server = Server()

    name = name.lower()
    name = name.split()
    if len(name) == 1:
        name.append('%')
    conn = pymysql.connect(user='anonymous', passwd='',
      host=server.host, port=server.port)
    cursor = conn.cursor()
    cursor.execute(
      '''select schema_name from information_schema.schemata
      where schema_name like %s;''', '\\_'.join(name) + '\\_core\\_%')
    result = list(cursor)
    cursor.close()
    # Each row from the query is a 1-tuple containing the database name;
    # :detuple:magic
    result, = zip(*result)
    # Check that there is only a single species.
    species = set(map(lambda n: n.split('_')[1], result))
    if len(species) != 1:
        return None
    def key(dbname):
        """Helper that extracts the sorting key from the database name
        so that we can sort them properly."""

        import re

        m = re.match(r'.*_(\d+)_(\d+)([a-z]?)$', dbname)
        g = m.groups()
        return (int(g[0]), int(g[1]), g[2])
    return sorted(result, key=key)

def Binomial(name, server=None):
    """Return a database for the most recent assembly of the species.

    *name* should be two words: 'Genus species'; or, one word: 'Genus'.
    Using just one word (genus only) works only if there is exactly
    one species in that genus in the Ensembl database (which generally
    seems to be the case).

    Example: Binomial('Homo sapiens')
    Binomial('mus') # for lazy mouse geneticists
    """

    return database(database=BinomialToEnsemblDBNameList(name, server)[-1],
      server=server)

class database:
    def __init__(self, database, server=None):
        """*database* should specify the name of an Ensembl database
        (usually containing a complete versioned organism assembly).

	*server* is an instance of the Server class that specifies
	the host and port.  If it's None, the default, a default
	server is used.

        Example: database('mus_musculus_core_56_37i')
        """

        if server is None:
            server = Server()

        self.database = database
        self.host = server.host
        self.port = server.port
        self.conn = pymysql.connect(user='anonymous', passwd='', db=database,
          host=self.host, port=self.port)
        # Actually not very interesting.  But we collect it anyway.
        self.server_info = self.conn.get_server_info()
        # Per-instance coord system cache.
        self._coord_system_id_cache = {}

    def __repr__(self):
        """'Friendly' printable representation."""

        return (
          'ensembl.database(database=%r, server=Server(host=%r, port=%d))' %
          (self.database, self.host, self.port))

    # In each instance this caches the result of _coord_system_id.
    _coord_system_id_cache = None
    def _coord_system_id(self, name):
        """Internal.  Determine the coord_system_id for the coord system
        *name*.  Results of database queries are cached.  *name* is
        typically 'chromosome', 'supercontig', 'contig', 'clone'.
        """

        if name in self._coord_system_id_cache:
          return self._coord_system_id_cache[name]

        cursor = self.conn.cursor()
        cursor.execute(
          '''select coord_system_id from coord_system where
            name = %s and attrib like '%%default_version%%';''', name)
        result = list(cursor)
        cursor.close()
        if len(result) == 0:
            # If *name* is valid, then coord_system table could be bogus.
            raise IntegrityError(
              'no coord_system_id for %r found' % name)
        if len(result) > 1:
            raise IntegrityError(
              'multiple coord_system_id for %r found (inspect coord_system table?)' % name)
        id = result[0][0]
        self._coord_system_id_cache[name] = id
        return id

    def version(self):
        """Ensembl version."""

        # See release notes:
        # http://www.ensembl.org/info/website/news.html?id=63

        cursor = self.conn.cursor()
        cursor.execute("""
          select meta_value from meta where meta_key='schema_version';""")
        result, = list(cursor)
        cursor.close()
        version = result[0]
        return version

    def chromosome_coord_system_id(self):
	"""Internal.  Return the coord_system_id that corresponds
	to a chromosome.
        """

        return self._coord_system_id('chromosome')

    def sequence_level_coord_system_id(self):
        """Internal.  Determine and cache the coord_system_id that
        corresponds to sequence level data.  Return it.
        """

        cursor = self.conn.cursor()
        cursor.execute(
          '''select coord_system_id from coord_system where
               attrib like '%sequence_level%';''')
        result = list(cursor)
        cursor.close()
        if len(result) == 0:
            raise IntegrityError('no coord_system_id for sequence_level found (inspect coord_system table?)')
        if len(result) > 1:
            raise IntegrityError('multiple coord_system_id for sequence_level found (inspect coord_system table?)')
        def cached(it=result[0][0]):
            return it
        self.sequence_level_coord_system_id = cached
        return cached()

    def fetch_region(self, coord_system, name):
        """Return a (sequence for a) generic genome region given its
        name.  *coord_system* specifies the name of the coord system
	(typically, 'chromosome', 'clone', and so on).  The *name*
	corresponds to the Ensembl name found in the seq_region
	table.
        """

        name = str(name)

        coord_system_id = self._coord_system_id(coord_system)
        cursor = self.conn.cursor()
        cursor.execute(
          '''select seq_region_id,length from seq_region where
          coord_system_id = %s and name = %s;''',
          (coord_system_id, name))
        result = list(cursor)
        cursor.close()
        if len(result) == 0:
            raise Error('no %r named %r found' % (coord_system, name))
        if len(result) > 1:
            raise IntegrityError('more than one %r named %r found' %
              (coord_system, name))
        seq_region_id,length = result[0]
        return SequenceById(self, seq_region_id)

    def fetch_gene_stable_id(self, stable_id):
        """Return a Gene object for the gene with stable id
        *stable_id*."""

        # Note: schema entirely reversely engineered.  No documentation.
        # Note: schema changed for ensembl release 67 or earlier
        # (stable_id has moved into the 'gene' table).
        cursor = self.conn.cursor()
        cursor.execute(
          '''select gene_id from gene where stable_id =
          %s;''', stable_id)
        result = list(cursor)
        cursor.close()
        if len(result) == 0:
            raise Error('no gene with stable id %r found' % (stable_id))
        if len(result) > 1:
            raise IntegrityError('more than one gene with stable id %r found' %
              stable_id)
        geneid, = result[0]
        return Gene(self, geneid)

    def fetch_gene_name(self, name):
        """Return a Gene object for the gene named *name*.

        (see :meth:`gene_id_from_name` for technicalities.)
        """

        ids = self.gene_id_from_name(name)
        if not ids:
            raise Error("No gene named %r." % name)
        if len(ids) > 1:
            raise Error("Too many (%d) genes named %r." % (len(ids), name))
        geneid, = ids
        return Gene(self, geneid)

    def fetch_transcripts_gene_id(self, geneid):
        """Return transcripts as a dictionary keyed by ensembl
        transcript stable ids (generally, these ids start "ENST"); each
        value is a Transcript instance.
        """

        # Slide 18 of [EBI 2006-11-13] very useful.

        # As of release 67 (quite probably, way earlier)
        # transcript_stable_id is in the transcript table.

        cursor = self.conn.cursor()
        if 0: cursor.execute('''select
            transcript.transcript_id,
            transcript_stable_id.stable_id,
            exon.exon_id,
            exon.seq_region_id,
              exon.seq_region_start,exon.seq_region_end,
              exon.seq_region_strand
          from exon join exon_transcript join gene join transcript
            left join transcript_stable_id
            on transcript.transcript_id = transcript_stable_id.transcript_id
          where
            transcript.transcript_id = exon_transcript.transcript_id and
            gene.gene_id = transcript.gene_id and
            exon_transcript.exon_id = exon.exon_id and
            gene.gene_id = %s
          order by transcript_id,seq_region_start;''', geneid)
        cursor.execute('''select
            transcript.transcript_id,
            transcript.stable_id,
            exon.exon_id,
            exon.seq_region_id,
              exon.seq_region_start,exon.seq_region_end,
              exon.seq_region_strand
          from exon join exon_transcript join gene join transcript
          where
            transcript.transcript_id = exon_transcript.transcript_id and
            gene.gene_id = transcript.gene_id and
            exon_transcript.exon_id = exon.exon_id and
            gene.gene_id = %s
          order by transcript_id,seq_region_start;''', geneid)
        def transcript_id(x): 
            """Extract the transcript id and its stable id, from a
            result row.
            """
            return x[:2]

        result = dict((tid[1],Transcript(database=self,
            dbid=tid[0], stable_id=tid[1],
            flat_exons=[row[-5:] for row in rows]))
          for tid,rows in itertools.groupby(cursor, transcript_id))
        cursor.close()
        return result
    
    def chromosome(self, name):
        """Return a (sequence for a) chromosome given its name.  The
        *name* corresponds to the Ensembl name found in the seq_region
        table."""

        name = str(name)

        return self.fetch_region('chromosome', name)

    def chromosomes(self):
        """Return a list of chromosomes.  Each chromosome is returned as
        a ('name', seq_region_id, length) triple; the *name* can be
        passed to the *chromosome* method.
        """
        
        cursor = self.conn.cursor()
        cursor.execute(
          '''select name, seq_region_id, length from seq_region
          where coord_system_id = %s;''',
          self.chromosome_coord_system_id())
        result = list(cursor)
        cursor.close()
        return [Chromosome(*x) for x in result]

    def gene_id_from_name(self, name):
        """Return a collection (usually a *set*) of all the Ensembl IDs
        for genes matching *name*.
        
        For those interested in the schema technicalities:
        the *xref.display_label* and *external_synonym.synonym*
        columns are searched.  Actually the synonym column isn't
        searched yet.
        """

        # Similar queries can be made for "transcript" and
        # "translation"; the Ensembl perl code exploits this.

        cursor = self.conn.cursor()
        cursor.execute(
          '''select o.ensembl_id from gene g join object_xref o on
          o.ensembl_id = g.gene_id join xref x on o.xref_id = x.xref_id where
          x.display_label = %s and o.ensembl_object_type = 'Gene' and
          g.is_current = 1;''', name)
        result = list(cursor)
        cursor.close()
        # :detuple:magic
        result, = zip(*result)
        return set(result)

    def fetch_dna(self, seq_region_id):
        """Primitive for returning DNA.  Returns the base sequence for
        the specified seq_region_id, which must be in the *dna* table.
        """

        seq_region_id = int(seq_region_id)

        cursor = self.conn.cursor()
        cursor.execute(
          '''select sequence from dna where seq_region_id = %s''',
          seq_region_id)
        result = list(cursor)
        cursor.close()
        if len(result) > 1:
            raise IntegrityError(
              'Multiple entries found for seq_region_id %d in dna table.' %
                seq_region_id)
        if len(result) == 0:
            return None
        return result[0][0]

def _assembly_find_contig(assembly, start):
    """Find the assembly record from which to start assembling
    for sequence data at index *start*.  *assembly* should be a
    list of tuples as returned by
    :meth:`sequence_level_assembly` (and in particular, should
    be in sorted order).
    """

    # Sadly, neither index nor bisect take a "key" argument.  So we
    # do a linear search here. :todo: improve. :todo: consider doing
    # database query(!).
    for i, c in enumerate(assembly):
        if c[1] >= start:
            return i
    return None

class Chromosome:
    def __init__(self, name, seq_region_id, length):
        self.name = name
        self.seq_region_id = seq_region_id
        self.length = length

    def __repr__(self):
        return "Chromosome(name=%(name)r, length=%(length)d)" % self.__dict__

class SequenceById:
    def __init__(self, database, seq_region_id, forward=True):
        self.database = database
        self.seq_region_id = seq_region_id
        self.debug = False
        self.forward = forward

    def Debug(self, debug=__import__('sys').stderr):
        self.debug = debug
        return self

    def sequence_level_assembly(self):
        """Returns a list of contig parts that assemble together to make
        the sequence; sorted.
        """
        cursor = self.database.conn.cursor()
        cursor.execute('''select
          a.asm_start,a.asm_end,a.cmp_seq_region_id,a.cmp_start,a.cmp_end,a.ori
          from assembly a join seq_region s on a.cmp_seq_region_id =
          s.seq_region_id where s.coord_system_id = %s and
          a.asm_seq_region_id = %s order by a.asm_start;''',
          (self.database.sequence_level_coord_system_id(),
          self.seq_region_id))
        ol = tuple(cursor)
        cursor.close()
        def cached(it=ol):
            return it
        self.sequence_level_assembly = cached
        return cached()

    def __repr__(self):
        return 'SequenceById(%s, %d, forward=%r)' % (
          self.database, self.seq_region_id, self.forward)

    def reversedComplement(self):
        """Return a Sequence that complements this one.  The sequences
        share the same coordinate system.  If *s* is this sequence and
        *r* = s.reverseComplement(), then the Slice r[x:y] is the
        reverse complement of the Slice s[x:y].
        """

        return SequenceById(self.database, self.seq_region_id, 
          forward=not self.forward)

    def __getitem__(self, key):
        # When a Sequence is indexed, it retrieves a (DNA) Slice according
        # to the Ensembl numbering, so the start of a sequence is generally
        # at index 1.

        # General warning: Python and drj use base,limit for ranges
        # where limit is not included (and generally start at 0);
        # ensembl generally use start,end for ranges and include end
        # (and generally start at 1).
        assembly = self.sequence_level_assembly()
        # The bounds of the entire assembly (in the chromosome
        # coordinate system).
        asmbase = assembly[0][0]
        asmlimit = assembly[-1][1] + 1
        if self.debug:
            print >> self.debug, "contigs", len(assembly), \
              "base", asmbase, "limit", asmlimit
        if type(key) != slice:
            try:
                good = asmbase <= key < asmlimit
            except:
                raise TypeError('integer required')
            if not good:
                raise IndexError('index out of range')
            key = slice(key, key+1)
        if key.step not in (None, 1):
            raise ValueError('slice step must be 1')
        key = _clamp_slice(key, asmbase, asmlimit)
        keybase = key.start
        keylimit = key.stop

        if keybase >= keylimit:
            # An empty slice.
            return Slice(chromosome=repr(self))
        return Slice(chromosome=self, base=keybase, limit=keylimit,
          strand=['reverse','forward'][bool(self.forward)])

    def asFASTA(self, out, linelength=79):
        """To the file object *out*, write the seuqnce in FASTA format.
        Lines will be limited to *linelength* characters (each line will
        have an additional newline).
        """

        return util.asFASTA(self[:], out,
          title=repr(self), linelength=linelength)

class Slice(object):
    """A data model for a slice of sequence data."""

    # Lazily loaded.
    _seq = None

    debug = None

    def __init__(self, **keys):
        """It is conventional to supply a *chromosome* key whose value
        is the string representation of the chromosome on which this
        Slice lies.
        """
        self.__dict__.update(keys)

    def __len__(self):
        return self.limit - self.base

    def __repr__(self):
        if not hasattr(self, 'base'):
            return "Slice()"
        return "%(chromosome)r[%(base)d:%(limit)d]" % self.__dict__

    def __str__(self):
        """(When converted to a string) return the sequence data."""
        return self.seq()

    def genes(self):
        """List of genes."""

        seq_region = self.chromosome
        cursor = self.chromosome.database.conn.cursor()
        cursor.execute("""select gene_id from gene
          where
            seq_region_id = %s and
            seq_region_start >= %s and
            seq_region_end < %s;""",
          (self.chromosome.seq_region_id, self.base, self.limit))
        result = [Gene(self.chromosome.database, dbid) for dbid, in cursor]
        return result

    def seq(self):
        """The sequence data."""

        if self._seq is not None:
            return self._seq

        result = ''

        # We could probably get all of the sequences in one query:
        # select blah from sequence join assembly where sequence.id ==
        # blsh and assembly start < END and assembly.end > START.

        # General warning: Python and drj use base,limit for ranges
        # where limit is not included (and generally start at 0);
        # ensembl generally use start,end for ranges and include end
        # (and generally start at 1).

        assembly = self.chromosome.sequence_level_assembly()

        base = self.base
        # We need to find which contig to start from, then iterate over
        # the contigs, starting there.
        i = _assembly_find_contig(assembly, base)
        for i,acontig in list(enumerate(assembly))[i:]:
            (asm_start, asm_end,
             cmp_seq_region_id, cmp_start, cmp_end, ori) = acontig
            if self.debug:
                print >> self.debug, "asm_start", asm_start, \
                  "asm_end", asm_end, "seq_region_id", cmp_seq_region_id, \
                  "cmp_start", cmp_start, "cmp_end", cmp_end, "ori", ori
                  
            if asm_start >= self.limit:
                break
            # A gap in the assembly?
            if base < asm_start:
                n = asm_start - base
                # Sanity check for very long gaps.
                if n > 9999999:
                    assert 0, "huge n: %d" % n
                # :todo: gaps represented by 'N'; make parameter?
                result += 'N'*n
                base += n
            assert asm_start <= base
            contig_sequence = self.chromosome.database.fetch_dna(
              cmp_seq_region_id)

            # Indexes (into the contig's DNA sequence) use the Ensembl
            # convention of starting at 1.  Here we convert to Python
            # convention.
            extract = contig_sequence[cmp_start-1:cmp_end]
            if ori == -1:
                extract = _complement(extract[::-1])
            # Slightly paranoid: check that the extract is as long as we
            # expect it to be (in other words, cmp_start and cmp_end are
            # within range).
            if self.debug:
                print >> self.debug, len(extract), cmp_end, cmp_start, cmp_end-cmp_start
            assert len(extract) == cmp_end - cmp_start + 1
            # chop off right-hand end
            if self.limit <= asm_end:
                d = self.limit - asm_start # note: negative
                extract = extract[:d]
            # chop off left-hand end
            if base > asm_start:
                d = base - asm_start
                extract = extract[d:]
            result += extract
            base += len(extract)
        # Should we check for a gap at the right-hand end?
        if self.strand == 'reverse':
            result = _complement(result[::-1])
        self._seq = result
        return result

class SubSequence(object):
    def __init__(self, sequence, start, end, ori=1):
        """A sequence that is a subsequence of another."""
        self.sequence = sequence
        self.start = start
        self.end = end
        if ori not in (-1, 1):
            raise Error('ori must be 1 or -1.')
        self.ori = ori

    def __getitem__(self, key):
        if type(key) == slice:
            key = _clamp_slice(key, 1, self.end - self.start + 2)
            if self.ori == 1:
                start = key.start + self.start - 1
                stop = key.stop + self.start - 1
                return self.sequence[start:stop:key.step]
            else:
                # When the underlying contig is stored reversed, then the
                # sequence will be reversed twice.  :todo: we can avoid
                # that by implementing a protocol saying whether we want
                # the reversed strand or not.
                start = self.end - key.stop + 2
                stop = self.end - key.start + 2
                result = self.sequence.reversedComplement()[start:stop:key.step]
                return result
        # :todo: this may not work
        # key is a integer, we assume
        if self.ori == 1:
            key = key + self.start - 1
            return self.sequence[key]
        else:
            key = self.end - key + 1
            return _complement(self.sequence[key])

class Gene(SubSequence):
    def __init__(self, database, dbid):
        self.database = database
        self.dbid = dbid
        cursor = database.conn.cursor()
        cursor.execute('''select
          seq_region_id,seq_region_start,seq_region_end,seq_region_strand
          from gene where gene_id = %s;''', dbid)
        result, = list(cursor)
        seq_id,start,end,ori = result
        seq = SequenceById(database, seq_id)
        super(Gene, self).__init__(seq, start, end, ori)

    def __repr__(self):
        return "Gene(%r, %d)" % (self.database, self.dbid)

    def Transcripts(self):
        return self.database.fetch_transcripts_gene_id(self.dbid)

class Exon(SubSequence):
    def __init__(self, database, dbid,
      seq_region_id, seq_region_start, seq_region_end, seq_region_strand):
        self.database = database
        self.dbid = dbid
        seq = SequenceById(database, seq_region_id)
        super(Exon, self).__init__(seq, seq_region_start,
          seq_region_end, seq_region_strand)

    def __repr__(self):
        return "Exon(%r, %d)" % (self.database, self.dbid)

class Transcript:
    def __init__(self, database, dbid, stable_id, flat_exons):
        """*flat_exons* is a list of exon descriptions.  Each exon is
        described by a 5-tuple of (exon_id, seq_region_id,
        seq_region_start, seq_region_end, seq_region_strand).
        """

        self.database = database
        self.dbid = dbid
        self.stable_id = stable_id
        self.flat_exons = flat_exons

    def __repr__(self):
       return ("Transcript(database=%(database)r, stable_id=%(stable_id)r)" %
         self.__dict__)
    
    def seq(self):
        """Return the transcript (DNA) sequence; the exons joined
        together.  Note: includes untranslated regions (at each end)."""

        result = ''
        # Assume flat_exons are in sorted order (which the original DB
        # query arranges).
        for ex in self.flat_exons:
            d = dict(zip(['dbid', 'seq_region_id', 'seq_region_start',
              'seq_region_end', 'seq_region_strand'],
              ex))
            exon = Exon(database=self.database, **d)
            result += exon[:].seq()
        return result

    def translated_seq(self):
        """Return the DNA sequence corresponding to the translated
	sequence.  This is a concatenation of the exon sequences,
	but with the untranslated regions removed.  Note that it
	is mRNA that is translated, but this function returns the
	corresponding DNA.
        """

        # (By inspection of the ENST00000380152 Transcript, I deduce)
        # the translation discards an initial portion and a terminal
        # portion.  (considering the *translation* database table) The
        # first base is base seq_start of start_exon_id (counting from
        # 1, in the usual ensembl convention); the last base if base
        # seq_end of end_exon_id (again counting from 1).  Refer to 
        # [EBI 2006-11-13] Slide 18.

        cursor = self.database.conn.cursor()
        cursor.execute('''select
            start_exon_id, seq_start, end_exon_id, seq_end
          from translation
          where transcript_id = %s;''', self.dbid)
        result, = list(cursor)
        start_exon_id, seq_start, end_exon_id, seq_end = result

        # *result* is None until we find the start exon.
        result = None
        for ex in self.flat_exons:
            d = dict(zip(['dbid', 'seq_region_id', 'seq_region_start',
              'seq_region_end', 'seq_region_strand'],
              ex))
            # Slightly special case for single exon transcript.  :todo:
            # could in fact pull this out of the loop.
            if d['dbid'] == start_exon_id == end_exon_id:
                exon = Exon(database=self.database, **d)
                result = exon[:].seq()
                result = result[:seq_end]
                result = result[seq_start-1:]
                return result
            if d['dbid'] == end_exon_id:
                # Must have already processed start exon.
                assert result is not None
                exon = Exon(database=self.database, **d)
                seq = exon[:].seq()
                seq = seq[:seq_end]
                result += seq
                return result
            if result is not None:
                # Interior exon.
                exon = Exon(database=self.database, **d)
                result += exon[:].seq()
            if d['dbid'] == start_exon_id:
                exon = Exon(database=self.database, **d)
                result = exon[:].seq()
                result = result[seq_start-1:]
        assert 0, "unreachable"

class External:
    def __init__(self, **k):
        self.__dict__.update(k)

    def __repr__(self):
        return """External(%r)""" % self.__dict__

def externals(thing):
    """*thing* is an object with a *database* and a *dbid* member.
    Returns a dict of External() objects."""

    database = thing.database
    dbid = thing.dbid
    cursor = database.conn.cursor()
    cursor.execute("""select
        external_db.db_name,xref.dbprimary_acc,xref.display_label
      from object_xref join xref on object_xref.xref_id = xref.xref_id
      join external_db on external_db.external_db_id = xref.external_db_id
      where ensembl_id = %s;""", dbid)
    result = list(cursor)
    cursor.close()
    return dict(
      (row[0],External(db_name=row[0],
        dbprimary_acc=row[1], display_label=row[2])) for row in result)


def _complement(s):
    """Return the complement of a DNA sequence in string form."""
    import re

    # See http://en.wikipedia.org/wiki/Nucleic_acid_notation for a
    # complete list.

    return re.sub('.',
        lambda m: dict(G='C', C='G', A='T', T='A', N='N')[m.group(0)], s)

def _clamp_slice(s, base, limit):
    """Convert (Python) slice *s* so that *s.start* and *s.stop* are within the
    bounds *base*,*limit*.  Return a fresh slice object.
    
    This is a bit more fiddly than it needs to be because [:] gives
    `slice(0,2147483647)` (which is madness), but [::] gives
    `slice(None,None)` (which is quite sensible).
    """

    b = s.start or base
    l = min(limit, s.stop or limit)
    return slice(b, l, s.step)

def debug():
    """Stuff previously found to be useful for debugging."""

    # Really a test of the improved pymysql than this module.
    if 0:
        homo=Binomial('homo')
        cursor=homo.conn.cursor()
        # This result happens to have just the right length to trigger a
        # bug in pymysql when the string length is encoded using a short
        # (0xFC)
        cursor.execute('''select sequence from dna where
          seq_region_id = 15754''')
        c = mus.conn.cursor()
        c.execute("""select seq_region_id,length from seq_region where
        coord_system_id = %s and name = 'MT';""", 1)
        print list(c)
        c = mus.conn.cursor()
        c.execute("""select seq_region_id,length from seq_region where
        coord_system_id = %d and name = %s;""", (1, 'MT'))
        print list(c)

def main(argv=None):
    if argv is None:
        import sys
        argv = sys.argv

    arg = argv[1:]
    if len(arg) > 0 and arg[0] == 'debug':
        debug()

if __name__ == '__main__':
    main()
