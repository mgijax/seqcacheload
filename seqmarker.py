#!/usr/local/bin/python
#
# seqmarker.py  
##########################################################################k
#
#  Purpose: This script creates the bcp file for SEQ_Marker_Cache
#           which caches sequence/marker pairs including:
#               1) all marker statuses
#	        2) all marker types
#               3) organisms: mouse, rat, human, dog, chimp
#               4) all sequence statuses except deleted
#	    This script also determines the representative genomic, 
#		transcript and protein sequence for each marker (if it can)
#	    A record in SEQ_Marker_Cache represents a uniq sequence, marker
#           reference relationship; there may be multiple references for 
#           a sequence to marker association, therefore multiple records
#           per sequence/marker pair
#
#  Usage:
#	seqmarker.py
#
#  Env Vars: Uses environment variables to determine Server and Database
#	  (DSQUERY and MGD).
#
#  Inputs: 1) mgd database
#          2) Configuration
#
#  Outputs: 1) log file
#           2) bcp file
#
#  Exit Codes:
#
#  History
#
#  02/07/2008	sc
#	- TR 8490 new genomic rep sequence algorithm
#
#  01/18/2007	lec
#	- include withdrawn markers as they may have sequence associations
#
#  02/02/2007	lec
#	- TR 8072; exclude deleted sequences from representative algorithm
#
#  06/28/2006	lec
#	- add dog & chimpanzee (TR 7508)
#
#  10/12/2005	lec
#	- add primary acc id
#
#  09/08/2005	lec
#	- PIRSF (TR 5972)
#	- add human/rat
#
#  10/23/2003	lec
#	- new (TR 3404, JSAM)
#

import sys
import os
import string
import db
import mgi_utils
import loadlib

# database errors
DB_ERROR = 'A database error occured: '
DB_CONNECT_ERROR = 'Connection to the database failed: '

# record delimiter
NL = '\n'

# column delimiter
DL = os.environ['COLDELIM']

# table for which we are creating bcp file
table = os.environ['TABLE']

# output directory
datadir = os.environ['CACHEDATADIR']

# date with which to record stamp database records
loaddate = loadlib.loaddate

# name of bcp file descriptor
outBCP = None

# representative sequence qualifier lookup by term
# looks like {qualifier:qualKey, ...}
qualByTermLookup = {}

# representative sequence qualifier lookup by term key
# looks like {qualKey:qualifier, ...}
qualByTermKeyLookup = {}

# marker lookup by genomic sequence key (to see if seq assoc with other markers)
# looks like {seqKey:[mkrKey1, ..., mkrKeyn], ...}
mkrsByGenomicSeqKeyLookup = {}

# represents all genomic seqs for the current marker by provider - see indexes
# each vega/genbank dictionary looks like {_Marker_key:[seqKey1, seqKeyN], ...}
# each ensembl/ncbi dictionary looks like (_Marker_key:{}, ...} where {} is  
# a db.sql result set
allgenomic = [{}, {}, {}, {}]

# indexes of allgenomic
VEGA = 0
ENSEMBL = 1
NCBI = 2
gGENBANK = 3

# represents the longest transcript for the current marker by provider 
# each dictionary looks like {_Marker_key:_Sequence_key, ...}
alltranscript = [{}, {}, {}, {}, {}, {}, {}]

# indexes of alltranscript:
# 0=RefSeq NR
# 1=GenBank RNA, not EST
# 2=Refseq XM
# 3=DFCI
# 4=DoTS
# 5=NIA
# 6=GenBank RNA, EST

#  represents the longest protein for the current marker by provider
# each dictionary looks like: {_Marker_key:_Sequence_key, ...}
allpolypeptide = [{}, {}, {}, {}]
# indexes of allpolypeptide:
# 0=SwissProt
# 1=RefSeq NP
# 2=TrEMBL
# 3=RefSeq XP

# each dictionary looks like {_Marker_key:_Sequence_key, ...}
# the set of genomic representative sequences for markers
genomic = {}

# the set of transcript representative sequences for markers
transcript = {}

# the set of protein representative sequences for markers
polypeptide = {}

# Purpose: Initialize db.sql settings, lookups, and file descriptor
# Returns: Nothing
# Assumes: Nothing
# Effects: connects to and queries a database
# Throws: Nothing

def init ():
    global qualByTermLookup, qualByTermKeyLookup, mkrsByGenomicSeqKeyLookup
    global outBCP
    
    db.useOneConnection(1)
    db.set_sqlLogFunction(db.sqlLogAll)

    print 'Initializing ...%s' % (mgi_utils.date())
    #
    # load representative sequence qualifer lookups
    #
    results = db.sql('select _Term_key, ' + \
        'term from VOC_Term_RepQualifier_View', 'auto')

    for r in results:
       qualByTermLookup[r['term']] = r['_Term_key']
       qualByTermKeyLookup[r['_Term_key']] = r['term']

    # query with which to load
    # genomic sequences associated with markers lookup
    # for VEGA Gene Model(85), Ensembl Gene Model (60),
    # NCBI Gene Model(59), GenBank DNA (9)

    # get the set of all GenBank DNA sequences
    db.sql('select a.accID as seqID, a._Object_key as _Sequence_key ' + \
        'into #gbDNA ' + \
        'from ACC_Accession a, SEQ_Sequence s ' + \
        'where a._LogicalDB_key = 9  ' + \
        'and a._MGIType_key = 19  ' + \
        'and a._Object_key = s._Sequence_key  ' + \
        'and s._SequenceType_key = 316347', None)
    # get the set of all Ensembl, NCBI, VEGA gene models
    db.sql('select a.accID as seqID, a._Object_key as _Sequence_key ' + \
	'into #gm ' + \
	'from ACC_Accession a ' + \
        'where a._LogicalDB_key in (59, 60, 85) ' + \
        'and a._MGIType_key = 19', None)
    # union the set
    db.sql('select seqID, _Sequence_key ' + \
	'into #allSeqs ' + \
	'from #gbDNA ' + \
	'union ' + \
	'select seqID, _Sequence_key ' + \
	'from #gm', None)
    db.sql('create nonclustered index idx_1 on ' + \
        '#allSeqs (seqID)', None)
    # get markers for these sequences
    results = db.sql('select s._Sequence_key, a._Object_key as _Marker_key ' + \
	'from #allSeqs s, ACC_Accession a ' + \
	'where a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key in (59, 60, 85, 9) ' + \
	'and s.seqID = a.accid ' + \
	'order by _Sequence_key ', 'auto')
    print 'length of results: %s ' % len(results)

    # load the lookup
    prevSeqKey = ''
    for r in results:
	seqKey = r['_Sequence_key']
	markerKey = r['_Marker_key']
        if seqKey != prevSeqKey:
	    mkrsByGenomicSeqKeyLookup[seqKey] = [markerKey]
        else:
	    mkrsByGenomicSeqKeyLookup[seqKey].append(markerKey)
	    #print 'appending marker %s to sequence %s' % (markerKey, seqKey)
	prevSeqKey = seqKey
    # create file descriptor for bcp file
    outBCP = open('%s/%s.bcp' % (datadir, table), 'w')
    return

# Purpose: Determines representative genomic, transcript, and protein
#          for the given marker
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def determineRepresentative(prevMarker):
    global genomic, transcript, polypeptide

    #
    # Determine Representative Genomic Sequence
    #

    # the current choice for genomic rep
    currGenomicRep = ''

    # If there is only one VEGA sequence for this marker and it is not
    # associated with any other markers, it is representative
    vegaKey = ''
    if allgenomic[VEGA].has_key(prevMarker) and \
	    len(allgenomic[VEGA][prevMarker]) == 1:
        vegaKey = allgenomic[VEGA][prevMarker][0]
	#print 'VEGA: %s numMarkers: %s' % (vegaKey, len(mkrsByGenomicSeqKeyLookup[vegaKey]))
        # it is representative if associated with only one marker
        if len(mkrsByGenomicSeqKeyLookup[vegaKey]) == 1:
            currGenomicRep = vegaKey

    # If we didnt' get a VEGA, look at Ensembl and NCBI
    if currGenomicRep == '':
        ensemblKey = ''
        ensemblLength = 0
        ncbiKey = ''
        ncbiLength = 0
        # get Ensembl id if only one Ensembl sequence for this marker
        # and it is not associated with any other markers
        if allgenomic[ENSEMBL].has_key(prevMarker) and \
		len(allgenomic[ENSEMBL][prevMarker]) == 1:
            ensemblResult = allgenomic[ENSEMBL][prevMarker][0]
            ensemblKey =  ensemblResult['_Sequence_key']
	    #print 'Ensembl: %s numMarkers: %s' % (ensemblKey, len(mkrsByGenomicSeqKeyLookup[ensemblKey]))

            if len(mkrsByGenomicSeqKeyLookup[ensemblKey]) == 1:
                ensemblLength = int(ensemblResult['length'])
            else:
                ensemblKey = ''

        # get NCBI id if only one NCBI sequence for this marker
        # and it is not associated with any other markers
        if allgenomic[NCBI].has_key(prevMarker) and \
		len(allgenomic[NCBI][prevMarker]) == 1:
            ncbiResult = allgenomic[NCBI][prevMarker][0]
            ncbiKey = ncbiResult['_Sequence_key']
	    #print 'NCBI: %s numMarkers: %s' % (ncbiKey, len(mkrsByGenomicSeqKeyLookup[ncbiKey]))

            if len(mkrsByGenomicSeqKeyLookup[ncbiKey]) == 1:
                ncbiLength = int(ncbiResult['length'])
            else:
                ncbiKey = ''

        # if there is an Ensembl AND an NCBI , take the longest
        if ensemblKey != '' and ncbiKey != '':
            # ensembl is rep if longest
            if ensemblLength > ncbiLength:
                currGenomicRep = ensemblKey
            # ncbi is rep if 1) lengths are equal 2) ncbi longest
            else:
                currGenomicRep = ncbiKey
        # take whichever is defined, if any
        elif ensemblKey != '':
            currGenomicRep = ensemblKey
        elif ncbiKey != '':
            currGenomicRep = ncbiKey

    # if we didn't get an Ensembl or an NCBI, look at GenBank
    if currGenomicRep == '':
        genbankKey = ''
        if allgenomic[gGENBANK].has_key(prevMarker) and \
		len(allgenomic[gGENBANK][prevMarker]) == 1:
            genbankKey = allgenomic[gGENBANK][prevMarker][0]
  	    #print 'GenBank: %s numMarkers: %s' % (genbankKey, len(mkrsByGenomicSeqKeyLookup[genbankKey]))

            # it is representative if associated with only one marker
            if len(mkrsByGenomicSeqKeyLookup[genbankKey]) == 1:
		currGenomicRep = genbankKey
    if currGenomicRep != '':
	genomic[prevMarker] = currGenomicRep

    #
    # Determine Representative Transcript Sequence
    #
    for i in range(len(alltranscript)):
        if alltranscript[i].has_key(prevMarker):
            transcript[prevMarker] = []
            transcript[prevMarker].append(alltranscript[i][prevMarker])
            break
    #
    # Determine Representative Protein Sequence
    #
    for i in range(len(allpolypeptide)):
        if allpolypeptide[i].has_key(prevMarker):
            polypeptide[prevMarker] = allpolypeptide[i][prevMarker]
            break
    return

# Purpose: formats and writes out record to bcp file
# Returns: Nothing
# Assumes: Nothing
# Effects: writes a record to a file
# Throws: Nothing

def writeRecord(r):
    outBCP.write(mgi_utils.prvalue(r['_Sequence_key']) + DL + \
	mgi_utils.prvalue(r['_Marker_key']) + DL + \
	mgi_utils.prvalue(r['_Organism_key']) + DL + \
	mgi_utils.prvalue(r['_Refs_key']) + DL)

    printedQualifier = 0
    if genomic.has_key(r['_Marker_key']):
	if genomic[r['_Marker_key']] == r['_Sequence_key']:
	    outBCP.write(mgi_utils.prvalue(qualByTermLookup['genomic']) + DL)
	    printedQualifier = 1

    if transcript.has_key(r['_Marker_key']):
       if r['_Sequence_key'] in transcript[r['_Marker_key']]:
	    outBCP.write(mgi_utils.prvalue(qualByTermLookup['transcript']) + DL)
	    printedQualifier = 1

    if polypeptide.has_key(r['_Marker_key']):
        if polypeptide[r['_Marker_key']] == r['_Sequence_key']:
	    outBCP.write(mgi_utils.prvalue(qualByTermLookup['polypeptide']) + DL)
	    printedQualifier = 1

    if not printedQualifier:
        outBCP.write(mgi_utils.prvalue(qualByTermLookup['Not Specified']) + DL)

    outBCP.write(mgi_utils.prvalue(r['_SequenceProvider_key']) + DL + \
	mgi_utils.prvalue(r['_SequenceType_key']) + DL + \
	mgi_utils.prvalue(r['_LogicalDB_key']) + DL + \
	r['accID'] + DL + \
	r['mdate'] + DL + \
        mgi_utils.prvalue(r['_User_key']) + DL + \
	mgi_utils.prvalue(r['_User_key']) + DL + \
        loaddate + DL + loaddate + NL)

# Purpose: Iterates through result set of sequence marker pairs
#          determining the representative sequence qualifier for each
#          (genomic, transcript, protein, or none of the above), writing
#          pairs out to bcp file
# Returns: Nothing
# Assumes: Nothing
# Effects: queries a database, writes records to a file
# Throws: Nothing

def createBCP():
    global allgenomic, alltranscript, allpolypeptide
    global outBCP

    print 'Processing ...%s' % (mgi_utils.date())
    #	
    # select only mouse, human, rat, dog & chimpanzee markers
    # with ANY marker status 

    db.sql('select _Marker_key, _Organism_key, _Marker_Type_key ' + \
	'into #markers from MRK_Marker ' + \
	'where _Organism_key in (1, 2, 40, 10, 13)', None)
    db.sql('create nonclustered index idx_key on ' + \
	'#markers (_Marker_key)', None)

    # select all non-MGI accession ids for markers 

    db.sql('select m._Marker_key, m._Organism_key, m._Marker_Type_key, ' + \
	'a._LogicalDB_key, a.accID, r._Refs_key, a._ModifiedBy_key, ' + \
	'mdate = convert(char(10), a.modification_date, 101) ' + \
	'into #markerAccs ' + \
	'from #markers m, ACC_Accession a, ACC_AccessionReference r ' + \
	'where m._Marker_key = a._Object_key ' + \
	'and a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key != 1 ' + \
	'and a._Accession_key = r._Accession_key', None)

    db.sql('create nonclustered index idx1 on #markerAccs ' + \
	'(_LogicalDB_key, accID)', None)

    # select all sequence annotations

    db.sql('select _Sequence_key = s._Object_key, ' + \
	'm._Marker_key, m._Organism_key, ' + \
	'm._Marker_Type_key, ' + \
	'm._LogicalDB_key, ' + \
	'm._Refs_key, m._ModifiedBy_key as _User_key, ' + \
	'm.mdate, m.accID ' + \
	'into #preallannot ' + \
	'from #markerAccs m, ACC_Accession s ' + \
	'where m.accID = s.accID ' + \
	'and m._LogicalDB_key = s._LogicalDB_key ' + \
	'and s._MGIType_key = 19 ', None)

    db.sql('create nonclustered index idx1 on ' + \
	'#preallannot (_Sequence_key)', None)

    # get the sequence provider and sequence type

    db.sql('select m._Sequence_key, m._Marker_key, m._Organism_key, ' + \
	'm._Marker_Type_key, ss._SequenceProvider_key, ' + \
	'ss._SequenceType_key, ' + \
	'm._LogicalDB_key, m._Refs_key, m._User_key, m.mdate, m.accID ' + \
	'into #allannot ' + \
	'from #preallannot m, SEQ_Sequence ss ' + \
	'where m._Sequence_key = ss._Sequence_key', None)

    db.sql('create nonclustered index idx1 on ' + \
	'#allannot (_Sequence_key)', None)

    # grab sequence's primary accID

    db.sql('select a._Sequence_key, a._Marker_key, a._Organism_key, ' + \
	'a._Marker_Type_key, a._SequenceProvider_key, ' + \
	'a._SequenceType_key, a._LogicalDB_key, ' + \
	'a._Refs_key, a._User_key, a.mdate, ac.accID ' + \
	'into #allseqannot ' + \
	'from #allannot a, ACC_Accession ac ' + \
	'where a._Sequence_key = ac._Object_key ' + \
	'and ac._MGIType_key = 19 ' + \
	'and ac.preferred = 1', None)

    db.sql('create nonclustered index idx1 on ' + \
	'#allseqannot (_Sequence_key, _Marker_key, _Refs_key)', None)

    # select records, grouping by sequence, marker and reference
    db.sql('select _Sequence_key, _Marker_key, _Organism_key, ' + \
	'_Marker_Type_key, _SequenceProvider_key, _SequenceType_key, ' + \
	'_LogicalDB_key, _Refs_key, _User_key, mdate = max(mdate), accID ' + \
	'into #finalannot ' + \
	'from #allseqannot ' + \
	'group by _Sequence_key, _Marker_key, _Refs_key', None)
    db.sql('create nonclustered index idx1 on #finalannot ' + \
	'(_Sequence_key, _Marker_key, _Refs_key, _User_key, mdate)', None)
    db.sql('create nonclustered index idx2 on #finalannot ' + \
	'(_Sequence_key, _Marker_key, _Marker_Type_key, accID)', None)
    db.sql('create nonclustered index idx3 on #finalannot ' + \
	'(_Marker_key)', None)

    db.sql('select distinct _Sequence_key, _Marker_key, ' + \
	'_Marker_Type_key, accID ' + \
	'into #deriveQuality ' + \
	'from #finalannot order by _Marker_key', None)
    db.sql('create nonclustered index idx1 on ' + \
	'#deriveQuality (_Sequence_key)', None)
    db.sql('create nonclustered index idx2 on ' + \
	'#deriveQuality (_Marker_key)', None)

    # do not include deleted sequences
    results = db.sql('select q._Sequence_key, q._Marker_key, ' + \
	'q._Marker_Type_key, q.accID, s._SequenceProvider_key, ' + \
	's._SequenceType_key, s.length ' + \
	'from #deriveQuality q, SEQ_Sequence s ' + \
	'where q._Sequence_key = s._Sequence_key ' + \
	'and s._SequenceStatus_key != 316343 ' + \
	'order by q._Marker_key, s._SequenceProvider_key', 'auto')

    # process derived representative values
    prevMarker = ''

    for r in results:
	m = r['_Marker_key']
	s = r['_Sequence_key']
	a = r['accID']

	if r['length'] is None:
	    seqlength = 0
	else:
	    seqlength = int(r['length'])

	providerKey = r['_SequenceProvider_key']
	seqTypeKey = r['_SequenceType_key']

        # lengths for transcript and polypeptide, we do genomic differently
	if prevMarker != m:
	    tlengths = [-1,-1,-1,-1,-1,-1,-1]
	    plengths = [-1,-1,-1,-1]

	    if prevMarker != '':
		determineRepresentative(prevMarker)

	# representative genomic
	# a genomic sequence can be representative only if 
	# 1) it is not also associated with another marker 
	#    e.g. if OTTMUSG00000012345 is associated with Mkr1 and Mkr2
	#    it cannot be a representative for either.
	# 2) the marker is not alsociated with another sequence of the same
	#    type .e.g. if Mkr1 is associated with OTTMUSG00000012345 (VEGA)
	#      and Mkr1 is also associated with OTTMUSG00000098765 then
	#      neither VEGA sequence can be the representative for Mkr1
	# Within the above constraints pick:
	# 1) VEGA
	# 2) NCBI/Ensembl coordinate. Pick longest if both, 
        #    if both same length pick NCBI
	# 3) GenBank DNA
	#

	# VEGA
	# for vega we don't care how  long it is
	if providerKey == 1865333:
	    if allgenomic[VEGA].has_key(m):
		allgenomic[VEGA][m].append(s)
	    else:
		allgenomic[VEGA][m] = [s]

	# Ensembl
	# if there is both Ensembl and NCBI we need length
	# so store the entire results in the dictionary
	elif (providerKey in [615429]):
            if allgenomic[ENSEMBL].has_key(m):
                allgenomic[ENSEMBL][m].append(r)
            else:
                allgenomic[ENSEMBL][m] = [r]

	# NCBI
	elif providerKey == 706915:
            if allgenomic[NCBI].has_key(m):
                allgenomic[NCBI][m].append(r)
            else:
                allgenomic[NCBI][m] = [r]

	# any GenBank; DNA
	# these are all GenBank provider terms by division
	# e.g. "GenBank/EMBL/DDBJ:Rodent" or "GenBank/EMBL/DDBJ:GSS"
	elif providerKey in \
	      [316380,316376,316379,316375,316377,316374,316373,316378,492451] \
	      and seqTypeKey == 316347:
            if allgenomic[gGENBANK].has_key(m):
                allgenomic[gGENBANK][m].append(s)
            else:
                allgenomic[gGENBANK][m] = [s]

	#
	# representative transcript
	#
	# longest NM_ or NR_ RefSeq
	# longest non-EST GenBank
	# longest XM_ or XR_ RefSeq
	# longest DFCI, DoTS, NIA Mouse Gene Index,
	# longest EST GenBank
	#

	# RefSeq
	if providerKey == 316372 and (string.find(a, 'NM_') > -1 or \
	      string.find(a, 'NR_') > -1):
	    if seqlength > tlengths[0]:
		alltranscript[0][m] = s
		tlengths[0] = seqlength

	# GenBank but not EST; RNA
	elif providerKey in [316380,316379,316375,316377,316374,316373,316378,492451] and seqTypeKey == 316346:
	    if seqlength > tlengths[1]:
		alltranscript[1][m] = s
		tlengths[1] = seqlength

	# RefSeq
	elif providerKey == 316372 \
	      and (string.find(a, 'XM_') > -1 \
	      or string.find(a, 'XR_') > -1):
	    if seqlength > tlengths[2]:
		alltranscript[2][m] = s
		tlengths[2] = seqlength

	# DFCI
	elif providerKey == 316381:
	    if seqlength > tlengths[3]:
		alltranscript[3][m] = s
		tlengths[3] = seqlength

	# DoTS
	elif providerKey == 316382:
	    if seqlength > tlengths[4]:
		alltranscript[4][m] = s
		tlengths[4] = seqlength

	# NIA
	elif providerKey == 316383:
	    if seqlength > tlengths[5]:
		alltranscript[5][m] = s
		tlengths[5] = seqlength

	# GenBank EST; RNA
	elif providerKey == 316376 and seqTypeKey == 316346:
	    if seqlength > tlengths[6]:
		alltranscript[6][m] = s
		tlengths[6] = seqlength

	#
	# representative polypeptide
	#
	# longest SWISS-PROT
	# longest NP_ RefSeq
	# longest TrEMBL
	# longest XP_ RefSeq
	#

	# SwissProt
	if providerKey == 316384 and seqlength > plengths[0]:
	    allpolypeptide[0][m] = s
	    plengths[0] = seqlength

	# RefSeq
	elif providerKey == 316372 and string.find(a, 'NP_') > -1 \
	      and seqlength > plengths[1]:
	    allpolypeptide[1][m] = s
	    plengths[1] = seqlength

	# TrEMBL
	elif providerKey == 316385 and seqlength > plengths[2]:
	    allpolypeptide[2][m] = s
	    plengths[2] = seqlength

	# RefSeq
	elif providerKey == 316372 and string.find(a, 'XP_') > -1 \
	      and seqlength > plengths[3]:
	    allpolypeptide[3][m] = s
	    plengths[3] = seqlength

	prevMarker = m

    # last record
    determineRepresentative(prevMarker)

    print 'Writing bcp file ...%s' % (mgi_utils.date())
    results = db.sql('select distinct _Sequence_key, _Marker_key, ' + \
	'_Organism_key, _SequenceProvider_key, _SequenceType_key, ' + \
	'_LogicalDB_key, _Refs_key, ' + \
	'_User_key, mdate, accID ' + \
	'from #finalannot', 'auto')
    
    # results are ordered by  _Sequence_key, _Marker_key, _Refs_key
    for r in results:
	writeRecord(r)


# Purpose: Perform cleanup steps for the script.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def finalize():
    global outBCP

    db.useOneConnection(0)
    outBCP.close()
    return

#
# Main Routine
#
try:
    init()
    createBCP()
    finalize()
except db.connection_exc, message:
    error = '%s%s' % (DB_CONNECT_ERROR, message)
    sys.stderr.write(message)
    sys.exit(message)
except db.error, message:
    error = '%s%s' % (DB_ERROR, message)
    sys.stderr.write(message)
    sys.exit(message)

