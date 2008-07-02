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
# each dictionary looks like (_Marker_key:{}, ...} 
# where {} is a db.sql result set
allgenomic = [{}, {}, {}, {}]

# indexes of allgenomic
VEGA = 0
ENSEMBL = 1
NCBI = 2
gGENBANK = 3

# indexes of genomic sequence dictionaries
SEQKEY=0
UNIQ=1
LENGTH=2

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

    # get the set of all preferred GenBank DNA sequences
    db.sql('select a.accID as seqID, a._Object_key as _Sequence_key ' + \
        'into #gbDNA ' + \
        'from ACC_Accession a, SEQ_Sequence s ' + \
        'where a._LogicalDB_key = 9  ' + \
        'and a._MGIType_key = 19  ' + \
	'and a.preferred = 1 ' + \
        'and a._Object_key = s._Sequence_key  ' + \
        'and s._SequenceType_key = 316347', None)
    # get the set of all Ensembl, NCBI, VEGA gene models
    db.sql('select a.accID as seqID, a._Object_key as _Sequence_key ' + \
	'into #gm ' + \
	'from ACC_Accession a ' + \
        'where a._LogicalDB_key in (59, 60, 85) ' + \
	'and a.preferred = 1 ' + \
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

def determineRepresentative(marker):
    global genomic, transcript, polypeptide
    print 'determineRep for marker: %s' % marker
    #
    # Determine Representative Transcript Sequence
    #
    for i in range(len(alltranscript)):
        if alltranscript[i].has_key(marker):
            transcript[marker] = []
            transcript[marker].append(alltranscript[i][marker])
            break
    #
    # Determine Representative Protein Sequence
    #
    for i in range(len(allpolypeptide)):
        if allpolypeptide[i].has_key(marker):
            polypeptide[marker] = allpolypeptide[i][marker]
            break

    #
    # Determine Representative Genomic Sequence
    #
    # see algorithm here: 
    # http://prodwww.informatics.jax.org/wiki/index.php/sw%
    #          3ARepresentative_sequence_algorithm
    # 

    ##-------------------------------------------------------------
    # Determine attributes for each provider and provider sequence
    # for this marker
    ##-------------------------------------------------------------

    ##--------------------------------------
    # The attributes
    ##--------------------------------------

    # * = VEGA|Ensembl|NCBI|GenBank
    # True if this marker has a * sequence
    # NOTE: we need these has* variables to determine multiple and not
    # uniq because, for example,  a marker can still have a VEGA Gene Model
    # if vegaIsSgl = False and vegaHasUniq = False 
    hasVega = False
    hasEnsembl = False
    hasNCBI = False
    hasGenBank = False

    # True if this marker has only one * id
    vegaIsSgl = False
    ensemblIsSgl = False
    ncbiIsSgl = False
    genbankIsSgl = False

    # True if this marker has a unique * sequence (associated with only
    # this marker)
    vegaHasUniq = False
    ensemblHasUniq = False
    ncbiHasUniq = False
    genbankHasUniq = False

    # list of dictionaries of * seqs for this marker
    # looks like [{UNIQ:True/False, SEQKEY:key, LENGTH:length}, ...]
    vegaSeqs = []
    ensemblSeqs = []
    ncbiSeqs = []
    genbankSeqs = []

    ##--------------------------------------
    # Get VEGA attributes 
    ##--------------------------------------    
    if allgenomic[VEGA].has_key(marker):

        hasVega = True

	# if this marker has only one VEGA sequence flag it as single
	if len(allgenomic[VEGA][marker]) == 1:
	    vegaIsSgl = True

	# get seqKey, seqLength, and uniqueness for each VEGA sequence 
	for result in allgenomic[VEGA][marker]:
	    seqDict =  {}
	    seqKey = result['_Sequence_key']
	    length = result['length']
	    seqDict[SEQKEY] = seqKey
	    seqDict[LENGTH] = length

	    value = False
	    if mkrsByGenomicSeqKeyLookup.has_key(seqKey) and \
		    len(mkrsByGenomicSeqKeyLookup[seqKey]) == 1:
		value = True
		vegaHasUniq = True
	    seqDict[UNIQ] = value
	    vegaSeqs.append(seqDict)
	print 'vegaseqs: %s' % vegaSeqs
    ##--------------------------------------
    # Get Ensembl attributes
    ##--------------------------------------
    if allgenomic[ENSEMBL].has_key(marker):
	
   	hasEnsembl = True
	# if this marker has only one ENSEMBL sequence flag it as single
	if len(allgenomic[ENSEMBL][marker]) == 1:
	    ensemblIsSgl = True
	# get seqKey, seqLength, and uniqueness for each ENSEMBL sequence
	for result in allgenomic[ENSEMBL][marker]:
	    seqDict =  {}
	    seqKey = result['_Sequence_key']
	    length = result['length']
	    seqDict[SEQKEY] = seqKey
	    seqDict[LENGTH] = length

	    value = False
	    if mkrsByGenomicSeqKeyLookup.has_key(seqKey) and \
		    len(mkrsByGenomicSeqKeyLookup[seqKey]) == 1:
		value = True
		ensemblHasUniq = True
	    seqDict[UNIQ] = value
	    ensemblSeqs.append(seqDict)
	print 'ensemblseqs: %s' % ensemblSeqs
	#print 'marker: %s hasEnsembl = %s' % (marker, hasEnsembl)
	#print 'marker: %s ensemblIsSgl = %s' % (marker, ensemblIsSgl)
	#print 'marker %s ensemblHasUniq = %s' % (marker, ensemblHasUniq)
    ##--------------------------------------
    # Get NCBI attributes
    ##--------------------------------------
    if allgenomic[NCBI].has_key(marker):
        
	hasNCBI = True

	# if this marker has only one NCBI sequence flag it as single
	if len(allgenomic[NCBI][marker]) == 1:
	    ncbiIsSgl = True

	# get seqKey, seqLength, and uniqueness for each NCBI sequence
	for result in allgenomic[NCBI][marker]:
	    seqDict =  {}
	    seqKey = result['_Sequence_key']
	    length = result['length']
	    seqDict[SEQKEY] = seqKey
	    seqDict[LENGTH] = length

	    value = False
	    if mkrsByGenomicSeqKeyLookup.has_key(seqKey) and \
		    len(mkrsByGenomicSeqKeyLookup[seqKey]) == 1:
		value = True
		ncbiHasUniq = True
	    seqDict[UNIQ] = value
	    ncbiSeqs.append(seqDict)
	print 'ncbiseqs: %s' % ncbiSeqs
    ##--------------------------------------
    # Get GenBank attributes
    ##--------------------------------------
    if allgenomic[gGENBANK].has_key(marker):
        hasGenBank = True 

	# if this marker has only one GENBANK sequence flag it as single
	if len(allgenomic[gGENBANK][marker]) == 1:
	    genbankIsSgl = True

	# get seqKey, seqLength, and uniqueness for each GENBANK sequence
	for result in allgenomic[gGENBANK][marker]:
	    seqDict =  {}
	    seqKey = result['_Sequence_key']
	    length = result['length']
	    seqDict[SEQKEY] = seqKey
	    seqDict[LENGTH] = length

	    value = False
	    if mkrsByGenomicSeqKeyLookup.has_key(seqKey) and \
		    len(mkrsByGenomicSeqKeyLookup[seqKey]) == 1:
		value = True
		genbankHasUniq = True
	    seqDict[UNIQ] = value
	    genbankSeqs.append(seqDict)
	print 'genbankseqs: %s' % genbankSeqs
    ##-----------------------------------------------------------------
    # Determine representative genomic for marker using provider
    # and provider sequence attributes
    ##-----------------------------------------------------------------

    genomicRep = ''
   
    hasSglUniqEnsembl = False
    hasSglUniqNCBI = False

    # marker has no GMs
    if not hasVega and not hasEnsembl and not hasNCBI: 
	# is there a longest uniq
	(s,l) = determineSeq(genbankSeqs, True, True)
	# no sequence found that match parameters if s == 0
	if s != 0:
	    genomicRep = s
	    print 'CASE 1'
	# if no longest uniq get longest
	if genomicRep == '':
	    (s,l) = determineSeq(genbankSeqs, True, False)
	    if s != 0:
		genomicRep = s
		print 'CASE 2'
	# else NO REPRESENTATIVE SEQUENCE value of genomicRep is still ''

    # marker has at least one GM, therefore WILL HAVE REP SEQUENCE
    else: 
	if vegaIsSgl and vegaHasUniq:
	    genomicRep = vegaSeqs[0][SEQKEY]
	    print 'CASE 3'
	else:    # no single uniq Vega exists
	    # check single uniq Ensembl and NCBI
	    if ensemblIsSgl and ensemblHasUniq:
		hasSglUniqEnsembl = True
	    if ncbiIsSgl and ncbiHasUniq:
		hasSglUniqNCBI = True
	    if hasSglUniqEnsembl and hasSglUniqNCBI:
		ensemblLen =   ensemblSeqs[0][LENGTH]
		ncbiLen = ncbiSeqs[0][LENGTH]
		value = determineShortest(ensemblLen, ncbiLen)
		# if ensembl and ncbi length equal or ncbi shorter pick ncbi
		if value == -1 or value == 1:
		    genomicRep = ncbiSeqs[0][SEQKEY]
		    print 'CASE 4'
		else:
		    genomicRep =  ensemblSeqs[0][SEQKEY]
		    print 'CASE 5'
	    elif hasSglUniqEnsembl:
		genomicRep = ensemblSeqs[0][SEQKEY]
		print 'CASE 6'
	    elif hasSglUniqNCBI:
		genomicRep = ncbiSeqs[0][SEQKEY]
	  	print 'CASE 7'
	    # only multiples (uniq or not) or single not-uniq left
	    else:
		# check uniq (must be multiple)
                if vegaHasUniq or ensemblHasUniq or ncbiHasUniq:
		    if vegaHasUniq:
			# pick shortest uniq, if tie pick one
			(s,l) = determineSeq(vegaSeqs, False, True)
			genomicRep = s
			print 'CASE 8'
		    elif ensemblHasUniq and ncbiHasUniq:
			# pick shortest uniq, if tie pick one
			(s_e, l_e) = determineSeq(ensemblSeqs, False, True)
			(s_n, l_n) = determineSeq(ncbiSeqs, False, True)
			value = determineShortest(l_e, l_n)
			if value == -1 or value == 1:
			    genomicRep = s_n
			    print 'CASE 9'
			else:
			    genomicRep = s_e
			    print 'CASE 10'
		    elif ensemblHasUniq:
			(s_e, l_e) = determineSeq(ensemblSeqs, False, True)
			genomicRep = s_e
			print 'CASE 11'
		    elif ncbiHasUniq:
			(s_n, l_n) = determineSeq(ncbiSeqs, False, True)
			genomicRep = s_n
			print 'CASE 12'
		# no uniques, only single or multiple non-uniq left
		else:
		    # check for vega sgl
		    if vegaIsSgl:
			genomicRep = vegaSeqs[0][SEQKEY]
			print 'CASE 13'
		    else:
			# check for ensembl and ncbi sgl
			if ensemblIsSgl or ncbiIsSgl:
			    if ensemblIsSgl and ncbiIsSgl:
				value = determineShortest( \
				    ensemblSeqs[0][LENGTH], \
				    ncbiSeqs[0][LENGTH])
				if value == -1 or value == 1:
				    genomicRep = ncbiSeqs[0][SEQKEY]
				    print 'CASE 14'
				else:
				    genomicRep =  ensemblSeqs[0][SEQKEY]
				    print 'CASE 15'
			    elif ensemblIsSgl:
				genomicRep = ensemblSeqs[0][SEQKEY]
				print 'CASE 16'
			    elif ncbiIsSgl:
				genomicRep = ncbiSeqs[0][SEQKEY]
				print 'CASE 17'
			# no singles, must be multiple non-uniq
			else:
			    if hasVega:
                                # pick shortest
				(s,l) = determineSeq(vegaSeqs, False, False)
				genomicRep = s
				print 'CASE 18'
                            elif hasEnsembl and hasNCBI:
                                # pick shortest, NCBI if tie
				(s_e, l_e) = determineSeq( \
				    ensemblSeqs, False, False)
				(s_n, l_n) = determineSeq( \
				    ncbiSeqs, False, False)
				value = determineShortest(l_e, l_n)
				if value == -1 or value == 1:
				    genomicRep = s_n
				    print 'CASE 19'
				else:
				    genomicRep =  s_e
				    print 'CASE 20'
                            elif hasEnsembl:
                                # pick shortest
				(s_e, l_e) = determineSeq( \
				    ensemblSeqs, False, False)
				genomicRep = s_e
				print 'CASE 21'
			    # must be an NCBI
                            else:
                                # pick shortest NCBI
				(s_n, l_n) = determineSeq( \
                                    ncbiSeqs, False, False)
				genomicRep = s_n
				print 'CASE 22'
    # if we found a genomicRep for this marker add it to the dictionary
    print 'genomicRep: %s' % genomicRep
    if genomicRep != '':
        genomic[marker] = genomicRep
    else:
	print 'CASE 23'
    return

# Purpose: Find the longest or shortest sequence given a dictionary like:
#          [{UNIQ:True/False, SEQKEY:key, LENGTH:length}, ...]
# Returns: tuple (seqKey, length) where sequence key is the longest/shortest 
#          uniq/notuniq depending on value of 'getLongest' and 'useUniq', 
#          else both members of the tuple are 0 
# Assumes: Nothing
# Throws: Nothing
#
def determineSeq(seqList, 	# list of dictionaries
                getLongest,     # boolean, determine longest if True, else
                                # shortest
		useUniq): 	# boolean, consider uniq only if True
    # current choice considering only uniq sequences
    # based on value of getLongest

    # a non-numeric default
    currUniqLen = ''
    currUniqSeqKey = 0
    # current choice considering all sequences based on value of getLongest
    currAllLen = ''
    currAllSeqKey = 0
    #print 'getLongest: %s, useUniq: %s' % (getLongest, useUniq)
    for seq in seqList:
        # if we are looking for the longest sequence
	if getLongest == True:
	     # current choice from the uniq set only
	    if seq[UNIQ] == True:
		l = determineLongest(currUniqLen, seq[LENGTH])
		# if seq[LENGTH] is longest or equal
		if l == 1 or l == -1:
		    currUniqLen = seq[LENGTH]
		    currUniqSeqKey = seq[SEQKEY]
	    # current choice from the full set 
	    l = determineLongest(currAllLen, seq[LENGTH])
	    # if seq[LENGTH] is longest or equal
	    if l == 1 or l == -1:
		currAllLen = seq[LENGTH]
		currAllSeqKey = seq[SEQKEY]
	# if we are looking for the shortest sequence
	else:
            if seq[UNIQ] == True:
                l = determineShortest(currUniqLen, seq[LENGTH])
		#print 'shortest of %s and %s is %s' % (currUniqLen, seq[LENGTH], l)
                # if seq[LENGTH] is shortest or equal
                if l == 1  or l == -1:
                    currUniqLen = seq[LENGTH]
                    currUniqSeqKey = seq[SEQKEY]
		#print 'currUniqLen: %s, currUniqSeqKey: %s' % (currUniqLen, currUniqSeqKey) 
            # current choice from the full set
            l = determineShortest(currAllLen, seq[LENGTH])
            # if seq[LENGTH] is shortest or equal
            if l == 1  or l == -1:
                currAllLen = seq[LENGTH]
                currAllSeqKey = seq[SEQKEY]
   
    if useUniq == True:
	return (currUniqSeqKey, currUniqLen)
    else:
	return (currAllSeqKey, currAllLen)

# Purpose: determine the longest length
# Returns: 0 if len1 longest, 1 if len2, -1 for tie
# Assumes: Nothing
# Throws: Nothing
def determineLongest (len1, len2): # integer sequence length
    # first comparison for a given seq set, one value will be the default of '' 
    if len1 == '':
	return 1
    elif len2 == '':
	return 0
    # 2nd - n comparisons both will be integers
    if len1 == len2:
	return -1
    elif  len1 > len2:
	return 0
    else:
	return 1

# Purpose: determine the shortest length
# Returns:  0 if len1 shortest, 1 if len2, -1 for tie
# Assumes: Nothing
# Throws: Nothing
def determineShortest (len1, len2): # integer sequence length
    # first comparison for a given seq set, one value will be the default of ''
    if len1 == '':
        return 1
    elif len2 == '':
        return 0
    # 2nd - n comparisons both will be integers
    if len1 == len2:
        return -1
    elif  len1 < len2:
        return 0
    else:
        return 1

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
    #db.sql('set rowcount 10000', None)
    db.sql('select _Marker_key, _Organism_key, _Marker_Type_key ' + \
	'into #markers from MRK_Marker ' + \
	'where _Organism_key in (1, 2, 40, 10, 13)', None)
	#'and _Marker_key in (6005, 6385, 6644)', None)
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

	# VEGA
	if providerKey == 1865333:
	    if allgenomic[VEGA].has_key(m):
		allgenomic[VEGA][m].append(r)
	    else:
		allgenomic[VEGA][m] = [r]

	# Ensembl
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
                allgenomic[gGENBANK][m].append(r)
            else:
                allgenomic[gGENBANK][m] = [r]

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

