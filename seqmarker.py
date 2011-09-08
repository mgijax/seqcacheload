#!/usr/local/bin/python
#
# seqmarker.py  
#####################################################################
#
#  Purpose: This script creates the bcp file for SEQ_Marker_Cache
#           which caches sequence/marker pairs including:
#
#	    A record in SEQ_Marker_Cache represents a uniq sequence, marker
#           reference relationship; there may be multiple references for 
#           a sequence to marker association, therefore multiple records
#           per sequence/marker pair
#
#  See for details: 
#  http://prodwww.informatics.jax.org/wiki/index.php/sw:Seqcacheload#2._Sequence_Marker_Cache_Load
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
# TO DO: 
# 1. string.lower term comparisons, 
# 2. string.strip on incoming equivalency mappings
# 3. update history
# 4. update comments
# 5. remove debug
#
#  History
#
# 09/08/11	sc
#	- TR10308/biotype conflict revision
#
# 04/04/2011	lec
#	- TR10658/add _Cache_key
#
#  03/01/2010   sc
#	- TR9774; update rep tran and prot sequence algorithm for 
#	 markers with VEGA or Ensembl representative Genomic sequences
#
#  02/18/2010	lec
#	- TR 9239; add rawbiotype, _BiotypeConflict_key, _Marker_Type_key
#       - add method to generate biotype lookup
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

#
# CONSTANTS
#

# database errors
DB_ERROR = 'A database error occured: '
DB_CONNECT_ERROR = 'Connection to the database failed: '

# record delimiter
NL = '\n'

# NEW - sequence to sequence qualifiers
# A transcript is transcribed from a genomic sequence
TRANSCRIBED_FROM_KEY = 5445464

# A protein is translated from a transcript sequence
TRANSLATED_FROM_KEY = 5445465

# column delimiter
DL = os.environ['COLDELIM']

# table for which we are creating bcp file
table = os.environ['TABLE']

# output directory
datadir = os.environ['CACHEDATADIR']

# when selecting the representative genomic sequence 
# prints case number, markerKey, four sets of genomic sequences and the
# representative sequence selected
debug = os.environ['SEQMARKER_DEBUG']

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

# biotype lookup by genomic sequence key + marker key
# {seqKey:mkrKey:[conflictType, rawBiotype], ...}
biotypeLookup = {}

# biotype default vocabulary = Not Applicable (_Vocab_key = 76)
biotypeDefaultConflict = '5420769'

# {markerKey:[GeneModel1, ...GeneModelN} 
markerToGMDict = {}

# represents all genomic seqs for the current marker by provider - see indexes
# each dictionary looks like (_Marker_key:{}, ...} 
# where {} is a db.sql result set
allgenomic = [{}, {}, {}, {}]

# indexes of allgenomic
VEGA = 0
ENSEMBL = 1
NCBI = 2
gGENBANK = 3

# genomic sequence provider terms, these are used when the 
# provider is VEGA or Ensembl, to determine the rep transcript and
# protein associated with the gene model.
genbank_prov = 'GENBANK'
ncbi_prov = 'NCBI'
ensembl_prov = 'ENSEMBL'
vega_prov = 'VEGA'

# Lookups from SEQ_Sequence_Assoc to determine relationships 
# between VEGA and Ensembl genomic, transcript, and protein sequences

# Looks like {gKey1:{tKey1:length, tKey2:length, ...},...,gKeyn:{tKeyn:length, ...}, ...}
transcriptLookupByGenomicKey = {}

# Looks like {pKey1:{tKey1:length, tKey2:length, ...},...,pKeyn:{tKeyn:length, ...}, ...}
transcriptLookupByProteinKey = {}

# Looks like {gKey1:{pKey1:length, pKey2:length, ...},...,gKeyn:{pKeyn:length, ...}, ...}
proteinLookupByGenomicKey = {}

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

# next max(_Cache_key)
nextMaxKey = 0

class GeneModel:
    # A representation of an gene model
    # as it applies to determining the biotype conflict
    def __init__(self):

        self.sequenceKey = None
        self.ldbKey = None
	self.rawBiotype = None
        self.equivalentBiotypeSet = None

# Purpose: Initialize db.sql settings, lookups, and file descriptor
# Returns: Nothing
# Assumes: Nothing
# Effects: connects to and queries a database
# Throws: Nothing

def init ():
    global qualByTermLookup, qualByTermKeyLookup, mkrsByGenomicSeqKeyLookup
    global proteinLookupByGenomicKey, transcriptLookupByGenomicKey
    global transcriptLookupByProteinKey, outBCP
    
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
    db.sql('create nonclustered index idx_1 on #gm (seqID)', None)
    db.sql('create nonclustered index idx_2 on #gm (_Sequence_key)', None)

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
    results = db.sql('select s._Sequence_key, ' + \
	'a._Object_key as _Marker_key ' + \
	'from #allSeqs s, ACC_Accession a ' + \
	'where a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key in (59, 60, 85, 9) ' + \
	'and s.seqID = a.accid ' + \
	'order by _Sequence_key ', 'auto')

    # load genomic sequences associated with markers lookup
    prevSeqKey = ''
    for r in results:
	seqKey = r['_Sequence_key']
	markerKey = r['_Marker_key']
        if seqKey != prevSeqKey:
	    mkrsByGenomicSeqKeyLookup[seqKey] = [markerKey]
        else:
	    mkrsByGenomicSeqKeyLookup[seqKey].append(markerKey)
	prevSeqKey = seqKey

    #
    # Load lookups determine relationships between VEGA and Ensembl 
    # genomic, transcript, and protein sequences
    #
    # Load transcriptLookupByGenomicKey 
    db.sql('select sa._Sequence_key_1 as transcriptKey, ' + \
	'ss.length as transcriptLength, ' + \
	'sa._Sequence_key_2 as genomicKey ' + \
	'into #transGen ' + \
	'from SEQ_Sequence_Assoc sa, SEQ_Sequence ss ' + \
	'where sa._Qualifier_key = %s ' % TRANSCRIBED_FROM_KEY + \
	'and sa._Sequence_key_1 = ss._Sequence_key', None)
    db.sql('create nonclustered index idx1 on ' + \
	'#transGen(transcriptKey)', None)
    results = db.sql('select * from #transGen ' + \
	'order by genomicKey', 'auto')
    prevGKey = ''
    for r in results:
	gKey = r['genomicKey']
	tKey = r['transcriptKey']
	tLength = r['transcriptLength']
	if gKey != prevGKey:
	    transcriptLookupByGenomicKey[gKey] = {}
	transcriptLookupByGenomicKey[gKey][tKey] = tLength
        prevGKey = gKey

    # Load proteinLookupByGenomicKey
    results = db.sql('select tg.genomicKey, ' + \
	'sa._Sequence_key_1 as proteinKey, ' + \
	'ss.length as proteinLength ' + \
	'from #transGen tg, SEQ_Sequence_Assoc sa, ' + \
	'SEQ_Sequence ss ' + \
	'where sa._Qualifier_key = %s ' % TRANSLATED_FROM_KEY + \
	'and tg.transcriptKey =  sa._Sequence_key_2 ' + \
	'and sa._Sequence_key_1 = ss._Sequence_key ' + \
	'order by tg.genomicKey', 'auto')

    prevGKey = ''
    for r in results:
        gKey = r['genomicKey']
        pKey = r['proteinKey']
	pLength = r['proteinLength']
        if gKey != prevGKey:
            proteinLookupByGenomicKey[gKey] = {}
	proteinLookupByGenomicKey[gKey][pKey] = pLength
        prevGKey = gKey
    # Load transcriptLookupByProteinKey
    # there should be only one transcript for a protein, but one never knows
    results = db.sql('select sa._Sequence_key_1 as proteinKey, ' + \
                'sa._Sequence_key_2 as transcriptKey, ' + \
		'ss.length as transcriptLength ' + \
                'from SEQ_Sequence_Assoc sa, ' + \
		'SEQ_Sequence ss ' + \
                'where sa._Qualifier_key = %s ' % TRANSLATED_FROM_KEY + \
		'and sa._Sequence_key_2 = ss._Sequence_key ' + \
                'order by sa._Sequence_key_1', 'auto')
    prevPKey = ''
    for r in results:
        pKey = r['proteinKey']
        tKey = r['transcriptKey']
	tLength = r['transcriptLength']
        if pKey != prevPKey:
            transcriptLookupByProteinKey[pKey] = {}
	transcriptLookupByProteinKey[pKey][tKey] = tLength
        prevPKey = pKey

    # generate biotype lookup
    generateBiotypeLookups()

    #
    # create file descriptor for bcp file
    #
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
    if debug == 'true':
	print 'determineRep for marker: %s' % marker
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
	if debug == 'true':
	    print 'vegaSeqs: %s' % vegaSeqs
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
	if debug == 'true':
	    print 'ensemblseqs: %s' % ensemblSeqs
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
	if debug == 'true':
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
	if debug == 'true':
	    print 'genbankseqs: %s' % genbankSeqs
    ##-----------------------------------------------------------------
    # Determine representative genomic for marker using provider
    # and provider sequence attributes
    ##-----------------------------------------------------------------

    genomicRepKey = 0
    genomicRepProvider = ''
    hasSglUniqEnsembl = False
    hasSglUniqNCBI = False

    # marker has no GMs
    if not hasVega and not hasEnsembl and not hasNCBI: 
	# is there a longest uniq
	(s,l) = determineSeq(genbankSeqs, True, True)
	# no sequence found that match parameters if s == 0
	if s != 0:
	    genomicRepKey = s
	    genomicRepProvider = genbank_prov
	    if debug == 'true':
		print 'CASE 1'
	# if no longest uniq get longest
	if genomicRepKey == 0:
	    (s,l) = determineSeq(genbankSeqs, True, False)
	    if s != 0:
		genomicRepKey = s
		genomicRepProvider = genbank_prov
		if debug == 'true':
		    print 'CASE 2'
	# else NO REPRESENTATIVE SEQUENCE value of genomicRepKey is still 0

    # marker has at least one GM, therefore WILL HAVE REP SEQUENCE
    else: 
	if vegaIsSgl and vegaHasUniq:
	    genomicRepKey = vegaSeqs[0][SEQKEY]
	    genomicRepProvider = vega_prov
	    if debug == 'true':
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
		    genomicRepKey = ncbiSeqs[0][SEQKEY]
		    genomicRepProvider = ncbi_prov
		    if debug == 'true':
			print 'CASE 4'
		else:
		    genomicRepKey =  ensemblSeqs[0][SEQKEY]
		    genomicRepProvider = ensembl_prov
		    if debug == 'true':
			print 'CASE 5'
	    elif hasSglUniqEnsembl:
		genomicRepKey = ensemblSeqs[0][SEQKEY]
		genomicRepProvider = ensembl_prov
		if debug == 'true':
		    print 'CASE 6'
	    elif hasSglUniqNCBI:
		genomicRepKey = ncbiSeqs[0][SEQKEY]
		genomicRepProvider = ncbi_prov
		if debug == 'true':
		    print 'CASE 7'
	    # only multiples (uniq or not) or single not-uniq left
	    else:
		# check uniq (must be multiple)
                if vegaHasUniq or ensemblHasUniq or ncbiHasUniq:
		    if vegaHasUniq:
			# pick shortest uniq, if tie pick one
			(s,l) = determineSeq(vegaSeqs, False, True)
			genomicRepKey = s
			genomicRepProvider = vega_prov
			if debug == 'true':
			    print 'CASE 8'
		    elif ensemblHasUniq and ncbiHasUniq:
			# pick shortest uniq, if tie pick one
			(s_e, l_e) = determineSeq(ensemblSeqs, False, True)
			(s_n, l_n) = determineSeq(ncbiSeqs, False, True)
			value = determineShortest(l_e, l_n)
			if value == -1 or value == 1:
			    genomicRepKey = s_n
			    genomicRepProvider = ncbi_prov
			    if debug == 'true':
				print 'CASE 9'
			else:
			    genomicRepKey = s_e
			    genomicRepProvider = ensembl_prov	
			    if debug == 'true':
				print 'CASE 10'
		    elif ensemblHasUniq:
			(s_e, l_e) = determineSeq(ensemblSeqs, False, True)
			genomicRepKey = s_e
			genomicRepProvider = ensembl_prov
			if debug == 'true':
			    print 'CASE 11'
		    elif ncbiHasUniq:
			(s_n, l_n) = determineSeq(ncbiSeqs, False, True)
			genomicRepKey = s_n
			genomicRepProvider = ncbi_prov
			if debug == 'true':
			    print 'CASE 12'
		# no uniques, only single or multiple non-uniq left
		else:
		    # check for vega sgl
		    if vegaIsSgl:
			genomicRepKey = vegaSeqs[0][SEQKEY]
			genomicRepProvider = vega_prov
			if debug == 'true':
			    print 'CASE 13'
		    else:
			# check for ensembl and ncbi sgl
			if ensemblIsSgl or ncbiIsSgl:
			    if ensemblIsSgl and ncbiIsSgl:
				value = determineShortest( \
				    ensemblSeqs[0][LENGTH], \
				    ncbiSeqs[0][LENGTH])
				if value == -1 or value == 1:
				    genomicRepKey = ncbiSeqs[0][SEQKEY]
				    genomicRepProvider = ncbi_prov
				    if debug == 'true':
					print 'CASE 14'
				else:
				    genomicRepKey =  ensemblSeqs[0][SEQKEY]
				    genomicRepProvider = ensembl_prov
			 	    if debug == 'true':
					print 'CASE 15'
			    elif ensemblIsSgl:
				genomicRepKey = ensemblSeqs[0][SEQKEY]
				genomicRepProvider = ensembl_prov
			 	if debug == 'true':
				    print 'CASE 16'
			    elif ncbiIsSgl:
				genomicRepKey = ncbiSeqs[0][SEQKEY]
				genomicRepProvider = ncbi_prov
				if debug == 'true':
				    print 'CASE 17'
			# no singles, must be multiple non-uniq
			else:
			    if hasVega:
                                # pick shortest
				(s,l) = determineSeq(vegaSeqs, False, False)
				genomicRepKey = s
				genomicRepProvider = vega_prov
				if debug == 'true':
				    print 'CASE 18'
                            elif hasEnsembl and hasNCBI:
                                # pick shortest, NCBI if tie
				(s_e, l_e) = determineSeq( \
				    ensemblSeqs, False, False)
				(s_n, l_n) = determineSeq( \
				    ncbiSeqs, False, False)
				value = determineShortest(l_e, l_n)
				if value == -1 or value == 1:
				    genomicRepKey = s_n
				    genomicRepProvider = ncbi_prov
				    if debug == 'true':
					print 'CASE 19'
				else:
				    genomicRepKey =  s_e
				    genomicRepProvider = ensembl_prov
				    if debug == 'true':
					print 'CASE 20'
                            elif hasEnsembl:
                                # pick shortest
				(s_e, l_e) = determineSeq( \
				    ensemblSeqs, False, False)
				genomicRepKey = s_e
				genomicRepProvider = ensembl_prov
				if debug == 'true':
				    print 'CASE 21'
			    # must be an NCBI
                            else:
                                # pick shortest NCBI
				(s_n, l_n) = determineSeq( \
                                    ncbiSeqs, False, False)
				genomicRepKey = s_n
				genomicRepProvider = ncbi_prov
				if debug == 'true':
				    print 'CASE 22'
    # if we found a genomicRepKey for this marker add it to the dictionary
    if debug == 'true':
	print 'genomicRepKey: %s' % genomicRepKey
	print 'genomicRepProvider: %s' % genomicRepProvider
    if genomicRepKey != 0:
        genomic[marker] = genomicRepKey
    else:
	if debug == 'true':
	    print 'CASE 23'

    #
    # Determine Representative Protein and Transcript Sequences
    #
    if genomicRepProvider == ensembl_prov or \
	    genomicRepProvider == vega_prov:
	determineVegaEnsProtTransRep(marker, genomicRepKey)
    else: # not VEGA or Ensembl
	determineNonVegaEnsProtRep(marker)
        determineNonVegaEnsTransRep(marker)
    return

def determineVegaEnsProtTransRep(marker, genomicRepKey):
    # Purpose: Determine the representative protein and transcript
    #     for 'marker'. When the rep genomic is VEGA or Ensembl
    #     the protein and transcript must be from same provider
    #     if they exist
    # Returns: nothing
    # Assumes: nothing
    # Effects: nothing
    # Throws: nothing

    global allpolypeptide, polypeptide, alltranscript, transcript
    global transcriptLookupByGenomicKey, transcriptLookupByProteinKey

    protRepKey = 0 	# default
    transRepKey = 0 	# default

    #
    # we determine the reprentative protein first per requirements -
    # see TR9774
    #
    if not proteinLookupByGenomicKey.has_key(genomicRepKey):
	# no prots for this genomic, get rep protein in the usual way
	determineNonVegaEnsProtRep(marker)

        # now get Vega/Ensembl transcript(s) for the genomicRepKey
	if transcriptLookupByGenomicKey.has_key(genomicRepKey):
	    # get the list of transcriptIds mapped to their length
	    #transDict- {t1:length, t2:length, ...}
	    transDict = transcriptLookupByGenomicKey[genomicRepKey]
	    # length of current longest transcript
            currentLongestTransLen = 0
            # Now determine the longest transcript
            for tKey in transDict.keys():
                tLength = transDict[tKey]
                if tLength > currentLongestTransLen:
                    currentLongestTransLen = tLength
                    transRepKey = tKey
	    if transRepKey != 0:
		transcript[marker] = transRepKey
	    else:
		print "This shouldn't happen 1"
		sys.exit("This shouldn't happen 1")
	else: # no trans for the genomic, get rep trans in the usual way
	    determineNonVegaEnsTransRep(marker)

    else: # there are proteins for the genomicRepKey, determine longest
	protDict = proteinLookupByGenomicKey[genomicRepKey]
	#print 'Proteins for genomicRepKey %s: %s' % (genomicRepKey, protDict)
	# length of current longest polypeptide
	currentLongestProtLen = 0
	protRepKey = 0
	# determine the longest polypeptide
	for pKey in protDict.keys():
	    pLength = protDict[pKey]
	    #print 'pKey: %s pLength %s type: %s currentLongestProtLen %s type: %s' % (pKey, pLength, type(pLength), currentLongestProtLen, type(currentLongestProtLen))
	    if int(pLength) > int(currentLongestProtLen):
		currentLongestProtLen = pLength
		protRepKey = pKey
	#print 'chosen LongestProtLen: %s protRepKey: %s' % (currentLongestProtLen, protRepKey)
	if protRepKey != 0:
	    polypeptide[marker] = protRepKey
	else: 
	    print "This shouldn't happen 2"
	    sys.exit("This shouldn't happen 2")
	# now get Vega/Ensembl transcript(s) for the protRepKey
	if transcriptLookupByProteinKey.has_key(protRepKey):
	    transDict = transcriptLookupByProteinKey[protRepKey]
	    # length of current longest transcript
	    currentLongestTransLen = 0
	    transRepKey = 0
	    # determine the longest transcript
	    for tKey in transDict.keys():
		tLength = transDict[tKey]
	        if tLength > currentLongestTransLen:
		    currentLongestTransLen = tLength
		    transRepKey = tKey
	    if transRepKey  != 0:
		transcript[marker] = transRepKey
	else:   # no VEGA or Ensembl protein i.e.
		# we have a protein w/o a transcript
	    print "This shouldn't happen 3"
            sys.exit("This shouldn't happen 3")
    return

def determineNonVegaEnsProtRep(marker):
    # Purpose: determine non-VEGA, non-Ensembl rep protein
    # Returns: nothing
    # Assumes: nothing
    # Effects: nothing
    # Throws: nothing

    global allpolypeptide, polypeptide
    for i in range(len(allpolypeptide)):
	if allpolypeptide[i].has_key(marker):
	    polypeptide[marker] = allpolypeptide[i][marker]
	    return
    return

def determineNonVegaEnsTransRep(marker):
    # Purpose: determine non-VEGA, non-Ensembl rep transcript
    # Returns: nothing
    # Assumes: nothing
    # Effects: nothing
    # Throws: nothing

    global alltranscript, transcript
    for i in range(len(alltranscript)):
	if alltranscript[i].has_key(marker):
	    #transcript[marker] = []
	    # why is this different then loading polypeptide i.e. 
	    # why append to list
	    #transcript[marker].append(alltranscript[i][marker])
	    # changed to be the same as Prot
	    transcript[marker] = alltranscript[i][marker]
	    return
    return

def determineSeq(seqList, 	# list of dictionaries
                getLongest,     # boolean, determine longest if True, else
                                # shortest
		useUniq): 	# boolean, consider uniq only if True
    # Purpose: Find the longest or shortest sequence given a dictionary like:
    #          [{UNIQ:True/False, SEQKEY:key, LENGTH:length}, ...]
    # Returns: tuple (seqKey, length) where sequence key is the longest/shortest
    #          uniq/notuniq depending on value of 'getLongest' and 'useUniq',
    #          else both members of the tuple are 0
    # Assumes: Nothing
    # Throws: Nothing
    #

    # current choice considering only uniq sequences
    # based on value of getLongest

    # a non-numeric default
    currUniqLen = ''
    currUniqSeqKey = 0
    # current choice considering all sequences based on value of getLongest
    currAllLen = ''
    currAllSeqKey = 0
    if debug == 'true':
	print 'getLongest: %s, useUniq: %s' % (getLongest, useUniq)
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
                # if seq[LENGTH] is shortest or equal
                if l == 1  or l == -1:
                    currUniqLen = seq[LENGTH]
                    currUniqSeqKey = seq[SEQKEY]
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

def determineLongest (len1, len2): # integer sequence length
    # Purpose: determine the longest length
    # Returns: 0 if len1 longest, 1 if len2, -1 for tie
    # Assumes: Nothing
    # Throws: Nothing

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

def determineShortest (len1, len2): # integer sequence length
    # Purpose: determine the shortest length
    # Returns:  0 if len1 shortest, 1 if len2, -1 for tie
    # Assumes: Nothing
    # Throws: Nothing

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


def generateBiotypeLookups():
    # Purpose: create lookups for use in determining biotype conflicts
    # for those markers associated with Gene Models. Create global lookup
    # 'biotypeLookup' for use in writing conflict attributes to bcp file
    # Returns:  Nothing
    # Assumes: Nothing
    # Throws: Nothing

    global biotypeLookup

    # do all comparisons in lower case
    NC_RNA_CONFIG_TERM=string.lower(string.strip(os.environ['NC_RNA_TERM']))
    ALL_FEATURES_CONFIG_TERM=string.lower(string.strip(os.environ['ALL_FEATURES_TERM']))

    # non-coding RNA gene feature types and its descendents
    ncRNAdescSet = set([])

    # all feature types descendents
    allFeatureTypesDescSet = set([])

    # map mcv feature type terms from the VOC_Term table to their keys
    mcvTermToKeyDict = {}

    # markers mapped to their MGI feature type. 9/8/11 - still currently only 1
    # {markerKey:[featureType1, ...featureTypeN}
    featureTypesDict = {}

    # provider raw biotypes mapped to their set of equivalent terms
    # equivalency dicts look like {rawTerm:[listOfEquivalentTerms], ...}
    NCBIEquivDict = {}
    EnsEquivDict = {}
    VEGAEquivDict = {}

    # conflict types (see _Vocab_key = 76)
    yesConflict = 5420767
    noConflict = 5420769

    print 'Initializing Biotype Lookups ... %s' % (mgi_utils.date())

    #
    # create set of non-coding RNA gene and its descendent terms
    #
    nonCodingRNAGeneTermKey = 6238162
    
    results = db.sql('''select t.term as descTerm, c._AncestorObject_key,
                c._DescendentObject_key
            from DAG_Closure c, VOC_Term t
            where c._DAG_key = 9
                and c._MGIType_key = 13
                and _AncestorObject_key = %s
                and c._DescendentObject_key = t._Term_key
            order by  c._AncestorObject_key, c._DescendentObject_key''' % \
                nonCodingRNAGeneTermKey, 'auto')
    # add the term itself
    ncRNAdescSet.add(nonCodingRNAGeneTermKey)
    for r in results:
        ncRNAdescSet.add(r['_DescendentObject_key'])

    #
    # create set of all feature types descendent terms
    #
    allFeatureTypesTermKey = 6238159

    results = db.sql('''select t.term as descTerm, c._AncestorObject_key,
                c._DescendentObject_key
            from DAG_Closure c, VOC_Term t
            where c._DAG_key = 9
                and c._MGIType_key = 13
                and _AncestorObject_key = %s
                and c._DescendentObject_key = t._Term_key
            order by  c._AncestorObject_key, c._DescendentObject_key''' % \
                allFeatureTypesTermKey, 'auto')
    for r in results:
        allFeatureTypesDescSet.add(r['_DescendentObject_key'])

    #
    # map  all feature type terms to their keys
    #
    results = db.sql('''select _Term_key, term
		from VOC_Term
		where _Vocab_key = 79''', 'auto')
    for r in results:
	mcvTermToKeyDict[string.lower(r['term'])] = r['_Term_key']

    #
    # create NCBI, Ensembl and VEGA equivalency Lookups 
    #
    # raw term/equivalency sets split on ','
    # raw term split from equivalent terms using ':'
    # equivalent terms split on '|'
    # VEGA example:
    # ..,KNOWN_polymorphic:protein coding gene,pseudogene:pseudogenic region|pseudogene|pseudogenic gene segment,.

    print 'Initializing NCBI raw biotype to equivalency mapping ... %s' % (mgi_utils.date())
    ncbiEquiv = string.lower(os.environ['NCBI_EQUIV'])
    mappingList = string.split(ncbiEquiv, ',')

    for m in mappingList:
        rawList = string.split(m, ':')
        raw = string.strip(rawList[0])
        equivList = string.split(rawList[1], '|')
	equivKeySet = set()
	for e in equivList:
	    e = string.strip(e)
	    if e == ALL_FEATURES_CONFIG_TERM:
		equivKeySet = equivKeySet.union(allFeatureTypesDescSet)
	    elif e == NC_RNA_CONFIG_TERM:
		equivKeySet = equivKeySet.union(ncRNAdescSet)
	    elif mcvTermToKeyDict.has_key(e):
                equivKeySet.add(mcvTermToKeyDict[e])
	    else:
		sys.exit('NCBI equivalency term does not resolve: %s' % e)
        NCBIEquivDict[raw] = equivKeySet

    print 'Initializing Ensembl raw biotype to equivalency mapping ... %s' % (mgi_utils.date())
    ensEquiv = string.lower(os.environ['ENSEMBL_EQUIV'])
    mappingList = string.split(ensEquiv, ',')

    for m in mappingList:
        rawList = string.split(m, ':')
        raw = string.strip(rawList[0])
        equivList = string.split(rawList[1], '|')
	equivKeySet = set()
	for e in equivList:
	    e = string.strip(e)
            if mcvTermToKeyDict.has_key(e):
                equivKeySet.add(mcvTermToKeyDict[e])
	    elif e == NC_RNA_CONFIG_TERM:
                equivKeySet = equivKeySet.union(ncRNAdescSet)
            else:
                sys.exit('Ensembl equivalency term does not resolve: %s' % e)
        EnsEquivDict[raw] = equivKeySet

    print 'Initializing VEGA raw biotype to equivalency mapping ... %s' % (mgi_utils.date())
    vegaEquiv = string.lower(os.environ['VEGA_EQUIV'])
    mappingList = string.split(vegaEquiv, ',')

    for m in mappingList:
	rawList = string.split(m, ':')
	raw = string.strip(rawList[0])
	equivList = string.split(rawList[1], '|')
	equivKeySet = set()
        for e in equivList:
	    e = string.strip(e)
            if mcvTermToKeyDict.has_key(e):
                equivKeySet.add(mcvTermToKeyDict[e])
	    elif e == NC_RNA_CONFIG_TERM:
		equivKeySet = equivKeySet.union(ncRNAdescSet)
            else:
                sys.exit('VEGA equivalency term does not resolve: %s' % e)
        VEGAEquivDict[raw] = equivKeySet

    #
    #   for each Marker associated with a NCBI (59), Ensembl (60), VEGA (85) 
    #       gene model sequence:
    #     get the gene model raw biotype and map it to its set of equivalent
    #     mcv terms
    #
    print 'Initializing  gene model lookup by marker key ... %s' % (mgi_utils.date())
    results = db.sql('''select s._Sequence_key, a._Object_key as _Marker_key,
		a._LogicalDB_key, g.rawBiotype
	 from #gm s, ACC_Accession a, SEQ_GeneModel g
	 where a._MGIType_key = 2
	 and a._LogicalDB_key in (59, 60, 85)
	 and s.seqID = a.accid
	 and s._Sequence_key = g._Sequence_key
	 order by s._Sequence_key
	 ''', 'auto')
	 #and m._Marker_key = 332286

    for r in results:
	markerKey = r['_Marker_key']
	ldbKey = r['_LogicalDB_key']
	sequenceKey = r['_Sequence_key']
	rawBiotype = r['rawBiotype']
	if rawBiotype == None:
            rawBiotype = 'null'
	currentEquivSet = set()

	# equivalencies are in lower case, so compare with lower biotype
	lowerRawBiotype = string.lower(rawBiotype)
	if ldbKey == 59 and NCBIEquivDict.has_key(lowerRawBiotype):
	    currentEquivSet =  NCBIEquivDict[lowerRawBiotype]
	elif ldbKey ==  60 and EnsEquivDict.has_key(lowerRawBiotype):
	    currentEquivSet = EnsEquivDict[lowerRawBiotype]
	elif ldbKey == 85 and VEGAEquivDict.has_key(lowerRawBiotype):
	    currentEquivSet = VEGAEquivDict[lowerRawBiotype]
	else:
	    print 'No equivalency set for sequenceKey: %s, ldbKey: %s, rawBiotype: %s' % (
		sequenceKey, ldbKey, rawBiotype)
	    continue

	# create a GeneModel object; map marker key to the GM object 
	gm = GeneModel()
	gm.sequenceKey = sequenceKey
	gm.ldbKey = ldbKey
	gm.rawBiotype = rawBiotype
	gm.equivalentBiotypeSet = currentEquivSet
	if not markerToGMDict.has_key(markerKey):
	     markerToGMDict[markerKey] = []
	markerToGMDict[markerKey].append(gm)

    #
    #   for each Marker in MRK_MCV_Cache 
    #        get the direct term(s) and map them to the marker
    #
    print 'Initializing MGI Marker feature type lookup ... %s' % (mgi_utils.date())
    results = db.sql('''
	select _Marker_key, _MCVTerm_key
	from MRK_MCV_Cache
	where qualifier = 'D' ''', 'auto')

    for r in results:
	key = r['_Marker_key']
        value = r['_MCVTerm_key']
        if not featureTypesDict.has_key(key):
	    featureTypesDict[key] = []
	featureTypesDict[key].append(value)

    # now iterate through the MRK_MCV_Cache markers
    print 'Iterating through all markers in MRK_MCV_Cache to determine conflicts ... %s' % (mgi_utils.date())
    markerKeys = featureTypesDict.keys()
    markerKeys.sort()

    for markerKey in markerKeys:
	# default conflict type
	conflictType = noConflict

	# get the list of MGI feature types for the current marker
	mkrFeatureTypeSet = set(featureTypesDict[markerKey])
	equivalencyList = []
	first = True
	gmIntersectSet = set()

	if markerToGMDict.has_key(markerKey):
	    for gm in markerToGMDict[markerKey]:
		s = gm.equivalentBiotypeSet
		if first == True:
		    gmIntersectSet = s
		    first = False
		else:
		    gmIntersectSet = gmIntersectSet.intersection(s)
	else:
	     # No gene models, we can move on to the next marker
	     continue

	# there are gene models so check for conflict
	# if the gene model set is empty that means conflicts btwn gene models
	if len(gmIntersectSet) == 0:
	    conflictType = yesConflict
	# otherwise the gene models agree, see if they agree with the marker
	else:
	    finalIntersectSet = mkrFeatureTypeSet.intersection(gmIntersectSet)
	    if len(finalIntersectSet) != 1:
		conflictType = yesConflict

	# now re-iterate thru the marker/sequences
	# and set the conflict key and raw biotype
	# all sequences for a given marker get the same conflict key value
	if markerToGMDict.has_key(markerKey):
	    for gm in  markerToGMDict[markerKey]:
		rawBiotype = gm.rawBiotype
		sequenceKey = gm.sequenceKey
		key = '%s:%s' % (markerKey, sequenceKey)
		value = [conflictType, rawBiotype]
		if not biotypeLookup.has_key(key):
		    biotypeLookup[key] = value
    return

def writeRecord(r):
    # Purpose: formats and writes out record to bcp file
    # Returns: Nothing
    # Assumes: Nothing
    # Effects: writes a record to a file
    # Throws: Nothing

    global nextMaxKey

    nextMaxKey = nextMaxKey + 1

    outBCP.write(
	str(nextMaxKey) + DL + \
	mgi_utils.prvalue(r['_Sequence_key']) + DL + \
	mgi_utils.prvalue(r['_Marker_key']) + DL + \
	mgi_utils.prvalue(r['_Organism_key']) + DL + \
	mgi_utils.prvalue(r['_Refs_key']) + DL)

    printedQualifier = 0
    if genomic.has_key(r['_Marker_key']):
	if genomic[r['_Marker_key']] == r['_Sequence_key']:
	    outBCP.write(mgi_utils.prvalue(qualByTermLookup['genomic']) + DL)
	    printedQualifier = 1

    if transcript.has_key(r['_Marker_key']):
	# used to be transcript[markerKey] used to be a list of one now 
	# just seqKey
	if r['_Sequence_key'] == transcript[r['_Marker_key']]:
	    outBCP.write(mgi_utils.prvalue(qualByTermLookup['transcript']) + DL)
	    printedQualifier = 1

    if polypeptide.has_key(r['_Marker_key']):
        if polypeptide[r['_Marker_key']] == r['_Sequence_key']:
	    outBCP.write(mgi_utils.prvalue(qualByTermLookup['polypeptide']) + \
		DL)
	    printedQualifier = 1

    if not printedQualifier:
        outBCP.write(mgi_utils.prvalue(qualByTermLookup['Not Specified']) + DL)

    #
    # get the biotype information from biotypeLookup
    # use defaults if there is no biotype record for this marker/sequence
    #
    biotypeKey = '%s:%s' % (r['_Marker_key'], r['_Sequence_key'])
    if biotypeLookup.has_key(biotypeKey):
	biotypeConflict = biotypeLookup[biotypeKey][0]
	biotypeRaw = biotypeLookup[biotypeKey][1]
    else:
        biotypeConflict = biotypeDefaultConflict
        biotypeRaw = None

    outBCP.write(mgi_utils.prvalue(r['_SequenceProvider_key']) + DL + \
	mgi_utils.prvalue(r['_SequenceType_key']) + DL + \
	mgi_utils.prvalue(r['_LogicalDB_key']) + DL + \
	mgi_utils.prvalue(r['_Marker_Type_key']) + DL + \
	mgi_utils.prvalue(biotypeConflict) + DL + \
	r['accID'] + DL + \
	mgi_utils.prvalue(biotypeRaw) + DL + \
	r['mdate'] + DL + \
        mgi_utils.prvalue(r['_User_key']) + DL + \
	mgi_utils.prvalue(r['_User_key']) + DL + \
        loaddate + DL + loaddate + NL)
    return

def createBCP():
    # Purpose: Iterates through result set of sequence marker pairs
    #          determining the representative sequence qualifier for each
    #          (genomic, transcript, protein, or none of the above), writing
    #          pairs out to bcp file
    # Returns: Nothing
    # Assumes: Nothing
    # Effects: queries a database, writes records to a file
    # Throws: Nothing

    global allgenomic, alltranscript, allpolypeptide
    global outBCP

    print 'Processing ...%s' % (mgi_utils.date())
    #	
    # select only mouse, human, rat, dog, chimpanzee & cattle markers
    # with ANY marker status 
    #db.sql('set rowcount 10000', None)
    db.sql('select _Marker_key, _Organism_key, _Marker_Type_key ' + \
	'into #markers from MRK_Marker ' + \
	'where _Organism_key in (1, 2, 40, 10, 13, 11) ', None)
	#'where _Organism_key = 1 ' + \
	#' and _Marker_key in (12179)', None)
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
    db.sql('create nonclustered index idx2 on #markerAccs ' + \
	'(accID)', None)
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
	'_Organism_key, _Marker_Type_key, _SequenceProvider_key, ' + \
	'_SequenceType_key, _LogicalDB_key, _Refs_key, ' + \
	'_User_key, mdate, accID ' + \
	'from #finalannot', 'auto')
    
    # results are ordered by  _Sequence_key, _Marker_key, _Refs_key
    for r in results:
	writeRecord(r)
    return

def finalize():
    # Purpose: Perform cleanup steps for the script.
    # Returns: Nothing
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing

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

