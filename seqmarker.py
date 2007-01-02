#!/usr/local/bin/python

'''
#
# Purpose:
#
# Create bcp file for SEQ_Marker_Cache
#
# Uses environment variables to determine Server and Database
# (DSQUERY and MGD).
#
# Usage:
#	seqmarker.py
#
# History
#
# 02/02/2007	lec
#	- TR 8072; exclude deleted sequences from representative algorithm
#
# 06/28/2006	lec
#	- add dog & chimpanzee (TR 7508)
#
# 10/12/2005	lec
#	- add primary acc id
#
# 09/08/2005	lec
#	- PIRSF (TR 5972)
#	- add human/rat
#
# 10/23/2003	lec
#	- new (TR 3404, JSAM)
#
'''

import sys
import os
import string
import db
import mgi_utils
import loadlib

NL = '\n'
DL = os.environ['COLDELIM']
table = os.environ['TABLE']
datadir = os.environ['CACHEDATADIR']
loaddate = loadlib.loaddate

outBCP = None
qualifiers1 = {}
qualifiers2 = {}
genomic = {}
transcript = {}
polypeptide = {}

def writeRecord(r):

    outBCP.write(mgi_utils.prvalue(r['sequenceKey']) + DL + \
	       	mgi_utils.prvalue(r['markerKey']) + DL + \
	       	mgi_utils.prvalue(r['organismKey']) + DL + \
	       	mgi_utils.prvalue(r['refsKey']) + DL)

    printedQualifier = 0
    if genomic.has_key(r['markerKey']):
	if genomic[r['markerKey']] == r['sequenceKey']:
	    outBCP.write(mgi_utils.prvalue(qualifiers1['genomic']) + DL)
	    printedQualifier = 1

    if transcript.has_key(r['markerKey']):
       if r['sequenceKey'] in transcript[r['markerKey']]:
	    outBCP.write(mgi_utils.prvalue(qualifiers1['transcript']) + DL)
	    printedQualifier = 1

    if polypeptide.has_key(r['markerKey']):
        if polypeptide[r['markerKey']] == r['sequenceKey']:
	    outBCP.write(mgi_utils.prvalue(qualifiers1['polypeptide']) + DL)
	    printedQualifier = 1

    if not printedQualifier:
        outBCP.write(mgi_utils.prvalue(qualifiers1['Not Specified']) + DL)

    outBCP.write(mgi_utils.prvalue(r['providerKey']) + DL + \
	mgi_utils.prvalue(r['typeKey']) + DL + \
	mgi_utils.prvalue(r['logicalDBKey']) + DL + \
	r['accID'] + DL + \
	r['mdate'] + DL + \
        mgi_utils.prvalue(r['userKey']) + DL + mgi_utils.prvalue(r['userKey']) + DL + \
        loaddate + DL + loaddate + NL)

def createBCP():
	global qualifiers1, qualifiers2
	global genomic, transcript, polypeptide
	global outBCP

	outBCP = open('%s/%s.bcp' % (datadir, table), 'w')

	results = db.sql('select _Term_key, term from VOC_Term_RepQualifier_View', 'auto')
	for r in results:
	   qualifiers1[r['term']] = r['_Term_key']
	   qualifiers2[r['_Term_key']] = r['term']

	# only mouse, human, rat, dog & chimpanzee markers

	db.sql('select _Marker_key, _Organism_key, _Marker_Type_key into #markers from MRK_Marker ' + \
		'where _Organism_key in (1, 2, 40, 10, 13) and _Marker_Status_key in (1,3)', None)
#		'where _Organism_key in (1, 2, 40, 10, 13) and _Marker_Status_key in (1,3) and _Marker_key = 31999', None)
	db.sql('create nonclustered index idx_key on #markers (_Marker_key)', None)

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

	db.sql('create nonclustered index idx1 on #markerAccs (_LogicalDB_key, accID)', None)

	# select all sequence annotations

	db.sql('select sequenceKey = s._Object_key, ' + \
		'markerKey = m._Marker_key, organismKey = m._Organism_key, markerType = m._Marker_Type_key, ' + \
		'logicalDBKey = m._LogicalDB_key, ' + \
		'refsKey = m._Refs_key, userKey = m._ModifiedBy_key, m.mdate, m.accID ' + \
		'into #preallannot ' + \
		'from #markerAccs m, ACC_Accession s ' + \
		'where m.accID = s.accID ' + \
		'and m._LogicalDB_key = s._LogicalDB_key ' + \
		'and s._MGIType_key = 19 ', None)

	db.sql('create nonclustered index idx1 on #preallannot (sequenceKey)', None)

	db.sql('select m.sequenceKey, m.markerKey, m.organismKey, m.markerType, ' + \
		'providerKey = ss._SequenceProvider_key, typeKey = ss._SequenceType_key, ' + \
		'm.logicalDBKey, m.refsKey, m.userKey, m.mdate, m.accID ' + \
		'into #allannot ' + \
		'from #preallannot m, SEQ_Sequence ss ' + \
		'where m.sequenceKey = ss._Sequence_key', None)

	db.sql('create nonclustered index idx1 on #allannot (sequenceKey)', None)

	# select annotations to all sequences; grab sequence's primary accID

	db.sql('select a.sequenceKey, a.markerKey, a.organismKey, a.markerType, ' + \
		'a.providerKey, a.typeKey, a.logicalDBKey, a.refsKey, ' + \
		'a.userKey, a.mdate, ac.accID ' + \
		'into #allseqannot ' + \
		'from #allannot a, ACC_Accession ac ' + \
		'where a.sequenceKey = ac._Object_key ' + \
		'and ac._MGIType_key = 19 ' + \
		'and ac.preferred = 1', None)

	db.sql('create nonclustered index idx1 on #allseqannot (sequenceKey, markerKey, refsKey)', None)

	# select records, grouping by sequence, marker and reference

	db.sql('select sequenceKey, markerKey, organismKey, markerType, ' + \
		'providerKey, typeKey, logicalDBKey, refsKey, userKey, mdate = max(mdate), accID ' + 
		'into #finalannot ' + \
		'from #allseqannot group by sequenceKey, markerKey, refsKey', None)
	db.sql('create nonclustered index idx1 on #finalannot (sequenceKey, markerKey, refsKey, userKey, mdate)', None)
	db.sql('create nonclustered index idx2 on #finalannot (sequenceKey, markerKey, markerType, accID)', None)
	db.sql('create nonclustered index idx3 on #finalannot (markerKey)', None)

	db.sql('select distinct sequenceKey, markerKey, markerType, accID into #deriveQuality ' + \
		'from #finalannot order by markerKey', None)
	db.sql('create nonclustered index idx1 on #deriveQuality (sequenceKey)', None)
	db.sql('create nonclustered index idx2 on #deriveQuality (markerKey)', None)

	allgenomic = [{}, {}]
	alltranscript = [{}, {}, {}, {}, {}, {}, {}]
	allpolypeptide = [{}, {}, {}, {}]
        multitranscript = []

	prevMarker = ''

	# process derived representative values; do not include deleted sequences

	results = db.sql('select q.sequenceKey, q.markerKey, q.markerType, q.accID, ' + \
		's._SequenceProvider_key, s._SequenceType_key, s.length ' + \
		'from #deriveQuality q, SEQ_Sequence s ' + \
		'where q.sequenceKey = s._Sequence_key ' + \
		'and s._SequenceStatus_key != 316343 ' + \
		'order by q.markerKey, s._SequenceProvider_key', 'auto')

	for r in results:

	    m = r['markerKey']
	    s = r['sequenceKey']
	    a = r['accID']
	    markerType = r['markerType']

	    if r['length'] is None:
		seqlength = 0
	    else:
		seqlength = int(r['length'])

	    providerKey = r['_SequenceProvider_key']
	    seqTypeKey = r['_SequenceType_key']

	    if prevMarker != m:

		glengths = [-1,-1,-1]
		tlengths = [-1,-1,-1,-1,-1,-1,-1]
		plengths = [-1,-1,-1,-1]

		if prevMarker != '':

		    # determine the one and only representative w/in each group

		    for i in range(2):
			if allgenomic[i].has_key(prevMarker):
			    genomic[prevMarker] = allgenomic[i][prevMarker]
			    break

		    for i in range(7):
			if alltranscript[i].has_key(prevMarker):
			    transcript[prevMarker] = []
			    transcript[prevMarker].append(alltranscript[i][prevMarker])
			    break

		    for i in range(4):
			if allpolypeptide[i].has_key(prevMarker):
			    polypeptide[prevMarker] = allpolypeptide[i][prevMarker]
			    break

		    # for markers that can have mulitple representative transcripts
		    if len(multitranscript) > 0:
			transcript[prevMarker] = multitranscript
                        multitranscript = []

	    #
	    # representative genomic
	    #
	    # longest VEGA
	    # longest NCBI/Ensembl coordinate OR longest GenBank DNA (tie goes to NCBI)
	    #

	    # VEGA
	    if providerKey == 1865333:
		if seqlength > glengths[2]:
		    allgenomic[0][m] = s
	            glengths[2] = seqlength

	    # NCBI, Ensembl
	    elif (providerKey in [706915,615429]) and seqlength > glengths[0]:
		allgenomic[0][m] = s
	        glengths[0] = seqlength

	    # NCBI
	    elif providerKey == 706915 and seqlength == glengths[0]:
		allgenomic[0][m] = s
	        glengths[0] = seqlength

	    # any GenBank; DNA
	    elif providerKey in [316380,316376,316379,316375,316377,316374,316373,316378,492451] and seqTypeKey == 316347 and seqlength > glengths[1]:
		allgenomic[1][m] = s
	        glengths[1] = seqlength

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
	    if providerKey == 316372 and (string.find(a, 'NM_') > -1 or string.find(a, 'NR_') > -1):
		if seqlength > tlengths[0]:
		    alltranscript[0][m] = s
	            tlengths[0] = seqlength

	    # GenBank but not EST; RNA
	    elif providerKey in [316380,316379,316375,316377,316374,316373,316378,492451] and seqTypeKey == 316346:
		if seqlength > tlengths[1]:
		    alltranscript[1][m] = s
	            tlengths[1] = seqlength

	    # RefSeq
	    elif providerKey == 316372 and (string.find(a, 'XM_') > -1 or string.find(a, 'XR_') > -1):
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
	    elif providerKey == 316372 and string.find(a, 'NP_') > -1 and seqlength > plengths[1]:
		allpolypeptide[1][m] = s
	        plengths[1] = seqlength

	    # TrEMBL
	    elif providerKey == 316385 and seqlength > plengths[2]:
		allpolypeptide[2][m] = s
	        plengths[2] = seqlength

	    # RefSeq
	    elif providerKey == 316372 and string.find(a, 'XP_') > -1 and seqlength > plengths[3]:
		allpolypeptide[3][m] = s
	        plengths[3] = seqlength

	    prevMarker = m

	# last record

        # determine the one and only representative w/in each group

        for i in range(2):
            if allgenomic[i].has_key(prevMarker):
                genomic[prevMarker] = allgenomic[i][prevMarker]
                break

        for i in range(7):
            if alltranscript[i].has_key(prevMarker):
		transcript[prevMarker] = []
		transcript[prevMarker].append(alltranscript[i][prevMarker])
                break

        for i in range(4):
            if allpolypeptide[i].has_key(prevMarker):
                polypeptide[prevMarker] = allpolypeptide[i][prevMarker]
                break

	# for markers that can have mulitple representative transcripts
        if len(multitranscript) > 0:
	    transcript[prevMarker] = multitranscript

	print 'qualifier end...%s' % (mgi_utils.date())

	print 'final results begin...%s' % (mgi_utils.date())
	results = db.sql('select distinct sequenceKey, markerKey, organismKey, providerKey, typeKey, logicalDBKey, refsKey, userKey, mdate, accID ' + \
		'from #finalannot', 'auto')
	print 'final results end...%s' % (mgi_utils.date())

	for r in results:
	    writeRecord(r)

	outBCP.close()

#
# Main Routine
#

db.useOneConnection(1)
db.set_sqlLogFunction(db.sqlLogAll)
createBCP()
db.useOneConnection(0)
