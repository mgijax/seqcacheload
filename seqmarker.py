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
DL = '|'
table = os.environ['TABLE']
loaddate = loadlib.loaddate

outBCP = None
qualifiers1 = {}
qualifiers2 = {}
seqProviders = {}
seqTypes = {}
genomic = {}
transcript = {}
polypeptide = {}

def writeRecord(r):

    outBCP.write(mgi_utils.prvalue(r['sequenceKey']) + DL + \
	       	mgi_utils.prvalue(r['markerKey']) + DL + \
	       	mgi_utils.prvalue(r['refsKey']) + DL)

    printedQualifier = 0
    if genomic.has_key(r['markerKey']):
	if genomic[r['markerKey']] == r['sequenceKey']:
	    outBCP.write(mgi_utils.prvalue(qualifiers1['genomic']) + DL)
	    printedQualifier = 1

    if transcript.has_key(r['markerKey']):
       if transcript[r['markerKey']] == r['sequenceKey']:
	    outBCP.write(mgi_utils.prvalue(qualifiers1['transcript']) + DL)
	    printedQualifier = 1

    if polypeptide.has_key(r['markerKey']):
        if polypeptide[r['markerKey']] == r['sequenceKey']:
	    outBCP.write(mgi_utils.prvalue(qualifiers1['polypeptide']) + DL)
	    printedQualifier = 1

    if not printedQualifier:
        outBCP.write(mgi_utils.prvalue(qualifiers1['Not Specified']) + DL)

    outBCP.write(r['mdate'] + DL + \
        mgi_utils.prvalue(r['userKey']) + DL + mgi_utils.prvalue(r['userKey']) + DL + \
        loaddate + DL + loaddate + NL)

def createBCP():
	global qualifiers1, qualifiers2, seqProviders, seqTypes
	global genomic, transcript, polypeptide
	global outBCP

	db.useOneConnection(1)

	outBCP = open('%s.bcp' % (table), 'w')

	results = db.sql('select _Term_key, term from VOC_Term_RepQualifier_View', 'auto')
	for r in results:
	   qualifiers1[r['term']] = r['_Term_key']
	   qualifiers2[r['_Term_key']] = r['term']

	# sequence providers
	results = db.sql('select _Term_key, term from VOC_Term_SequenceProvider_View', 'auto')
	for r in results:
	    seqProviders[r['_Term_key']] = r['term']

	# sequence types
	results = db.sql('select _Term_key, term from VOC_Term_SequenceType_View', 'auto')
	for r in results:
	    seqTypes[r['_Term_key']] = r['term']

	# only mouse markers

	db.sql('select _Marker_key into #mouse from MRK_Marker ' + \
		'where _Organism_key = 1 and _Marker_Status_key in (1,3)', None)
#		'where _Organism_key = 1 and _Marker_Status_key in (1,3) and _Marker_key = 31999', None)
	db.sql('create nonclustered index idx_key on #mouse (_Marker_key)', None)

	# select all non-MGI accession ids for mouse markers 

	db.sql('select m._Marker_key, a._LogicalDB_key, a.accID, r._Refs_key, a._ModifiedBy_key, ' + \
		'mdate = convert(char(10), a.modification_date, 101) ' + \
		'into #mouseAccs ' + \
		'from #mouse m, ACC_Accession a, ACC_AccessionReference r ' + \
		'where m._Marker_key = a._Object_key ' + \
		'and a._MGIType_key = 2 ' + \
		'and a._LogicalDB_key != 1 ' + \
		'and a._Accession_key = r._Accession_key', None)

	db.sql('create nonclustered index idx1 on #mouseAccs (_LogicalDB_key, accID)', None)

	# select all mouse annotations

	db.sql('select sequenceKey = s._Object_key, markerKey = m._Marker_key, ' + \
		'refsKey = m._Refs_key, userKey = m._ModifiedBy_key, m.mdate, m.accID ' + \
		'into #allannot ' + \
		'from #mouseAccs m, ACC_Accession s ' + \
		'where m.accID = s.accID ' + \
		'and m._LogicalDB_key = s._LogicalDB_key ' + \
		'and s._MGIType_key = 19 ', None)

	db.sql('create nonclustered index idx1 on #allannot (sequenceKey)', None)

	# select annotations to all sequences

	db.sql('select a.sequenceKey, a.markerKey, a.refsKey, a.userKey, a.mdate, a.accID ' + \
		'into #allseqannot ' + \
		'from #allannot a, SEQ_Sequence ss ' + \
		'where a.sequenceKey = ss._Sequence_key ', None)

	db.sql('create nonclustered index idx1 on #allseqannot (sequenceKey, markerKey, refsKey)', None)

	# select records, grouping by sequence, marker and reference

	db.sql('select sequenceKey, markerKey, refsKey, userKey, mdate = max(mdate), accID ' + 
		'into #finalannot ' + \
		'from #allseqannot group by sequenceKey, markerKey, refsKey', None)
	db.sql('create nonclustered index idx1 on #finalannot (sequenceKey, markerKey, refsKey, userKey, mdate)', None)
	db.sql('create nonclustered index idx2 on #finalannot (markerKey)', None)

	db.sql('select distinct sequenceKey, markerKey, accID into #deriveQuality ' + \
		'from #finalannot order by markerKey', None)
	db.sql('create nonclustered index idx1 on #deriveQuality (sequenceKey)', None)
	db.sql('create nonclustered index idx2 on #deriveQuality (markerKey)', None)

	#
	# manually curated sequences...
	#

	results = db.sql('select _Sequence_key, _Marker_key, _Qualifier_key from MRK_CuratedRepSequence', 'auto')

	allgenomic = [{}, {}, {}]
	alltranscript = [{}, {}, {}, {}, {}, {}, {}, {}]
	allpolypeptide = [{}, {}, {}, {}, {}]
	prevMarker = ''

	# process manually curated representative values
	# bucket 0 = highest priority sequence (manually curated)
	# bucket 1 = next highest priority sequence
	# etc.

	for r in results:
	    m = r['_Marker_key']
	    s = r['_Sequence_key']
	    q = r['_Qualifier_key']

	    if qualifiers2[q] == 'genomic':
		allgenomic[0][m] = s
	    elif qualifiers2[q] == 'transcript':
		alltranscript[0][m] = s
	    elif qualifiers2[q] == 'polypeptide':
		allpolypeptide[0][m] = s

	# process derived representative values

	results = db.sql('select q.sequenceKey, q.markerKey, q.accID, ' + \
		's._SequenceProvider_key, s._SequenceType_key, s.length ' + \
		'from #deriveQuality q, SEQ_Sequence s ' + \
		'where q.sequenceKey = s._Sequence_key ' + \
		'order by q.markerKey, s._SequenceProvider_key', 'auto')

	for r in results:

	    m = r['markerKey']
	    s = r['sequenceKey']
	    a = r['accID']

	    if r['length'] is None:
		seqlength = 0
	    else:
		seqlength = int(r['length'])

	    provider = seqProviders[r['_SequenceProvider_key']]
	    sType = seqTypes[r['_SequenceType_key']]

	    if prevMarker != m:
	        glengths = [-1,-1,-1]
		tlengths = [-1,-1,-1,-1,-1,-1,-1,-1]
		plengths = [-1,-1,-1,-1,-1]

		if prevMarker != '':

		    # determine the one and only representative w/in each group
		    for i in range(3):
			if allgenomic[i].has_key(prevMarker):
			    genomic[prevMarker] = allgenomic[i][prevMarker]
			    break

		    for i in range(8):
			if alltranscript[i].has_key(prevMarker):
			    transcript[prevMarker] = alltranscript[i][prevMarker]
			    break

		    for i in range(5):
			if allpolypeptide[i].has_key(prevMarker):
			    polypeptide[prevMarker] = allpolypeptide[i][prevMarker]
			    break

	    #
	    # representative genomic
	    #
	    # longest NCBI/Ensembl coordinate OR longest GenBank DNA
	    # tie goes to NCBI
	    #

	    if (provider == 'NCBI Gene Model' or provider == 'Ensembl Gene Model') and seqlength > glengths[1]:
		allgenomic[1][m] = s
	        glengths[1] = seqlength
	    elif provider == 'NCBI Gene Model' and seqlength == glengths[1]:
		allgenomic[1][m] = s
	        glengths[1] = seqlength
	    elif string.find(provider, 'GenBank') > -1 and sType == 'DNA' and seqlength > glengths[2]:
		allgenomic[2][m] = s
	        glengths[2] = seqlength

	    #
	    # representative transcript
	    #
	    # longest NM_ or NR_ RefSeq
	    # longest non-EST GenBank
	    # longest XM_ or XR_ RefSeq
	    # longest TIGR, DoTS, NIA Mouse Gene Index,
	    # longest EST GenBank
	    #

	    elif provider == 'RefSeq' and (string.find(a, 'NM_') > -1 or string.find(a, 'NR_') > -1) and seqlength > tlengths[1]:
		alltranscript[1][m] = s
	        tlengths[1] = seqlength
	    elif string.find(provider, 'GenBank') > -1 and provider != 'GenBank/EMBL/DDBJ:EST' and sType == 'RNA' and seqlength > tlengths[2]:
		alltranscript[2][m] = s
	        tlengths[2] = seqlength
	    elif provider == 'RefSeq' and (string.find(a, 'XM_') > -1 or string.find(a, 'XR_') > -1) and seqlength > tlengths[3]:
		alltranscript[3][m] = s
	        tlengths[3] = seqlength
	    elif provider == 'TIGR Mouse Gene Index' and seqlength > tlengths[4]:
		alltranscript[4][m] = s
	        tlengths[4] = seqlength
	    elif provider == 'DoTS' and seqlength > tlengths[5]:
		alltranscript[5][m] = s
	        tlengths[5] = seqlength
	    elif provider == 'NIA Mouse Gene Index' and seqlength > tlengths[6]:
		alltranscript[6][m] = s
	        tlengths[6] = seqlength
	    elif provider == 'GenBank/EMBL/DDBJ:EST' and sType == 'RNA' and seqlength > tlengths[7]:
		alltranscript[7][m] = s
	        tlengths[7] = seqlength

	    #
	    # representative polypeptide
	    #
	    # longest SWISS-PROT
	    # longest NP_ RefSeq
	    # longest TrEMBL
	    # longest XP_ RefSeq
	    #

	    if provider == 'SWISS-PROT' and seqlength > plengths[1]:
		allpolypeptide[1][m] = s
	        plengths[1] = seqlength
	    elif provider == 'RefSeq' and string.find(a, 'NP_') > -1 and seqlength > plengths[2]:
		allpolypeptide[2][m] = s
	        plengths[2] = seqlength
	    elif provider == 'TrEMBL' and seqlength > plengths[3]:
		allpolypeptide[3][m] = s
	        plengths[3] = seqlength
	    elif provider == 'RefSeq' and string.find(a, 'XP_') > -1 and seqlength > plengths[4]:
		allpolypeptide[4][m] = s
	        plengths[4] = seqlength

	    prevMarker = m

	# last record
        # determine the one and only representative w/in each group
        for i in range(3):
            if allgenomic[i].has_key(prevMarker):
                genomic[prevMarker] = allgenomic[i][prevMarker]
                break

        for i in range(8):
            if alltranscript[i].has_key(prevMarker):
                transcript[prevMarker] = alltranscript[i][prevMarker]
                break

        for i in range(5):
            if allpolypeptide[i].has_key(prevMarker):
                polypeptide[prevMarker] = allpolypeptide[i][prevMarker]
                break

	print 'qualifier end...%s' % (mgi_utils.date())

	print 'final results begin...%s' % (mgi_utils.date())
	results = db.sql('select distinct sequenceKey, markerKey, refsKey, userKey, mdate from #finalannot', 'auto')
	print 'final results end...%s' % (mgi_utils.date())

	for r in results:
	    writeRecord(r)

	outBCP.close()
	db.useOneConnection(0)

#
# Main Routine
#

db.set_sqlLogFunction(db.sqlLogAll)
createBCP()
