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
import db
import mgi_utils
import loadlib

NL = '\n'
DL = '|'
table = os.environ['TABLE']
userKey = 0
loaddate = loadlib.loaddate

def createBCP():

	print 'Creating %s.bcp...%s' % (table, mgi_utils.date())

	outBCP = open('%s.bcp' % (table), 'w')

	qresults = db.sql('select _Term_key from VOC_Term_RepQualifier_View where term = "Not Specified"', 'auto')
	for r in qresults:
	   qualityKey = r['_Term_key']

	cmds = []

	# only mouse markers

	cmds.append('select _Marker_key into #mouse from MRK_Marker where _Organism_key = 1 and _Marker_Status_key in (1,3)')
	cmds.append('create nonclustered index idx_key on #mouse (_Marker_key)')

	# select annotations to dummy sequences

	cmds.append('select m.accID, sequenceKey = s._Object_key, ' + \
		'markerKey = m._Object_key, refsKey = ar._Refs_key, ' + \
		'mdate = convert(char(10), m.modification_date, 101) ' + \
		'into #dummy ' + \
		'from #mouse mm, ACC_Accession m, ACC_Accession s, ACC_AccessionReference ar, ' + \
		'SEQ_Sequence ss, VOC_Term t ' + \
		'where mm._Marker_key = m._Object_key ' + \
		'and m._MGIType_key = 2 ' + \
		'and m.accID = s.accID ' + \
		'and m._LogicalDB_key = s._LogicalDB_key ' + \
		'and s._MGIType_key = 19 ' + \
		'and m._Accession_key = ar._Accession_key ' + \
		'and s._Object_key = ss._Sequence_key ' + \
		'and ss._SequenceStatus_key = t._Term_key ' + \
		'and t.term = "Not Loaded"')

	cmds.append('create nonclustered index idx_acc on #dummy (accID)')
	cmds.append('create nonclustered index idx_seq on #dummy (sequenceKey)')
	cmds.append('create nonclustered index idx_mrk on #dummy (markerKey)')
	cmds.append('create nonclustered index idx_ref on #dummy (refsKey)')
	cmds.append('create nonclustered index idx_mdt on #dummy (mdate)')

	cmds.append('select m.accID, sequenceKey = s._Object_key, ' + \
		'markerKey = m._Object_key, refsKey = ar._Refs_key, ' + \
		'mdate = convert(char(10), m.modification_date, 101) ' + \
		'into #sequences ' + \
		'from #mouse mm, ACC_Accession m, ACC_Accession s, ACC_AccessionReference ar, ' + \
		'SEQ_Sequence ss, VOC_Term t ' + \
		'where mm._Marker_key = m._Object_key ' + \
		'and m._MGIType_key = 2 ' + \
		'and m.accID = s.accID ' + \
		'and m._LogicalDB_key = s._LogicalDB_key ' + \
		'and s._MGIType_key = 19 ' + \
		'and m._Accession_key = ar._Accession_key ' + \
		'and s._Object_key = ss._Sequence_key ' + \
		'and ss._SequenceStatus_key = t._Term_key ' + \
		'and t.term != "Not Loaded"')

	cmds.append('create nonclustered index idx_acc on #sequences (accID)')
	cmds.append('create nonclustered index idx_seq on #sequences (sequenceKey)')
	cmds.append('create nonclustered index idx_mrk on #sequences (markerKey)')
	cmds.append('create nonclustered index idx_ref on #sequences (refsKey)')
	cmds.append('create nonclustered index idx_mdt on #sequences (mdate)')

	#  don't cache marker/dummy annotations if the real thing exists

	cmds.append('select distinct sequenceKey, markerKey, refsKey, mdate = max(mdate) from #sequences ' + \
		'group by sequenceKey, markerKey, refsKey ' + \
		'union ' + \
		'select distinct sequenceKey, markerKey, refsKey, mdate = max(mdate) from #dummy d ' + \
		'where not exists (select 1 from #sequences s where d.accID = s.accID) ' + \
		'group by sequenceKey, markerKey, refsKey')

	results = db.sql(cmds, 'auto')

	for r in results[-1]:

		outBCP.write(mgi_utils.prvalue(r['sequenceKey']) + DL + \
		       	mgi_utils.prvalue(r['markerKey']) + DL + \
		       	mgi_utils.prvalue(r['refsKey']) + DL + \
		       	mgi_utils.prvalue(r['qualityKey']) + DL + \
			r['mdate'] + DL + \
			str(userKey) + DL + str(userKey) + DL + \
			loaddate + DL + loaddate + NL)

	outBCP.close()

#
# Main Routine
#

userKey = loadlib.verifyUser(os.environ['DBUSER'], 1, None)
print '%s' % mgi_utils.date()
createBCP()
print '%s' % mgi_utils.date()

