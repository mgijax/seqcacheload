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

#	qresults = db.sql('select _Term_key from VOC_Term_RepQualifier_View where term = "Not Specified"', 'auto')
#	for r in qresults:
#	   qualityKey = r['_Term_key']

	cmds = []

	# only mouse markers

	cmds.append('select _Marker_key into #mouse from MRK_Marker where _Organism_key = 1 and _Marker_Status_key in (1,3)')
	cmds.append('create nonclustered index idx_key on #mouse (_Marker_key)')

	# select all non-MGI accession ids for mouse markers 

	cmds.append('select m._Marker_key, a._LogicalDB_key, a.accID, r._Refs_key, ' + \
		'mdate = convert(char(10), a.modification_date, 101) ' + \
		'into #mouseAccs ' + \
		'from #mouse m, ACC_Accession a, ACC_AccessionReference r ' + \
		'where m._Marker_key = a._Object_key ' + \
		'and a._MGIType_key = 2 ' + \
		'and a._LogicalDB_key != 1 ' + \
		'and a._Accession_key = r._Accession_key')

	cmds.append('create nonclustered index idx1 on #mouseAccs (_Marker_key)')
	cmds.append('create nonclustered index idx2 on #mouseAccs (_LogicalDB_key)')
	cmds.append('create nonclustered index idx3 on #mouseAccs (accID)')
	cmds.append('create nonclustered index idx4 on #mouseAccs (mdate)')

	# select all mouse annotations

	cmds.append('select sequenceKey = s._Object_key, markerKey = m._Marker_key, refsKey = m._Refs_key, m.mdate ' + \
		'into #allannot ' + \
		'from #mouseAccs m, ACC_Accession s ' + \
		'where m.accID = s.accID ' + \
		'and m._LogicalDB_key = s._LogicalDB_key ' + \
		'and s._MGIType_key = 19 ')

	cmds.append('create nonclustered index idx1 on #allannot (sequenceKey)')

	# select annotations to dummy sequences

	cmds.append('select a.sequenceKey, a.markerKey, a.refsKey, a.mdate ' + \
		'into #dummyannot ' + \
		'from #allannot a, SEQ_Sequence ss ' + \
		'where a.sequenceKey = ss._Sequence_key ' + \
		'and ss._SequenceStatus_key = 316345')

	cmds.append('create nonclustered index idx1 on #dummyannot (sequenceKey)')
	cmds.append('create nonclustered index idx2 on #dummyannot (markerKey)')
	cmds.append('create nonclustered index idx3 on #dummyannot (refsKey)')
	cmds.append('create nonclustered index idx4 on #dummyannot (mdate)')

	# select annotations to nondummy sequences

	cmds.append('select a.sequenceKey, a.markerKey, a.refsKey, a.mdate ' + \
		'into #nondummyannot ' + \
		'from #allannot a, SEQ_Sequence ss ' + \
		'where a.sequenceKey = ss._Sequence_key ' + \
		'and ss._SequenceStatus_key != 316345')

	cmds.append('create nonclustered index idx1 on #nondummyannot (sequenceKey)')
	cmds.append('create nonclustered index idx2 on #nondummyannot (markerKey)')
	cmds.append('create nonclustered index idx3 on #nondummyannot (refsKey)')
	cmds.append('create nonclustered index idx4 on #nondummyannot (mdate)')

	#  don't cache marker/dummy annotations if the real thing exists

	cmds.append('delete #dummyannot from #dummyannot d, #nondummyannot n ' + \
		'where d.sequenceKey = n.sequenceKey and d.markerKey = n.markerKey')

	# select records, grouping by sequence, marker and reference

	cmds.append('select sequenceKey, markerKey, refsKey, mdate = max(mdate) ' + 
		'from #dummyannot group by sequenceKey, markerKey, refsKey')

	cmds.append('select sequenceKey, markerKey, refsKey, mdate = max(mdate) ' + 
		'from #nondummyannot group by sequenceKey, markerKey, refsKey')

	results = db.sql(cmds, 'auto')

	for r in results[-1]:

#		       	mgi_utils.prvalue(r['qualityKey']) + DL + \
		outBCP.write(mgi_utils.prvalue(r['sequenceKey']) + DL + \
		       	mgi_utils.prvalue(r['markerKey']) + DL + \
		       	mgi_utils.prvalue(r['refsKey']) + DL + \
			r['mdate'] + DL + \
			str(userKey) + DL + str(userKey) + DL + \
			loaddate + DL + loaddate + NL)

	for r in results[-2]:

#		       	mgi_utils.prvalue(r['qualityKey']) + DL + \
		outBCP.write(mgi_utils.prvalue(r['sequenceKey']) + DL + \
		       	mgi_utils.prvalue(r['markerKey']) + DL + \
		       	mgi_utils.prvalue(r['refsKey']) + DL + \
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

