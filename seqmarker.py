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

	db.useOneConnection(1)

	print 'Creating %s.bcp...%s' % (table, mgi_utils.date())

	outBCP = open('%s.bcp' % (table), 'w')

	qresults = db.sql('select _Term_key from VOC_Term_RepQualifier_View where term = "Not Specified"', 'auto')
	for r in qresults:
	   qualityKey = r['_Term_key']

	# only mouse markers

	print 'markers begin...%s' % (mgi_utils.date())
	cmds = []
	cmds.append('select _Marker_key into #mouse from MRK_Marker where _Organism_key = 1 and _Marker_Status_key in (1,3)')
	cmds.append('create nonclustered index idx_key on #mouse (_Marker_key)')
	db.sql(cmds, None)
	print 'markers end...%s' % (mgi_utils.date())

	# select all non-MGI accession ids for mouse markers 

	print 'accessions begin...%s' % (mgi_utils.date())
	cmds = []
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
	db.sql(cmds, None)
	print 'accessions end...%s' % (mgi_utils.date())

	# select all mouse annotations

	print 'all annotations begin...%s' % (mgi_utils.date())
	cmds = []
	cmds.append('select sequenceKey = s._Object_key, markerKey = m._Marker_key, refsKey = m._Refs_key, m.mdate ' + \
		'into #allannot ' + \
		'from #mouseAccs m, ACC_Accession s ' + \
		'where m.accID = s.accID ' + \
		'and m._LogicalDB_key = s._LogicalDB_key ' + \
		'and s._MGIType_key = 19 ')

	cmds.append('create nonclustered index idx1 on #allannot (sequenceKey)')
	db.sql(cmds, None)
	print 'all annotations end...%s' % (mgi_utils.date())

	# select annotations to nondummy sequences

	print 'non-dummy begin...%s' % (mgi_utils.date())
	cmds = []
	cmds.append('select a.sequenceKey, a.markerKey, a.refsKey, a.mdate ' + \
		'into #nondummyannot ' + \
		'from #allannot a, SEQ_Sequence ss ' + \
		'where a.sequenceKey = ss._Sequence_key ' + \
		'and ss._SequenceStatus_key != 316345')

	cmds.append('create nonclustered index idx1 on #nondummyannot (sequenceKey)')
	cmds.append('create nonclustered index idx2 on #nondummyannot (markerKey)')
	cmds.append('create nonclustered index idx3 on #nondummyannot (refsKey)')
	cmds.append('create nonclustered index idx4 on #nondummyannot (mdate)')
	db.sql(cmds, None)
	print 'non-dummy end...%s' % (mgi_utils.date())

	# select annotations to dummy sequences

	print 'dummy begin...%s' % (mgi_utils.date())
	cmds = []
	cmds.append('select a.sequenceKey, a.markerKey, a.refsKey, a.mdate ' + \
		'into #dummyannot1 ' + \
		'from #allannot a, SEQ_Sequence ss ' + \
		'where a.sequenceKey = ss._Sequence_key ' + \
		'and ss._SequenceStatus_key = 316345')

	cmds.append('create nonclustered index idx1 on #dummyannot1 (sequenceKey)')
	cmds.append('create nonclustered index idx2 on #dummyannot1 (markerKey)')
	cmds.append('create nonclustered index idx3 on #dummyannot1 (refsKey)')
	cmds.append('create nonclustered index idx4 on #dummyannot1 (mdate)')
	db.sql(cmds, None)
	print 'dummy end...%s' % (mgi_utils.date())

	# select annotations to dummy sequences which are not also in non-dummy sequences

	print 'dummy2 begin...%s' % (mgi_utils.date())
	cmds = []
	cmds.append('select d.sequenceKey, d.markerKey, d.refsKey, d.mdate ' + \
		'into #dummyannot ' + \
		'from #dummyannot1 d ' + \
		'where not exists (select 1 from #nondummyannot n ' + \
		'where d.sequenceKey = n.sequenceKey ' + \
		'and d.markerKey = n.markerKey)')

	cmds.append('create nonclustered index idx1 on #dummyannot (sequenceKey)')
	cmds.append('create nonclustered index idx2 on #dummyannot (markerKey)')
	cmds.append('create nonclustered index idx3 on #dummyannot (refsKey)')
	cmds.append('create nonclustered index idx4 on #dummyannot (mdate)')
	db.sql(cmds, None)
	print 'dummy2 end...%s' % (mgi_utils.date())

	# select records, grouping by sequence, marker and reference

	print 'final results begin...%s' % (mgi_utils.date())
	cmds = []
	cmds.append('select sequenceKey, markerKey, refsKey, mdate = max(mdate) ' + 
		'from #dummyannot group by sequenceKey, markerKey, refsKey')

	cmds.append('select sequenceKey, markerKey, refsKey, mdate = max(mdate) ' + 
		'from #nondummyannot group by sequenceKey, markerKey, refsKey')
	results = db.sql(cmds, 'auto')
	print 'final results end...%s' % (mgi_utils.date())

	for r in results[-1]:

		outBCP.write(mgi_utils.prvalue(r['sequenceKey']) + DL + \
		       	mgi_utils.prvalue(r['markerKey']) + DL + \
		       	mgi_utils.prvalue(r['refsKey']) + DL + \
		       	mgi_utils.prvalue(qualityKey) + DL + \
			r['mdate'] + DL + \
			str(userKey) + DL + str(userKey) + DL + \
			loaddate + DL + loaddate + NL)

	for r in results[-2]:

		outBCP.write(mgi_utils.prvalue(r['sequenceKey']) + DL + \
		       	mgi_utils.prvalue(r['markerKey']) + DL + \
		       	mgi_utils.prvalue(r['refsKey']) + DL + \
		       	mgi_utils.prvalue(qualityKey) + DL + \
			r['mdate'] + DL + \
			str(userKey) + DL + str(userKey) + DL + \
			loaddate + DL + loaddate + NL)

	outBCP.close()
	db.useOneConnection(0)

#
# Main Routine
#

userKey = loadlib.verifyUser(os.environ['DBUSER'], 1, None)
print '%s' % mgi_utils.date()
createBCP()
print '%s' % mgi_utils.date()

