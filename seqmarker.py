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
#	seqmarker.py [markerkey]
#
# If markerkey is provided, then only create the bcp file for that marker.
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

	cmds = []
	cmds.append('select sequenceKey = s._Object_key, ' + \
		'markerKey = m._Object_key, refsKey = ar._Refs_key, ' + \
		'mdate = convert(char(10), m.modification_date, 101) ' + \
		'into #sequences ' + \
		'from ACC_Accession s, ACC_Accession m, ACC_AccessionReference ar ' + \
		'where s._MGIType_key = 19 ' + \
		'and s.accID = m.accID ' + \
		'and m._MGIType_key = 2 ' + \
		'and s._LogicalDB_key = m._LogicalDB_key ' + \
		'and m._Accession_key = ar._Accession_key')

	cmds.append('create nonclustered index idx_seq on #sequences (sequenceKey)')
	cmds.append('create nonclustered index idx_mrk on #sequences (markerKey)')
	cmds.append('create nonclustered index idx_ref on #sequences (refsKey)')
	cmds.append('create nonclustered index idx_mdt on #sequences (mdate)')

	cmds.append('select distinct sequenceKey, markerKey, refsKey, mdate from #sequences')

	results = db.sql(cmds, 'auto')

	for r in results[-1]:

		outBCP.write(mgi_utils.prvalue(r['sequenceKey']) + DL + \
		       	mgi_utils.prvalue(r['markerKey']) + DL + \
		       	mgi_utils.prvalue(r['refsKey']) + DL + \
			mdate + DL + \
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

