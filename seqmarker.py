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

def createBCP(markerKey):

	print 'Creating %s.bcp...%s' % (table, mgi_utils.date())

	outBCP = open('%s.bcp' % (table), 'w')

	cmd = 'select distinct sequenceKey = s._Object_key, ' + \
		'markerKey = m._Object_key, refsKey = ar._Refs_key,
		'mdate = convert(char(10), m.modification_date, 101) ' + \
		'from ACC_Accession s, ACC_Accession m, ACC_AccessionReference ar ' + \
		'where s._MGIType_key = 19 ' + \
		'and s.accID = m.accID ' + \
		'and m._MGIType_key = 2 ' + \
		'and s._LogicalDB_key = m._LogicalDB_key ' + \
		'and m._Accession_key = ar._Accession_key'

	if markerKey is not None:
		cmd = cmd + ' and m._Object_key = %s\n' % markerKey

	results = db.sql(cmd, 'auto')

	for r in results:

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

if len(sys.argv) == 2:
	markerKey = sys.argv[1]
else:
	markerKey = None

print '%s' % mgi_utils.date()
createBCP(markerKey)
print '%s' % mgi_utils.date()

