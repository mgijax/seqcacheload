#!/usr/local/bin/python

'''
#
# Purpose:
#
# Create bcp file for SEQ_MarkerCache
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
#	- sequence include GenBank (9), SwissProt (13), TrEMBL (41), RefSeq (27)
#	  TIGR (35), DoTs (36)
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

	cmd = 'select distinct sequenceKey = s._Object_key, markerKey = m._Object_key ' + \
		'from SEQ_Sequence_Acc_View s, MRK_Acc_View m ' + \
		'where s.accID = m.accID ' + \
		'and m._LogicalDB_Key in (9, 13, 27, 35, 36, 41)'

	if markerKey is not None:
		cmd = cmd + 'and m._Object_key = %s\n' % markerKey

	results = db.sql(cmd, 'auto')

	for r in results:

		outBCP.write(mgi_utils.prvalue(r['sequenceKey']) + DL + \
		       	mgi_utils.prvalue(r['markerKey']) + DL + \
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

# Log all SQL commands
#db.set_sqlLogFunction(db.sqlLogAll)

createBCP(markerKey)

print '%s' % mgi_utils.date()

