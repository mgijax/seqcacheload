#!/usr/local/bin/python

'''
#
# Purpose:
#
# Create bcp file for SEQ_Probe_Cache
#
# Uses environment variables to determine Server and Database
# (DSQUERY and MGD).
#
# Usage:
#	seqprobe.py [probekey]
#
# If probekey is provided, then only create the bcp file for that probe.
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

def createBCP(probeKey):

	print 'Creating %s.bcp...%s' % (table, mgi_utils.date())

	outBCP = open('%s.bcp' % (table), 'w')

	cmd = 'select distinct sequenceKey = s._Object_key, probeKey = p._Object_key, ' + \
		'refsKey = ar._Refs_key ' + \
		'from SEQ_Sequence_Acc_View s, PRB_Acc_View p, ACC_AccessionReference ar ' + \
		'where s.accID = p.accID ' + \
		'and s._LogicalDB_key = p._LogicalDB_key ' + \
		'and p._Accession_key = ar._Accession_key'

	if probeKey is not None:
		cmd = cmd + 'and p._Object_key = %s\n' % probeKey

	results = db.sql(cmd, 'auto')

	for r in results:

		outBCP.write(mgi_utils.prvalue(r['sequenceKey']) + DL + \
		       	mgi_utils.prvalue(r['probeKey']) + DL + \
		       	mgi_utils.prvalue(r['refsKey']) + DL + \
			str(userKey) + DL + str(userKey) + DL + \
			loaddate + DL + loaddate + NL)

	outBCP.close()

#
# Main Routine
#

userKey = loadlib.verifyUser(os.environ['DBUSER'], 1, None)

if len(sys.argv) == 2:
	probeKey = sys.argv[1]
else:
	probeKey = None

print '%s' % mgi_utils.date()
createBCP(probeKey)
print '%s' % mgi_utils.date()

