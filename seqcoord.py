#!/usr/local/bin/python

'''
#
# Purpose:
#
# Create bcp file for SEQ_Coord_Cache
#
# Uses environment variables to determine Server and Database
# (DSQUERY and MGD).
#
# Usage:
#	seqcoord.py [coordkey]
#
# If coordkey is provided, then only create the bcp file for that coord.
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

        cmd = 'select mcf.startBP, mcf.endBP, mcf.strand, c.chromosome, sequenceKey = a._Object_key ' + \
              'from MAP_Coord_Feature mcf, MAP_Coordinate mc, VOC_Term t, MRK_Chromosome c, MRK_Acc_View a ' + \
	      'where mcf._Map_key = mc._Map_key ' + \
	      'and mc._MapType_key = t._Term_key ' + \
	      'and t.term = "Assembly" ' + \
	      'and mc._Chromosome_key = c._Chromosome_key ' + \
	      'and mcf._MGIType_key = 2 ' + \
	      'and mcf._Object_key = a._Object_key ' + \
	      'and a._LogicalDB_key = ?'

	results = db.sql(cmd, 'auto')

	for r in results:

		outBCP.write(str(r['sequenceKey']) + DL + \
			r['chromosome'] + DL + \
			str(r['startBP']) + DL + \
			str(r['endBP']) + DL + \
			str(r['strand']) + DL + \
			str(userKey) + DL + str(userKey) + DL + \
			loaddate + DL + loaddate + NL)

	outBCP.close()

#
# Main Routine
#

userKey = loadlib.verifyUser(os.environ['DBUSER'], 1, None)

if len(sys.argv) == 2:
	coordKey = sys.argv[1]
else:
	coordKey = None

print '%s' % mgi_utils.date()

# Log all SQL commands
#db.set_sqlLogFunction(db.sqlLogAll)

createBCP()

print '%s' % mgi_utils.date()

