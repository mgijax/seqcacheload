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
# 07/07/2004	lec
#	- Assembly (TR 5395)
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

        cmds.append('select mc._Map_key, mc.seqRetrievalParam, ' + \
	      'mcf._Object_key, mcf.startCoordinate, mcf.endCoordinate, mcf.strand, c.chromosome ' + \
	      'into #sequences ' + \
              'from MAP_Coordinate mc, VOC_Term t, MAP_Coord_Feature mcf, MRK_Chromosome c ' + \
	      'where mc._MapType_key = t._Term_key ' + \
	      'and t.term = "Assembly" ' + \
	      'and mc._MGIType_key = type for chromosome ' + \
	      'and mc._Object_key = c._Chromosome_key ' + \
	      'and mc._Map_key = mcf._Map_key ' + \
	      'and mcf._MGIType_key = 19')

	cmds.append('create index idx_object on #sequences(_Object_key)')

	cmds.append('select m.*, s.version, s.description, provider = t.term ' + \
              'from #sequences s, ACC_Accession a, SEQ_Sequence s, VOC_Term t ' + \
	      'where m._Object_key = a._Object_key ' + \
	      'and a._MGIType_key = 19 ' + \
	      'and a._Object_key = s._Sequence_key ' + \
	      'and s._SequenceProvider_key = t._Term_key')

	results = db.sql(cmds, 'auto')

	for r in results[-1]:

		outBCP.write(str(r['_Map_key']) + DL + \
			str(r['_Object_key']) + DL + \
			r['chromosome'] + DL + \
			str(r['startCoordinate']) + DL + \
			str(r['endCoordinate']) + DL + \
			str(r['strand']) + DL + \
			str(r['provider']) + DL + \
			str(r['version']) + DL + \
			str(r['description']) + DL + \
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

