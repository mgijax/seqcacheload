#!/usr/local/bin/python

'''
#
# Purpose:
#
# Create bcp file for SEQ_Description_Cache
#
# Uses environment variables to determine Server and Database
# (DSQUERY and MGD).
#
# Usage:
#	seqdescription.py
#
# History
#
# 03/30/2004	lec
#	- JSAM
#
'''

import sys
import os
import db
import mgi_utils
import loadlib

NL = '\n'
DL = os.environ['COLDELIM']
table = os.environ['TABLE']
datadir = os.environ['CACHEDATADIR']
userKey = 0
loaddate = loadlib.loaddate

def createBCP():

	print 'Creating %s.bcp...%s' % (table, mgi_utils.date())

	outBCP = open('%s/%s.bcp' % (datadir, table), 'w')

	cmds = []

	cmds.append('select distinct _Sequence_key into #sequences from SEQ_Marker_Cache')
	cmds.append('create nonclustered index idx_seq on #sequences (_Sequence_key)')
	cmds.append('select s._Sequence_key, s.description ' + \
		'from #sequences c, SEQ_Sequence s ' + \
		'where c._Sequence_key = s._Sequence_key ' + \
		'and s.description is not null')

	results = db.sql(cmds, 'auto')

	for r in results[-1]:

		outBCP.write(mgi_utils.prvalue(r['_Sequence_key']) + DL + \
		       	mgi_utils.prvalue(r['description']) + DL + \
			str(userKey) + DL + str(userKey) + DL + \
			loaddate + DL + loaddate + NL)

	outBCP.close()

#
# Main Routine
#

db.useOneConnection(1)
db.set_sqlLogFunction(db.sqlLogAll)
userKey = loadlib.verifyUser(os.environ['MGI_DBUSER'], 1, None)
print '%s' % mgi_utils.date()
createBCP()
print '%s' % mgi_utils.date()
db.useOneConnection(0)

