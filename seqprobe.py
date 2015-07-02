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
# 11/23/2004	lec
#	- added createExcluded() for TR 6118 (GXD Gray data load)
#
# 10/23/2003	lec
#	- new (TR 3404, JSAM)
#
'''

import sys
import os
import mgi_utils
import loadlib
import db

db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

NL = '\n'
DL = os.environ['COLDELIM']
table = os.environ['TABLE']
datadir = os.environ['CACHEDATADIR']
loaddate = loadlib.loaddate

def createExcluded():

    excludeNote = 'The source of the material used to create this cDNA probe was different than that used to create the GenBank sequence record.'

    print 'excluded begin...%s' % (mgi_utils.date())
    cmds = []
    cmds.append('select _Probe_key INTO TEMPORARY TABLE excluded from PRB_Notes ' + \
	'where note like \'The source of the material used to create this cDNA probe was different%\'')
    cmds.append('create index idx1 on excluded(_Probe_key)')
    db.sql(cmds, None)
    print 'excluded end...%s' % (mgi_utils.date())

def createBCP():

	outBCP = open('%s/%s.bcp' % (datadir, table), 'w')

	print 'sequences1 begin...%s' % (mgi_utils.date())
	cmds = []
	cmds.append('select s._Object_key as sequenceKey, p._Object_key as probeKey, p._Accession_key ' + \
		'INTO TEMPORARY TABLE sequences1 ' + \
		'from ACC_Accession s, ACC_Accession p ' + \
		'where s._MGIType_key = 19 ' + \
		'and s.accID = p.accID ' + \
		'and p._MGIType_key = 3 ' + \
		'and s._LogicalDB_key = p._LogicalDB_key ')
	cmds.append('create index idx2 on sequences1 (sequenceKey)')
	cmds.append('create index idx3 on sequences1 (probeKey)')
	cmds.append('create index idx4 on sequences1 (_Accession_key)')
	db.sql(cmds, None)
	print 'sequences1 end...%s' % (mgi_utils.date())

	print 'deletion begin...%s' % (mgi_utils.date())
	db.sql('delete from sequences1 using excluded e where sequences1.probeKey = e._Probe_key', None)
	print 'deletion end...%s' % (mgi_utils.date())
	db.commit()

	print 'sequences2 begin...%s' % (mgi_utils.date())
	cmds = []
	cmds.append('select s.sequenceKey, s.probeKey, ar._Refs_key as refskey, ' + \
		'ar._ModifiedBy_key as userKey, ar.modification_date as mdate ' + \
		'INTO TEMPORARY TABLE sequences2 ' + \
		'from sequences1 s, ACC_AccessionReference ar ' + \
		'where s._Accession_key = ar._Accession_key')
	cmds.append('create index idx5 on sequences2 (sequenceKey, probeKey, refsKey, userKey, mdate)')
	cmds.append('create index idx6 on sequences2 (userKey)')
	cmds.append('create index idx7 on sequences2 (mdate)')
	db.sql(cmds, None)
	print 'sequences2 end...%s' % (mgi_utils.date())

	print 'final begin...%s' % (mgi_utils.date())
	results = db.sql('select distinct sequenceKey, probeKey, refsKey, max(userKey) as userKey, max(mdate) as mdate from sequences2 ' + \
		'group by sequenceKey, probeKey, refsKey', 'auto')
	print 'final end...%s' % (mgi_utils.date())

	for r in results:
		outBCP.write(mgi_utils.prvalue(r['sequenceKey']) + DL + \
		       	mgi_utils.prvalue(r['probeKey']) + DL + \
		       	mgi_utils.prvalue(r['refsKey']) + DL + \
			r['mdate'] + DL + \
        		mgi_utils.prvalue(r['userKey']) + DL + mgi_utils.prvalue(r['userKey']) + DL + \
			loaddate + DL + loaddate + NL)
	outBCP.close()

#
# Main Routine
#

db.useOneConnection(1)
print '%s' % mgi_utils.date()
createExcluded()
createBCP()
db.useOneConnection(0)
print '%s' % mgi_utils.date()

