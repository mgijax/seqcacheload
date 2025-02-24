
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
#	seqcoord.py
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
import mgi_utils
import loadlib
import db

NL = '\n'
DL = os.environ['COLDELIM']
table = os.environ['TABLE']
datadir = os.environ['CACHEDATADIR']
userKey = 0
loaddate = loadlib.loaddate

def createBCP():

        print('Creating %s.bcp...%s' % (table, mgi_utils.date()))

        outBCP = open('%s/%s.bcp' % (datadir, table), 'w')

        cmd = '''
            select distinct mc._Map_key, mc.version, t2.abbreviation as mapUnits,
            mcf._Object_key, mcf.startCoordinate, mcf.endCoordinate, mcf.strand,
            c.chromosome, t3.term as provider 
            from MAP_Coordinate mc, MAP_Coord_Feature mcf,
            MRK_Chromosome c, SEQ_Sequence s, VOC_Term t1, VOC_Term t2, VOC_Term t3 
            where mc._MapType_key = t1._Term_key 
            and t1.term = 'Assembly' 
            and mc._Units_key = t2._Term_key 
            and mc._MGIType_key = 27 
            and mc._Object_key = c._Chromosome_key 
            and mc._Map_key = mcf._Map_key 
            and mcf._MGIType_key = 19 
            and mcf._Object_key = s._Sequence_key 
            and s._SequenceProvider_key = t3._Term_key
            '''

        results = db.sql(cmd, 'auto')

        for r in results:

                outBCP.write(str(r['_Map_key']) + DL + \
                        str(r['_Object_key']) + DL + \
                        r['chromosome'] + DL + \
                        str(r['startCoordinate']) + DL + \
                        str(r['endCoordinate']) + DL + \
                        str(r['strand']) + DL + \
                        str(r['mapUnits']) + DL + \
                        str(r['provider']) + DL + \
                        str(r['version']) + DL + \
                        str(userKey) + DL + str(userKey) + DL + \
                        loaddate + DL + loaddate + NL)

        outBCP.close()

#
# Main Routine
#

userKey = loadlib.verifyUser(os.environ['MGD_DBUSER'], 1, None)

db.useOneConnection(1)
print('%s' % mgi_utils.date())
createBCP()
print('%s' % mgi_utils.date())
db.useOneConnection(0)
