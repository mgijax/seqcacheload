#!/usr/local/bin/python
#
# seqbiotype.py
##########################################################################k
#
#  Purpose: This script updates the SEQ_Marker_Cache
#	    rawbiotype, _BiotypeConflict_key
#	    values by _Marker_key
#
#	for each Marker in SEQ_Marker_Cache
#	  that contains a NCBI/Ensembl/VEGA sequence in SEQ_GeneModel:
#
#	  determine if the MGD marker types is in conflict with the biotype
#
#	  if MGD marker type != "pseudogene" then translate to "gene" and verify
#
#	  if MGD marker type is equal to all biotypes, then there is no conflict (Not Applicable)
#		else there is a conflict
#
#  Usage:
#	seqbiotype.py
#
#  Env Vars: Uses environment variables to determine Server and Database
#	  (DSQUERY and MGD).
#
#  Inputs: 1) mgd database
#          2) Configuration
#
#  Outputs: 1) log file
#
#  Exit Codes:
#
#  History
#
#  02/18/2010	lec
#	- TR 9239; update rawbiotype, _BiotypeConflict_key from SEQ_Marker_Cache
#

import sys
import os
import string
import db
import mgi_utils
import loadlib

#
# from configuration file
#
user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']

# database errors
DB_ERROR = 'A database error occured: '
DB_CONNECT_ERROR = 'Connection to the database failed: '

# sql update script
updateSQL = '''
update SEQ_Marker_Cache 
set _BiotypeConflict_key = %s, rawbiotype = "%s"
where _Marker_key = %s
and _Sequence_key = %s
'''

# list of marker/sequence, marker type
markerList = {}

# conflict types
yesconflictKey = 5420767
noconflictKey = 5420769

# 'gene', 'pseudogene' : MRK_Types
geneTypeKey = 1
pseudoTypeKey = 7

# Purpose: Initialize db.sql settings, lookups, and file descriptor
# Returns: Nothing
# Assumes: Nothing
# Effects: connects to and queries a database
# Throws: Nothing

def init ():

    global markerList
    
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)
    db.set_sqlLogFunction(db.sqlLogAll)

    print 'Initializing ...%s' % (mgi_utils.date())

    #
    # select all SEQ_Marker_Cache sequence/markers that contain
    # sequences in SEQ_GeneModel
    # where logicaldb in NCBI (59), Ensembl (60), VEGA (85)
    #

    results = db.sql('''
	  select distinct s._Sequence_key, s._Marker_key, s._Marker_Type_key, 
			  m._GMMarker_Type_key, m.rawBiotype
	  from SEQ_Marker_Cache s, SEQ_GeneModel m
	  where s._Sequence_key = m._Sequence_key
	  and _LogicalDB_key in (59, 60, 85)
	  ''', 'auto')
	  #and s._Marker_key = 192012

    for r in results:
	key = r['_Marker_key']
	value = r
	if not markerList.has_key(key):
	    markerList[key] = []
        markerList[key].append(value)

    return

# Purpose: 
# Returns: Nothing
# Assumes: Nothing
# Effects: 
# Throws: Nothing

def createSQL():

    execSQL = ''

    # by marker

    for m in markerList:

        typeList = []

	# for each sequence...

        for s in markerList[m]:

	    markerKey = s['_Marker_key']
	    mgdtypeKey = s['_Marker_Type_key']
	    gmtypeKey = s['_GMMarker_Type_key']

	    # if mgdtypeKey != pseudogene, then set to "gene"

	    if mgdtypeKey != pseudoTypeKey:
	        mgdtype = geneTypeKey

	    # store each type once

	    if typeList.count(mgdtypeKey) == 0:
	        typeList.append(mgdtypeKey)

	    if typeList.count(gmtypeKey) == 0:
	        typeList.append(gmtypeKey)

	# if more than one marker type appears in the list, 
	# then there is a conflict

	if len(typeList) > 1:
	    conflictKey = yesconflictKey
        else:
	    conflictKey = noconflictKey

	# now re-iterate thru the marker/sequences
	# and set the conflict key and raw biotype
	# all sequences for a given marker get the same
	# conflict key value

        for s in markerList[m]:
	    rawbiotype = s['rawBiotype']
	    sequenceKey = s['_Sequence_key']
	    execSQL = updateSQL % (conflictKey, rawbiotype, markerKey, sequenceKey)
	    print execSQL
            db.sql(execSQL, None)

    return

# Purpose: Perform cleanup steps for the script.
# Returns: Nothing
# Assumes: Nothing
# Effects: Nothing
# Throws: Nothing

def finalize():

    db.useOneConnection(0)
    return

#
# Main Routine
#
try:
    init()
    createSQL()
    finalize()
except db.connection_exc, message:
    error = '%s%s' % (DB_CONNECT_ERROR, message)
    sys.stderr.write(message)
    sys.exit(message)
except db.error, message:
    error = '%s%s' % (DB_ERROR, message)
    sys.stderr.write(message)
    sys.exit(message)

