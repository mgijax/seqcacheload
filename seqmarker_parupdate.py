#
# seqmarker_parupdate.py  
#####################################################################
#
#  Purpose: This script updates seq_marker_cache - removing the bogus
#       entries for par marker pairs (x and y) and sequences (GMs) which
#       both have the same NCBI ID associated with them.
#
#  See for details: 
#
#  Usage:
#	seqmarker_parupdate.py
#
#  Env Vars: Uses environment variables to determine Server and Database
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
# 03/24/2023    sc
#       FL2b project PAR Epic
#

import sys
import os
import mgi_utils
import loadlib
import db

db.setTrace()

#
# CONSTANTS
#

# 
debug = os.environ['SEQMARKER_DEBUG']

# column delimiter
DL = os.environ['COLDELIM']

# record delimiter
NL = '\n'

# database errors
DB_ERROR = 'A database error occured: '
DB_CONNECT_ERROR = 'Connection to the database failed: '


# Purpose: Initialize db.sql settings, lookups, and file descriptor
# Returns: Nothing
# Assumes: Nothing
# Effects: connects to and queries a database
# Throws: Nothing

def doDeletes():
    
    db.useOneConnection(1)

    #
    # using the map_gm_coord_cache_view find the seqmarker cache entries that
    # need to be updated
    #
    print('Initializing ...%s' % (mgi_utils.date()))

    # join back to the accession table via ID to get marker information for NCBI markers on 'XY', i
    # get the last char on the symbol as symbolChromosome
    db.sql(''' select a._object_key as _marker_key, m.symbol as markerSymbol, 
            UPPER(RIGHT(m.symbol, 1)) as symbolChromosome, m.chromosome as geneticChromosome, gm.*
        into temporary table ncbi
        from map_gm_coord_feature_view gm, acc_accession a, mrk_marker m
        where gm.seqID = a.accID
        and a._mgitype_key = 2
        and a._logicaldb_key in (55)
        and m.chromosome = 'XY' -- genetic
        and gm.genomicChromosome in ('X', 'Y')
        and a._object_key = m._marker_key''', None)

    # These need to be removed from seq_marker_cache by _marker_key/_sequence_key 
    results = db.sql('''select _marker_key, _sequence_key 
        from ncbi
        where (genomicChromosome = 'Y' and symbolChromosome = 'X')
        or (genomicChromosome = 'X' and symbolChromosome = 'Y')''', 'auto')

    # Now update seq_marker_cache
    print('Updating seq_marker_cache')
    print('len(results): %s' % len(results))
    for r in results:
        markerKey = r['_marker_key']
        sequenceKey = r['_sequence_key']
        cmd = '''delete from seq_marker_cache
            where _marker_key = %s
            and _sequence_key = %s''' % (markerKey, sequenceKey)
        print(cmd)
        db.sql(cmd, None)
        db.commit()

    return

#
# Main Routine
#
try:
    doDeletes()
    db.useOneConnection(0)
except db.connection_exc as message:
    error = '%s%s' % (DB_CONNECT_ERROR, message)
    sys.stderr.write(message)
    sys.exit(message)
except db.error as message:
    error = '%s%s' % (DB_ERROR, message)
    sys.stderr.write(message)
    sys.exit(message)

