#!/usr/local/bin/python

# $Name$

#
# Program: seqdummy.py
#
# Original Author: Lori Corbani
#
# Purpose:
#
#	To load new Sequences into:
#
#	SEQ_Sequence
#	SEQ_Source_Assoc
#	ACC_Accession
#
# Requirements Satisfied by This Program:
#
# Usage:
#	seqdummy.py
#
# Envvars:
#
# Outputs:
#
#       BCP files:
#
#       SEQ_Sequence.bcp                master Sequence records
#	SEQ_Source_Assoc.bcp		Sequence/Source Association records
#       ACC_Accession.bcp               Accession records
#
#       Diagnostics file of all input parameters and SQL commands
#       Error file
#
# Exit Codes:
#
# Assumes:
#
#	That no one else is adding such records to the database.
#
# Bugs:
#
# Implementation:
#

import sys
import os
import string
import getopt
import accessionlib
import db
import mgi_utils
import loadlib
import sourceloadlib

#globals

TAB = '\t'		# tab
CRT = '\n'		# carriage return/newline

seqFile = ''          # file descriptor
sourceFile = ''		# file descriptor
accFile = ''            # file descriptor

seqTable = 'SEQ_Sequence'
sourceTable = 'SEQ_Source_Assoc'
accTable = 'ACC_Accession'

seqFileName = seqTable + '.bcp'
sourceFileName = sourceTable + '.bcp'
accFileName = accTable + '.bcp'

seqKey = 0              # SEQ_Sequence._Sequence_key
assocKey = 0		# SEQ_Source_Assoc._Assoc_key
accKey = 0              # ACC_Accession._Accession_key
userKey = 0		# MGI_User._User_key

mgiTypeKey = 19		# Sequence
statusKey = 316345	# "Not Loaded" Sequence Status
mouseSourceKey = 47395
nonmouseSourceKey = 48166
createdBy = 'mgd_dbo'
notLoaded = 'Not Loaded'

typeDict = {}
qualityDict = {}
providerDict = {}

loaddate = loadlib.loaddate

# Purpose: prints error message and exits
# Returns: nothing
# Assumes: nothing
# Effects: exits with exit status
# Throws: nothing

def exit(
    status,          # numeric exit status (integer)
    message = None   # exit message (string)
    ):

    if message is not None:
        sys.stderr.write('\n' + str(message) + '\n')
 
    try:
	seqFile.close()
	sourceFile.close()
	accFile.close()
    except:
        pass

    db.useOneConnection(0)
    sys.exit(status)
 
# Purpose: process command line options
# Returns: nothing
# Assumes: nothing
# Effects: initializes global variables
#          exits if files cannot be opened
# Throws: nothing

def init():
    global seqFile, sourceFile, accFile
    global typeDict, qualityDict, providerDict
 
    db.useOneConnection(1)
 
    try:
        seqFile = open(seqFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % seqFileName)

    try:
        sourceFile = open(sourceFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % sourceFileName)

    try:
        accFile = open(accFileName, 'w')
    except:
        exit(1, 'Could not open file %s\n' % accFileName)

    # Log all SQL
    db.set_sqlLogFunction(db.sqlLogAll)

    results = db.sql('select _Term_key, term from VOC_Term_SequenceType_View', 'auto')
    for r in results:
        typeDict[r['term']] = r['_Term_key']

    results = db.sql('select _Term_key, term from VOC_Term_SequenceQuality_View', 'auto')
    for r in results:
        qualityDict[r['term']] = r['_Term_key']

    results = db.sql('select _Term_key, term from VOC_Term_SequenceProvider_View', 'auto')
    for r in results:
        providerDict[r['term']] = r['_Term_key']

    return

# Purpose:  sets global primary key variables
# Returns:  nothing
# Assumes:  nothing
# Effects:  sets global primary key variables
# Throws:   nothing

def setPrimaryKeys():

    global seqKey, assocKey, accKey, userKey

    results = db.sql('select maxKey = max(_Sequence_key) + 1 from %s' % (seqTable), 'auto')
    seqKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Assoc_key) + 1 from %s' % (sourceTable), 'auto')
    assocKey = results[0]['maxKey']

    results = db.sql('select maxKey = max(_Accession_key) + 1 from %s' % (accTable), 'auto')
    accKey = results[0]['maxKey']

    userKey = loadlib.verifyUser(createdBy, 0, None)

# Purpose:  processes data
# Returns:  nothing
# Assumes:  nothing
# Effects:  verifies and processes each line in the input file
# Throws:   nothing

def process():

    global seqKey, assocKey, accKey

    # generate table of all mouse molecular segments Acc IDs whose GenBank SeqIDs
    # are not represented as Sequence objects.

    db.sql('select distinct a.accID, a._LogicalDB_key ' + \
	'into #probeaccs1 ' + \
	'from ACC_Accession a, PRB_Probe p, PRB_Source ps ' + \
	'where a._MGIType_key = 3 ' + \
	'and a._LogicalDB_key = 9 ' + \
	'and a._Object_key = p._Probe_key ' + \
	'and p._Source_key = ps._Source_key ' + \
	'and ps._Organism_key = 1 ' + \
	'and not exists (select 1 from ACC_Accession s ' + \
	'where s._MGIType_key = 19 and s._LogicalDB_key = a._LogicalDB_key and s.accID = a.accID)', None)

    # generate table of all mouse marker Acc IDs whose GenBank, SWISSProt, RefSeq,
    # TIGR, DoTS, TrEMBL IDs are not represented as Sequence objects.

    db.sql('select distinct a.accID, a._LogicalDB_key ' + \
	'into #markeraccs1 ' + \
	'from ACC_Accession a, MRK_Marker m ' + \
	'where a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key in (9,13,27,35,36,41,53) ' + \
	'and a._Object_key = m._Marker_key ' + \
	'and m._Organism_key = 1 ' + \
	'and not exists (select 1 from ACC_Accession s ' + \
	'where s._MGIType_key = 19 and s._LogicalDB_key = a._LogicalDB_key and s.accID = a.accID)', None)

    # generate table of all non-mouse molecular segments Acc IDs whose GenBank SeqIDs
    # are not represented as Sequence objects.

    db.sql('select distinct a.accID, a._LogicalDB_key ' + \
	'into #probeaccs2 ' + \
	'from ACC_Accession a, PRB_Probe p, PRB_Source s ' + \
	'where a._MGIType_key = 3 ' + \
	'and a._LogicalDB_key = 9 ' + \
	'and a._Object_key = p._Probe_key ' + \
	'and p._Source_key = s._Source_key ' + \
	'and s._Organism_key != 1 ' + \
	'and not exists (select 1 from ACC_Accession s ' + \
	'where s._MGIType_key = 19 and s._LogicalDB_key = a._LogicalDB_key and s.accID = a.accID)', None)

    # generate table of all non-mouse marker Acc IDs whose GenBank, SWISSProt, RefSeq,
    # TIGR, DoTS, TrEMBL IDs are not represented as Sequence objects.

    db.sql('select distinct a.accID, a._LogicalDB_key ' + \
	'into #markeraccs2 ' + \
	'from ACC_Accession a, MRK_Marker m ' + \
	'where a._MGIType_key = 2 ' + \
	'and a._LogicalDB_key in (9,13,27,35,36,41,53) ' + \
	'and a._Object_key = m._Marker_key ' + \
	'and m._Organism_key != 1 ' + \
	'and not exists (select 1 from ACC_Accession s ' + \
	'where s._MGIType_key = 19 and s._LogicalDB_key = a._LogicalDB_key and s.accID = a.accID)', None)

    # union these 4 sets together to form one unique set

    db.sql('select accID, _LogicalDB_key, isMouse = 1 ' + \
	'into #allaccs ' + \
	'from #probeaccs1 ' + \
	'union ' + \
	'select accID, _LogicalDB_key, isMouse = 1 ' + \
	'from #markeraccs1 ' + \
	'union ' + \
	'select accID, _LogicalDB_key, isMouse = 0 ' + \
	'from #probeaccs2 ' + \
	'union ' + \
	'select accID, _LogicalDB_key, isMouse = 0 ' + \
	'from #markeraccs2', None)

    results = db.sql('select * from #allaccs', 'auto')
    for r in results:

	accID = r['accID']
	logicalDB = r['_LogicalDB_key']

	if r['isMouse'] == 1:
	    sourceKey = mouseSourceKey
        else:
	    sourceKey = nonmouseSourceKey

	virtual = 1

        # change values for specific cases

        if logicalDB == 9:
            typeKey = typeDict["Not Loaded"]
            qualityKey = qualityDict["Not Loaded"]
            providerKey = providerDict["GenBank/EMBL/DDBJ"]
            virtual = 0

        elif logicalDB == 27:     # RefSeq
            typeKey = typeDict["RNA"]
            qualityKey = qualityDict["High"]
            providerKey = providerDict["RefSeq"]

        elif logicalDB == 35:     # TIGR
            typeKey = typeDict["RNA"]
            qualityKey = qualityDict["Low"]
            providerKey = providerDict["TIGR Mouse Gene Index"]

        elif logicalDB == 36:     # DoTS
            typeKey = typeDict["RNA"]
            qualityKey = qualityDict["Low"]
            providerKey = providerDict["DoTS"]

        elif logicalDB == 13:     # SwissProt
            typeKey = typeDict["Polypeptide"]
            qualityKey = qualityDict["High"]
            providerKey = providerDict["SWISS-PROT"]

        elif logicalDB == 41:     # TrEMBL
            typeKey = typeDict["Polypeptide"]
            qualityKey = qualityDict["Low"]
            providerKey = providerDict["TrEMBL"]

        elif logicalDB == 53:     # NIA Mouse Gene Index
            typeKey = typeDict["RNA"]
            qualityKey = qualityDict["Low"]
            providerKey = providerDict["NIA Mouse Gene Index"]

        seqFile.write('%d|%d|%d|%d|%d|||||%d|%s|%s|%s|%s|%s|%s|%s|%s||%s|%s|%s|%s|%s|%s\n' \
            % (seqKey, typeKey, qualityKey, statusKey, providerKey, virtual, \
               notLoaded, notLoaded, notLoaded, notLoaded, notLoaded, notLoaded, notLoaded, notLoaded, \
	       loaddate, loaddate, userKey, userKey, loaddate, loaddate))

        sourceFile.write('%d|%d|%d|%s|%s|%s|%s\n' % (assocKey, seqKey, sourceKey, userKey, userKey, loaddate, loaddate))

	prefixPart, numericPart = accessionlib.split_accnum(accID)
        accFile.write('%s|%s|%s|%s|%s|%d|%d|0|1|%s|%s|%s|%s\n' \
                % (accKey, accID, prefixPart, numericPart, logicalDB, seqKey, mgiTypeKey, userKey, userKey, loaddate, loaddate))

        seqKey = seqKey + 1
	assocKey = assocKey + 1
	accKey = accKey + 1

#
# Main
#

init()
setPrimaryKeys()
process()
exit(0)

