
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
#	SEQ_Sequence_Raw
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
import accessionlib
import mgi_utils
import loadlib
import db

db.setTrace()

#globals

NL = '\n'
DL = os.environ['COLDELIM']
datadir = os.environ['CACHEDATADIR']

seqFile = ''          	# file descriptor
rawFile = ''		# file descriptor
sourceFile = ''		# file descriptor
accFile = ''            # file descriptor

seqTable = 'SEQ_Sequence'
rawTable = 'SEQ_Sequence_Raw'
sourceTable = 'SEQ_Source_Assoc'
accTable = 'ACC_Accession'

seqFileName = datadir + '/' + seqTable + '.bcp'
rawFileName = datadir + '/' + rawTable + '.bcp'
sourceFileName = datadir + '/' + sourceTable + '.bcp'
accFileName = datadir + '/' + accTable + '.bcp'

seqKey = 0              # SEQ_Sequence._Sequence_key
assocKey = 0		# SEQ_Source_Assoc._Assoc_key
accKey = 0              # ACC_Accession._Accession_key
userKey = 0		# MGI_User._User_key

mgiTypeKey = 19		# Sequence
statusKey = 316345	# "Not Loaded" Sequence Status
mouseSourceKey = 47395
nonmouseSourceKey = 48166
notLoaded = 'Not Loaded'

loaddate = loadlib.loaddate


def init():
    """
    Initialize database connection
    Open output files
    """
    global seqFile, rawFile, sourceFile, accFile
 
    db.useOneConnection(1)
 
    seqFile = open(seqFileName, 'w')

    rawFile = open(rawFileName, 'w')

    sourceFile = open(sourceFileName, 'w')

    accFile = open(accFileName, 'w')


def setPrimaryKeys():
    """
    Assign global primary key variables
        using max keys from database
    """

    global seqKey, assocKey, accKey, userKey

    results = db.sql("select max(_Sequence_key) + 1 as maxKey from %s" % (seqTable), "auto")
    seqKey = results[0]["maxKey"]

    results = db.sql('''select nextval('seq_source_assoc_seq') as nextKey''', 'auto')
    assocKey = results[0]["nextKey"]

    results = db.sql("select max(_Accession_key) + 1 as maxKey from %s" % (accTable), "auto")
    accKey = results[0]["maxKey"]

    userKey = loadlib.verifyUser(os.environ['MGD_DBUSER'], 1, None)


def process():
    """
    Query database to determine if dummy sequences need to
        be created
    Generates appropriate BCP files
    """

    global seqKey, assocKey, accKey

    # generate table of all mouse molecular segments Acc IDs whose GenBank SeqIDs
    # are not represented as Sequence objects.

    db.sql("""select a.accID, a._LogicalDB_key, ps._Organism_key 
        INTO TEMPORARY TABLE probeaccs1 
        from ACC_Accession a, PRB_Probe p, PRB_Source ps 
        where a._MGIType_key = 3 
        and a._LogicalDB_key = 9 
        and a._Object_key = p._Probe_key 
        and p._Source_key = ps._Source_key 
        and ps._Organism_key = 1 
        and not exists (select 1 from ACC_Accession s 
            where s._MGIType_key = 19 
                and s._LogicalDB_key = a._LogicalDB_key 
                and lower(s.accID) = lower(a.accID)
        )""", None)

    # generate table of all mouse marker Acc IDs whose GenBank, SWISSProt, RefSeq,
    # TrEMBL IDs are not represented as Sequence objects.

    db.sql("""select a.accID, a._LogicalDB_key, m._Organism_key 
        INTO TEMPORARY TABLE markeraccs1 
        from ACC_Accession a, MRK_Marker m 
        where a._MGIType_key = 2 
        and a._LogicalDB_key in (9,13,27,41) 
        and a._Object_key = m._Marker_key 
        and m._Organism_key = 1 
        and m._Marker_Status_key in (1,2)
        and not exists (select 1 from ACC_Accession s 
            where s._MGIType_key = 19 
                and s._LogicalDB_key = a._LogicalDB_key 
                and lower(s.accID) = lower(a.accID)
        )""", None)

    # generate table of all non-mouse molecular segments Acc IDs whose GenBank SeqIDs
    # are not represented as Sequence objects.

    db.sql("""select a.accID, a._LogicalDB_key, s._Organism_key 
        INTO TEMPORARY TABLE probeaccs2 
        from ACC_Accession a, PRB_Probe p, PRB_Source s 
        where a._MGIType_key = 3 
        and a._LogicalDB_key = 9 
        and a._Object_key = p._Probe_key 
        and p._Source_key = s._Source_key 
        and s._Organism_key != 1 
        and not exists (select 1 from ACC_Accession s 
            where s._MGIType_key = 19 
                and s._LogicalDB_key = a._LogicalDB_key 
                and lower(s.accID) = lower(a.accID)
        )""", None)

    # generate table of all non-mouse marker Acc IDs whose GenBank, SWISSProt, RefSeq,
    # TrEMBL IDs are not represented as Sequence objects.

    db.sql("""select a.accID, a._LogicalDB_key, m._Organism_key 
        INTO TEMPORARY TABLE markeraccs2 
        from ACC_Accession a, MRK_Marker m 
        where a._MGIType_key = 2 
        and a._LogicalDB_key in (9,13,27,41) 
        and a._Object_key = m._Marker_key 
        and m._Organism_key != 1 
        and not exists (select 1 from ACC_Accession s 
            where s._MGIType_key = 19 
                and s._LogicalDB_key = a._LogicalDB_key 
                and lower(s.accID) = lower(a.accID)
        )""", None)

    # union these 4 sets together to form one unique set

    db.sql('select accID, _LogicalDB_key, _Organism_key ' + \
        'INTO TEMPORARY TABLE allaccs ' + \
        'from probeaccs1 ' + \
        'union ' + \
        'select accID, _LogicalDB_key, _Organism_key ' + \
        'from markeraccs1 ' + \
        'union ' + \
        'select accID, _LogicalDB_key, _Organism_key ' + \
        'from probeaccs2 ' + \
        'union ' + \
        'select accID, _LogicalDB_key, _Organism_key ' + \
        'from markeraccs2', None)

    results = db.sql('select * from allaccs', 'auto')
    for r in results:

        accID = r['accID']
        logicalDB = r['_LogicalDB_key']
        organism = r['_Organism_key']

        if organism == 1:
            sourceKey = mouseSourceKey
        else:
            sourceKey = nonmouseSourceKey

        virtual = 1

        # change values for specific cases

        # types:  316347 (DNA), 316346 (RNA), 316348 (polypeptide), 316349 (not loaded)
        # quality:  316338 (high), 316339 (medium), 316340 (low), 316341 (not loaded)
        # provider: 316380 (GenBank/EMBL/DDBJ), 316372 (RefSeq)
        #           316384 (SwissProt), 316385 (TrEMBL)

        if logicalDB == 9:	# GenBank
            typeKey = 316349
            qualityKey = 316341
            providerKey = 316380
            virtual = 0

        elif logicalDB == 27:     # RefSeq
            typeKey = 316349
            qualityKey = 316338
            providerKey = 316372

        elif logicalDB == 13:     # SwissProt
            typeKey = 316348
            qualityKey = 316338
            providerKey = 316384

        elif logicalDB == 41:     # TrEMBL
            typeKey = 316348
            qualityKey = 316340
            providerKey = 316385

        seqFile.write(mgi_utils.prvalue(seqKey) + DL + \
                mgi_utils.prvalue(typeKey) + DL + \
                mgi_utils.prvalue(qualityKey) + DL + \
                mgi_utils.prvalue(statusKey) + DL + \
                mgi_utils.prvalue(providerKey) + DL + \
                mgi_utils.prvalue(organism) + DL + \
                DL + DL + DL + DL + \
                mgi_utils.prvalue(virtual) + DL + \
                DL + \
                loaddate + DL + loaddate + DL + \
                str(userKey) + DL + str(userKey) + DL + \
                loaddate + DL + loaddate + NL)

        rawFile.write(mgi_utils.prvalue(seqKey) + DL + \
                notLoaded + DL + \
                notLoaded + DL + \
                notLoaded + DL + \
                notLoaded + DL + \
                notLoaded + DL + \
                notLoaded + DL + \
                notLoaded + DL + \
                notLoaded + DL + \
                str(userKey) + DL + str(userKey) + DL + \
                loaddate + DL + loaddate + NL)

        sourceFile.write(mgi_utils.prvalue(assocKey) + DL + \
                mgi_utils.prvalue(seqKey) + DL + \
                mgi_utils.prvalue(sourceKey) + DL + \
                str(userKey) + DL + str(userKey) + DL + \
                loaddate + DL + loaddate + NL)

        prefixPart, numericPart = accessionlib.split_accnum(accID)
        accFile.write(mgi_utils.prvalue(accKey) + DL + \
                mgi_utils.prvalue(accID) + DL + \
                mgi_utils.prvalue(prefixPart) + DL + \
                mgi_utils.prvalue(numericPart) + DL + \
                mgi_utils.prvalue(logicalDB) + DL + \
                mgi_utils.prvalue(seqKey) + DL + \
                mgi_utils.prvalue(mgiTypeKey) + DL + \
                '0' + DL + \
                '1' + DL + \
                str(userKey) + DL + str(userKey) + DL + \
                loaddate + DL + loaddate + NL)

        seqKey = seqKey + 1
        assocKey = assocKey + 1
        accKey = accKey + 1


if __name__ == "__main__":
    try:
        init()
        setPrimaryKeys()
        process()

    finally:
        # always close output files
        seqFile.close()
        rawFile.close()
        sourceFile.close()
        accFile.close()
