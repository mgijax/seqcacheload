#!/bin/csh -fx

#
# Usage:  seqcoord.csh
#
# History
#
# lec	10/23/2003
#

cd `dirname ${0` && source Configuration

cd ${SEQCACHEBCPDIR}

setenv LOG	${SEQCACHEBCPDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

setenv TABLE	SEQ_Coord_Cache

date >>& ${LOG}

# Create the bcp file

../seqcoord.py >>& ${LOG}

# Allow bcp into database and truncate tables

${DBUTILSBINDIR}/turnonbulkcopy.csh ${DBSERVER} ${DBNAME} >>& ${LOG}
${SCHEMADIR}/table/${TABLE}_truncate.object >>& ${LOG}

# Drop indexes
${SCHEMADIR}/index/${TABLE}_drop.object >>& ${LOG}

# BCP new data into tables
cat ${DBPASSWORDFILE} | bcp ${DBNAME}..${TABLE} in ${TABLE}.bcp -c -t\| -S${DBSERVER} -U${DBUSER} >>& ${LOG}

# Create indexes
${SCHEMADIR}/index/${TABLE}_create.object >>& ${LOG}

date >>& ${LOG}
