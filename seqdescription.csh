#!/bin/csh -fx

#
# Usage:  seqdescription.csh
#
# History
#
# 03/30/2004	lec
#	- JSAM
#

cd `dirname $0` && source ./Configuration

cd ${SEQCACHEBCPDIR}

setenv LOG	${SEQCACHEBCPDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

setenv TABLE	SEQ_Description_Cache
setenv COLDELIM	"&=&"

date >>& ${LOG}

# Create the bcp file

../seqdescription.py >>& ${LOG}

# Allow bcp into database and truncate tables

${DBUTILSBINDIR}/turnonbulkcopy.csh ${DBSERVER} ${DBNAME} >>& ${LOG}
${SCHEMADIR}/table/${TABLE}_truncate.object >>& ${LOG}

# Drop indexes
${SCHEMADIR}/index/${TABLE}_drop.object >>& ${LOG}

# BCP new data into tables
cat ${DBPASSWORDFILE} | bcp ${DBNAME}..${TABLE} in ${TABLE}.bcp -c -t"${COLDELIM}" -S${DBSERVER} -U${DBUSER} >>& ${LOG}

# Create indexes
${SCHEMADIR}/index/${TABLE}_create.object >>& ${LOG}

date >>& ${LOG}
