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

cd ${CACHEDATADIR}

setenv LOG      ${CACHELOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

setenv TABLE	SEQ_Description_Cache
setenv COLDELIM	"&=&"

date | tee -a ${LOG}

# Create the bcp file

${CACHEINSTALLDIR}/seqdescription.py | tee -a ${LOG}

if ( -z ${TABLE}.bcp ) then
echo 'BCP Files are empty' | tee -a ${LOG}
exit 0
endif

# Allow bcp into database and truncate tables

${DBUTILSBINDIR}/turnonbulkcopy.csh ${DBSERVER} ${DBNAME} | tee -a ${LOG}
${SCHEMADIR}/table/${TABLE}_truncate.object | tee -a ${LOG}

# Drop indexes
${SCHEMADIR}/index/${TABLE}_drop.object | tee -a ${LOG}

# BCP new data into tables
cat ${DBPASSWORDFILE} | bcp ${DBNAME}..${TABLE} in ${TABLE}.bcp -c -t"${COLDELIM}" -S${DBSERVER} -U${DBUSER} | tee -a ${LOG}

# Create indexes
${SCHEMADIR}/index/${TABLE}_create.object | tee -a ${LOG}

date | tee -a ${LOG}
