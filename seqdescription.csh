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

setenv LOG      ${CACHELOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

setenv TABLE	SEQ_Description_Cache

date | tee -a ${LOG}

# Create the bcp file

./seqdescription.py | tee -a ${LOG}

if ( -z ${TABLE}.bcp ) then
echo 'BCP Files are empty' | tee -a ${LOG}
exit 0
endif

# truncate table

${MGD_DBSCHEMADIR}/table/${TABLE}_truncate.object | tee -a ${LOG}

# Drop indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_drop.object | tee -a ${LOG}

# BCP new data into tables
cat ${MGD_DBPASSWORDFILE} | bcp ${MGD_DBNAME}..${TABLE} in ${CACHEDATADIR}/${TABLE}.bcp -c -t"${FIELDDELIM}" -e ${CACHEDATADIR}/${TABLE}.bcp.error -S${MGD_DBSERVER} -U${MGD_DBUSER} | tee -a ${LOG}

# Create indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_create.object | tee -a ${LOG}

date | tee -a ${LOG}
