#!/bin/csh -f

#
# Usage:  seqdescription.csh
#
# History
#
# 03/30/2004	lec
#	- JSAM
#

cd `dirname $0` && source ./Configuration

setenv LOG      ${CACHELOGSDIR}/`basename $0 .csh`.log
rm -rf ${LOG}
touch ${LOG}

setenv TABLE	SEQ_Description_Cache

setenv COLDELIM '\t'

date | tee -a ${LOG}

# Create the bcp file

./seqdescription.py | tee -a ${LOG}

if ( -z ${TABLE}.bcp ) then
echo 'BCP Files are empty' | tee -a ${LOG}
exit 0
endif

# truncate table

${SCHEMADIR}/table/${TABLE}_truncate.object | tee -a ${LOG}

# Drop indexes
${SCHEMADIR}/index/${TABLE}_drop.object | tee -a ${LOG}

# BCP new data into tables
${BCP_CMD} ${TABLE} ${CACHEDATADIR} ${TABLE}.bcp ${COLDELIM} ${LINEDELIM} ${PG_DB_SCHEMA} | tee -a ${LOG}

# Create indexes
${SCHEMADIR}/index/${TABLE}_create.object | tee -a ${LOG}

date | tee -a ${LOG}
