#!/bin/csh -f

#
# Usage:  seqmarker.csh
#
# History
#
# lec	02/18/2010
#	- TR9239; rawbiotype, _BiotypeConflict_key
#
# lec	10/23/2003
#

cd `dirname $0` && source ./Configuration

setenv LOG      ${CACHELOGSDIR}/`basename $0 .csh`.log
rm -rf ${LOG}
touch ${LOG}

setenv TABLE	SEQ_Marker_Cache

date | tee -a ${LOG}

# Create the bcp file

${PYTHON} ./seqmarker.py | tee -a ${LOG}

date | tee -a ${LOG}

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
