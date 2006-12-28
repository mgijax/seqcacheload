#!/bin/csh -fx

#
# Usage:  seqmarker.csh
#
# History
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

./seqmarker.py | tee -a ${LOG}

date | tee -a ${LOG}

if ( -z ${TABLE}.bcp ) then
echo 'BCP Files are empty' | tee -a ${LOG}
exit 0
endif

# truncate table

${MGD_DBSCHEMADIR}/table/${TABLE}_truncate.object | tee -a ${LOG}

# Drop indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_drop.object | tee -a ${LOG}

# BCP new data into tables
${MGI_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} ${TABLE} ${CACHEDATADIR} ${TABLE}.bcp ${COLDELIM} ${LINEDELIM} | tee -a ${LOG}

# Create indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_create.object | tee -a ${LOG}

date | tee -a ${LOG}
