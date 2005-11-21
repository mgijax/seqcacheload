#!/bin/csh -fx

#
# Usage:  seqcoord.csh
#
# History
#
# lec	10/23/2003
#

cd `dirname $0` && source ./Configuration

setenv LOG	${CACHELOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

setenv TABLE	SEQ_Coord_Cache

date | tee -a ${LOG}

# Create the bcp file

./seqcoord.py | tee -a ${LOG}

if ( -z ${TABLE}.bcp ) then
echo 'BCP Files are empty' | tee -a ${LOG}
exit 0
endif

# Truncate table

${MGD_DBSCHEMADIR}/table/${TABLE}_truncate.object | tee -a ${LOG}

# Drop indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_drop.object | tee -a ${LOG}

# BCP new data into tables
cat ${MGD_DBPASSWORDFILE} | bcp ${MGD_DBNAME}..${TABLE} in ${CACHEDATADIR}/${TABLE}.bcp -c -t"${FIELDDELIM}" -S${MGD_DBSERVER} -U${MGD_DBUSER} | tee -a ${LOG}

# Create indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_create.object | tee -a ${LOG}

date | tee -a ${LOG}
