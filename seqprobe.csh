#!/bin/csh -fx

#
# Usage:  seqprobe.csh
#
# History
#
# lec	10/23/2003
#

cd `dirname $0` && source ./Configuration

cd ${CACHEBCPDIR}

setenv LOG	${CACHEBCPDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

setenv TABLE	SEQ_Probe_Cache

date | tee -a ${LOG}

# Create the bcp file

../seqprobe.py | tee -a ${LOG}

if ( -z ${TABLE}.bcp ) then
echo 'BCP Files are empty' >>& $LOG
exit 0
endif

# Allow bcp into database and truncate tables

${DBUTILSBINDIR}/turnonbulkcopy.csh ${DBSERVER} ${DBNAME} | tee -a ${LOG}
${SCHEMADIR}/table/${TABLE}_truncate.object | tee -a ${LOG}

# Drop indexes
${SCHEMADIR}/index/${TABLE}_drop.object | tee -a ${LOG}

# BCP new data into tables
cat ${DBPASSWORDFILE} | bcp ${DBNAME}..${TABLE} in ${TABLE}.bcp -c -t\| -S${DBSERVER} -U${DBUSER} | tee -a ${LOG}

# Create indexes
${SCHEMADIR}/index/${TABLE}_create.object | tee -a ${LOG}

${DBUTILSBINDIR}/updateStatistics.csh ${DBSERVER} ${DBNAME} ${TABLE} | tee -a ${LOG}

date | tee -a ${LOG}
