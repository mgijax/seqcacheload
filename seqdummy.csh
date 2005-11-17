#!/bin/csh -fx

#
# Usage:  seqdummy.csh
#
# History
#
# lec	11/08/2005
#	- TR 7094/MGI 3.5
#
# lec	10/27/2005
#

cd `dirname $0` && source ./Configuration

cd ${CACHEDATADIR}

setenv LOG      ${CACHELOGSDIR}/`basename $0`.log
rm -rf ${LOG}
touch ${LOG}

date | tee -a ${LOG}

# Create the bcp file

${CACHEINSTALLDIR}/seqdummy.py | tee -a ${LOG}

if ( -z SEQ_Sequence.bcp ) then
echo 'BCP Files are empty' | tee -a ${LOG}
exit 0
endif

date | tee -a ${LOG}

set a=`wc -l SEQ_Sequence.bcp`
set b=`echo $a | cut -f1 -d " "`

if ( $b > 2000 ) then
    ${SCHEMADIR}/index/SEQ_Sequence_drop.object | tee -a ${LOG}
endif

# Drop index and triggers

${SCHEMADIR}/trigger/SEQ_Sequence_drop.object | tee -a ${LOG}
${SCHEMADIR}/trigger/SEQ_Source_Assoc_drop.object | tee -a ${LOG}
${SCHEMADIR}/trigger/ACC_Accession_drop.object | tee -a ${LOG}

# BCP new data into tables
cat ${DBPASSWORDFILE} | bcp ${DBNAME}..SEQ_Sequence in SEQ_Sequence.bcp -c -t\| -S${DBSERVER} -U${DBUSER} | tee -a ${LOG}
cat ${DBPASSWORDFILE} | bcp ${DBNAME}..SEQ_Source_Assoc in SEQ_Source_Assoc.bcp -c -t\| -S${DBSERVER} -U${DBUSER} | tee -a ${LOG}
cat ${DBPASSWORDFILE} | bcp ${DBNAME}..ACC_Accession in ACC_Accession.bcp -c -t\| -S${DBSERVER} -U${DBUSER} | tee -a ${LOG}

# Re-create index and triggers

if ( $b > 2000 ) then
    ${SCHEMADIR}/index/SEQ_Sequence_create.object | tee -a ${LOG}
endif

${SCHEMADIR}/trigger/SEQ_Sequence_create.object | tee -a ${LOG}
${SCHEMADIR}/trigger/SEQ_Source_Assoc_create.object | tee -a ${LOG}
${SCHEMADIR}/trigger/ACC_Accession_create.object | tee -a ${LOG}

date | tee -a ${LOG}
