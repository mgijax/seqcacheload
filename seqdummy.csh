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

# Drop index and triggers

#${MGD_DBSCHEMADIR}/index/SEQ_Sequence_drop.object | tee -a ${LOG}
${MGD_DBSCHEMADIR}/trigger/SEQ_Sequence_drop.object | tee -a ${LOG}
${MGD_DBSCHEMADIR}/trigger/SEQ_Source_Assoc_drop.object | tee -a ${LOG}
${MGD_DBSCHEMADIR}/trigger/ACC_Accession_drop.object | tee -a ${LOG}

# BCP new data into tables
cat ${MGD_DBPASSWORDFILE} | bcp ${MGD_DBNAME}..SEQ_Sequence in SEQ_Sequence.bcp -c -t\| -S${MGD_DBSERVER} -U${MGD_DBUSER} | tee -a ${LOG}
cat ${MGD_DBPASSWORDFILE} | bcp ${MGD_DBNAME}..SEQ_Sequence_Raw in SEQ_Sequence_Raw.bcp -c -t\| -S${MGD_DBSERVER} -U${MGD_DBUSER} | tee -a ${LOG}
cat ${MGD_DBPASSWORDFILE} | bcp ${MGD_DBNAME}..SEQ_Source_Assoc in SEQ_Source_Assoc.bcp -c -t\| -S${MGD_DBSERVER} -U${MGD_DBUSER} | tee -a ${LOG}
cat ${MGD_DBPASSWORDFILE} | bcp ${MGD_DBNAME}..ACC_Accession in ACC_Accession.bcp -c -t\| -S${MGD_DBSERVER} -U${MGD_DBUSER} | tee -a ${LOG}

# Re-create index and triggers

#${MGD_DBSCHEMADIR}/index/SEQ_Sequence_create.object | tee -a ${LOG}
${MGD_DBSCHEMADIR}/trigger/SEQ_Sequence_create.object | tee -a ${LOG}
${MGD_DBSCHEMADIR}/trigger/SEQ_Source_Assoc_create.object | tee -a ${LOG}
${MGD_DBSCHEMADIR}/trigger/ACC_Accession_create.object | tee -a ${LOG}

date | tee -a ${LOG}