#!/bin/csh -f

#
# Usage:  seqdummy.csh
#
# History
#
# lec	03/10/2011
#	- trigger SEQ_Source_Assoc has been removed from the system
#
# lec	11/08/2005
#	- TR 7094/MGI 3.5
#
# lec	10/27/2005
#

cd `dirname $0` && source ./Configuration

setenv LOG      ${CACHELOGSDIR}/`basename $0 .csh`.log
rm -rf ${LOG}
touch ${LOG}

date | tee -a ${LOG}

# Create the bcp file

./seqdummy.py | tee -a ${LOG}

if ( -z SEQ_Sequence.bcp ) then
echo 'BCP Files are empty' | tee -a ${LOG}
exit 0
endif

if ( -z SEQ_Sequence.bcp ) then
echo 'BCP Files are empty' | tee -a ${LOG}
exit 0
endif

date | tee -a ${LOG}

set a=`wc -l ${CACHEDATADIR}/SEQ_Sequence.bcp`
set b=`echo $a | cut -f1 -d " "`

# Drop index and triggers

if ( $b > 3000 ) then
    ${MGD_DBSCHEMADIR}/index/SEQ_Sequence_drop.object | tee -a ${LOG}
endif

# Drop index and triggers

${MGD_DBSCHEMADIR}/trigger/SEQ_Sequence_drop.object | tee -a ${LOG}
${MGD_DBSCHEMADIR}/trigger/ACC_Accession_drop.object | tee -a ${LOG}

# BCP new data into tables
${MGI_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} SEQ_Sequence ${CACHEDATADIR} SEQ_Sequence.bcp ${COLDELIM} ${LINEDELIM} | tee -a ${LOG}
${MGI_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} SEQ_Sequence_Raw ${CACHEDATADIR} SEQ_Sequence_Raw.bcp ${COLDELIM} ${LINEDELIM} | tee -a ${LOG}
${MGI_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} SEQ_Source_Assoc ${CACHEDATADIR} SEQ_Source_Assoc.bcp ${COLDELIM} ${LINEDELIM} | tee -a ${LOG}
${MGI_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} ACC_Accession ${CACHEDATADIR} ACC_Accession.bcp ${COLDELIM} ${LINEDELIM} | tee -a ${LOG}

# Re-create index and triggers

if ( $b > 3000 ) then
    ${MGD_DBSCHEMADIR}/index/SEQ_Sequence_create.object | tee -a ${LOG}
endif

${MGD_DBSCHEMADIR}/trigger/SEQ_Sequence_create.object | tee -a ${LOG}
${MGD_DBSCHEMADIR}/trigger/SEQ_Source_Assoc_create.object | tee -a ${LOG}
${MGD_DBSCHEMADIR}/trigger/ACC_Accession_create.object | tee -a ${LOG}

date | tee -a ${LOG}
