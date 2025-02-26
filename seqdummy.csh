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

${PYTHON} ./seqdummy.py | tee -a ${LOG}

if ( ! -s SEQ_Sequence.bcp ) then
echo 'BCP Files are empty : done' | tee -a ${LOG}
exit 0
endif

echo 'BCP Files are not empty : continue' | tee -a ${LOG}
date | tee -a ${LOG}

set a=`wc -l ${CACHEDATADIR}/SEQ_Sequence.bcp`
set b=`echo $a | cut -f1 -d " "`

# Drop index and triggers

if ( $b > 3000 ) then
    ${SCHEMADIR}/index/SEQ_Sequence_drop.object | tee -a ${LOG}
endif

# Drop index and triggers

# BCP new data into tables
${BCP_CMD} SEQ_Sequence ${CACHEDATADIR} SEQ_Sequence.bcp ${COLDELIM} ${LINEDELIM} ${PG_DB_SCHEMA} | tee -a ${LOG}
${BCP_CMD} SEQ_Sequence_Raw ${CACHEDATADIR} SEQ_Sequence_Raw.bcp ${COLDELIM} ${LINEDELIM} ${PG_DB_SCHEMA} | tee -a ${LOG}
${BCP_CMD} SEQ_Source_Assoc ${CACHEDATADIR} SEQ_Source_Assoc.bcp ${COLDELIM} ${LINEDELIM} ${PG_DB_SCHEMA} | tee -a ${LOG}
${BCP_CMD} ACC_Accession ${CACHEDATADIR} ACC_Accession.bcp ${COLDELIM} ${LINEDELIM} ${PG_DB_SCHEMA} | tee -a ${LOG}

# update serialization on seq_source_assoc
cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 | tee -a ${LOG}
select setval('seq_source_assoc_seq', (select max(_Assoc_key) from SEQ_Source_Assoc));
EOSQL

# Re-create index and triggers

if ( $b > 3000 ) then
    ${SCHEMADIR}/index/SEQ_Sequence_create.object | tee -a ${LOG}
endif

date | tee -a ${LOG}
