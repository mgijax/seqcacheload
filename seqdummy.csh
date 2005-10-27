#!/bin/csh -fx

#
# Usage:  seqdummy.csh
#
# History
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

date | tee -a ${LOG}

# Allow bcp into database

${DBUTILSBINDIR}/turnonbulkcopy.csh ${DBSERVER} ${DBNAME} | tee -a ${LOG}

# BCP new data into tables
cat ${DBPASSWORDFILE} | bcp ${DBNAME}..SEQ_Sequence in SEQ_Sequence.bcp -c -t\| -S${DBSERVER} -U${DBUSER} | tee -a ${LOG}
cat ${DBPASSWORDFILE} | bcp ${DBNAME}..SEQ_Source_Assoc in SEQ_Source_Assoc.bcp -c -t\| -S${DBSERVER} -U${DBUSER} | tee -a ${LOG}
cat ${DBPASSWORDFILE} | bcp ${DBNAME}..ACC_Accession in ACC_Accession.bcp -c -t\| -S${DBSERVER} -U${DBUSER} | tee -a ${LOG}

date | tee -a ${LOG}
