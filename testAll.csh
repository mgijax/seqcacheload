#!/bin/csh -f

echo "Running all seq cache loads"
foreach load ( seqcoord.csh seqmarker.csh seqdescription.csh seqprobe.csh seqdummy.csh )
    setenv DB_TYPE sybase
    ./$load
    setenv DB_TYPE postgres
    ./$load
end

echo "Performing test"
python ${MGD_DBUTILS}/bin/comparePostgresTable.py seq_coord_cache seq_marker_cache seq_description_cache seq_probe_cache
python ${MGD_DBUTILS}/bin/comparePostgresTable.py -c acc_accession seq_sequence_raw seq_sequence seq_source_assoc

echo "Tests successful"
