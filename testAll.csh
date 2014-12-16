#!/bin/csh -f

echo "Running all seq cache loads"
foreach load ( seqcoord.csh seqmarcker.csh seqdescription.csh )
    setenv DB_TYPE sybase
    ./$load
    setenv DB_TYPE postgres
    ./$load
end

echo "Performing test"
python ${MGD_DBUTILS}/bin/comparePostgresTable.py seq_coord_cache seq_marker_cache seq_description_cache

echo "Tests successful"
