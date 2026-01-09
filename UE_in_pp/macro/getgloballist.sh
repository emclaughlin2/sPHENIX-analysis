#!/bin/bash

DST_GLOBAL_LIST="run28_jet30_dst_global.list"

# Run SQL queries and redirect results to files
psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 28 AND dsttype = 'DST_GLOBAL' AND filename LIKE '%pythia8_Jet30-0000000028%' ORDER BY segment;" > "$DST_GLOBAL_LIST"