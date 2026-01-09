#!/bin/bash

DST_TRKR_SEED_LIST="run28_jet10_dst_trkr_seed.list"
DST_TRACKS_LIST="run28_jet10_dst_tracks.list"

# Run SQL queries and redirect results to files
psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 28 AND dsttype = 'DST_TRACKSEEDS' AND filename LIKE '%pythia8_Jet10-0000000028%' ORDER BY segment;" > "$DST_TRKR_SEED_LIST"
psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 28 AND dsttype = 'DST_TRACKS' AND filename LIKE '%pythia8_Jet10-0000000028%' ORDER BY segment;" > "$DST_TRACKS_LIST"
