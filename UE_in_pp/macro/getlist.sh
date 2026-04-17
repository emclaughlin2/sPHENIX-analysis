#!/bin/bash

DST_CALO_LIST="run28_herwig_jet50_dst_calo_cluster.list"
DST_TRUTH_LIST="run28_herwig_jet50_g4hits.list"
DST_GLOBAL_LIST="run28_herwig_jet50_dst_mbd_epd.list"
DST_TRUTH_JET_LIST="run28_herwig_jet50_dst_truth_jet.list"
DST_TRACKS_LIST="run28_herwig_jet50_dst_tracks.list"

# Run SQL queries and redirect results to files
psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 28 AND dsttype = 'DST_CALO_CLUSTER' AND filename LIKE '%Herwig_Jet50-0000000028%' ORDER BY segment;" > "$DST_CALO_LIST"
psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 28 AND dsttype = 'G4Hits' AND filename LIKE '%Herwig_Jet50-0000000028%' ORDER BY segment;" > "$DST_TRUTH_LIST"
psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 28 AND dsttype = 'DST_MBD_EPD' AND filename LIKE '%Herwig_Jet50-0000000028%' ORDER BY segment;" > "$DST_GLOBAL_LIST"
psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 28 AND dsttype = 'DST_TRUTH_JET' AND filename LIKE '%Herwig_Jet50-0000000028%' ORDER BY segment;" > "$DST_TRUTH_JET_LIST"
#psql -h sphnxdbmaster.sdcc.bnl.gov -d FileCatalog -t -A -q -F" " --command "SELECT filename FROM datasets WHERE runnumber = 22 AND dsttype = 'DST_TRACKS' AND filename LIKE '%Herwig_Jet50-0000000022%' ORDER BY segment;" > "$DST_TRACKS_LIST"
