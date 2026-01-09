#!/bin/bash

lines=3
args="-l $lines --numeric-suffixes=0 --suffix-length=4 --additional-suffix=.list"
#split $args run28_jet10_dst_tracks.list "../mc_list_files/run28_jet10_dst_tracks_"
#split $args run28_jet30_dst_tracks.list "../mc_list_files/run28_jet30_dst_tracks_"
split $args run28_jet10_dst_trkr_seed.list "../mc_list_files/run28_jet10_dst_trkr_seed_"
split $args run28_jet30_dst_trkr_seed.list "../mc_list_files/run28_jet30_dst_trkr_seed_"