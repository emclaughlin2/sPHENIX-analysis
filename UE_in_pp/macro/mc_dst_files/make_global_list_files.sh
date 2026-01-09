#!/bin/bash

lines=3
args="-l $lines --numeric-suffixes=0 --suffix-length=4 --additional-suffix=.list"
split $args run28_jet10_dst_global.list "../mc_list_files/run28_jet10_dst_global_"
split $args run28_jet30_dst_global.list "../mc_list_files/run28_jet30_dst_global_"