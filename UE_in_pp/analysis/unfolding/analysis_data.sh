#!/bin/bash
# file name: firstcondor.sh

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.536
export HOME=/sphenix/u/egm2153
export MYINSTALL=$HOME/install
source $OPT_SPHENIX/bin/setup_local.sh $MYINSTALL
export SPHENIX=$MYINSTALL

export LD_LIBRARY_PATH=/sphenix/user/phnxbld/workarea/sPHENIX_SL7.3/almalinux-9.2/gcc-14.2.0/need_root_version/root-6.32.06/RooUnfold-3.0.5/build:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=/sphenix/user/phnxbld/workarea/sPHENIX_SL7.3/almalinux-9.2/gcc-14.2.0/need_root_version/root-6.32.06/RooUnfold-3.0.5/build:$ROOT_INCLUDE_PATH

root.exe -q -b -l 'analysis_data.C('$1',"'$2'")'
