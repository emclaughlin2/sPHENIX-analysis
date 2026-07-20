#!/bin/bash
# file name: firstcondor.sh

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.536
export HOME=/sphenix/u/egm2153
export MYINSTALL=$HOME/install
source $OPT_SPHENIX/bin/setup_local.sh $MYINSTALL
export SPHENIX=$MYINSTALL
root.exe -q -b -l 'analysis_sim_primary_particle_effects.C("'$1'",'$2','$3','$4',"'$5'")'
