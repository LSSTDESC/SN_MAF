#!/bin/bash

location=$1
script_loc=$2

source ${PWD}/sn_maf/setups/setup_metric.sh ${location} ${script_loc}

# now update python path
echo "updating python path"
export SN_UTILS_DIR=${PWD}/sn_utils
export SALT2_DIR=${PWD}/sn_utils/SALT2_Files
export PYTHONPATH=${thedir}:$PYTHONPATH
