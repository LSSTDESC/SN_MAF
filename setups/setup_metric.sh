#!/bin/bash

#setup lsst_sims

#source /global/common/software/lsst/cori-haswell-gcc/stack/setup_w_2018_13-sims_2_7_0.sh
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_sims/sims_2_8_0/loadLSST.bash
setup lsst_sims

export PYTHONPATH=${PWD}/SN_MAF/SN_Metrics:$PYTHONPATH
export PYTHONPATH=${PWD}/SN_MAF/SN_Stackers:$PYTHONPATH
export PYTHONPATH=${PWD}/SN_MAF/SN_Simulations:$PYTHONPATH
export PYTHONPATH=${PWD}/SN_MAF/SN_Tools:$PYTHONPATH

