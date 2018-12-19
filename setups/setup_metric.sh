#!/bin/bash

location=$1

declare -A arr
arr['NERSC']='/global/common/software/lsst/cori-haswell-gcc/stack/setup_w_2018_13-sims_2_7_0.sh'
arr['CCIN2P3']='/cvmfs/sw.lsst.eu/linux-x86_64/lsst_sims/sims_2_8_0/loadLSST.bash'
arr['mylaptop']='/cvmfs/sw.lsst.eu/linux-x86_64/lsst_sims/sims_2_8_0/loadLSST.bash'

#echo ${arr[@]}

if [ -z ${location} ]
then
    echo 'This script requires a parameter' ${location}
    echo 'Possible values are: ' ${!arr[@]}
    echo 'Please rerun this script with one of these values'
    return
fi


if [[ ! ${arr[$location]} ]]
then
    echo 'Problem here. We do not know where you would like to run'
    echo 'Maybe you should edit and modify this setup script'
else
    thescript=${arr[$location]}
    source ${thescript}
    setup lsst_sims

    export PYTHONPATH=${PWD}/SN_MAF/SN_Metrics:$PYTHONPATH
    export PYTHONPATH=${PWD}/SN_MAF/SN_Stackers:$PYTHONPATH
    export PYTHONPATH=${PWD}/SN_MAF/SN_Simulations:$PYTHONPATH
    export PYTHONPATH=${PWD}/SN_MAF/SN_Tools:$PYTHONPATH
    
fi
