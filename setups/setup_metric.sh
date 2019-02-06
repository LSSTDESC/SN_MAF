#!/bin/bash

location=$1
script_loc=$2

declare -A arr
#arr['NERSC']='/global/common/software/lsst/cori-haswell-gcc/stack/setup_w_2018_13-sims_2_7_0.sh'
arr['NERSC']='/global/common/software/lsst/cori-haswell-gcc/stack/setup_w_2018_19-sims_2_8_0.sh'
arr['CCIN2P3']='/cvmfs/sw.lsst.eu/linux-x86_64/lsst_sims/sims_2_8_0/loadLSST.bash'
arr['mylaptop']='/cvmfs/sw.lsst.eu/linux-x86_64/lsst_sims/sims_2_8_0/loadLSST.bash'

if [ ! -z ${script_loc} ]
then
    arr[$1]=$2
fi
#echo ${arr[@]}

if [ -z ${location} ]
then
    echo 'This script requires at least one  parameter' ${location}
    echo 'Possible values are: ' ${!arr[@]}
    echo 'Please rerun this script with one of these values or use:'
    echo 'python SN_MAF/setups/setup_metric.sh MYENV full_path_to_release_stack'
    echo 'where the stack release should include the lsst_sim package'
    return
fi


if [[ ! ${arr[$location]} ]]
then
    echo 'Problem here. We do not know where you would like to run'
    echo 'and which setup (stack) script you would like to use' 
    echo 'Try the following command'
    echo 'python sn_maf/setups/setup_metric.sh MYENV full_path_to_release_stack'
    echo 'where the stack release should include the lsst_sim package'
else
    thescript=${arr[$location]}
    source ${thescript}
    setup lsst_sims
    export PYTHONPATH=${PWD}:$PYTHONPATH
    #export PYTHONPATH=${PWD}/sn_maf/sn_metrics:$PYTHONPATH
    #export PYTHONPATH=${PWD}/sn_maf/sn_stackers:$PYTHONPATH
    #export PYTHONPATH=${PWD}/sn_maf/sn_simulations:$PYTHONPATH
    #export PYTHONPATH=${PWD}/sn_maf/sn_tools:$PYTHONPATH
    #export PYTHONPATH=${PWD}/sn_utils:$PYTHONPATH

    #checking whether hdf5 is accessible localy or not
    thedir=${PWD}/lib/python3.6/site-packages/
    lib='h5py'
    echo $thedir
    if [ -d ${thedir}$lib ]
    then
	echo $lib 'already installed -> updating PYTHONPATH'
    else
	echo $lib 'not installed -> installing with pip'
	pip install --prefix=${PWD} ${lib}==2.7.1
    fi
    export PYTHONPATH=${thedir}:$PYTHONPATH
fi
