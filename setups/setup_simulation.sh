#!/bin/bash

location=$1
script_loc=$2

source ${PWD}/sn_maf/setups/setup_metric.sh ${location} ${script_loc}


thedir=${PWD}/lib/python3.6/site-packages/
if [ ! -d ${PWD}/sn_catalog_Simulations ]
then
    if [ ! -d  "sn_catalog_simulations" ]; then
	echo "Seems that the simulation packages you need are not installed."
	echo "Will do it for you..."
	where='https://github.com/lsstdesc'
	echo "...taking packages from "${where}

	for pack in sn_catalog_simulations sn_utils
	do
	    echo 'Cloning' $pack
	    git clone $where/$pack.git
	    echo "Moving to the dev branch"
	    cd $pack
	    git branch
	    git checkout dev
	    git branch
	    cd ..
	done
    fi
    fi
# now update python path
echo "updating python path"
#export PYTHONPATH=${PWD}/sn_catalog_simulations:$PYTHONPATH
#export PYTHONPATH=${PWD}/sn_catalog_simulations/sn_simulation:$PYTHONPATH
#export PYTHONPATH=${PWD}/sn_catalog_simulations/sim_sncosmo:$PYTHONPATH
#export PYTHONPATH=${PWD}/sn_catalog_simulations/sim_snsim:$PYTHONPATH
#export PYTHONPATH=${PWD}/sn_catalog_simulations/sim_snana:$PYTHONPATH
#export PYTHONPATH=${PWD}/sn_catalog_simulations/sim_snfast:$PYTHONPATH

#export PYTHONPATH=${PWD}/sn_utils/utils:$PYTHONPATH
#export PYTHONPATH=${PWD}/sn_utils:$PYTHONPATH
export SN_UTILS_DIR=${PWD}/sn_utils
export SALT2_DIR=${PWD}/sn_utils/SALT2_Files
export PYTHONPATH=${thedir}:$PYTHONPATH
