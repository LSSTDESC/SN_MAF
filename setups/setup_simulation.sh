#!/bin/bash

location=$1
script_loc=$2

source ${PWD}/sn_maf/setups/setup_metric.sh ${location} ${script_loc}


thedir=${PWD}/lib/python3.6/site-packages/

pack="sn_catalog_simulations"
if [ ! -d  $pack ]; then
	echo "Seems that the $pack package you need is not installed."
	echo "Will do it for you..."
	where='https://github.com/lsstdesc'
	echo "...taking packages from "${where}
	
	echo 'Cloning' $pack
	git clone $where/$pack.git
	echo "Moving to the dev_stable branch"
	cd $pack
	git branch
	git checkout dev_stable
	git branch
	cd ..
fi

# now update python path
echo "updating python path"
export SN_UTILS_DIR=${PWD}/sn_utils
export SALT2_DIR=${PWD}/sn_utils/SALT2_Files
export PYTHONPATH=${thedir}:$PYTHONPATH
