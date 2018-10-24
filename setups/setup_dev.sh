#!/bin/bash

#setup lsst_sims

#source /global/common/software/lsst/cori-haswell-gcc/stack/setup_w_2018_13-sims_2_7_0.sh
source /cvmfs/sw.lsst.eu/linux-x86_64/lsst_sims/sims_2_8_0/loadLSST.bash
setup lsst_sims

export PYTHONPATH=${PWD}/SN_MAF/SN_Metrics:$PYTHONPATH
export PYTHONPATH=${PWD}/SN_MAF/SN_Stackers:$PYTHONPATH

thedir=${PWD}/lib/python3.6/site-packages/
if [ ! -d ${PWD}/SN_Catalog_Simulations ]
then
    echo "Seems that the simulation packages you need are not installed."
    echo "Will do it for you."
    where='https://github.com/lsstdesc'

    for pack in SN_Catalog_Simulations SN_Utils
    do
	echo 'Cloning' $pack
	git clone $where/$pack.git
	cd $pack
	git branch
	git checkout dev
	git branch
	cd ..
    done
   
    #checking whether hdf5 is accessible localy or not
    lib='h5py'
    echo $thedir
    if [ -d ${thedir}$lib ]
    then
	echo $lib 'already installed -> updating PYTHONPATH'
    else
	echo $lib 'not installed -> installing with pip'
	pip install --prefix=${PWD} ${lib}==2.7.1
    fi
    fi
# now update python path
echo "updating python path"
export PYTHONPATH=${PWD}/SN_Catalog_Simulations/SN_Simulation:$PYTHONPATH
export PYTHONPATH=${PWD}/SN_Catalog_Simulations/Sim_SNCosmo:$PYTHONPATH
export PYTHONPATH=${PWD}/SN_Catalog_Simulations/Sim_SNSim:$PYTHONPATH
export PYTHONPATH=${PWD}/SN_Catalog_Simulations/Sim_SNAna:$PYTHONPATH

export PYTHONPATH=${PWD}/SN_Utils/Utils:$PYTHONPATH
export SN_UTILS_DIR=${PWD}/SN_Utils
export SALT2_DIR=${PWD}/SN_Utils/SALT2_Files
export PYTHONPATH=${thedir}:$PYTHONPATH  
