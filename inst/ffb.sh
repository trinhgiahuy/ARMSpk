#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile FFB
if [ ! -f $ROOTDIR/FFB/bin/les3x.mpi ]; then
	cd $ROOTDIR/FFB/src
	if [ ! -f ./metis-5.1.0/bin/graphchk ]; then
		wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
		tar xzf metis-5.1.0.tar.gz
		cd ./metis-5.1.0/
		make config cc=icc prefix=`pwd`
		make install
		cd $ROOTDIR/FFB/src
	fi
	if [ ! -f ./REVOCAP_Refiner-1.1.01/lib/x86_64-linux-intel/libRcapRefiner.a ]; then
		tar xzf $ROOTDIR/dep/REVOCAP_Refiner-1.1.01.tgz
		cd ./REVOCAP_Refiner-1.1.01
		rm ./MakefileConfig.in; ln -s ./MakefileConfig.LinuxIntelCompiler ./MakefileConfig.in
		sed -i '/#include <sstream>/a #include <stdlib.h>' ./RevocapIO/kmbHecmwIO_V3.cpp
		make
		cd $ROOTDIR/FFB/src
	fi
	cp ./make_setting.intel ./make_setting
	sed -i -e 's/^DEFINE += -DNO_METIS/#DEFINE += -DNO_METIS/g' -e "s#\$(HOME)/opt_intel/metis5#$ROOTDIR/FFB/src/metis-5.1.0#g" ./make_setting
	sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt_intel/REVOCAP_Refiner#$ROOTDIR/FFB/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/x86_64-linux-intel #" ./make_setting
	make
	cd $ROOTDIR
fi

