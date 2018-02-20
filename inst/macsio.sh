#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile MACSio
if [ ! -f $ROOTDIR/MACSio/macsio/macsio ]; then
	cd $ROOTDIR/MACSio/
	if [ ! -f $ROOTDIR/dep/json-cwx/lib/libjson-cwx.a ]; then
		cd $ROOTDIR/dep/json-cwx/json-cwx
		./autogen.sh
		./configure --prefix=`pwd`/../
		make install
		cd $ROOTDIR/MACSio/
	fi
	if [ ! -f $ROOTDIR/dep/silo-4.10.2/bin/silofile ]; then
		cd $ROOTDIR/dep/
		wget https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz
		tar xzf silo-4.10.2.tar.gz
		cd silo-4.10.2/
		./configure --prefix=`pwd`
		make install
		cd $ROOTDIR/MACSio/
	fi
	mkdir -p build; cd build
	cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
	make
	make install
	cd $ROOTDIR
fi

