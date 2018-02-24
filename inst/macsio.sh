#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort
alias ar=`which xiar`
alias ld=`which xild`

BM="MACSio"
if [ ! -f $ROOTDIR/$BM/macsio/macsio ]; then
	cd $ROOTDIR/$BM/
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am < $ROOTDIR/patches/*1-${BM}*.patch; fi
	if [ ! -f $ROOTDIR/dep/json-cwx/lib/libjson-cwx.a ]; then
		cd $ROOTDIR/dep/json-cwx/json-cwx
		./autogen.sh
		./configure --prefix=`pwd`/../ CC=icc CFLAGS="-O2 -ipo -xHost"
		make
		make install
		cd $ROOTDIR/$BM/
	fi
	if [ ! -f $ROOTDIR/dep/silo-4.10.2/bin/silofile ]; then
		cd $ROOTDIR/dep/
		wget https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz
		tar xzf silo-4.10.2.tar.gz
		cd silo-4.10.2/
		./configure --prefix=`pwd` CC=icc CFLAGS="-O2 -ipo -xHost" CXX=icpc CXXFLAGS="-O2 -ipo -xHost" FC=ifort FCFLAGS="-O2 -ipo -xHost" F77=ifort FFLAGS="-O2 -ipo -xHost"
		make install
		cd $ROOTDIR/$BM/
	fi
	rm -rf build; mkdir -p build; cd build
	cmake -DCMAKE_C_COMPILER=`which mpicc` -DCMAKE_C_FLAGS="-O3 -ipo -xHost" -DCMAKE_CXX_COMPILER=`which mpicxx` -DCMAKE_CXX_FLAGS="-O3 -ipo -xHost" -DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
	make
	make install
	cd $ROOTDIR
fi

