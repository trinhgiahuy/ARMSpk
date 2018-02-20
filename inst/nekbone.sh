#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile Nekbone
if [ ! -f $ROOTDIR/Nekbone/test/nek_mgrid/nekbone ]; then
	cd $ROOTDIR/Nekbone/test/nek_mgrid
	sed -i -e 's/lp = 10)/lp = 576)/' -e 's/lelt=100/lelt=1024/' SIZE
	./makenek NotUsedCasename $ROOTDIR/Nekbone/src
	cd $ROOTDIR
fi

