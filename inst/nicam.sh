#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile NICAM
if [ ! -f $ROOTDIR/NICAM/bin/nhm_driver ]; then
	cd $ROOTDIR/NICAM/src
	export NICAM_SYS=Linux64-intel-impi
	make
	cd '../test/case/jablonowski'
	make
	cd $ROOTDIR
fi

