#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile CoMD
if [ ! -f $ROOTDIR/CoMD/bin/CoMD-openmp-mpi ]; then
	cd $ROOTDIR/CoMD/src-openmp/
	cp Makefile.vanilla Makefile
	make
	cd $ROOTDIR
fi

