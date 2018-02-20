#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile mVMC
if [ ! -f $ROOTDIR/MVMC/src/vmc.out ]; then
	cd $ROOTDIR/MVMC/src
	sed -i -e 's/-openmp/-fopenmp/g' -e 's/-opt-prefetch=3/-qopt-prefetch=3/g' -e "s#L/usr/local/intel/composer_xe_2013/mkl#L$MKLROOT#g"  Makefile_intel
	make intel
	cd $ROOTDIR
fi

