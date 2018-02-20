#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile XSBench
if [ ! -f $ROOTDIR/XSBench/src/XSBench ]; then
	cd $ROOTDIR/XSBench/src
	sed -i -e 's/^COMPILER.*= gnu/COMPILER = intel/' -e 's/^MPI.* = no/MPI = yes/' -e 's/-openmp/-fopenmp/' Makefile
	make
	cd $ROOTDIR
fi

