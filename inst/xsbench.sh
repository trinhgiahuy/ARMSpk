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

BM="XSBench"
if [ ! -f $ROOTDIR/$BM/src/XSBench ]; then
	cd $ROOTDIR/$BM/
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	sed -i -e 's/^COMPILER.*= gnu/COMPILER = intel/' -e 's/^MPI.* = no/MPI = yes/' -e 's/-openmp/-fopenmp/' Makefile
	make
	cd $ROOTDIR
fi

