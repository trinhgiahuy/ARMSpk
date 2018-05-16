#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source `cat $ROOTDIR/conf/intel.cfg` intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort
alias ar=`which xiar`
alias ld=`which xild`

BM="NTChem"
VERSION="fcafcc4fec195d8a81c19affd1a3b83f7bab4285"
if [ ! -f $ROOTDIR/$BM/bin/rimp2.exe ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	TOP_DIR=`pwd`
	TYPE=intel
	./config_mine
	mkdir -p $ROOTDIR/$BM/bin
	make CC=mpicc CXX=mpicxx F77C=mpif77 F90C=mpif90
	cd $ROOTDIR
fi

