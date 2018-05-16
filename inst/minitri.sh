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

BM="MiniTri"
VERSION="9771c71f3d25023fc50bc6e84a905d6d50e81151"
if [ ! -f $ROOTDIR/$BM/miniTri/linearAlgebra/MPI/miniTri.exe ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/miniTri/linearAlgebra/MPI
	make
	cd $ROOTDIR/$BM/miniTri/linearAlgebra/openmp
	make
	# get an valid input
	if [ ! -f $ROOTDIR/$BM/bcsstk30.mtx ]; then
		cd $ROOTDIR/$BM/
		wget ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcsstruc5/bcsstk30.mtx.gz
		gunzip bcsstk30.mtx.gz
	fi
	cd $ROOTDIR
fi

