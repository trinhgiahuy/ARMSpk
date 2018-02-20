#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile miniTri
if [ ! -f $ROOTDIR/MiniTri/miniTri/linearAlgebra/MPI/miniTri.exe ]; then
	cd $ROOTDIR/MiniTri/miniTri/linearAlgebra/MPI
	make
	cd $ROOTDIR/MiniTri/miniTri/linearAlgebra/openmp
	sed -i -e 's/= g++/= icpc/' Makefile
	sed -i -r '/Time to compute miniTri/ s#^(.*)$#//\1#' miniTri.cc
	make
	# get an valid input
	if [ ! -f $ROOTDIR/MiniTri/bcsstk30.mtx ]; then
		cd $ROOTDIR/MiniTri
		wget ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcsstruc5/bcsstk30.mtx.gz
		gunzip bcsstk30.mtx.gz
	fi
	cd $ROOTDIR
fi

