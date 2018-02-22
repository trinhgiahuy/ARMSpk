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

BM="NTChem"
if [ ! -f $ROOTDIR/$BM/bin/rimp2.exe ]; then
	cd $ROOTDIR/$BM/
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am < $ROOTDIR/patches/*1-${BM}*.patch; fi
	TOP_DIR=`pwd`
	TYPE=intel
	cp platforms/config_mine.${TYPE} ./config_mine
	sed -i -e 's/-openmp/-fopenmp/g' ./config/linux64_mpif90_omp_intel_proto.makeconfig.in
	./config_mine
	mkdir -p $ROOTDIR/$BM/bin
	make CC=mpicc CXX=mpicxx F77C=mpif77 F90C=mpif90
	cd $ROOTDIR
fi

