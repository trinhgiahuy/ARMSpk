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

BM="AMG"
if [ ! -f $ROOTDIR/$BM/test/amg ]; then
	cd $ROOTDIR/$BM/
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am < $ROOTDIR/patches/*1-${BM}*.patch; fi
	# avx512 on KNL/KNM causes errors in AMG exec
	if [[ $HOSTNAME = *"lyon"* ]] || [[ $HOSTNAME = *"mill"* ]]; then
		sed -i -e 's/xHost/xCORE-AVX2/g' ./Makefile.include
	fi
	make
	cd $ROOTDIR
fi

