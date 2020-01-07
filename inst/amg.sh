#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort
alias ar=`which xiar`
alias ld=`which xild`
export ADVISOR_2018_DIR=${ADVISOR_2019_DIR}

BM="AMG"
VERSION="295de9693eaabf6f7330ac3a35fd9bd4ad030522"
if [ ! -f $ROOTDIR/$BM/test/amg ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	# avx512 on KNL/KNM causes errors in AMG exec
	if [[ $HOSTNAME = *"${IKNLHOST}"* ]] || [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
		sed -i -e 's/xHost/xCORE-AVX2/g' ./Makefile.include
	fi
	make
	cd $ROOTDIR
fi

