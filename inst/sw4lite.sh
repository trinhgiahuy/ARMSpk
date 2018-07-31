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

BM="SW4lite"
VERSION="5ab8063ecdc94bdb59a5e65396c85bd54f9e0916"
if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	HHOST="${XEONHOST}"
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	HHOST="${IKNLHOST}"
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	HHOST="${IKNMHOST}"
else
	echo "Unsupported host"
	exit
fi
if [ ! -f $ROOTDIR/$BM/optimize_mp_${HHOST}/sw4lite ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	sed -i -e "s/^HOSTNAME := /HOSTNAME := ${HHOST} #/g" Makefile
	sed -i -e "s/quadknl/${HHOST}/g" Makefile
	sed -i -e "s#/opt/intel/compilers_and_libraries_2017/linux#`dirname $MKLROOT`#g"  Makefile
	if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
		sed -i -e "s/-xmic-avx512/#NOKNL-xmic-avx512/g" Makefile
	fi
	make
	cd $ROOTDIR
fi

