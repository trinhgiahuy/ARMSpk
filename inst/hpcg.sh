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

BM="HPCG"
VERSION="5422fecd0a009a8731d0bd96b957d443297a53bc"
if [ ! -f $ROOTDIR/$BM/bin/build/xhpcg ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	mkdir -p build; cd build
	if [[ $HOSTNAME = *"kiev"* ]]; then
		CONF="IMPI_IOMP_AVX2"
	elif [[ $HOSTNAME = *"lyon"* ]]; then
		CONF="IMPI_IOMP_KNL"
	elif [[ $HOSTNAME = *"mill"* ]]; then
		CONF="IMPI_IOMP_KNL"
	else
		echo "Unsupported host"
		exit
	fi
	../configure $CONF
	make
	cd $ROOTDIR
fi

