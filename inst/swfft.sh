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

BM="SWFFT"  # fortran version is 5-10% faster in my tests
VERSION="d0ef31454577740fbb87618cc35789b7ef838238"
if [ ! -f $ROOTDIR/$BM/build.openmp/TestFDfft ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	if [ ! -f $ROOTDIR/$BM/fftw/bin/fftw-wisdom ]; then
		wget http://fftw.org/fftw-3.3.4.tar.gz
		tar xzf fftw-3.3.4.tar.gz
		cd ./fftw-3.3.4/
		./configure --prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=icc
		if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
			make -j CFLAGS="-O3 -ipo -xHost -xCORE-AVX2 -fp-model fast=2 -no-prec-div -qoverride-limits"
		else
			make -j CFLAGS="-O3 -ipo -xHost -xCORE-AVX512 -fp-model fast=2 -no-prec-div -qoverride-limits"
		fi
		make install
		cd $ROOTDIR/$BM/
	fi
	export oldPATH=$PATH
	export PATH=$ROOTDIR/$BM/fftw/bin:$oldPATH
	make -f GNUmakefile.openmp
	export PATH=$oldPATH
	unset oldPATH
	cd $ROOTDIR
fi

