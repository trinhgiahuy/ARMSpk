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

BM="SWFFT"  # fortran version is 5-10% faster in my tests
VERSION="d0ef31454577740fbb87618cc35789b7ef838238"
if [ ! -f $ROOTDIR/$BM/build.xeon/TestDfft ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	if [ ! -f $ROOTDIR/fftw-3.3.4/bin/fftw-wisdom ]; then
		wget http://fftw.org/fftw-3.3.4.tar.gz
		tar xzf fftw-3.3.4.tar.gz
		cd ./fftw-3.3.4/
		./configure --prefix=`pwd`/../fftw-xmic --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=icc
		make -j CFLAGS="-O3 -ipo -xHost -xMIC-AVX512 -fp-model fast=2 -no-prec-div -qoverride-limits"
		make install
		make distclean
		./configure --prefix=`pwd`/../fftw-xeon --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=icc
		make -j CFLAGS="-O3 -ipo -xHost -xCORE-AVX2 -fp-model fast=2 -no-prec-div -qoverride-limits"
		make install
		cd $ROOTDIR/$BM/
	fi
	oldPATH=$PATH
	export PATH=$oldPATH:`pwd`/fftw-xeon/bin
	make -f GNUmakefile.openmp
	mv build.openmp build.xeon
	make -f GNUmakefile.openmp clean
	export PATH=$oldPATH:`pwd`/fftw-xmic/bin
	make -f GNUmakefile.openmp
	mv build.openmp build.xmic
	make -f GNUmakefile.openmp clean
	export PATH=$oldPATH
	cd $ROOTDIR
fi

