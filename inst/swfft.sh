#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile SWFFT
if [ ! -f $ROOTDIR/SWFFT/build.xeon/TestDfft ]; then
	cd $ROOTDIR/SWFFT
	if [ ! -f $ROOTDIR/fftw-3.3.4/bin/fftw-wisdom ]; then
		wget http://fftw.org/fftw-3.3.4.tar.gz
		tar xzf fftw-3.3.4.tar.gz
		cd ./fftw-3.3.4/
		./configure --prefix=`pwd`/../fftw-xmic --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=icc
		make -j CFLAGS="-O3 -xMIC-AVX512 -fp-model fast=2 -no-prec-div -qoverride-limits"
		make install
		make distclean
		./configure --prefix=`pwd`/../fftw-xeon --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=icc
		make -j CFLAGS="-O3 -xCORE-AVX2 -fp-model fast=2 -no-prec-div -qoverride-limits"
		make install
		cd $ROOTDIR/SWFFT
	fi
	# fortran version is 5-10% faster in my tests
	sed -i -e 's/^default: nativec/default: fortran/' GNUmakefile
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

