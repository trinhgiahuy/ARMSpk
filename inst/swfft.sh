#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
if [ -z $1 ]; then
	source $ROOTDIR/conf/intel.cfg
	source $INTEL_PACKAGE intel64 > /dev/null 2>&1
	export I_MPI_CC=icc
	export I_MPI_CXX=icpc
	export I_MPI_F77=ifort
	export I_MPI_F90=ifort
	alias ar=`which xiar`
	alias ld=`which xild`
	export ADVISOR_2018_DIR=${ADVISOR_2019_DIR}

	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load openmpi@3.1.6%intel@19.0.1.144
	export OMPI_CC=$I_MPI_CC
	export OMPI_CXX=$I_MPI_CXX
	export OMPI_F77=$I_MPI_F77
	export OMPI_FC=$I_MPI_F90
else
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
	spack load openmpi@3.1.6%gcc@8.4.0
	export OMPI_CC=gcc
	export OMPI_CXX=g++
	export OMPI_F77=gfortran
	export OMPI_FC=gfortran
fi

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
		if [ -z $1 ]; then
			./configure --prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=icc
			if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
				make -j CFLAGS="-O3 -ipo -xHost -xCORE-AVX2 -fp-model fast=2 -no-prec-div -qoverride-limits"
			else
				make -j CFLAGS="-O3 -ipo -xHost -xCORE-AVX512 -fp-model fast=2 -no-prec-div -qoverride-limits"
			fi
		else
			./configure --prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=gcc
			make -j CFLAGS="-O3 -march=native"
		fi
		make install
		cd $ROOTDIR/$BM/
	fi
	export oldPATH=$PATH
	export PATH=$ROOTDIR/$BM/fftw/bin:$oldPATH
	if [ -z $1 ]; then
		sed -i -e 's/-L${ADVISOR/-lmpi_cxx -static -static-intel -qopenmp-link=static -L${ADVISOR/' ./GNUmakefile
	else
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -lmpi_cxx -static#g' ./GNUmakefile
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -lmpi_cxx -static#g' ./GNUmakefile.openmp
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	fi
	make -f GNUmakefile.openmp
	export PATH=$oldPATH
	unset oldPATH
	cd $ROOTDIR
fi

