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
elif [[ "$1" = *"gnu"* ]]; then
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
	spack load openmpi@3.1.6%gcc@8.4.0
	export OMPI_CC=gcc
	export OMPI_CXX=g++
	export OMPI_F77=gfortran
	export OMPI_FC=gfortran
elif [[ "$1" = *"fuji"* ]]; then
	sleep 0
elif [[ "$1" = *"gem5"* ]]; then
	module load FujitsuCompiler/202007
else
	echo 'wrong compiler'
	exit 1
fi

BM="HPCG"
VERSION="5422fecd0a009a8731d0bd96b957d443297a53bc"
if [ ! -f $ROOTDIR/$BM/build/bin/xhpcg ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	if [ -z $1 ]; then
		git apply --check $ROOTDIR/patches/*1-${BM}*.patch
		if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	else
		git apply --check $ROOTDIR/patches/*2-${BM}*.patch
		if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*2-${BM}*.patch; fi
	fi
	mkdir -p build; cd build
	if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
		CONF="IMPI_IOMP_AVX2"
	elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
		CONF="IMPI_IOMP_KNL"
	elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
		CONF="IMPI_IOMP_KNL"
	fi
	if [ -z $1 ]; then
		../configure $CONF
		sed -i -e 's/mpiicpc/mpicxx/' -e 's/= -L${ADVISOR/= -static -static-intel -qopenmp-link=static -L${ADVISOR/' ./setup/Make.IMPI_IOMP_AVX2
	elif [[ "$1" = *"gnu"* ]]; then
		../configure MPI_GCC_OMP
		sed -i -e 's/-O3/-O3 -march=native -static/g' ./setup/Make.MPI_GCC_OMP
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		../configure MPI_GCC_OMP
		sed -i -e 's/^CXX .*=.*/CXX = mpiFCCpx/g' -e 's/-O3/-Nclang -Ofast -ffj-ocl -mllvm -polly -flto/g' ./setup/Make.MPI_GCC_OMP
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"gem5"* ]]; then
		../configure MPI_GCC_OMP
		sed -i -e 's/^CXX .*=.*/CXX = FCCpx/g' -e 's/-O3/-Nclang -Ofast -ffj-no-largepage -ffj-ocl -mllvm -polly -flto/g' -e 's/^HPCG_OPTS .*=.*/HPCG_OPTS = -DHPCG_NO_MPI/g' ./setup/Make.MPI_GCC_OMP
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	fi
	make
	cd $ROOTDIR
fi

