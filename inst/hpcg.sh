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

source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load openmpi@3.1.6%intel@19.0.1.144
export OMPI_CC=$I_MPI_CC
export OMPI_CXX=$I_MPI_CXX
export OMPI_F77=$I_MPI_F77
export OMPI_FC=$I_MPI_F90

BM="HPCG"
VERSION="5422fecd0a009a8731d0bd96b957d443297a53bc"
if [ ! -f $ROOTDIR/$BM/build/bin/xhpcg ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	mkdir -p build; cd build
	if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
		CONF="IMPI_IOMP_AVX2"
	elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
		CONF="IMPI_IOMP_KNL"
	elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
		CONF="IMPI_IOMP_KNL"
	else
		echo "Unsupported host"
		exit
	fi
	../configure $CONF
	sed -i -e 's/mpiicpc/mpicxx/' -e 's/= -L${ADVISOR/= -static -static-intel -qopenmp-link=static -L${ADVISOR/' ./setup/Make.IMPI_IOMP_AVX2
	make
	cd $ROOTDIR
fi

