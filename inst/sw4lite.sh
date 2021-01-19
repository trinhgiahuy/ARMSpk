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
	sed -i -e "s/^HOSTNAME := /HOSTNAME := ${HHOST} #/g" ./Makefile
	sed -i -e "s/quadknl/${HHOST}/g" ./Makefile
	sed -i -e "s#/opt/intel/compilers_and_libraries_2017/linux#`dirname $MKLROOT`#g"  ./Makefile
	if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
		sed -i -e "s/-xmic-avx512/#NOKNL-xmic-avx512/g" ./Makefile
	fi
	sed -i -e 's/mpifort/mpif90/' -e 's/mpiifort/mpif90/' -e 's/mpiicpc/mpicxx/' -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' -e "s#MKL_PATH = .*#MKL_PATH = $MKLROOT/lib/intel64#" -e 's#-lmkl_intel_lp64 -lmkl_core -lmkl_sequential#-Wl,--start-group ${MKL_PATH}/libmkl_intel_lp64.a ${MKL_PATH}/libmkl_sequential.a ${MKL_PATH}/libmkl_core.a -Wl,--end-group#' -e 's/-lintlc/-lirc/' ./Makefile
	make
	cd $ROOTDIR
fi

