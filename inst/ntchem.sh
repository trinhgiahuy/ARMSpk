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

BM="NTChem"
VERSION="fcafcc4fec195d8a81c19affd1a3b83f7bab4285"
if [ ! -f $ROOTDIR/$BM/bin/rimp2.exe ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./src/mp2/GNUmakefile
	TOP_DIR=`pwd`
	TYPE=intel
	./config_mine
	mkdir -p $ROOTDIR/$BM/bin
	make CC=mpicc CXX=mpicxx F77C=mpif77 F90C=mpif90
	unset TOP_DIR
	unset TYPE
	cd $ROOTDIR
fi

