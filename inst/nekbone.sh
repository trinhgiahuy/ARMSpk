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

BM="Nekbone"
VERSION="d681db43f45f5437e5258b3663b5a92c078cfb57"
if [ ! -f $ROOTDIR/$BM/test/nek_mgrid/nekbone ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/test/nek_mgrid
	sed -i -e 's/lp = 10)/lp = 576)/' -e 's/lelt=100/lelt=1024/' SIZE
	sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' $ROOTDIR/$BM/src/makefile.template
	./makenek NotUsedCasename $ROOTDIR/$BM/src
	cd $ROOTDIR
fi

