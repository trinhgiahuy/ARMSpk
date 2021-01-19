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

BM="MVMC"
VERSION="7c58766b180ccb1035e4c220208b64ace3c49cf2"
if [ ! -f $ROOTDIR/$BM/src/vmc.out ] || [ "x`ls -s $ROOTDIR/$BM/src/vmc.out | awk '{print $1}'`" = "x0" ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	sed -i -e 's/blacs_intelmpi/blacs_openmpi/' -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' -e "s#L/usr/local/intel/composer_xe_2013/mkl#L$MKLROOT#g" ./Makefile_intel
	for x in `/bin/grep -r 'init_by_array(' . | cut -d':' -f1 | sort -u`; do sed -i -e 's/init_by_array(/xxxinit_by_array(/' $x; done
	make intel
	cd $ROOTDIR
fi

