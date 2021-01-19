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

BM="HPL"
if [ ! -f $ROOTDIR/$BM/bin/Linux_Intel64/xhpl ]; then
	mkdir -p $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/
	wget http://www.netlib.org/benchmark/hpl/hpl-2.2.tar.gz
	tar xzf ./hpl-2.2.tar.gz -C $ROOTDIR/$BM --strip-components 1
	patch -p1 < $ROOTDIR/patches/*1-${BM}*.patch
	sed -i -e 's/mpiicc/mpicc/' -e 's/ -L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" ./Make.Linux_Intel64
	make arch=Linux_Intel64
	cd $ROOTDIR
fi

