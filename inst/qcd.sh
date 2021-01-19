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

source $ROOTDIR/conf/qcd.sh

BM="QCD"
VERSION="07277047b170529caa5fcd164afd814e70286ce4"
if [ ! -f $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2_111 ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	sed -i -e 's/-shared-intel -mcmodel=medium/-static -static-intel -qopenmp-link=static/' -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./make.ifort.inc
	for TEST in $TESTCONF; do
		PX="`echo $TEST | cut -d '|' -f3`"
		PY="`echo $TEST | cut -d '|' -f4`"
		PZ="`echo $TEST | cut -d '|' -f5`"
		if [ ! -f $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2_${PX}${PY}${PZ} ]; then
			make clean >> /dev/null 2>&1
			make MAKE_INC=make.ifort.inc CLASS=2 PX=${PX} PY=${PY} PZ=${PZ}
			mv $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2 $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2_${PX}${PY}${PZ}
		fi
	done
	cd $ROOTDIR
fi

