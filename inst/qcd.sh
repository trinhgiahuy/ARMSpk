#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source `cat $ROOTDIR/conf/intel.cfg` intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort
alias ar=`which xiar`
alias ld=`which xild`

BM="QCD"
VERSION="07277047b170529caa5fcd164afd814e70286ce4"
if [ ! -f $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class1 ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	make MAKE_INC=make.ifort.inc CLASS=1 PX=1 PY=1 PZ=1
	make MAKE_INC=make.ifort.inc CLASS=2 PX=1 PY=1 PZ=1
	make MAKE_INC=make.ifort.inc CLASS=3 PX=1 PY=1 PZ=1
	cd $ROOTDIR
fi

