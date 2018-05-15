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

BM="NICAM"
VERSION="3f758000ffce6ee95a27fb6099f654ecdc5e3add"
if [ ! -f $ROOTDIR/$BM/bin/nhm_driver ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	export NICAM_SYS=Linux64-intel-impi
	make ENABLE_OPENMP=1
	cd '../test/case/jablonowski'
	make ENABLE_OPENMP=1
	cd $ROOTDIR
fi

