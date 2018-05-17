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

if [ ! -f $ROOTDIR/dep/modylas-mini-1.0.0.tar.gz ]; then
	echo "ERR: Cannot find modylas-mini-1.0.0.tar.gz"
	echo "Please download from: http://hpci-aplfs.aics.riken.jp/fiber/modylas-mini.html and place modylas-mini-1.0.0.tar.gz in ./dep subfolder"
	exit
fi

BM="MODYLAS" # (req. license agreement on website)
if [ ! -f $ROOTDIR/$BM/src/modylas_mini ]; then
	mkdir -p $ROOTDIR/$BM/
	tar xzf $ROOTDIR/dep/modylas-mini-1.0.0.tar.gz -C $ROOTDIR/$BM --strip-components 1
	cd $ROOTDIR/$BM/
	patch -p1 < $ROOTDIR/patches/*1-${BM}*.patch
	cd $ROOTDIR/$BM/src
	rm make_setting; ln -s make_setting.intel make_setting
	make
	cd $ROOTDIR
fi

