#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort
alias ar=`which xiar`
alias ld=`which xild`

BM="SW4lite"
if [ ! -f $ROOTDIR/$BM/optimize_mp_kiev/sw4lite ]; then
	cd $ROOTDIR/$BM/
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am < $ROOTDIR/patches/*1-${BM}*.patch; fi
	sed -i -e "s/^HOSTNAME := /HOSTNAME := kiev #/g" Makefile
	sed -i -e "s/quadknl/kiev/g" -e "s#/opt/intel/compilers_and_libraries_2017/linux#`dirname $MKLROOT`#g"  Makefile
	sed -i -e "s/-xmic-avx512/#NOKNL-xmic-avx512/g" Makefile
	make
	sed -i -e "s/kiev/mill/g" Makefile
	sed -i -e "s/#NOKNL-xmic-avx512/-xmic-avx512/g" Makefile
	make
	cd $ROOTDIR
fi
