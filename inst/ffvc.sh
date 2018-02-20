#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile FFVC
if [ ! -f $ROOTDIR/FFVC/bin/ffvc_mini ]; then
	cd $ROOTDIR/FFVC/src
	sed -i -e 's/-openmp/-fopenmp/' make_setting.intel
	rm make_setting; ln -s make_setting.intel make_setting
	make
	cd $ROOTDIR
fi

