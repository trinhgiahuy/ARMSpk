#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile CCS QCD
if [ ! -f $ROOTDIR/QCD/src/ccs_qcd_solver_bench_class1 ]; then
	cd $ROOTDIR/QCD/src
	sed -i -e 's/-openmp/-fopenmp/' make.ifort.inc
	make MAKE_INC=make.ifort.inc CLASS=1 PX=1 PY=1 PZ=1
	make MAKE_INC=make.ifort.inc CLASS=2 PX=1 PY=1 PZ=1
	make MAKE_INC=make.ifort.inc CLASS=3 PX=1 PY=1 PZ=1
	cd $ROOTDIR
fi

