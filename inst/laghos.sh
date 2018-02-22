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

BM="Laghos"
if [ ! -f $ROOTDIR/$BM/laghos ]; then
	cd $ROOTDIR/$BM/
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am < $ROOTDIR/patches/*1-${BM}*.patch; fi
	if [ ! -f ./hypre-2.10.0b/src/hypre/lib/libHYPRE.a ]; then
		wget https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.10.0b.tar.gz
		tar xzf hypre-2.10.0b.tar.gz
		cd ./hypre-2.10.0b/src
		./configure --disable-fortran -with-openmp
		sed -i -e 's/ -openmp/ -fopenmp/g' ./config/Makefile.config
		make -j
		cd $ROOTDIR/$BM/
	fi
	if [ ! -f ./metis-4.0.3/graphchk ]; then
		wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
		tar xzf metis-4.0.3.tar.gz
		cd ./metis-4.0.3/
		make
		cd $ROOTDIR/$BM/
	fi
	if [ ! -f $ROOTDIR/dep/mfem/libmfem.a ]; then
		cd $ROOTDIR/dep/mfem/
		git checkout laghos-v1.0
		sed -i -e 's#^OPTIM_FLAGS = -O3#OPTIM_FLAGS = -O3 -fopenmp#' config/defaults.mk
		sed -i -e "s#@MFEM_DIR@/../hypre#@MFEM_DIR@/../../$BM/hypre#" config/defaults.mk
		sed -i -e "s#@MFEM_DIR@/../metis-4.0#@MFEM_DIR@/../../$BM/metis-4.0.3#" config/defaults.mk
		make parallel -j
		cd $ROOTDIR/$BM/
	fi
	sed -i -e 's#^CXXFLAGS = \$(MFEM_CXXFLAGS)$#CXXFLAGS = \$(MFEM_CXXFLAGS) -fopenmp#' makefile
	sed -i -e 's#MFEM_DIR = ../mfem$#MFEM_DIR = ../dep/mfem#' makefile
	sed -i -e 's#LAGHOS_LIBS = $(MFEM_LIBS)$#LAGHOS_LIBS = $(MFEM_LIBS) -lirc -lsvml#' makefile
	make
	cd $ROOTDIR
fi

