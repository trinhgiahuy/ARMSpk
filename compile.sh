#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile AMG -> comes w/ 2 problems
if [ ! -f ./AMG/test/amg ]; then
	cd ./AMG/
	make
	make clean
	cd $ROOTDIR
fi

# compile CANDLE -> comes w/ 7 problems
if [ ! -f $HOME/anaconda2/bin/anaconda ]; then
	cd ./CANDLE/
	curl -o Anaconda2-4.3.1-Linux-x86_64.sh https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh
	chmod u+x ./Anaconda2-4.3.1-Linux-x86_64.sh
	./Anaconda2-4.3.1-Linux-x86_64.sh -b
	export PATH=$HOME/anaconda2/bin:$PATH
	conda install -y -c conda-forge tensorflow
	conda install -y -c anaconda hdf5=1.8.17
	conda install -y -c anaconda theano
	conda install -y -c conda-forge keras=2
	conda install -y -c conda-forge opencv
	conda install -y -c conda-forge tqdm
	conda update -y -c conda-forge numpy
	cd $ROOTDIR
fi

# compile CoMD
if [ ! -f ./CoMD/bin/CoMD-openmp-mpi ]; then
	cd ./CoMD/src-openmp/
	cp Makefile.vanilla Makefile
	make
	cd $ROOTDIR
fi

# compile Laghos
if [ ! -f ./Laghos/laghos ]; then
	cd ./Laghos/
	if [ ! -f hypre-2.10.0b/src/hypre/lib/libHYPRE.a ]; then
		wget https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.10.0b.tar.gz
		tar xzf hypre-2.10.0b.tar.gz
		cd ./hypre-2.10.0b/src
		./configure --disable-fortran
		make -j
		cd $ROOTDIR/Laghos/
	fi
	if [ ! -f metis-4.0.3/graphchk ]; then
		wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
		tar xzf metis-4.0.3.tar.gz
		cd ./metis-4.0.3/
		make
		cd $ROOTDIR/Laghos/
	fi
	if [ ! -f ../dep/mfem/libmfem.a ]; then
		cd ../dep/mfem/
		git checkout laghos-v1.0
		sed -i -e 's#@MFEM_DIR@/../hypre#@MFEM_DIR@/../../Laghos/hypre#' config/defaults.mk
		sed -i -e 's#@MFEM_DIR@/../metis-4.0$#@MFEM_DIR@/../../Laghos/metis-4.0.3#' config/defaults.mk
		make parallel -j
		cd $ROOTDIR/Laghos/
	fi
	sed -i -e 's#MFEM_DIR = ../mfem$#MFEM_DIR = ../dep/mfem#' makefile
	sed -i -e 's#LAGHOS_LIBS = $(MFEM_LIBS)$#LAGHOS_LIBS = $(MFEM_LIBS) -lirc -lsvml#' makefile
	make
	cd $ROOTDIR
fi

# compile MACSio
if [ ! -f ./MACSio/macsio/macsio ]; then
	cd ./MACSio/
	if [ ! -f ../dep/json-cwx/lib/libjson-cwx.a ]; then
		cd ../dep/json-cwx/json-cwx
		./autogen.sh
		./configure --prefix=`pwd`/../
		make install
		cd $ROOTDIR/MACSio/
	fi
	if [ ! -f ../dep/silo-4.10.2/bin/silofile ]; then
		cd ../dep/
		wget https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz
		tar xzf silo-4.10.2.tar.gz
		cd silo-4.10.2/
		./configure --prefix=`pwd`
		make install
		cd $ROOTDIR/MACSio/
	fi
	mkdir -p build; cd build
	cmake -DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
	make
	make install
	cd $ROOTDIR
fi

# compile miniAMR
if [ ! -f ./MiniAMR/ref/ma.x ]; then
	cd ./MiniAMR/ref/
	sed -i -e 's#= cc#= mpicc#' Makefile
	make
	cd $ROOTDIR
fi

# compile miniFE
if [ ! -f ./MiniFE/mkl/src/miniFE.x ]; then
	cd ./MiniFE/mkl/src
	make
	cd $ROOTDIR
	cd ./MiniFE/openmp-opt/src
	make
	cd $ROOTDIR
	cd ./MiniFE/openmp-opt-knl/src
	make
	cd $ROOTDIR
fi

