#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
if [ -z $1 ]; then
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
elif [[ "$1" = *"gnu"* ]]; then
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
	spack load openmpi@3.1.6%gcc@8.4.0
	export OMPI_CC=gcc
	export OMPI_CXX=g++
	export OMPI_F77=gfortran
	export OMPI_FC=gfortran
elif [[ "$1" = *"fuji"* ]]; then
	module load FujitsuCompiler/202007
fi

BM="Laghos"
VERSION="9a074521257434e0b9acff9e59ff10e3e881bc32"
if [ ! -f $ROOTDIR/$BM/laghos ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	if [ ! -f ./hypre-2.10.0b/src/hypre/lib/libHYPRE.a ]; then
		wget https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.10.0b.tar.gz
		tar xzf hypre-2.10.0b.tar.gz
		cd ./hypre-2.10.0b/src
		if [ -z $1 ]; then
			./configure --disable-fortran -with-openmp CC=mpicc CFLAGS="-O3 -ipo -xHost" CXX=mpicxx CXXFLAGS="-O3 -ipo -xHost" F77=mpif77 FFLAGS="-O3 -ipo -xHost"
		elif [[ "$1" = *"gnu"* ]]; then
			./configure --disable-fortran -with-openmp CC=mpicc CFLAGS="-O3 -march=native" CXX=mpicxx CXXFLAGS="-O3 -march=native" F77=mpif77 FFLAGS="-O3 -march=native"
		elif [[ "$1" = *"fuji"* ]]; then
			./configure --host=aarch64-unknown-linux-gnu --build=x84_64-unknown-linux-gnu --disable-fortran -with-openmp --without-MPI CC=fccpx CFLAGS="-O3" CXX=FCCpx CXXFLAGS="-O3" F77=frtpx FFLAGS="-O3"
		fi
		sed -i -e 's/ -openmp/ -fopenmp/g' ./config/Makefile.config
		make -j
		cd $ROOTDIR/$BM/
	fi
	if [ ! -f ./metis-4.0.3/graphchk ]; then
		wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
		tar xzf metis-4.0.3.tar.gz
		cd ./metis-4.0.3/
		if [ -z $1 ]; then
			sed -i -e 's/CC = cc/CC = icc/g' -e 's/OPTFLAGS = -O2\s*$/OPTFLAGS = -O2 -ipo -xHost/g' ./Makefile.in
		elif [[ "$1" = *"gnu"* ]]; then
			sed -i -e 's/CC = cc/CC = gcc/g' -e 's/OPTFLAGS = -O2\s*$/OPTFLAGS = -O2 -march=native/g' ./Makefile.in
		elif [[ "$1" = *"fuji"* ]]; then
			sed -i -e 's/CC = cc/CC = fccpx/g' ./Makefile.in
		fi
		make
		cd $ROOTDIR/$BM/
	fi
	if [ ! -f $ROOTDIR/dep/mfem/libmfem.a ]; then
		cd $ROOTDIR/dep/mfem/
		git checkout laghos-v1.0
		git apply --check $ROOTDIR/patches/*1-mfem*.patch
		if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-mfem*.patch; fi
		if [[ "$1" = *"gnu"* ]]; then
			sed -i -e 's/icpc/g++/g' -e 's/-ipo -xHost/-march=native/g' ./config/defaults.mk
		elif [[ "$1" = *"fuji"* ]]; then
			sed -i -e 's/icpc/FCCpx/g' -e 's/-ipo -xHost//g' ./config/defaults.mk
		fi
		if [[ "$1" != *"fuji"* ]]; then
			make config MFEM_USE_MPI=YES MPICXX=mpicxx MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
		else
			make config CMAKE_CXX_COMPILER=FCCpx MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
		fi
		cd $ROOTDIR/$BM/
	fi
	if [ -z $1 ]; then
		sed -i -e 's/= -L${ADVISOR/= -static -static-intel -qopenmp-link=static -L${ADVISOR/' ./makefile
		make
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -static#g' -e 's/-ipo -xHost/-march=native/g' -e 's/ -lirc -lsvml//g' ./makefile
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
		make
	elif [[ "$1" = *"fuji"* ]]; then
		cd serial/
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -Bstatic#g' -e 's/$(LAGHOS_LIBS) $(LDFLAGS)/$(LDFLAGS) $(LAGHOS_LIBS)/g' -e 's/-ipo -xHost//g' -e 's/ -lirc -lsvml//g' ./makefile
		sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' ./laghos_solver_s.hpp
		make
		cd -
		cp serial/laghos .
	fi
	cd $ROOTDIR
fi

