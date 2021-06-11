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
elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
	echo "ERR: does not compile on login node; please use compute node"; exit 1
elif lscpu | grep 'sve' >/dev/null 2>&1 && [[ "$1" = *"fuji"* ]]; then
	sleep 0
elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"gem5"* ]]; then
	echo "ERR: does not compile on login node; please use compute node"; exit 1
elif lscpu | grep 'sve' >/dev/null 2>&1 && [[ "$1" = *"gem5"* ]]; then
	module load FujitsuCompiler/202007
else
	echo 'wrong compiler'
	exit 1
fi

BM="Laghos"
VERSION="9a074521257434e0b9acff9e59ff10e3e881bc32"
if [ ! -f $ROOTDIR/$BM/laghos ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	if [ ! -f ./hypre-2.10.0b/src/hypre/lib/libHYPRE.a ]; then
		if ! [ -f ./hypre-2.10.0b.tar.gz ]; then wget https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.10.0b.tar.gz; fi
		tar xzf hypre-2.10.0b.tar.gz
		cd ./hypre-2.10.0b/src
		if [ -z $1 ]; then
			./configure --disable-fortran -with-openmp CC=mpicc CFLAGS="-O3 -ipo -xHost" CXX=mpicxx CXXFLAGS="-O3 -ipo -xHost" F77=mpif77 FFLAGS="-O3 -ipo -xHost"
		elif [[ "$1" = *"gnu"* ]]; then
			./configure --disable-fortran -with-openmp CC=mpicc CFLAGS="-O3 -march=native" CXX=mpicxx CXXFLAGS="-O3 -march=native" F77=mpif77 FFLAGS="-O3 -march=native"
		elif [[ "$1" = *"fuji"* ]]; then
			cd config/; rm -f config.guess config.sub; wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.guess'; wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.sub'; cd -
			./configure --disable-fortran -with-openmp CC=mpifcc CFLAGS="-Nclang -Ofast -ffj-ocl -mllvm -polly -flto" CXX=mpiFCC CXXFLAGS="-Nclang -Ofast -ffj-ocl -mllvm -polly -flto" F77=mpifrt FFLAGS="-Nclang -Ofast -ffj-ocl -mllvm -polly -flto"
		elif [[ "$1" = *"gem5"* ]]; then
			cd config/; rm -f config.guess config.sub; wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.guess'; wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.sub'; cd -
			./configure --disable-fortran -with-openmp --without-MPI CC=fcc CFLAGS="-Nclang -Ofast -ffj-no-largepage -ffj-ocl -mllvm -polly -flto" CXX=FCC CXXFLAGS="-Nclang -Ofast -ffj-no-largepage -ffj-ocl -mllvm -polly -flto" F77=frt FFLAGS="-Nclang -Ofast -ffj-no-largepage -ffj-ocl -mllvm -polly -flto"
		fi
		sed -i -e 's/ -openmp/ -fopenmp/g' ./config/Makefile.config
		make -j
		cd $ROOTDIR/$BM/
	fi
	if [ ! -f ./metis-4.0.3/graphchk ]; then
		if ! [ -f ./metis-4.0.3.tar.gz ]; then wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz; fi
		tar xzf ./metis-4.0.3.tar.gz
		cd ./metis-4.0.3/
		if [ -z $1 ]; then
			sed -i -e 's/CC = cc/CC = icc/g' -e 's/OPTFLAGS = -O2\s*$/OPTFLAGS = -O2 -ipo -xHost/g' ./Makefile.in
		elif [[ "$1" = *"gnu"* ]]; then
			sed -i -e 's/CC = cc/CC = gcc/g' -e 's/OPTFLAGS = -O2\s*$/OPTFLAGS = -O2 -march=native/g' ./Makefile.in
		elif [[ "$1" = *"fuji"* ]]; then
			sed -i -e 's/CC = cc/CC = fcc/g' -e 's/OPTFLAGS = -O2\s*$/OPTFLAGS = -Nclang -O2 -march=armv8.3-a+sve -mllvm -polly -flto/g' ./Makefile.in
		elif [[ "$1" = *"gem5"* ]]; then
			sed -i -e 's/CC = cc/CC = fcc/g' -e 's/OPTFLAGS = -O2\s*$/OPTFLAGS = -Nclang -O2 -march=armv8.3-a+sve -mllvm -polly -flto -ffj-no-largepage/g' ./Makefile.in
		fi
		make
		cd $ROOTDIR/$BM/
	fi
	if [ ! -f $ROOTDIR/dep/mfem/libmfem.a ]; then
		cd $ROOTDIR/dep/mfem/
		if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"laghos-v1.0"* ]]; then git checkout laghos-v1.0; fi
		git apply --check $ROOTDIR/patches/*1-mfem*.patch
		if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-mfem*.patch; fi
		if [[ "$1" = *"gnu"* ]]; then
			sed -i -e 's/icpc/g++/g' -e 's/-ipo -xHost/-march=native/g' ./config/defaults.mk
		elif [[ "$1" = *"fuji"* ]]; then
			sed -i -e 's/icpc/FCC/g' -e 's/-ipo -xHost/-Nclang -Ofast -ffj-ocl -mllvm -polly -flto/g' ./config/defaults.mk
		elif [[ "$1" = *"gem5"* ]]; then
			sed -i -e 's/icpc/FCC/g' -e 's/-ipo -xHost/-Nclang -Ofast -ffj-no-largepage -ffj-ocl -mllvm -polly -flto/g' ./config/defaults.mk
		fi
		if [ -z $1 ] || [[ "$1" = *"gnu"* ]]; then
			make config MFEM_USE_MPI=YES MPICXX=mpicxx MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
		elif [[ "$1" = *"fuji"* ]]; then
			make config MFEM_USE_MPI=YES MPICXX=mpiFCC MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
		elif [[ "$1" = *"gem5"* ]]; then
			make config CMAKE_CXX_COMPILER=FCC MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
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
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -flto#g' -e 's/-ipo -xHost/-Nclang -Ofast -ffj-ocl -mllvm -polly -flto/g' -e 's/ -lirc -lsvml//g' ./makefile
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
		make
	elif [[ "$1" = *"fuji"* ]]; then
		cd serial/
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -flto#g' -e 's/$(LAGHOS_LIBS) $(LDFLAGS)/$(LDFLAGS) $(LAGHOS_LIBS)/g' -e 's/-ipo -xHost/-Nclang -Ofast -ffj-no-largepage -ffj-ocl -mllvm -polly -flto/g' -e 's/ -lirc -lsvml//g' ./makefile
		sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' ./laghos_solver_s.hpp
		make
		cd -
		cp serial/laghos .
	fi
	cd $ROOTDIR
fi

