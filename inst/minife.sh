#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
if [ -z $1 ]; then
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
	source `echo $INTEL_PACKAGE | cut -d'/' -f-3`/mkl/bin/mklvars.sh intel64 > /dev/null 2>&1
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
	spack load openmpi@3.1.6%gcc@8.4.0
	export OMPI_CC=gcc
	export OMPI_CXX=g++
	export OMPI_F77=gfortran
	export OMPI_FC=gfortran
elif [[ "$1" = *"fuji"* ]]; then
	sleep 0
elif [[ "$1" = *"gem5"* ]]; then
	sleep 0; #module load FujitsuCompiler/202007
else
	echo 'wrong compiler'
	exit 1
fi

BM="MiniFE"
VERSION="daeddf3bfaf3b521a932245fad9871336b53c166"
if [ ! -f $ROOTDIR/$BM/mkl/src/miniFE.x ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	for SUB in mkl openmp-opt openmp-opt-knl; do
		cd $ROOTDIR/$BM/$SUB/src
		if [ -z $1 ]; then
			sed -i -e 's/mpiicpc/mpicxx/' -e 's/=-L${ADVISOR/=-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./Makefile
		elif [[ "$1" = *"gnu"* ]]; then
			sed -i -e 's/mpiicpc/mpicxx/' -e 's/-ipo -xHost/-march=native/g' -e 's/-ipo -xMIC-AVX512/-march=native/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's#-L${ADVISOR_2018_DIR}/lib64 -littnotify#-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -static#g' -e 's# -mkl # -m64 -I$(MKLROOT)/include #g' ./Makefile
			for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
		elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
			if [[ "$SUB" = *"mkl"* ]] || [[ "$SUB" = *"knl"* ]]; then continue; fi
			sed -i -e 's/mpicxx/mpiFCCpx/g' -e 's/mpicc/mpifccpx/g' -e 's/-ipo -xHost/-Nclang -Ofast -ffj-ocl -mllvm -polly -flto/g' -e 's/-ipo -xMIC-AVX512/-Nclang -Ofast -ffj-ocl -mllvm -polly -flto/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's#-L${ADVISOR_2018_DIR}/lib64 -littnotify# -flto#g' ./Makefile
			for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
		elif [[ "$1" = *"gem5"* ]]; then
			if [[ "$SUB" = *"mkl"* ]] || [[ "$SUB" = *"knl"* ]]; then continue; fi
			sed -i -e 's/mpicxx/FCCpx/g' -e 's/mpicc/fccpx/g' -e 's/-DHAVE_MPI/-DHAVE_NO_MPI/g' -e 's/-ipo -xHost/-Nclang -Ofast -ffj-no-largepage -ffj-ocl -mllvm -polly/g' -e 's/-ipo -xMIC-AVX512/-Nclang -Ofast -ffj-no-largepage -ffj-ocl -mllvm -polly/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's#-L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
			for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' -e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' -e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE; done
		fi
		make
		if [[ "$1" = *"fuji"* ]]; then cp miniFE.x ../../mkl/src/; fi
	done
	cd $ROOTDIR
fi

