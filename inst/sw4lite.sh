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
	export LD_LIBRARY_PATH=$ROOTDIR/dep/mpistub/lib:$LD_LIBRARY_PATH
else
	echo 'wrong compiler'
	exit 1
fi

BM="SW4lite"
VERSION="5ab8063ecdc94bdb59a5e65396c85bd54f9e0916"
if [[ $HOSTNAME = *"${XEONHOST}"* ]] || [[ "`hostname -s`" = *"fn01"* ]]; then
	HHOST="${XEONHOST}"
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	HHOST="${IKNLHOST}"
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	HHOST="${IKNMHOST}"
else
	echo "Unsupported host"
	exit
fi
if [ ! -f $ROOTDIR/$BM/optimize_mp_${HHOST}/sw4lite ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	sed -i -e "s/^HOSTNAME := /HOSTNAME := ${HHOST} #/g" ./Makefile
	sed -i -e "s/quadknl/${HHOST}/g" ./Makefile
	sed -i -e "s#/opt/intel/compilers_and_libraries_2017/linux#`dirname $MKLROOT`#g"  ./Makefile
	if [[ $HOSTNAME = *"${XEONHOST}"* ]] || [[ "`hostname -s`" = *"fn01"* ]]; then
		sed -i -e "s/-xmic-avx512/#NOKNL-xmic-avx512/g" ./Makefile
	fi
	if [ -z $1 ]; then
		sed -i -e 's/mpifort/mpif90/' -e 's/mpiifort/mpif90/' -e 's/mpiicpc/mpicxx/' -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' -e "s#MKL_PATH = .*#MKL_PATH = $MKLROOT/lib/intel64#" -e 's#-lmkl_intel_lp64 -lmkl_core -lmkl_sequential#-Wl,--start-group ${MKL_PATH}/libmkl_intel_lp64.a ${MKL_PATH}/libmkl_sequential.a ${MKL_PATH}/libmkl_core.a -Wl,--end-group#' -e 's/-lintlc/-lirc/' ./Makefile
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/mpifort/mpif90/' -e 's/mpiifort/mpif90/' -e 's/mpiicpc/mpicxx/' -e 's/= icc/= gcc/' -e 's/-ipo -xHost/-march=native/g' -e 's# -I${ADVISOR_2018_DIR}/include# -m64 -I${MKLROOT}/include#g' -e 's#EXTRA_LINK_FLAGS = .*#EXTRA_LINK_FLAGS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgfortran -lquadmath -lpthread -lm -ldl -static#' ./Makefile
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"fuji"* ]]; then
		sed -i -e 's/mpifort/mpifrtpx/' -e 's/mpiifort/mpifrtpx/' -e 's/mpiicpc/mpiFCCpx/' -e 's/= icc/= mpifccpx/' -e 's/-ipo -xHost/-Kfast/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s#EXTRA_LINK_FLAGS = .*#EXTRA_LINK_FLAGS = -SSL2BLAMP#g" ./Makefile
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/mpifort/frtpx/' -e 's/mpiifort/frtpx/' -e 's/mpiicpc/FCCpx/' -e 's/= icc/= fccpx/' -e 's/-ipo -xHost/-Kfast -Knolargepage/g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s#EXTRA_LINK_FLAGS = .*#EXTRA_LINK_FLAGS = -SSL2BLAMP -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi#g" ./Makefile
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' -e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' -e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE; done
	fi
	make
	cd $ROOTDIR
fi

