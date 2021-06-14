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
	sleep 0
elif [[ "$1" = *"gem5"* ]]; then
	export LD_LIBRARY_PATH=$ROOTDIR/dep/mpistub/lib:$LD_LIBRARY_PATH
else
	echo 'wrong compiler'
	exit 1
fi


source $ROOTDIR/conf/qcd.sh

BM="QCD"
VERSION="07277047b170529caa5fcd164afd814e70286ce4"
if [ ! -f $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2_111 ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	if [ -z $1 ]; then
		TYPE=ifort
		sed -i -e 's/-shared-intel -mcmodel=medium/-static -static-intel -qopenmp-link=static/' -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./make.${TYPE}.inc
	elif [[ "$1" = *"gnu"* ]]; then
		TYPE=gfortran
		sed -i -e 's/-march=core2 -msse3/-march=native -fno-inline-small-functions/' -e 's/LDFLAGS = /LDFLAGS = -static /' ./make.${TYPE}.inc
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		TYPE=fx10
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"gem5"* ]]; then
		TYPE=fx10
		sed -i -e 's/mpifrtpx/frtpx/g' -e 's/mpifccpx/fccpx/g' -e "s#INCLUDE =.*#INCLUDE = -I./ -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s#\$(FFLAGS).*#\$(FFLAGS) -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" -e "s!#LIBS = !LIBS = -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort!g" -e "s#-Kprefetch.*#-Kprefetch -I$ROOTDIR/dep/mpistub/include/mpistub#" ./make.${TYPE}.inc
		sed -i -e "s#\$(PFLAGS).*#\$(PFLAGS) -I$ROOTDIR/dep/mpistub/include/mpistub -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" ./ma_prof/src/Makefile
		sed -i -e '/use mpi/d' ./comlib.F90
		sed -i "0,/implicit none/s//implicit none\n  include 'mpif.h'/" ./comlib.F90
		sed -i -e '/use mpi/d' -e "/implicit none/a \  include 'mpif.h'" ./ma_prof/src/mod_maprof.F90
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' -e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' -e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE; done
	fi
	prevPX=0; prevPY=0; prevPZ=0
	for TEST in $TESTCONF; do
		PX="`echo $TEST | cut -d '|' -f3`"
		PY="`echo $TEST | cut -d '|' -f4`"
		PZ="`echo $TEST | cut -d '|' -f5`"
		if [ $prevPX -eq $PX ] && [ $prevPY -eq $PY ] && [ $prevPZ -eq $PZ ] ; then continue; fi
		if [[ "$1" = *"gem5"* ]] && [[ "$((${PX}*${PY}*${PZ}))" -ne 1 ]]; then break; fi
		if [ ! -f $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2_${PX}${PY}${PZ} ]; then
			make clean >> /dev/null 2>&1
			make MAKE_INC=make.${TYPE}.inc CLASS=2 PX=${PX} PY=${PY} PZ=${PZ}
			mv $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2 $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2_${PX}${PY}${PZ}
		fi
	done
	unset TYPE
	cd $ROOTDIR
fi

