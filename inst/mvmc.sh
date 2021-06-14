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
elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
	sleep 0
elif [[ "$1" = *"gem5"* ]]; then
	#module load FujitsuCompiler/202007
	export LD_LIBRARY_PATH=$ROOTDIR/dep/mpistub/lib:$LD_LIBRARY_PATH
else
	echo 'wrong compiler'
	exit 1
fi

BM="MVMC"
VERSION="7c58766b180ccb1035e4c220208b64ace3c49cf2"
if [ ! -f $ROOTDIR/$BM/src/vmc.out ] || [ "x`ls -s $ROOTDIR/$BM/src/vmc.out | awk '{print $1}'`" = "x0" ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	for x in `/bin/grep -r 'init_by_array(' . | cut -d':' -f1 | sort -u`; do sed -i -e 's/init_by_array(/xxxinit_by_array(/' $x; done
	if [ -z $1 ]; then
		sed -i -e 's/blacs_intelmpi/blacs_openmpi/' -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' -e "s#L/usr/local/intel/composer_xe_2013/mkl#L$MKLROOT#g" ./Makefile_intel
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's# -I${ADVISOR_2018_DIR}/include# -m64 -I${MKLROOT}/include#g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -static#g' -e "s#L/usr/local/intel/composer_xe_2013/mkl#L$MKLROOT#g" -e 's/blacs_intelmpi/blacs_openmpi/' -e 's/mkl_intel_thread/mkl_gnu_thread/' -e 's/-lpthread/-lgomp -lpthread -lm -ldl/' -e 's/ -qopt-prefetch=3 -nofor-main//' -e 's/ -vec-report//' ./Makefile_intel
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's/ ifort/ gfortran/' -e 's/-implicitnone/-fimplicit-none/' ./pfapack/Makefile_intel
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's/ icc/ gcc/' -e 's/-no-ansi-alias//' ./sfmt/Makefile_intel
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		cp ./Makefile_kei ./Makefile_intel
		cp ./pfapack/Makefile_kei ./pfapack/Makefile_intel
		cp ./sfmt/Makefile_kei ./sfmt/Makefile_intel
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"gem5"* ]]; then
		cp ./Makefile_kei ./Makefile_intel
		cp ./pfapack/Makefile_kei ./pfapack/Makefile_intel
		cp ./sfmt/Makefile_kei ./sfmt/Makefile_intel
		sed -i -e 's/^CC .*=.*/CC = fccpx/' -e 's/^FC .*=.*/FC = frtpx/' -e "s#CFLAGS = #CFLAGS = -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s#^LIB = # -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi#g" ./Makefile_intel
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' -e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' -e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE; done
	fi
	make intel
	cd $ROOTDIR
fi

