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
elif [[ "$1" = *"fuji"* ]]; then
	module load FujitsuCompiler/202007
	export LD_LIBRARY_PATH=$ROOTDIR/dep/mpistub/lib:$LD_LIBRARY_PATH
else
	echo 'wrong compiler'
	exit 1
fi

BM="SWFFT"  # fortran version is 5-10% faster in my tests
VERSION="d0ef31454577740fbb87618cc35789b7ef838238"
if [ ! -f $ROOTDIR/$BM/build.openmp/TestFDfft ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	if [ ! -f $ROOTDIR/$BM/fftw/bin/fftw-wisdom ]; then
		if [ ! -f fftw-3.3.4.tar.gz ]; then wget http://fftw.org/fftw-3.3.4.tar.gz; fi
		tar xzf fftw-3.3.4.tar.gz
		cd ./fftw-3.3.4/
		if [ -z $1 ]; then
			./configure --prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=icc
			if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
				make -j CFLAGS="-O3 -ipo -xHost -xCORE-AVX2 -fp-model fast=2 -no-prec-div -qoverride-limits"
			else
				make -j CFLAGS="-O3 -ipo -xHost -xCORE-AVX512 -fp-model fast=2 -no-prec-div -qoverride-limits"
			fi
		elif [[ "$1" = *"gnu"* ]]; then
			./configure --prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=gcc
			make -j CFLAGS="-O3 -march=native"
		elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
			./configure --host=aarch64-unknown-linux-gnu --build=x84_64-unknown-linux-gnu \
				--prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran CC=fccpx
			make -j CFLAGS="-O3 -Kfast"
		elif [[ "$1" = *"fuji"* ]]; then
			./configure --host=aarch64-unknown-linux-gnu --build=x84_64-unknown-linux-gnu \
				--prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran CC=fccpx
			make -j CFLAGS="-O3"
		fi
		make install
		cd $ROOTDIR/$BM/
	fi
	export oldPATH=$PATH
	export PATH=$ROOTDIR/$BM/fftw/bin:$oldPATH
	if [ -z $1 ]; then
		sed -i -e 's/-L${ADVISOR/-lmpi_cxx -static -static-intel -qopenmp-link=static -L${ADVISOR/' ./GNUmakefile
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -lmpi_cxx -static#g' ./GNUmakefile
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -lmpi_cxx -static#g' ./GNUmakefile.openmp
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		sed -i -e 's/?= mpicc/?= mpifccpx/g' -e 's/?= mpicxx/?= mpiFCCpx/g' -e 's/?= mpif90/?= mpifrtpx/g' -e 's/-ipo -xHost/-Kfast/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify##g" -e 's/DFFT_MPI_FLDFLAGS ?=.*/DFFT_MPI_FLDFLAGS ?= --linkstl=libfjc++/g' ./GNUmakefile
		sed -i -e 's/-ipo -xHost/-Kfast/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify##g" ./GNUmakefile.openmp
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"fuji"* ]]; then
		sed -i -e 's/?= mpicc/?= fccpx/g' -e 's/?= mpicxx/?= FCCpx/g' -e 's/?= mpif90/?= frtpx/g' -e 's/-ipo -xHost//g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" -e 's/DFFT_MPI_FLDFLAGS ?=.*/DFFT_MPI_FLDFLAGS ?= -lfjc++ -lfjc++abi -lfjdemgl/g' -e "s# -lm# -L$ROOTDIR/$BM -Bstatic -lm#g" ./GNUmakefile
		sed -i -e 's/-ipo -xHost//g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" -e "s# -lm# -L$ROOTDIR/$BM -Bstatic -lm#g" ./GNUmakefile.openmp
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' -e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' -e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE; done
		sed -i -e '/use mpi/d' -e "/implicit none/a \  include 'mpif.h'" ./FDistribution.f90
		sed -i -e '/use mpi/d' -e "/implicit none/a \  include 'mpif.h'" ./TestFDfft.f90
	fi
	make -f GNUmakefile.openmp
	export PATH=$oldPATH
	unset oldPATH
	cd $ROOTDIR
fi

