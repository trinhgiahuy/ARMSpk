#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="SWFFT"  # fortran version is 5-10% faster in my tests
VERSION="d0ef31454577740fbb87618cc35789b7ef838238"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/build.openmp/TestFDfft ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	if [ ! -f $ROOTDIR/$BM/fftw/bin/fftw-wisdom ]; then
		URL="http://fftw.org/fftw-3.3.4.tar.gz"; DEP=$(basename $URL)
		if [ ! -f $ROOTDIR/dep/${DEP} ]; then if ! wget ${URL} -O $ROOTDIR/dep/${DEP}; then echo "ERR: download failed for ${URL}"; exit 1; fi; fi; tar xzf $ROOTDIR/dep/${DEP}
		cd ./fftw-3.3.4/
		if [[ "$1" = *"intel"* ]]; then
			./configure --prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=icc
			make -j CFLAGS="-O3 -ipo -xHost -fp-model fast=2 -no-prec-div -qoverride-limits"
		elif [[ "$1" = *"gnu"* ]]; then
			if [ -n "$FJMPI" ]; then ./configure --prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran CC=gcc;
			else                     ./configure --prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=gcc; fi
			make -j CFLAGS="-O3 -march=native -fno-lto"
		elif [[ "$1" = *"fujitrad"* ]]; then
			./configure --host=aarch64-unknown-linux-gnu --build=x84_64-unknown-linux-gnu \
				--prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran CC=fcc
			make -j CFLAGS="-Kfast,ocl,largepage"
		elif [[ "$1" = *"fujiclang"* ]]; then
			./configure --host=aarch64-unknown-linux-gnu --build=x84_64-unknown-linux-gnu \
				--prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran CC=fcc
			make -j CFLAGS="-Nclang -Ofast -mcpu=a64fx+sve -ffj-ocl -ffj-largepage -fno-lto"
		elif [[ "$1" = *"gem5"* ]]; then
			./configure --host=aarch64-unknown-linux-gnu --build=x84_64-unknown-linux-gnu \
				--prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran CC=fcc
			make -j CFLAGS="-Nclang -Ofast -mcpu=a64fx+sve -ffj-ocl -ffj-no-largepage -fno-lto"
		elif [[ "$1" = *"llvm12"* ]]; then
			./configure --host=aarch64-unknown-linux-gnu --build=x84_64-unknown-linux-gnu \
				--prefix=`pwd`/../fftw --disable-mpi --enable-openmp --disable-fortran CC=clang
			make -j CFLAGS="-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -mllvm -polly -mllvm -polly-vectorizer=polly -fno-lto -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)"
		fi
		make install
		cd $ROOTDIR/$BM/
	fi
	export oldPATH=$PATH
	export PATH=$ROOTDIR/$BM/fftw/bin:$oldPATH
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/-L${ADVISOR/-lmpi_cxx -static -static-intel -qopenmp-link=static -L${ADVISOR/' ./GNUmakefile
	elif [[ "$1" = *"gnu"* ]]; then
		if [ -n "$FJMPI" ]; then sed -i -e 's/?= mpicc/?= mpifcc/g' -e 's/?= mpicxx/?= mpiFCC/g' -e 's/?= mpif90/?= mpifrt/g' ./GNUmakefile; fi
		sed -i -e 's/-ipo -xHost/-march=native -fno-lto/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -lmpi_cxx -fno-lto ${MAYBESTATIC}#g" ./GNUmakefile
		sed -i -e 's/-ipo -xHost/-march=native -fno-lto/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -lmpi_cxx -fno-lto ${MAYBESTATIC}#g" ./GNUmakefile.openmp
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's/?= mpicc/?= mpifcc/g' -e 's/?= mpicxx/?= mpiFCC/g' -e 's/?= mpif90/?= mpifrt/g' -e 's/-ipo -xHost -cpp/-Kfast,openmp,ocl,largepage,lto -cpp/g' -e 's/-ipo -xHost/-Kfast,openmp,ocl,largepage/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify##g" -e 's/DFFT_MPI_FLDFLAGS ?=.*/DFFT_MPI_FLDFLAGS ?= --linkstl=libfjc++/g' ./GNUmakefile
		sed -i -e 's/-ipo -xHost -fopenmp -cpp/-Kfast,openmp,ocl,largepage,lto -cpp/g' -e 's/-ipo -xHost/-Kfast,openmp,ocl,largepage/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify##g" ./GNUmakefile.openmp
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's/?= mpicc/?= mpifcc/g' -e 's/?= mpicxx/?= mpiFCC/g' -e 's/?= mpif90/?= mpifrt/g' -e 's/-ipo -xHost -cpp/-mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto -cpp/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify##g" ./GNUmakefile
		sed -i -e 's/-ipo -xHost -fopenmp -cpp/-Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto -cpp/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -fno-lto/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify##g" ./GNUmakefile.openmp
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/?= mpicc/?= fcc/g' -e 's/?= mpicxx/?= FCC/g' -e 's/?= mpif90/?= frt/g' -e 's/-ipo -xHost -cpp/-mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,nolargepage,nolto -cpp/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" ./GNUmakefile
		sed -i -e 's/-ipo -xHost -fopenmp -cpp/-Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,nolargepage,nolto -cpp/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" ./GNUmakefile.openmp
		sed -i -e '/use mpi/d' -e "/implicit none/a \  include 'mpif.h'" ./FDistribution.f90
		sed -i -e '/use mpi/d' -e "/implicit none/a \  include 'mpif.h'" ./TestFDfft.f90
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's/?= mpicc/?= mpifcc/g' -e 's/?= mpicxx/?= mpiFCC/g' -e 's/?= mpif90/?= mpifrt/g' -e 's/-ipo -xHost -cpp/-mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto -cpp/g' -e 's/-ipo -xHost/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g"  ./GNUmakefile
		sed -i -e 's/-ipo -xHost -fopenmp -cpp/-mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto -cpp/g' -e 's/-ipo -xHost/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./GNUmakefile.openmp
	fi
	make -f GNUmakefile.openmp
	export PATH=$oldPATH
	unset oldPATH
	cd $ROOTDIR
fi

