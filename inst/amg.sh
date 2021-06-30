#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="AMG"
VERSION="295de9693eaabf6f7330ac3a35fd9bd4ad030522"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/test/amg ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	# avx512 on KNL/KNM causes errors in AMG exec
	instrument_kernel "$1" $ROOTDIR/$BM/
	if [[ "$1" = *"intel"* ]]; then
		if [ -n ${IKNLHOST} ] || [ -n ${IKNMHOST} ]; then
			# avoid odd crashes with avx512 on KNL/KNM
			sed -i -e 's/xHost/xCORE-AVX2/g' ./Makefile.include
		fi
		sed -i -e 's/INCLUDE_LFLAGS = /INCLUDE_LFLAGS = -static -static-intel -qopenmp-link=static /' ./Makefile.include
	elif [[ "$1" = *"gnu"* ]]; then
		#XXX: fugaku gcc lacks something -> ar: amg_linklist.o: plugin needed to handle lto object
		if [ -n "$FJMPI" ]; then sed -i -e 's/^CC =.*/CC = mpifcc/g' ./Makefile.include; fi
		sed -i -e "s/-ipo -xHost/-march=native -fno-lto ${MAYBESTATIC}/g" ./Makefile.include
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile.include
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's/^CC =.*/CC = mpifcc/g' -e 's/-ipo -xHost/-Kfast,openmp,ocl,largepage/g' ./Makefile.include
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile.include
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's/^CC =.*/CC = mpifcc/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/g' ./Makefile.include
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile.include
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/define HYPRE_MPI_INT MPI_LONG_LONG.*/define HYPRE_MPI_INT MPI_LONG_LONG_INT/g' ./HYPRE.h
		sed -i -e 's/^CC =.*/CC = fcc/g' -e 's/ -DTIMER_USE_MPI//g' -e 's/-ipo -xHost/-DHYPRE_SEQUENTIAL=1 -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' ./Makefile.include
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile.include
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's/^CC =.*/CC = mpifcc/g' -e 's/-ipo -xHost/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin/g' ./Makefile.include
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./Makefile.include
	fi
	make
	cd $ROOTDIR
fi
