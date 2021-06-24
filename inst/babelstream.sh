#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

BM="BabelStream"
VERSION="d9b089a0f94e9423b0653ca7ca533bd04c8501cb"
if [ ! -f $ROOTDIR/$BM/omp-stream ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	if [[ "$1" = *"intel"* ]]; then
		make -f OpenMP.make COMPILER=INTEL TARGET=CPU EXTRA_FLAGS="-static -static-intel -qopenmp-link=static"
	elif [[ "$1" = *"gnu"* ]]; then
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
		make -f OpenMP.make COMPILER=GNU TARGET=CPU EXTRA_FLAGS="-static"
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's/^COMPILER_GNU =.*/COMPILER_GNU = FCC/' ./OpenMP.make
		make -f OpenMP.make COMPILER=GNU TARGET=CPU EXTRA_FLAGS="-Kfast,openmp,ocl,largepage,assume=memory_bandwidth"
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's/^COMPILER_CLANG =.*/COMPILER_CLANG = FCC/' -e 's/-fopenmp=libomp/-fopenmp/' ./OpenMP.make
		make -f OpenMP.make COMPILER=CLANG TARGET=CPU EXTRA_FLAGS="-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto"
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/^COMPILER_CLANG =.*/COMPILER_CLANG = FCC/' -e 's/-fopenmp=libomp/-fopenmp/' ./OpenMP.make
		make -f OpenMP.make COMPILER=CLANG TARGET=CPU EXTRA_FLAGS="-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto"
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's/-fopenmp=libomp/-fopenmp/' ./OpenMP.make
		make -f OpenMP.make COMPILER=CLANG TARGET=CPU EXTRA_FLAGS="-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)"
	fi
	cd $ROOTDIR
fi

