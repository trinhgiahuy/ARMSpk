#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="XSBench"
VERSION="4772cf0194e2ae6d6752c5cacb8cf063fbfef7d0"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/src/XSBench ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/src
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./Makefile
	elif [[ "$1" = *"gnu"* ]]; then
		if [ -n "$FJMPI" ]; then sed -i -e 's/CC = mpicc/CC = mpifcc/' ./Makefile; fi
		sed -i -e 's/= intel/= gnu/' -e 's/-flto/-flto -march=native/' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -flto ${MAYBESTATIC}#g" ./Makefile
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's/= intel/= gnu/' -e 's/CC = mpicc/CC = mpifcc/' -e 's/-flto/-Kfast,openmp,ocl,largepage/' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's/= intel/= gnu/' -e 's/CC = mpicc/CC = mpifcc/' -e 's/-flto/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/^MPI.*= yes/MPI = no/' -e 's/= intel/= gnu/' -e 's/CC = gcc/CC = fcc/' -e 's/-flto/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's/= intel/= gnu/' -e 's/CC = mpicc/CC = mpifcc/' -e 's/-flto/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin/' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./Makefile
	fi
	make
	cd $ROOTDIR
fi

