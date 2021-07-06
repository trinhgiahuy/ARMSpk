#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="HPCG"
VERSION="5422fecd0a009a8731d0bd96b957d443297a53bc"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/build/bin/xhpcg ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	if [[ "$1" = *"intel"* ]]; then
		git apply --check $ROOTDIR/patches/*1-${BM}*.patch
		if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	else
		git apply --check $ROOTDIR/patches/*2-${BM}*.patch
		if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*2-${BM}*.patch; fi
	fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	mkdir -p build; cd build
	if [ -n "${IKNLHOST}" ] || [ -n "${IKNMHOST}" ]; then
		CONF="IMPI_IOMP_KNL"
	elif [ -n "${XEONHOST}" ] && [[ "$1" = *"intel"* ]]; then
		CONF="IMPI_IOMP_AVX2"
	else
		CONF="MPI_GCC_OMP"
	fi
	if [[ "$1" = *"intel"* ]]; then
		../configure $CONF
		sed -i -e 's/mpiicpc/mpicxx/' -e 's/= -L${ADVISOR/= -static -static-intel -qopenmp-link=static -L${ADVISOR/' ./setup/Make.$CONF
	elif [[ "$1" = *"gnu"* ]]; then
		../configure $CONF
		if [ -n "$FJMPI" ]; then sed -i -e 's/^CXX .*=.*/CXX = mpiFCC/g' ./setup/Make.$CONF; fi
		sed -i -e "s/-O3/-O3 -march=native -flto ${MAYBESTATIC}/g" ./setup/Make.$CONF
		sed -i -e 's/shared(local_residual/shared(n, local_residual/g' ../src/ComputeResidual.cpp
	elif [[ "$1" = *"fujitrad"* ]]; then
		../configure $CONF
		sed -i -e 's/^CXX .*=.*/CXX = mpiFCC/g' -e 's/-O3/-Kfast,openmp,ocl,largepage/g' -e 's/-ftree-vectorizer-verbose=0//g' ./setup/Make.$CONF
	elif [[ "$1" = *"fujiclang"* ]]; then
		../configure $CONF
		sed -i -e 's/^CXX .*=.*/CXX = mpiFCC/g' -e 's/-O3/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/g' -e 's/-ftree-vectorizer-verbose=0//g' ./setup/Make.$CONF
	elif [[ "$1" = *"gem5"* ]]; then
		../configure $CONF
		sed -i -e 's/^CXX .*=.*/CXX = FCC/g' -e 's/-O3/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's/-ftree-vectorizer-verbose=0//g' -e 's/^HPCG_OPTS .*=.*/HPCG_OPTS = -DHPCG_NO_MPI/g' ./setup/Make.$CONF
	elif [[ "$1" = *"llvm12"* ]]; then
		../configure $CONF
		sed -i -e 's/^CXX .*=.*/CXX = mpiFCC/g' -e 's/-O3/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=full/g' -e "s#^LINKFLAGS .*= #LINKFLAGS = -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" -e 's/-ftree-vectorizer-verbose=0//g' ./setup/Make.$CONF
		sed -i -e 's/shared(local_residual/shared(n, local_residual/g' ../src/ComputeResidual.cpp
	fi
	make
	cd $ROOTDIR
fi

