#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="CoMD"
VERSION="3d48396b77ca8caa3124bc2391f9139c3ffb556c"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/bin/CoMD-openmp-mpi ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src-openmp/
	#cp Makefile.vanilla Makefile
	instrument_kernel "$1" $ROOTDIR/$BM/
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/OTHER_LIB =/OTHER_LIB = -static -static-intel -qopenmp-link=static/' ./Makefile
	elif [[ "$1" = *"gnu"* ]]; then
		if [ -n "$FJMPI" ]; then sed -i -e 's/^CC =.*/CC = mpifcc/' ./Makefile; fi
		sed -i -e "s/-ipo -xHost/-march=native -flto ${MAYBESTATIC}/g" ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's/^CC =.*/CC = mpifcc/' -e 's/-ipo -xHost/-Kfast,openmp,ocl,largepage/g' ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's/^CC =.*/CC = mpifcc/' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/g' ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/^DO_MPI =.*/DO_MPI = OFF/g' -e 's/^CC =.*/CC = fcc/' ./Makefile
		sed -i -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		sed -i -e '/.*include.*stdio\.h/i #ifndef _POSIX_C_SOURCE\n#define _POSIX_C_SOURCE 199309L\n#endif' -e '/.*include.*time\.h/i #include <sys\/types.h>\n#include <signal.h>' -e '/.*include.*mpi\.h/d' $ROOTDIR/$BM/src-openmp/CoMD.c
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's/^CC =.*/CC = mpifcc/' -e 's/-ipo -xHost/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -mllvm -polly -mllvm -polly-vectorizer=polly -flto=full/g' ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./Makefile
	fi
	make
	#missing mpi in fuji version caused change in binary name -> fix that
	if [[ "$1" = *"gem5"* ]]; then mv $ROOTDIR/$BM/bin/CoMD-openmp $ROOTDIR/$BM/bin/CoMD-openmp-mpi; fi
	cd $ROOTDIR
fi

