#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

BM="Nekbone"
VERSION="d681db43f45f5437e5258b3663b5a92c078cfb57"
if [ ! -f $ROOTDIR/$BM/test/nek_mgrid/nekbone ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/test/nek_mgrid
	sed -i -e 's/lp = 10)/lp = 576)/' -e 's/lelt=100/lelt=1024/' SIZE
	for FILE in `/usr/bin/grep 'log2(' -r ../../ | cut -d':' -f1 | sort -u`; do sed -i -e 's/log2(/log2XXX(/g' $FILE; done
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' $ROOTDIR/$BM/src/makefile.template
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's# -I${ADVISOR_2018_DIR}/include# -m64 -I${MKLROOT}/include#g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -static#g' -e 's# -mkl # -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -static #g' -e 's# -mkl-static# -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -static#' $ROOTDIR/$BM/src/makefile.template
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's# -mkl-static# -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl#g' ./makenek
		sed -i -e 's/-fdefault-real-8/-fdefault-real-8 -fdefault-double-8/g' $ROOTDIR/$BM/src/makenek.inc
	elif [[ "$1" = *"fujitrad"* ]]; then
		# fancy flags from https://arxiv.org/pdf/2009.11806.pdf
		sed -i -e 's/-ipo -xHost/-CcdRR8 -Cpp -Fixed -Kfast,openmp,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" $ROOTDIR/$BM/src/makefile.template
		sed -i -e 's/gfortran/frt/g' -e 's/-fdefault-real-8 -x f77-cpp-input/-CcdRR8 -Cpp -Fixed -Kfast,openmp,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' $ROOTDIR/$BM/src/makenek.inc
		sed -i -e 's/^CC=.*/CC=mpifcc/g' -e 's/^F77=.*/F77=mpifrt/g' -e 's/-ipo -xHost/-CcdRR8 -Cpp -Fixed -Kfast,openmp,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" ./makenek
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" $ROOTDIR/$BM/src/makefile.template
		sed -i -e 's/gfortran/frt/g' -e 's/-fdefault-real-8/-fdefault-real-8 -fdefault-double-8 -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' $ROOTDIR/$BM/src/makenek.inc
		sed -i -e 's/^CC=.*/CC=mpifcc/g' -e 's/^F77=.*/F77=mpifrt/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" ./makenek
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/^IFMPI.*/IFMPI:=false/' -e "s#-ipo -xHost#-CcdRR8 -Cpp -Fixed -Kfast,ocl,nolargepage,nolto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation#g" -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" $ROOTDIR/$BM/src/makefile.template
		sed -i -e 's/gfortran/frt/g' -e "s#-fdefault-real-8 -x f77-cpp-input#-CcdRR8 -Cpp -Fixed -Kfast,ocl,nolargepage,nolto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation#g" $ROOTDIR/$BM/src/makenek.inc
		sed -i -e 's/^CC=.*/CC=fcc/g' -e 's/^F77=.*/F77=frt/g' -e 's/^#IFMPI="f/IFMPI="f/g' -e "s#-ipo -xHost#-CcdRR8 -Cpp -Fixed -Kfast,ocl,nolargepage,nolto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation#g" -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" ./makenek
		sed -i -e '277,287s#/#!#g' $ROOTDIR/$BM/src/makenek.inc
		sed -i -e 's/^CFLAGS+=-DMPI/#CFLAGS+=-DMPI/g' $ROOTDIR/$BM/src/jl/Makefile
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's/-ipo -xHost/-CcdRR8 -Cpp -Fixed -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" $ROOTDIR/$BM/src/makefile.template
		sed -i -e 's/gfortran/frt/g' -e 's/-fdefault-real-8 -x f77-cpp-input/-CcdRR8 -Cpp -Fixed -Kfast,openmp,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' $ROOTDIR/$BM/src/makenek.inc
	fi
	./makenek NotUsedCasename $ROOTDIR/$BM/src
	cd $ROOTDIR
fi

