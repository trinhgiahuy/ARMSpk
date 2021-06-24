#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

BM="MiniFE"
VERSION="daeddf3bfaf3b521a932245fad9871336b53c166"
if [ ! -f $ROOTDIR/$BM/mkl/src/miniFE.x ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	for SUB in mkl openmp-opt openmp-opt-knl; do
		cd $ROOTDIR/$BM/$SUB/src
		if [[ "$1" = *"intel"* ]]; then
			sed -i -e 's/mpiicpc/mpicxx/' -e 's/=-L${ADVISOR/=-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./Makefile
		elif [[ "$1" = *"gnu"* ]]; then
			sed -i -e 's/mpiicpc/mpicxx/' -e 's/-ipo -x[a-zA-Z0-9\-]*/-march=native/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's#-L${ADVISOR_2018_DIR}/lib64 -littnotify#-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -static#g' -e 's# -mkl # -m64 -I$(MKLROOT)/include #g' ./Makefile
		elif [[ "$1" = *"fujitrad"* ]]; then
			if [[ "$SUB" = *"mkl"* ]] || [[ "$SUB" = *"knl"* ]]; then continue; fi
			sed -i -e 's/mpicxx/mpiFCC/g' -e 's/mpicc/mpifcc/g' -e 's/-ipo -x[a-zA-Z0-9\-]*/-Kfast,openmp,ocl,largepage/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's#-L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		elif [[ "$1" = *"fujiclang"* ]]; then
			if [[ "$SUB" = *"mkl"* ]] || [[ "$SUB" = *"knl"* ]]; then continue; fi
			sed -i -e 's/mpicxx/mpiFCC/g' -e 's/mpicc/mpifcc/g' -e 's/-ipo -x[a-zA-Z0-9\-]*/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's#-L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		elif [[ "$1" = *"gem5"* ]]; then
			if [[ "$SUB" = *"mkl"* ]] || [[ "$SUB" = *"knl"* ]]; then continue; fi
			sed -i -e 's/mpicxx/FCC/g' -e 's/mpicc/fcc/g' -e 's/-DHAVE_MPI/-DHAVE_NO_MPI/g' -e 's/-ipo -x[a-zA-Z0-9\-]*/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's#-L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		fi
		make
		if [[ "$1" = *"fuji"* ]] || [[ "$1" = *"gem5"* ]] || [[ "$1" = *"llvm12"* ]]; then cp miniFE.x ../../mkl/src/; fi
	done
	cd $ROOTDIR
fi

