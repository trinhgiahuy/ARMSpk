#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

source $ROOTDIR/conf/qcd.sh

BM="QCD"
VERSION="07277047b170529caa5fcd164afd814e70286ce4"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2_111 ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/src
	if [[ "$1" = *"intel"* ]]; then
		TYPE=ifort
		sed -i -e 's/-shared-intel -mcmodel=medium/-static -static-intel -qopenmp-link=static/' -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./make.${TYPE}.inc
	elif [[ "$1" = *"gnu"* ]]; then
		#XXX: maybe this helps here too... like magic, fucking broken maprof code
		sed -i '/call maprof_set_fp_ops(SEC_BICGSTAB/,+17d' ./ccs_qcd_solver_bench.F90
		TYPE=gfortran
		if [ -n "$FJMPI" ]; then sed -i -e 's/ mpicc/ mpifcc/g' -e 's/ mpif90/ mpifrt/g' ./make.${TYPE}.inc; fi
		sed -i -e 's/-march=core2 -msse3/-march=native -fno-lto -fno-inline-small-functions/' -e "s/LDFLAGS = /LDFLAGS = -fno-lto ${MAYBESTATIC} /" ./make.${TYPE}.inc
		sed -i -e 's/mcmodel=medium/mcmodel=large/g' ./make.${TYPE}.inc
	elif [[ "$1" = *"fujitrad"* ]]; then
		# crashing in some stupid yaml shit with fujitsu compilers
		sed -i '/call maprof_set_fp_ops(SEC_BICGSTAB/,+17d' ./ccs_qcd_solver_bench.F90
		TYPE=fx10
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make.${TYPE}.inc
		sed -i -e "s#-Kprefetch.*#-Kprefetch -Kfast,openmp,ocl,largepage,lto#" ./make.${TYPE}.inc
	elif [[ "$1" = *"fujiclang"* ]]; then
		# crashing in some stupid yaml shit with fujitsu compilers
		sed -i '/call maprof_set_fp_ops(SEC_BICGSTAB/,+17d' ./ccs_qcd_solver_bench.F90
		TYPE=fx10
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make.${TYPE}.inc
		sed -i -e "s#-Kprefetch.*#-Kprefetch -Nclang -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto#" ./make.${TYPE}.inc
	elif [[ "$1" = *"gem5"* ]]; then
		# crashing in some stupid yaml shit with fujitsu compilers
		sed -i '/call maprof_set_fp_ops(SEC_BICGSTAB/,+17d' ./ccs_qcd_solver_bench.F90
		TYPE=fx10
		sed -i -E 's/mpi(fcc|FCC|frt)px/\1/g' ./make.${TYPE}.inc
		sed -i -e "s#INCLUDE =.*#INCLUDE = -I./ -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s#\$(FFLAGS).*#\$(FFLAGS) -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" -e "s#-Kprefetch.*#-Kprefetch -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,nolargepage,nolto -I$ROOTDIR/dep/mpistub/include/mpistub#" ./make.${TYPE}.inc
		sed -i -e "s#^LIBS += -lmaprof_f#LIBS += -lmaprof_f -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" ./Makefile
		sed -i -e "s#\$(PFLAGS).*#\$(PFLAGS) -I$ROOTDIR/dep/mpistub/include/mpistub -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" ./ma_prof/src/Makefile
		sed -i -e '/use mpi/d' ./comlib.F90
		sed -i "0,/implicit none/s//implicit none\n  include 'mpif.h'/" ./comlib.F90
		sed -i -e '/use mpi/d' -e "/implicit none/a \  include 'mpif.h'" ./ma_prof/src/mod_maprof.F90
	elif [[ "$1" = *"llvm12"* ]]; then
		# crashing in some stupid yaml shit with fujitsu compilers
		sed -i '/call maprof_set_fp_ops(SEC_BICGSTAB/,+17d' ./ccs_qcd_solver_bench.F90
		TYPE=fx10
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make.${TYPE}.inc
		sed -i -e "s#-Kprefetch.*#-Kprefetch -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto#" ./make.${TYPE}.inc
	fi
	for TEST in $TESTCONF; do
		PX="`echo $TEST | cut -d '|' -f3`"
		PY="`echo $TEST | cut -d '|' -f4`"
		PZ="`echo $TEST | cut -d '|' -f5`"
		if [ ! -f $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2_${PX}${PY}${PZ} ]; then
			make clean >> /dev/null 2>&1
			make MAKE_INC=make.${TYPE}.inc CLASS=2 PX=${PX} PY=${PY} PZ=${PZ}
			mv $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2 $ROOTDIR/$BM/src/ccs_qcd_solver_bench_class2_${PX}${PY}${PZ}
		fi
		if [[ "$1" = *"gem5"* ]] && [ $PX -eq 1 ] && [ $PY -eq 1 ] && [ $PZ -eq 1 ]; then break; fi
	done
	unset TYPE
	cd $ROOTDIR
fi

