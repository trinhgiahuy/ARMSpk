#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

if [[ "$1" = *"gem5"* ]]; then
	echo "ERR: no point, because this one needs MPI"; exit 1
fi

BM="NICAM"
VERSION="3f758000ffce6ee95a27fb6099f654ecdc5e3add"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/bin/nhm_driver ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/src
	export NICAM_SYS=Linux64-intel-impi
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./Makefile
		sed -i -e 's/mpiifort/mpif90/' -e 's/mpiicc/mpicc/' -e 's/-shared-intel/-static -static-intel -qopenmp-link=static/g' -e 's/^LFLAGS = /LFLAGS = -static -static-intel -qopenmp-link=static /' ../sysdep/Makedef.${NICAM_SYS}
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -flto ${MAYBESTATIC}#g" ./Makefile
		cp ../sysdep/Makedef.Linux64-gnu-openmpi ../sysdep/Makedef.${NICAM_SYS}
		if [ -n "$FJMPI" ]; then sed -i -e 's/ mpif90/ mpifrt/g' -e 's/ mpicc/ mpifcc/g' -e 's/ -m64//g' ../sysdep/Makedef.${NICAM_SYS}; fi
		sed -i -e 's/-O3/-O3 -march=native -fallow-argument-mismatch -fno-lto/g' -e 's/ -pedantic-errors//' -e 's/std=f2003/std=f2008/' -e "s/^LFLAGS = /LFLAGS = -fno-lto ${MAYBESTATIC} /" ../sysdep/Makedef.${NICAM_SYS}
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		cp ../sysdep/Makedef.FX10 ../sysdep/Makedef.${NICAM_SYS}
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ../sysdep/Makedef.${NICAM_SYS}
		sed -i -e 's/^PERF_MONIT.*/PERF_MONIT = -Ntl_notrt -U_FIPP_ -U_FAPP_/g' -e 's/FAST.*= -K/FAST = -Kfast,openmp,ocl,largepage,lto -K/g' -e 's/^CFLAGS =.*/CFLAGS = -Kfast,openmp,ocl,largepage/g' ../sysdep/Makedef.${NICAM_SYS}
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		cp ../sysdep/Makedef.FX10 ../sysdep/Makedef.${NICAM_SYS}
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ../sysdep/Makedef.${NICAM_SYS}
		sed -i -e 's/^PERF_MONIT.*/PERF_MONIT = -Ntl_notrt -U_FIPP_ -U_FAPP_/g' -e 's/FAST.*= -K/FAST = -Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto -K/g' -e 's/^CFLAGS =.*/CFLAGS = -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -fno-lto/g' ../sysdep/Makedef.${NICAM_SYS}
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./Makefile
		cp ../sysdep/Makedef.FX10 ../sysdep/Makedef.${NICAM_SYS}
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ../sysdep/Makedef.${NICAM_SYS}
		sed -i -e 's/^PERF_MONIT.*/PERF_MONIT = -Ntl_notrt -U_FIPP_ -U_FAPP_/g' -e 's/FAST.*= -K/FAST = -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto -K/g' -e 's/^CFLAGS =.*/CFLAGS = -Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -fno-lto/g' ../sysdep/Makedef.${NICAM_SYS}
	fi
	# need this because its also used in runs, so overwrite with compiler specifics instead of changing them
	make ENABLE_OPENMP=1
	cd '../test/case/jablonowski'
	make ENABLE_OPENMP=1
	unset NICAM_SYS
	cd $ROOTDIR
fi

