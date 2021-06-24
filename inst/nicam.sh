#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

if [[ "$1" = *"gem5"* ]]; then
	echo "ERR: no point, because this one needs MPI"; exit 1
fi

BM="NICAM"
VERSION="3f758000ffce6ee95a27fb6099f654ecdc5e3add"
if [ ! -f $ROOTDIR/$BM/bin/nhm_driver ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/src
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./Makefile
		sed -i -e 's/mpiifort/mpif90/' -e 's/mpiicc/mpicc/' -e 's/-shared-intel/-static -static-intel -qopenmp-link=static/g' -e 's/^LFLAGS = /LFLAGS = -static -static-intel -qopenmp-link=static /' ../sysdep/Makedef.Linux64-intel-impi
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -static#g' ./Makefile
		cp ../sysdep/Makedef.Linux64-gnu-openmpi ../sysdep/Makedef.Linux64-intel-impi
		sed -i -e 's/-O3/-O3 -march=native/g' -e 's/ -pedantic-errors//' -e 's/std=f2003/std=f2008/' -e 's/^LFLAGS = /LFLAGS = -static /' ../sysdep/Makedef.Linux64-intel-impi
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		cp ../sysdep/Makedef.FX10 ../sysdep/Makedef.Linux64-intel-impi
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' -e 's/-Kfast/-Kfast,openmp,ocl,largepage,lto/g' ../sysdep/Makedef.Linux64-intel-impi
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		cp ../sysdep/Makedef.Linux64-gnu-openmpi ../sysdep/Makedef.Linux64-intel-impi
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' -e 's/FAST = -O3/FAST = -O3 -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' -e 's/CFLAGS = -O3/CFLAGS = -O3 -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/g' ../sysdep/Makedef.Linux64-intel-impi
	fi
	# need this because its also used in runs, so overwrite with compiler specifics instead of changing them
	export NICAM_SYS=Linux64-intel-impi
	make ENABLE_OPENMP=1
	cd '../test/case/jablonowski'
	make ENABLE_OPENMP=1
	unset NICAM_SYS
	cd $ROOTDIR
fi

