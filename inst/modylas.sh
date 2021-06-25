#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

if [[ "$1" = *"fuji"* ]] || [[ "$1" = *"gem5"* ]]; then
	echo "WRN: DOES NOT compile in -Nclang mode"
fi

if [ ! -f $ROOTDIR/dep/modylas-mini-1.0.0.tar.gz ]; then
	echo "ERR: Cannot find modylas-mini-1.0.0.tar.gz"
	echo "Please download from: http://hpci-aplfs.aics.riken.jp/fiber/modylas-mini.html and place modylas-mini-1.0.0.tar.gz in ./dep subfolder"
	exit
fi

BM="MODYLAS" # (req. license agreement on website)
if [ ! -f $ROOTDIR/$BM/src/modylas_mini ]; then
	mkdir -p $ROOTDIR/$BM/
	tar xzf $ROOTDIR/dep/modylas-mini-1.0.0.tar.gz -C $ROOTDIR/$BM --strip-components 1
	cd $ROOTDIR/$BM/
	patch -p1 --forward < $ROOTDIR/patches/*1-${BM}*.patch
	instrument_kernel "$1" $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/src
	if [[ "$1" = *"intel"* ]]; then
		rm make_setting; ln -s make_setting.intel make_setting
		sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./make_setting
	elif [[ "$1" = *"gnu"* ]]; then
		rm make_setting; ln -s make_setting.gcc make_setting
		sed -i -e 's/-O3/-O3 -march=native\nLIBS += -static/g' ./make_setting
	elif [[ "$1" = *"fujitrad"* ]] || [[ "$1" = *"fujiclang"* ]]; then
		# crashing in some stupid yaml shit with fujitsu compilers
		sed -i -e 's/FFLAGS += -DPROF_MAPROF/FFLAGS += -DNO_PROF_MAPROF/g' ./Makefile
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/-Kfast/-Kfast,openmp,ocl,largepage/g' -e '/^FFLAGS = /a CFLAGS = -Kfast,openmp,ocl,largepage' ./make_setting
	elif [[ "$1" = *"gem5"* ]]; then
		# crashing in some stupid yaml shit with fujitsu compilers
		sed -i -e 's/FFLAGS += -DPROF_MAPROF/FFLAGS += -DNO_PROF_MAPROF/g' ./Makefile
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/ = mpi/ = /g' -e "s#^FFLAGS = #FFLAGS = -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s#mfunc=2#mfunc=2 -I$ROOTDIR/dep/mpistub/include/mpistub\nLIBS += -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" -e 's/-Kfast/-Kfast,openmp,ocl,nolargepage,nolto/g' -e '/^FFLAGS = /a CFLAGS = -Kfast,openmp,ocl,nolargepage' ./make_setting
		sed -i -e '/use mpi/d' -e "/implicit none/a \  include 'mpif.h'" ./ma_prof/src/mod_maprof.F90
	fi
	make
	cd $ROOTDIR
fi

