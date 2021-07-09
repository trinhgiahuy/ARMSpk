#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
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
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
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
		sed -i -e 's/FFLAGS += -DPROF_MAPROF/FFLAGS += -DNO_PROF_MAPROF/g' ./Makefile
		rm make_setting; ln -s make_setting.gcc make_setting
		if [ -n "$FJMPI" ]; then sed -i -e 's/ mpif90/ mpifrt/g' -e 's/ mpicc/ mpifcc/g' ./make_setting; fi
		sed -i -e "s/-O3/-O3 -march=native -fallow-argument-mismatch -fno-lto\nLIBS += -fno-lto ${MAYBESTATIC}/g" ./make_setting
	elif [[ "$1" = *"fujitrad"* ]]; then
		# crashing in some stupid yaml shit with fujitsu compilers
		sed -i -e 's/FFLAGS += -DPROF_MAPROF/FFLAGS += -DNO_PROF_MAPROF/g' ./Makefile
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/-Kfast/-Kfast,openmp,ocl,largepage/g' -e '/^FFLAGS = /a CFLAGS = -Kfast,openmp,ocl,largepage' ./make_setting
	elif [[ "$1" = *"fujiclang"* ]]; then
		# crashing in some stupid yaml shit with fujitsu compilers
		sed -i -e 's/FFLAGS += -DPROF_MAPROF/FFLAGS += -DNO_PROF_MAPROF/g' ./Makefile
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/-Kfast/-Nclang -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto/g' -e '/^FFLAGS = /a CFLAGS = -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage' ./make_setting
		#XXX: LD is frt, so clang lto does not work
	elif [[ "$1" = *"gem5"* ]]; then
		# crashing in some stupid yaml shit with fujitsu compilers
		sed -i -e 's/FFLAGS += -DPROF_MAPROF/FFLAGS += -DNO_PROF_MAPROF/g' ./Makefile
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/ = mpi/ = /g' -e "s#^FFLAGS = #FFLAGS = -I$ROOTDIR/dep/mpistub/include/mpistub #g" -e "s#mfunc=2#mfunc=2 -I$ROOTDIR/dep/mpistub/include/mpistub\nLIBS += -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" -e 's/-Kfast/-Nclang -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,nolargepage,nolto/g' -e "/^FFLAGS = /a CFLAGS = -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto -I$ROOTDIR/dep/mpistub/include/mpistub" ./make_setting
		sed -i -e '/use mpi/d' -e "/implicit none/a \  include 'mpif.h'" ./ma_prof/src/mod_maprof.F90
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's/FFLAGS += -DPROF_MAPROF/FFLAGS += -DNO_PROF_MAPROF/g' ./Makefile
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/-Kfast/-mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto/g' -e '/^FFLAGS = /a CFLAGS = -Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly' -e "s#^CFLAGS = #CFLAGS = -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./make_setting
		#XXX: LD is frt, so clang lto does not work
	fi
	make
	cd $ROOTDIR
fi

