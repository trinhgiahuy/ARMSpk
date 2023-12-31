#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="FFVC"
VERSION="890a3f9bb3a5cf358504063a1751383b7d46f86d"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/bin/ffvc_mini ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	instrument_kernel "$1" $ROOTDIR/$BM/
	if [[ "$1" = *"intel"* ]]; then
		rm make_setting; ln -s make_setting.intel make_setting
		sed -i -e 's/= -lifport/= -static -static-intel -qopenmp-link=static -lifport/' ./make_setting
	elif [[ "$1" = *"gnu"* ]]; then
		rm make_setting; ln -s make_setting.gcc make_setting
		if [ -n "$FJMPI" ]; then sed -i -e 's/^CXX .*=.*/CXX = mpiFCC/g' -e 's/^F90 .*=.*/F90 = mpifrt/g' ./make_setting; fi
		#XXX: RAGE.... segfaults w/o -g, i give up, this is getting beyond stupid
		sed -i -e '/= -DPROF_MAPROF$/d' -e 's/= -lgfortran/= /g' -e 's/-O3/-O3 -march=native -fallow-argument-mismatch -fallow-invalid-boz -fno-lto/g' ./make_setting
		if [ -z "${MAYBESTATIC}" ]; then
			sed -i -e "s/\$(LIBS).*/\$(LIBS) -lgfortran/g" ./FFV/Makefile
		else
			sed -i -e "s/\$(LIBS).*/\$(LIBS) ${MAYBESTATIC} -l:libgfortran.a -l:libquadmath.a/g" ./FFV/Makefile
		fi
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	elif [[ "$1" = *"fujitrad"* ]]; then
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e '/= -DPROF_MAPROF$/d' -e 's/Kfast/Kfast,openmp,ocl,largepage/g' ./make_setting
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	elif [[ "$1" = *"fujiclang"* ]]; then
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e '/= -DPROF_MAPROF$/d' -e 's/^CXXFLAGS.*/CXXFLAGS = -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage/g' -e 's/-Cpp -Kfast/-Cpp -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto/g' ./make_setting
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	elif [[ "$1" = *"gem5"* ]]; then
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e '/= -DPROF_MAPROF$/d' -e 's/^CXXFLAGS.*/CXXFLAGS = -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's/-Cpp -Kfast/-Cpp -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,nolargepage,nolto/g' -e 's/= mpi/= /g' -e "s#^LIBS.*#LIBS = -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi\nCXXFLAGS += -I$ROOTDIR/dep/mpistub/include/mpistub\nF90FLAGS += -I$ROOTDIR/dep/mpistub/include/mpistub#g" ./make_setting
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	elif [[ "$1" = *"llvm12"* ]]; then
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e '/= -DPROF_MAPROF$/d' -e 's/^CXXFLAGS.*/CXXFLAGS = -Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=full/g' -e 's/-Cpp -Kfast/-Cpp -Kfast,openmp,ocl,largepage/g' -e "s#--linkfortran#-L$(readlink -f $(dirname $(which mpifcc))/../lib64) $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjhpctag.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjlang08.o -lfj90i -lfj90fmt_sve -lfj90f -lfjsrcinfo -lfj90rt#g" -e "s#^LIBS.*#LIBS = -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./make_setting
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	fi
	make
	cd $ROOTDIR
fi

