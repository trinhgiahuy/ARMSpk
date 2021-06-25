#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

BM="FFVC"
VERSION="890a3f9bb3a5cf358504063a1751383b7d46f86d"
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
		#XXX: RAGE.... segfaults w/o -g, i give up, this is getting beyond stupid
		sed -i -e 's/= -lgfortran/= /g' -e 's/-O3/-O3 -g -march=native/g' ./make_setting
		sed -i -e 's/\$(LIBS).*/\$(LIBS) -static -l:libgfortran.a -l:libquadmath.a/g' ./FFV/Makefile
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	elif [[ "$1" = *"fujitrad"* ]]; then
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/Kfast/Kfast,openmp,ocl,largepage/g' ./make_setting
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	elif [[ "$1" = *"fujiclang"* ]]; then
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/^CXXFLAGS.*/CXXFLAGS = -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/g' -e 's/-Cpp -Kfast/-Cpp -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' ./make_setting
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	elif [[ "$1" = *"gem5"* ]]; then
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/^CXXFLAGS.*/CXXFLAGS = -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's/-Cpp -Kfast/-Cpp -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -Kfast,ocl,nolargepage,nolto/g' -e 's/= mpi/= /g' -e "s#^LIBS.*#LIBS = -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi\nCXXFLAGS += -I$ROOTDIR/dep/mpistub/include/mpistub\nF90FLAGS += -I$ROOTDIR/dep/mpistub/include/mpistub#g" ./make_setting
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	elif [[ "$1" = *"llvm12"* ]]; then
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/^CXXFLAGS.*/CXXFLAGS = -Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin/g' -e "s#--linkfortran#-L$(readlink -f $(dirname $(which mpifcc))/../lib64) $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjhpctag.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/lib64/fjlang08.o -lfj90i -lfj90fmt_sve -lfj90f -lfjsrcinfo -lfj90rt#g" -e "s#^LIBS.*#LIBS = -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./make_setting
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	fi
	make
	cd $ROOTDIR
fi

