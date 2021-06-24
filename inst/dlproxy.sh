#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

BM="DLproxy"
VERSION="5087c437452c6cc3dbcf0bbaf40648c3155cfca9"
if [ ! -f $ROOTDIR/$BM/benchmarks/conv_gemm/main ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/benchmarks/conv_gemm/
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/CXXFLAGS=/CXXFLAGS=-static -static-intel -qopenmp-link=static /' ./Makefile
		sed -i -e 's#-mkl#-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl#g' ./Makefile
		make CC=icc CXX=icpc compile
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/-ipo -xHost/-march=native -static/g' ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		sed -i -e 's#-I${MKLROOT}/include#-m64 -I${MKLROOT}/include#g' -e 's#-mkl#-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl#g' ./Makefile
		make CC=gcc CXX=g++ compile
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's/-ipo -xHost/-Kfast,openmp,ocl,largepage/g' ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		sed -i -e "s#-DUSE_MKL -I\${MKLROOT}/include#-m64 -I$(dirname `which fcc`)/../include#g" ./Makefile
		sed -i -e "s#-mkl#-SSL2BLAMP#g" ./Makefile
		make CC=fcc CXX=FCC compile
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/g' ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		sed -i -e "s#-DUSE_MKL -I\${MKLROOT}/include#-m64 -I$(dirname `which fcc`)/../include#g" ./Makefile
		sed -i -e "s#-mkl#-SSL2BLAMP#g" ./Makefile
		make CC=fcc CXX=FCC compile
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		sed -i -e "s#-DUSE_MKL -I\${MKLROOT}/include#-m64 -I$(dirname `which fcc`)/../include#g" ./Makefile
		sed -i -e "s#-mkl#-SSL2BLAMP#g" ./Makefile
		make CC=fcc CXX=FCC compile
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's/-ipo -xHost/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin/g' ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./Makefile
		sed -i -e "s#-DUSE_MKL -I\${MKLROOT}/include#-m64 -I$(dirname `which fcc`)/../include#g" ./Makefile
		sed -i -e "s#-mkl#-lfj90rt2 -lssl2mtexsve -lssl2mtsve -lfj90i -lfj90fmt_sve -lfj90f -lfjsrcinfo -lfj90rt -lfjprofcore -lfjprofomp#g" ./Makefile
		make CC=clang CXX=clang++ compile
	fi
	cd $ROOTDIR
fi

