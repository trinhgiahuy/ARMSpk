#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
if [ -z $1 ]; then
	source $ROOTDIR/conf/intel.cfg
	source $INTEL_PACKAGE intel64 > /dev/null 2>&1
	export I_MPI_CC=icc
	export I_MPI_CXX=icpc
	export I_MPI_F77=ifort
	export I_MPI_F90=ifort
	alias ar=`which xiar`
	alias ld=`which xild`
	export ADVISOR_2018_DIR=${ADVISOR_2019_DIR}

	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load openmpi@3.1.6%intel@19.0.1.144
	export OMPI_CC=$I_MPI_CC
	export OMPI_CXX=$I_MPI_CXX
	export OMPI_F77=$I_MPI_F77
	export OMPI_FC=$I_MPI_F90
elif [[ "$1" = *"gnu"* ]]; then
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
	spack load openmpi@3.1.6%gcc@8.4.0
	export OMPI_CC=gcc
	export OMPI_CXX=g++
	export OMPI_F77=gfortran
	export OMPI_FC=gfortran
elif [[ "$1" = *"fuji"* ]]; then
	sleep 0
elif [[ "$1" = *"gem5"* ]]; then
	module load FujitsuCompiler/202007
else
	echo 'wrong compiler'
	exit 1
fi

BM="DLproxy"
VERSION="5087c437452c6cc3dbcf0bbaf40648c3155cfca9"
if [ ! -f $ROOTDIR/$BM/benchmarks/conv_gemm/main ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/benchmarks/conv_gemm/
	if [ -z $1 ]; then
		sed -i -e 's/CXXFLAGS=/CXXFLAGS=-static -static-intel -qopenmp-link=static /' ./Makefile
		sed -i -e 's#-mkl#-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl#g' ./Makefile
		make CC=icc CXX=icpc compile
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/-ipo -xHost/-march=native -static/g' ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		sed -i -e 's#-I${MKLROOT}/include#-m64 -I${MKLROOT}/include#g' -e 's#-mkl#-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl#g' ./Makefile
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
		make CC=gcc CXX=g++ compile
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		sed -i -e 's/-ipo -xHost/-Nclang -Ofast -ffj-ocl -mllvm -polly -flto/g' ./Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -flto#g' ./Makefile
		sed -i -e "s#-DUSE_MKL -I\${MKLROOT}/include#-m64 -I$(dirname `which fccpx`)/../include#g" -e 's#-I${ADVISOR_2018_DIR}/include##g' -e 's#-L${ADVISOR_2018_DIR}/lib64 -littnotify#-flto#g' ./Makefile
		sed -i -e "s#-mkl#-SSL2BLAMP#g" ./Makefile
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
		make CC=fccpx CXX=FCCpx compile
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/-ipo -xHost/-Nclang -Ofast -ffj-no-largepage -ffj-ocl -mllvm -polly -flto/g' ./Makefile
		sed -i -e "s#-DUSE_MKL -I\${MKLROOT}/include#-m64 -I$(dirname `which fccpx`)/../include#g" -e 's#-I${ADVISOR_2018_DIR}/include##g' -e 's#-L${ADVISOR_2018_DIR}/lib64 -littnotify#-flto#g' ./Makefile
		sed -i -e "s#-mkl#-SSL2BLAMP#g" ./Makefile
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
		make CC=fccpx CXX=FCCpx compile
	fi
	cd $ROOTDIR
fi
