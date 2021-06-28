#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

BM="SW4lite"
VERSION="5ab8063ecdc94bdb59a5e65396c85bd54f9e0916"
HHOST="${XEONHOST}${IKNLHOST}${IKNMHOST}${FUJIHOST}"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/optimize_mp_${HHOST}/sw4lite ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	sed -i -e "s/^HOSTNAME := /HOSTNAME := ${HHOST} #/g" ./Makefile
	sed -i -e "s/quadknl/${HHOST}/g" ./Makefile
	sed -i -e "s#/opt/intel/compilers_and_libraries_2017/linux#`dirname $MKLROOT`#g"  ./Makefile
	if [ -z "${IKNLHOST}" ] && [ -z "${IKNMHOST}" ]; then
		sed -i -e "s/-xmic-avx512/#NOKNL-xmic-avx512/g" ./Makefile
	fi
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/mpifort/mpif90/' -e 's/mpiifort/mpif90/' -e 's/mpiicpc/mpicxx/' -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' -e "s#MKL_PATH = .*#MKL_PATH = $MKLROOT/lib/intel64#" -e 's#-lmkl_intel_lp64 -lmkl_core -lmkl_sequential#-Wl,--start-group ${MKL_PATH}/libmkl_intel_lp64.a ${MKL_PATH}/libmkl_sequential.a ${MKL_PATH}/libmkl_core.a -Wl,--end-group#' -e 's/-lintlc/-lirc/' ./Makefile
	elif [[ "$1" = *"gnu"* ]]; then
		if [ -n "$FJMPI" ]; then sed -i -e 's/mpifort/mpifrt/' -e 's/mpiifort/mpifrt/' -e 's/mpiicpc/mpiFCC/' -e 's/= icc/= mpifcc/' ./Makefile;
		else                     sed -i -e 's/mpifort/mpif90/' -e 's/mpiifort/mpif90/' -e 's/mpiicpc/mpicxx/' -e 's/= icc/= gcc/' ./Makefile; fi
		sed -i -e 's/-ipo -xHost/-march=native -flto/g' -e 's# -I${ADVISOR_2018_DIR}/include# -m64 -I${MKLROOT}/include#g' -e "s#EXTRA_LINK_FLAGS = .*#EXTRA_LINK_FLAGS = -Wl,--start-group \${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \${MKLROOT}/lib/intel64/libmkl_sequential.a \${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgfortran -lquadmath -lpthread -lm -ldl -flto ${MAYBESTATIC}#" ./Makefile
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's/mpifort/mpifrt/' -e 's/mpiifort/mpifrt/' -e 's/mpiicpc/mpiFCC/' -e 's/= icc/= mpifcc/' -e 's/-ipo -xHost/-Kfast,openmp,ocl,largepage,lto/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s#EXTRA_LINK_FLAGS = .*#EXTRA_LINK_FLAGS = -SSL2BLAMP#g" ./Makefile
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's/mpifort/mpifrt/' -e 's/mpiifort/mpifrt/' -e 's/mpiicpc/mpiFCC/' -e 's/= icc/= mpifcc/' -e 's/FORT_FLAGS = -ipo -xHost/FORT_FLAGS = -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s#EXTRA_LINK_FLAGS = .*#EXTRA_LINK_FLAGS = -SSL2BLAMP#g" ./Makefile
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/mpifort/frt/' -e 's/mpiifort/frt/' -e 's/mpiicpc/FCC/' -e 's/= icc/= fcc/' -e 's/FORT_FLAGS = -ipo -xHost/FORT_FLAGS = -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,nolargepage,nolto/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s#EXTRA_LINK_FLAGS = .*#EXTRA_LINK_FLAGS = -SSL2BLAMP -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi#g" ./Makefile
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's/mpifort/mpifrt/' -e 's/mpiifort/mpifrt/' -e 's/mpiicpc/mpiFCC/' -e 's/= icc/= mpifcc/' -e 's/FORT_FLAGS = -ipo -xHost/FORT_FLAGS = -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' -e 's/-ipo -xHost/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s#EXTRA_LINK_FLAGS = .*#EXTRA_LINK_FLAGS = -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib) $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjhpctag.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjlang08.o -lfj90rt2 -lssl2mtexsve -lssl2mtsve -lfj90i -lfj90fmt_sve -lfj90f -lfjsrcinfo -lfj90rt -lfjprofcore -lfjprofomp#g" ./Makefile
	fi
	make
	cd $ROOTDIR
fi

