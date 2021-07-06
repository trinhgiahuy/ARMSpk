#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="MiniTri"
VERSION="9771c71f3d25023fc50bc6e84a905d6d50e81151"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/miniTri/linearAlgebra/MPI/miniTri.exe ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	for SUB in MPI openmp; do
		cd $ROOTDIR/$BM/miniTri/linearAlgebra/$SUB
		if [[ "$1" = *"intel"* ]]; then
			sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./Makefile
		elif [[ "$1" = *"gnu"* ]]; then
			if [ -n "$FJMPI" ]; then sed -i -e 's/ mpicxx/ mpiFCC/g' ./Makefile; fi
			sed -i -e 's/ icpc/ g++/g' -e 's/-ipo -xHost/-march=native -fno-lto/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fno-lto ${MAYBESTATIC}#g" ./Makefile
		elif [[ "$1" = *"fujitrad"* ]]; then
			sed -i -e 's/ mpicxx/ mpiFCC/g' -e 's/ icpc/ FCC/g' -e 's/-ipo -xHost/-Kfast,ocl,largepage/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		elif [[ "$1" = *"fujiclang"* ]]; then
			sed -i -e 's/ mpicxx/ mpiFCC/g' -e 's/ icpc/ FCC/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -ffj-ocl -ffj-largepage -flto/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		elif [[ "$1" = *"gem5"* ]]; then
			if [[ "$SUB" = *"MPI"* ]]; then continue; fi
			sed -i -e 's/ icpc/ FCC/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		elif [[ "$1" = *"llvm12"* ]]; then
			sed -i -e 's/mpicxx/mpiFCC/g' -e 's/ icpc/ clang++/g' -e 's/$(CCC) $(INC/mpiFCC $(INC/g' -e 's/-ipo -xHost/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=full/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./Makefile
		fi
		make
	done
	# get an valid input
	URL="ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcsstruc5/bcsstk30.mtx.gz"; DEP=$(basename $URL)
	if [ ! -f $ROOTDIR/dep/${DEP} ]; then if ! wget ${URL} -O $ROOTDIR/dep/${DEP}; then echo "ERR: download failed for ${URL}"; exit 1; fi; fi
	cp $ROOTDIR/dep/${DEP} $ROOTDIR/$BM/; cd $ROOTDIR/$BM/; gunzip bcsstk30.mtx.gz
	cd $ROOTDIR
fi

