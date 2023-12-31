#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="MiniAMR"
VERSION="07297b3a2a46ebf08bc9be33d8d28ea21c0f4956"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/ref/ma.x ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	for SUB in ref openmp; do
		cd $ROOTDIR/$BM/$SUB
		sed -i -e 's/^LDFLAGS =.*/LDFLAGS = $(CFLAGS)/g' ./Makefile
		for x in `/bin/grep -r 'exchange(' . | cut -d':' -f1 | sort -u`; do sed -i -e 's/exchange(/xxxexchange(/' $x; done
		if [[ "$1" = *"intel"* ]]; then
			sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./Makefile
		elif [[ "$1" = *"gnu"* ]]; then
			if [ -n "$FJMPI" ]; then sed -i -e 's/^CC .*=.*/CC = mpifcc/' -e 's/^LD .*=.*/LD = mpifcc/' ./Makefile; fi
			sed -i -e 's/-ipo -xHost/-march=native -fno-lto -fcommon/g' -e 's/-qopenmp/-fopenmp/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fno-lto ${MAYBESTATIC}#g" ./Makefile
		elif [[ "$1" = *"fujitrad"* ]]; then
			sed -i -e 's/^CC .*=.*/CC = mpifcc/' -e 's/^LD .*=.*/LD = mpifcc/' -e 's/-ipo -xHost/-Kfast,openmp,ocl,largepage/g' -e 's/-qopenmp/-fopenmp/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify##g" ./Makefile
		elif [[ "$1" = *"fujiclang"* ]]; then
			sed -i -e 's/^CC .*=.*/CC = mpifcc/' -e 's/^LD .*=.*/LD = mpifcc/' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/g' -e 's/-qopenmp/-fopenmp/g' -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify##g" ./Makefile
		elif [[ "$1" = *"gem5"* ]]; then
			sed -i -e 's/^CC .*=.*/CC = fcc/' -e 's/^LD .*=.*/LD = fcc/' -e 's/ -lm / /g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's/-qopenmp/-fopenmp/g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi#g" ./Makefile
		elif [[ "$1" = *"llvm12"* ]]; then
			sed -i -e 's/^CC .*=.*/CC = mpifcc/' -e 's/^LD .*=.*/LD = mpifcc/' -e 's/-ipo -xHost/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=full -fcommon/g' -e 's/-qopenmp/-fopenmp/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./Makefile
			#XXX need -fcommon , see https://www.bountysource.com/issues/94002370-clang-11-generates-many-multiple-definitions-when-linking-mmg
		fi
		make
	done
	cd $ROOTDIR
fi

