#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="polybench"
source $ROOTDIR/conf/${BM}.sh
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if ! [ $(find $ROOTDIR/$BM/linear-algebra/blas/gemm -executable -type f | wc -l) -ge 1 ]; then
	mkdir -p $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/
	URL="https://downloads.sourceforge.net/project/polybench/polybench-c-4.2.1-beta.tar.gz"; DEP=$(basename $URL)
	if [ ! -f $ROOTDIR/dep/${DEP} ]; then if ! wget ${URL} -O $ROOTDIR/dep/${DEP}; then echo "ERR: download failed for ${URL}"; exit 1; fi; fi
	tar xzf $ROOTDIR/dep/${DEP} -C $ROOTDIR/$BM --strip-components 1
	if patch --dry-run -s -f -p1 < $ROOTDIR/patches/*1-${BM}*.patch; then
		patch -p1 --forward < $ROOTDIR/patches/*1-${BM}*.patch
	fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	if [[ "$1" = *"intel"* ]]; then
		#COMPILE="icc -O3 -xHost -static -static-intel -I${ADVISOR_2018_DIR}/include"
		COMPILE="icc -O3 -march=corei7 -static -static-intel -I${ADVISOR_2018_DIR}/include"
		LINK="-L${ADVISOR_2018_DIR}/lib64 -littnotify -L/usr/lib64 -lm"
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e '/^#define STARTSDE/i #ifndef  MYLITTLEMERMAID_H\n#define  MYLITTLEMERMAID_H\nstatic void __attribute__ ((noinline)) ariel_enable() { asm (""); return; }\nstatic void __attribute__ ((noinline)) ariel_fence()  { asm (""); return; }\n#endif //MYLITTLEMERMAID_H' -e 's/__SSC_MARK(0x111);/ariel_enable(); __SSC_MARK(0x111);/g' -e 's/__SSC_MARK(0x222);/__SSC_MARK(0x222); ariel_fence();/g' $FILE; done
	elif [[ "$1" = *"gnu"* ]]; then
		COMPILE="gcc -O3 -march=native -flto -funroll-loops -ffast-math -ftree-vectorize"
		LINK="-flto ${MAYBESTATIC} -L/usr/lib64 -lm"
	elif [[ "$1" = *"fujitrad"* ]]; then
		COMPILE="fcc -Kfast,openmp,ocl,largepage"
		LINK=""
	elif [[ "$1" = *"fujiclang"* ]]; then
		COMPILE="fcc -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto"
		LINK=""
	elif [[ "$1" = *"gem5"* ]]; then
		COMPILE="fcc -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto"
		LINK=""
	elif [[ "$1" = *"arm"* ]]; then
		COMPILE="armclang -Ofast -ffast-math -march=armv8.2-a+sve -mcpu=a64fx -mtune=a64fx -flto"
		LINK="-lm"
	elif [[ "$1" = *"llvm12"* ]]; then
		COMPILE="clang -Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=full"
		LINK="-fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)"
	fi
	for BMconf in ${BINARYS}; do
		BName="`echo ${BMconf} | cut -d '|' -f1`"
		DSize="`echo ${BMconf} | cut -d '|' -f2`"
		rm -f `basename ${BM}`
		${COMPILE} -I`dirname ${BName}` -I./utilities ./utilities/polybench.c ${BName}.c -D${DSize}_DATASET -DPOLYBENCH_TIME -o ${BName}_${DSize} ${LINK}
		#clang -O2 -fprofile-instr-generate -I`dirname ${BName}` -I./utilities ./utilities/polybench.c ${BName}.c -D${DSize}_DATASET -DPOLYBENCH_TIME -o ${BName}_${DSize} ${LINK} -lm
		#rm -f code-*.profraw code.profdata; LLVM_PROFILE_FILE="code-%p.profraw" numactl -l -C 24 ${BName}_${DSize}; rm ${BName}_${DSize}; llvm-profdata merge -output=code.profdata code-*.profraw
		#${COMPILE} -fprofile-instr-use=code.profdata -I`dirname ${BName}` -I./utilities ./utilities/polybench.c ${BName}.c -D${DSize}_DATASET -DPOLYBENCH_TIME -o ${BName}_${DSize} ${LINK}
	done
	find $ROOTDIR/$BM/ -executable -type f -exec echo {} \; -exec ldd {} \;
	cd $ROOTDIR
fi

