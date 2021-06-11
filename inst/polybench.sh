#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
if [ -z $1 ]; then
	source $ROOTDIR/conf/intel.cfg
	source $INTEL_PACKAGE intel64 > /dev/null 2>&1
	alias ar=`which xiar`
	alias ld=`which xild`
	export ADVISOR_2018_DIR=${ADVISOR_2019_DIR}
elif [[ "$1" = *"gnu"* ]]; then
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
elif [[ "$1" = *"fuji"* ]]; then
	sleep 0
elif [[ "$1" = *"gem5"* ]]; then
	module load FujitsuCompiler/201903
elif [[ "$1" = *"arm"* ]]; then
	module load /opt/arm/modulefiles/A64FX/RHEL/8/arm-linux-compiler-20.3/armpl/20.3.0
	#module load /opt/arm/modulefiles/A64FX/RHEL/8/gcc-9.3.0/armpl/20.3.0
fi

BM="polybench"
source $ROOTDIR/conf/${BM}.sh
if ! [ $(find $ROOTDIR/$BM/linear-algebra/blas/gemm -executable -type f | wc -l) -ge 1 ]; then
	mkdir -p $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/
	if [ ! -f ./polybench-c-4.2.1-beta.tar.gz ]; then wget https://downloads.sourceforge.net/project/polybench/polybench-c-4.2.1-beta.tar.gz; fi
	tar xzf ./polybench-c-4.2.1-beta.tar.gz -C $ROOTDIR/$BM --strip-components 1
	if patch --dry-run -s -f -p1 < $ROOTDIR/patches/*1-${BM}*.patch; then
		patch -p1 --forward < $ROOTDIR/patches/*1-${BM}*.patch
	fi
	if [ -z $1 ]; then
		#COMPILE="icc -O3 -xHost -static -static-intel -I${ADVISOR_2018_DIR}/include -I./utilities"
		COMPILE="icc -O3 -march=corei7 -static -static-intel -I${ADVISOR_2018_DIR}/include -I./utilities"
		LINK="-L${ADVISOR_2018_DIR}/lib64 -littnotify -L/usr/lib64 -lm"
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e '/^#define STARTSDE/i #ifndef  MYLITTLEMERMAID_H\n#define  MYLITTLEMERMAID_H\nstatic void __attribute__ ((noinline)) ariel_enable() { asm (""); return; }\nstatic void __attribute__ ((noinline)) ariel_fence()  { asm (""); return; }\n#endif //MYLITTLEMERMAID_H' -e 's/__SSC_MARK(0x111);/ariel_enable(); __SSC_MARK(0x111);/g' -e 's/__SSC_MARK(0x222);/__SSC_MARK(0x222); ariel_fence();/g' $FILE; done
	elif [[ "$1" = *"gnu"* ]]; then
		COMPILE="gcc -O3 -march=native -funroll-loops -ffast-math -ftree-vectorize -I./utilities"
		LINK="-static -L/usr/lib64 -lm"
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		COMPILE="fccpx -Nclang -Ofast -ffj-ocl -mllvm -polly -flto -I./utilities"
		LINK=""
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"gem5"* ]]; then
		COMPILE="fccpx -Nclang -Ofast -ffj-no-largepage -ffj-ocl -mllvm -polly -flto -I./utilities"
		LINK=""
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"arm"* ]]; then
		COMPILE="armclang -mcpu=a64fx -march=armv8.2-a+sve -Ofast -ffast-math -flto -I./utilities"
		#COMPILE="gcc -mcpu=native -march=armv8.2-a+sve -O3 -I./utilities"
		LINK="-lm"
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	fi
	for BMconf in ${BINARYS}; do
		BName="`echo ${BMconf} | cut -d '|' -f1`"
		DSize="`echo ${BMconf} | cut -d '|' -f2`"
		rm -f `basename ${BM}`
		${COMPILE} -I`dirname ${BName}` utilities/polybench.c ${BName}.c -D${DSize}_DATASET -DPOLYBENCH_TIME -o ${BName}_${DSize} ${LINK}
		#${COMPILE} -prof-gen -prof-dir $(dirname ${BName}) -I`dirname ${BName}` utilities/polybench.c ${BName}.c -D${DSize}_DATASET -DPOLYBENCH_TIME -o ${BName}_${DSize} ${LINK}
		#numactl -l -C 7 ${BName}_${DSize}; rm ${BName}_${DSize}
		#${COMPILE} -qopt-report -prof-use -prof-dir $(dirname ${BName}) -I`dirname ${BName}` utilities/polybench.c ${BName}.c -D${DSize}_DATASET -DPOLYBENCH_TIME -o ${BName}_${DSize} ${LINK}
	done
	find $ROOTDIR/$BM/ -executable -type f -exec echo {} \; -exec ldd {} \;
	cd $ROOTDIR
fi

