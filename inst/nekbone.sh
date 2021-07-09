#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="Nekbone"
VERSION="d681db43f45f5437e5258b3663b5a92c078cfb57"
source ${ROOTDIR}/conf/nekbone.sh
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/test/nek_mgrid/nekbone ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	subMPI="$( for C in ${TESTCONF}; do echo ${C} | cut -d '|' -f1; done | sort -g -u )"
	for NumMPI in ${subMPI}; do
		cp -r "$ROOTDIR/$BM/test/nek_mgrid" "$ROOTDIR/$BM/test/nek_mgrid_nummpi${NumMPI}"
		cd "$ROOTDIR/$BM/test/nek_mgrid_nummpi${NumMPI}"
		sed -i -e "s/lp = 10)/lp = ${NumMPI})/" -e "s/lelt=100/lelt=$((${ielN}/${NumMPI}))/" SIZE
		for FILE in `/usr/bin/grep 'log2(' -r ../../ | cut -d':' -f1 | sort -u`; do sed -i -e 's/log2(/log2XXX(/g' $FILE; done
		if [[ "$1" = *"intel"* ]]; then
			sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' $ROOTDIR/$BM/src/makefile.template
		elif [[ "$1" = *"gnu"* ]]; then
			if [ -n "$MKLROOT" ]; then
				sed -i -e 's/-ipo -xHost/-march=native -flto/g' -e 's# -I${ADVISOR_2018_DIR}/include# -m64 -I${MKLROOT}/include#g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -flto ${MAYBESTATIC}#g" -e "s# -mkl # -Wl,--start-group \${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \${MKLROOT}/lib/intel64/libmkl_gnu_thread.a \${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -flto ${MAYBESTATIC} #g" -e "s# -mkl-static# -Wl,--start-group \${MKLROOT}/lib/intel64/libmkl_gf_lp64.a \${MKLROOT}/lib/intel64/libmkl_gnu_thread.a \${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -flto ${MAYBESTATIC}#" $ROOTDIR/$BM/src/makefile.template
				sed -i -e 's/-ipo -xHost/-march=native/g' -e 's# -mkl-static# -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl#g' ./makenek
			elif [ -n "$FJBLAS" ]; then
				sed -i -e 's/-ipo -xHost/-march=native -flto/g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$(dirname `which fcc`)/../include#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -flto ${MAYBESTATIC}#g" -e "s# -mkl # -L$(readlink -f $(dirname $(which mpifcc))/../lib64) $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjhpctag.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjlang08.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjomp.o -lfj90rt2 -lssl2mtexsve -lssl2mtsve -lfj90i -lfj90fmt_sve -lfj90f -lfjsrcinfo -lfj90rt -lfjprofcore -lfjprofomp -lelf -lm -flto ${MAYBESTATIC} #g" -e "s# -mkl-static# -L$(readlink -f $(dirname $(which mpifcc))/../lib64) $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjhpctag.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjlang08.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjomp.o -lfj90rt2 -lssl2mtexsve -lssl2mtsve -lfj90i -lfj90fmt_sve -lfj90f -lfjsrcinfo -lfj90rt -lfjprofcore -lfjprofomp -lelf -lm -flto ${MAYBESTATIC} #g" $ROOTDIR/$BM/src/makefile.template
				sed -i -e 's/-ipo -xHost/-march=native/g' -e "s# -mkl-static# -L$(readlink -f $(dirname $(which mpifcc))/../lib64) $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjhpctag.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjlang08.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjomp.o -lfj90rt2 -lssl2mtexsve -lssl2mtsve -lfj90i -lfj90fmt_sve -lfj90f -lfjsrcinfo -lfj90rt -lfjprofcore -lfjprofomp -lelf -lm#g" ./makenek
			fi
			if [ -n "$FJMPI" ]; then sed -i -e 's/^CC=.*/CC=mpifcc/g' -e 's/^F77=.*/F77=mpifrt/g' ./makenek; fi
			sed -i -e 's/-fdefault-real-8/-fdefault-real-8 -fdefault-double-8 -fallow-argument-mismatch/g' $ROOTDIR/$BM/src/makenek.inc
		elif [[ "$1" = *"fujitrad"* ]]; then
			# fancy flags from https://arxiv.org/pdf/2009.11806.pdf
			sed -i -e 's/-ipo -xHost/-CcdRR8 -Cpp -Fixed -Kfast,openmp,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's# -mkl-static##g' -e "s# -mkl# -CcdRR8 -Cpp -Fixed -Kfast,openmp,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation -SSL2BLAMP#g" $ROOTDIR/$BM/src/makefile.template
			sed -i -e 's/gfortran/frt/g' -e 's/-fdefault-real-8 -x f77-cpp-input/-CcdRR8 -Cpp -Fixed -Kfast,openmp,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' $ROOTDIR/$BM/src/makenek.inc
			sed -i -e 's/^CC=.*/CC=mpifcc/g' -e 's/^F77=.*/F77=mpifrt/g' -e 's/-ipo -xHost/-CcdRR8 -Cpp -Fixed -Kfast,openmp,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" ./makenek
		elif [[ "$1" = *"fujiclang"* ]]; then
			sed -i -e 's/-ipo -xHost//g' -e 's/^cFL2.*=/cFL2 = -CcdRR8 -Cpp -Fixed -Nclang -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's# -mkl-static##g' -e "s# -mkl# -CcdRR8 -Cpp -Fixed -Nclang -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation -SSL2BLAMP#g" -e '/^cFL2/a MYCFL2 = -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -fno-lto $(G) $(PPS_C)' -e 's/CC).*-c $(cFL2/CC) -c $(MYCFL2/g' $ROOTDIR/$BM/src/makefile.template
			sed -i -e 's/gfortran/frt/g' -e 's/-fdefault-real-8 -x f77-cpp-input/-CcdRR8 -Cpp -Fixed -Kfast,openmp,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' $ROOTDIR/$BM/src/makenek.inc
			sed -i -e 's/^CC=.*/CC=mpifcc/g' -e 's/^F77=.*/F77=mpifrt/g' -e 's/-ipo -xHost/-CcdRR8 -Cpp -Nclang -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" ./makenek
		elif [[ "$1" = *"gem5"* ]]; then
			if [ "${NumMPI}" -ne 1 ]; then continue; fi
			sed -i -e 's/^IFMPI.*/IFMPI:=false/' -e 's/^cFL2.*=/cFL2 = -CcdRR8 -Cpp -Fixed -Nclang -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,nolargepage,nolto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's# -mkl-static##g' -e "s# -mkl# -CcdRR8 -Cpp -Fixed -Nclang -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,nolargepage,nolto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation -SSL2BLAMP#g" -e '/^cFL2/a MYCFL2 = -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto $(G) $(PPS_C)' -e 's/CC).*-c $(cFL2/CC) -c $(MYCFL2/g' $ROOTDIR/$BM/src/makefile.template
			sed -i -e 's/gfortran/frt/g' -e "s#-fdefault-real-8 -x f77-cpp-input#-CcdRR8 -Cpp -Fixed -Nclang -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,nolargepage,nolto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation#g" $ROOTDIR/$BM/src/makenek.inc
			sed -i -e 's/^CC=.*/CC=fcc/g' -e 's/^F77=.*/F77=frt/g' -e 's/^#IFMPI="f/IFMPI="f/g' -e "s#-ipo -xHost#-CcdRR8 -Cpp -Fixed -Nclang -mcpu=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,nolargepage,nolto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation#g" -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" ./makenek
			sed -i -e '277,287s#/#!#g' $ROOTDIR/$BM/src/makenek.inc
			sed -i -e 's/^CFLAGS+=-DMPI/#CFLAGS+=-DMPI/g' $ROOTDIR/$BM/src/jl/Makefile
		elif [[ "$1" = *"llvm12"* ]]; then
			sed -i -e 's/-ipo -xHost//g' -e 's/^cFL2.*=/cFL2 = -CcdRR8 -Cpp -Fixed -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's# -mkl-static##g' -e "s# -mkl# -CcdRR8 -Cpp -Fixed -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation -SSL2BLAMP#g" -e '/^cFL2/a MYCFL2 = -Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -fno-lto $(G) $(PPS_C)' -e 's/CC).*-c $(cFL2/CC) -c $(MYCFL2/g' $ROOTDIR/$BM/src/makefile.template
			sed -i -e 's/gfortran/frt/g' -e 's/-fdefault-real-8 -x f77-cpp-input/-CcdRR8 -Cpp -Fixed -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' $ROOTDIR/$BM/src/makenek.inc
			sed -i -e 's/^CC=.*/CC=mpifcc/g' -e 's/^F77=.*/F77=mpifrt/g' -e 's/-ipo -xHost/-CcdRR8 -Cpp -Fixed -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kopenmp -Kfast,ocl,largepage,lto -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" ./makenek
		fi
		./makenek NotUsedCasename $ROOTDIR/$BM/src
	done
	cd $ROOTDIR
fi

