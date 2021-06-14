#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
if [ -z $1 ]; then
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
	source `echo $INTEL_PACKAGE | cut -d'/' -f-3`/mkl/bin/mklvars.sh intel64 > /dev/null 2>&1
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
	spack load openmpi@3.1.6%gcc@8.4.0
	export OMPI_CC=gcc
	export OMPI_CXX=g++
	export OMPI_F77=gfortran
	export OMPI_FC=gfortran
elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
	sleep 0
elif [[ "$1" = *"gem5"* ]]; then
	#echo "ERR: cannot use this one either"; exit 1
	sleep 0; #module load FujitsuCompiler/202007
else
	echo 'wrong compiler'
	exit 1
fi

BM="Nekbone"
VERSION="d681db43f45f5437e5258b3663b5a92c078cfb57"
if [ ! -f $ROOTDIR/$BM/test/nek_mgrid/nekbone ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/test/nek_mgrid
	sed -i -e 's/lp = 10)/lp = 576)/' -e 's/lelt=100/lelt=1024/' SIZE
	if [ -z $1 ]; then
		sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' $ROOTDIR/$BM/src/makefile.template
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's# -I${ADVISOR_2018_DIR}/include# -m64 -I${MKLROOT}/include#g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -static#g' -e 's# -mkl # -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -static #g' -e 's# -mkl-static# -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -static#' $ROOTDIR/$BM/src/makefile.template
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's# -mkl-static# -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl#g' ./makenek
		sed -i -e 's/-fdefault-real-8/-fdefault-real-8 -fdefault-double-8/g' $ROOTDIR/$BM/src/makenek.inc
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r ../../ | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"fuji"* ]]; then
		# fancy flags from https://arxiv.org/pdf/2009.11806.pdf
		sed -i -e 's/-ipo -xHost/-CcdRR8 -Cpp -Fixed -O3 -Kfast -KA64FX -KSVE -KARMV8_3_A -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" $ROOTDIR/$BM/src/makefile.template
		sed -i -e 's/gfortran/frt/g' -e 's/-fdefault-real-8 -x f77-cpp-input/-CcdRR8 -Cpp -Fixed -O3 -Kfast -KA64FX -KSVE -KARMV8_3_A -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' $ROOTDIR/$BM/src/makenek.inc
		sed -i -e 's/^CC=.*/CC=mpifcc/g' -e 's/^F77=.*/F77=mpifrt/g' -e 's/-ipo -xHost/-CcdRR8 -Cpp -Fixed -O3 -Kfast -KA64FX -KSVE -KARMV8_3_A -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation/g' -e 's# -mkl-static##g' -e "s# -mkl# -SSL2BLAMP#g" ./makenek
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r ../../ | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
		for FILE in `/usr/bin/grep 'log2(' -r ../../ | cut -d':' -f1 | sort -u`; do sed -i -e 's/log2(/log2XXX(/g' $FILE; done
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/^IFMPI=.*/IFMPI=false/' -e "s#-ipo -xHost#-CcdRR8 -Cpp -Fixed -O3 -Kfast -KA64FX -KSVE -KARMV8_3_A -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation -Knolargepage -I$ROOTDIR/dep/mpistub/include/mpistub#g" $ROOTDIR/$BM/src/makefile.template
		sed -i -e 's/gfortran/frt/g' -e "s#-fdefault-real-8 -x f77-cpp-input#-CcdRR8 -Cpp -Fixed -O3 -Kfast -KA64FX -KSVE -KARMV8_3_A -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation -Knolargepage -I$ROOTDIR/dep/mpistub/include/mpistub#g" $ROOTDIR/$BM/src/makenek.inc
		sed -i -e 's/^CC=.*/CC=fccpx/g' -e 's/^F77=.*/F77=frtpx/g' -e "s#-ipo -xHost#-CcdRR8 -Cpp -Fixed -O3 -Kfast -KA64FX -KSVE -KARMV8_3_A -Kassume=noshortloop -Kassume=memory_bandwidth -Kassume=notime_saving_compilation -Knolargepage -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e 's# -mkl-static##g' -e "s# -mkl# -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -SSL2BLAMP#g" ./makenek
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r ../../ | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' -e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' -e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE; done
		for FILE in `/usr/bin/grep 'log2(' -r ../../ | cut -d':' -f1 | sort -u`; do sed -i -e 's/log2(/log2XXX(/g' $FILE; done
	fi
	./makenek NotUsedCasename $ROOTDIR/$BM/src
	cd $ROOTDIR
fi

