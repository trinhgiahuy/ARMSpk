#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="NTChem"
VERSION="fcafcc4fec195d8a81c19affd1a3b83f7bab4285"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/bin/rimp2.exe ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	if [[ "$1" = *"intel"* ]]; then
		TYPE=intel
		sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./src/mp2/GNUmakefile
	elif [[ "$1" = *"gnu"* ]]; then
		TYPE=gfortran
		if [ -n "$FJMPI" ]; then sed -i -e 's/ mpif77/ mpifrt/g' -e 's/ mpif90/ mpifrt/g' ./config/linux64_mpif90_omp_gfortran.makeconfig.in; fi
		sed -i -e 's/ -fc=gfortran/ -flto/' ./config/linux64_mpif90_omp_gfortran.makeconfig.in
		if [ -n "$MKLROOT" ]; then
			sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -Wl,--start-group \${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \${MKLROOT}/lib/intel64/libmkl_gnu_thread.a \${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -flto ${MAYBESTATIC}#g" ./src/mp2/GNUmakefile
			sed -i -e 's# -I${ADVISOR_2018_DIR}/include# -m64 -I${MKLROOT}/include#g' ./src/util_lib/GNUmakefile
		elif [ -n "$FJBLAS" ]; then
			sed -i -e 's# -m64##g' ./config/linux64_mpif90_omp_gfortran.makeconfig.in
			sed -i -e 's# -m64##g' ./src/mp2/GNUmakefile
			sed -i -e 's# -m64##g' ./src/util_lib/GNUmakefile
			sed -i -e "s# -I\${ADVISOR_2018_DIR}/include# -I$(dirname `which fcc`)/../include#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -L$(readlink -f $(dirname $(which mpifcc))/../lib64) $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjhpctag.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjlang08.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjomp.o -lfj90rt2 -lssl2mtexsve -lssl2mtsve -lfj90i -lfj90fmt_sve -lfj90f -lfjsrcinfo -lfj90rt -lfjprofcore -lfjprofomp -lelf -flto ${MAYBESTATIC}#g" ./src/mp2/GNUmakefile
			sed -i -e "s# -I\${ADVISOR_2018_DIR}/include# -I$(dirname `which fcc`)/../include#g" ./src/util_lib/GNUmakefile
		fi
		sed -e 's/-SSL2BLAMP//' -e 's/mpifrtpx_omp_k_fx10/mpif90_omp_gfortran/' platforms/config_mine.K > config_mine
		sed -i -e 's/int /long /' ./src/util_lib/ssc.c
		sed -i -e 's/integer(C_INT)/integer(kind=8)/' ./src/mp2/mp2_main_mpiomp.f90
		sed -i -e '/CHARACTER(/d' -e '/INTEGER :: MLeng/a\      CHARACTER(LEN=MLeng) :: Chara' ./src/util_lib/util_transchar.f90
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./config/linux64_mpifrtpx_omp_k_fx10.makeconfig.in
		sed -i -e "s# -lmpi_f90 -lmpi_f77##g" -e 's/-Kfast/-Kfast,openmp,ocl,largepage,lto/g' ./config/linux64_mpifrtpx_omp_k_fx10.makeconfig.in
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./src/mp2/GNUmakefile
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include##g" ./src/util_lib/GNUmakefile
		cp platforms/config_mine.K config_mine
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./config/linux64_mpifrtpx_omp_k_fx10.makeconfig.in
		sed -i -e "s# -lmpi_f90 -lmpi_f77##g" -e 's/INC) -Kfast/INC) -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -fno-lto/g' -e 's/INCMOD) -Kfast/INCMOD) -Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' ./config/linux64_mpifrtpx_omp_k_fx10.makeconfig.in
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./src/mp2/GNUmakefile
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include##g" ./src/util_lib/GNUmakefile
		cp platforms/config_mine.K config_mine
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./config/linux64_mpifrtpx_omp_k_fx10.makeconfig.in
		sed -i -e "s# -lmpi_f90 -lmpi_f77# -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" -e 's/INC) -Kfast/INC) -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's/INCMOD) -Kfast/INCMOD) -Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,nolargepage,nolto/g' -e "s#-Kfast#-I$ROOTDIR/dep/mpistub/include/mpistub -Kfast#g" ./config/linux64_mpifrtpx_omp_k_fx10.makeconfig.in
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./src/mp2/GNUmakefile
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" ./src/util_lib/GNUmakefile
		cp platforms/config_mine.K config_mine
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./config/linux64_mpifrtpx_omp_k_fx10.makeconfig.in
		sed -i -e "s# -lmpi_f90 -lmpi_f77##g" -e 's/INC) -Kfast/INC) -Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -fno-lto/g' -e 's/INCMOD) -Kfast/INCMOD) -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' ./config/linux64_mpifrtpx_omp_k_fx10.makeconfig.in
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./src/mp2/GNUmakefile
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include##g" ./src/util_lib/GNUmakefile
		cp platforms/config_mine.K config_mine
	fi
	TOP_DIR=`pwd`
	./config_mine
	mkdir -p $ROOTDIR/$BM/bin
	sed -i -e 's/^#include/include/g' ./GNUmakefile
	if [[ "$1" = *"intel"* ]]; then
		make CC=mpicc CXX=mpicxx F77C=mpif77 F90C=mpif90
	elif [[ "$1" = *"gnu"* ]]; then
		if [ -n "$FJMPI" ]; then make CC=mpifcc CXX=mpiFCC F77C=mpifrt F90C=mpifrt;
		else                     make CC=mpicc CXX=mpicxx F77C=mpif77 F90C=mpif90; fi
	elif [[ "$1" = *"fujitrad"* ]] || [[ "$1" = *"fujiclang"* ]]; then
		make CC=mpifcc CXX=mpiFCC F77C=mpifrt F90C=mpifrt
	elif [[ "$1" = *"gem5"* ]]; then
		make CC=fcc CXX=FCC F77C=frt F90C=frt LD=FCC
	elif [[ "$1" = *"llvm12"* ]]; then
		make CC=mpifcc CXX=mpiFCC F77C=mpifrt F90C=mpifrt LD=mpifrt
	fi
	unset TOP_DIR
	unset TYPE
	cd $ROOTDIR
fi

