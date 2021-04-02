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
elif [[ "$1" = *"fuji"* ]]; then
	echo "ERR: cannot use this one either, because peach's outdated SSL2 has no dpotrf_ and other fn"; exit 1
	module load FujitsuCompiler/202007
else
	echo 'wrong compiler'
	exit 1
fi

BM="NTChem"
VERSION="fcafcc4fec195d8a81c19affd1a3b83f7bab4285"
if [ ! -f $ROOTDIR/$BM/bin/rimp2.exe ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	if [ -z $1 ]; then
		TYPE=intel
		sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./src/mp2/GNUmakefile
	elif [[ "$1" = *"gnu"* ]]; then
		TYPE=gfortran
		sed -i -e 's/ -fc=gfortran//' ./config/linux64_mpif90_omp_gfortran.makeconfig.in
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -static#g' ./src/mp2/GNUmakefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include# -m64 -I${MKLROOT}/include#g' ./src/util_lib/GNUmakefile
		sed -e 's/-SSL2BLAMP//' -e 's/mpifrtpx_omp_k_fx10/mpif90_omp_gfortran/' platforms/config_mine.K > config_mine
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e 's/int /long /' $FILE; done
		sed -i -e 's/integer(C_INT)/integer(kind=8)/' ./src/mp2/mp2_main_mpiomp.f90
		sed -i -e '/CHARACTER(/d' -e '/INTEGER :: MLeng/a\      CHARACTER(LEN=MLeng) :: Chara' ./src/util_lib/util_transchar.f90
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		sed -i -e "s# -lmpi_f90 -lmpi_f77##g" ./config/linux64_mpifrtpx_omp_k_fx10.makeconfig.in
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include##g" -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./src/mp2/GNUmakefile
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include##g" ./src/util_lib/GNUmakefile
		cp platforms/config_mine.K config_mine
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e 's/int /long /' $FILE; done
	elif [[ "$1" = *"fuji"* ]]; then
		SSL2LIB=/opt/FJT/FJTMathlibs_201903/lib64					# new Fj version lacks ssl2
		ln -s $(dirname `which fccpx`)/../lib64/libfj90rt2.a $ROOTDIR/$BM/libfj90rt.a	# fix broken linker
		sed -i -e "s# -lmpi_f90 -lmpi_f77# -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" ./config/linux64_mpifrtpx_omp_k_fx10.makeconfig.in
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -Bstatic#g' ./src/mp2/GNUmakefile
		sed -i -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" ./src/util_lib/GNUmakefile
		sed -e "s#-SSL2BLAMP#\"-L$ROOTDIR/$BM/ -L$SSL2LIB -SSL2BLAMP\"#g" platforms/config_mine.K > config_mine
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' -e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' -e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE; done
	fi
	TOP_DIR=`pwd`
	./config_mine
	mkdir -p $ROOTDIR/$BM/bin
	if [ -z $1 ] || [[ "$1" = *"gnu"* ]]; then
		make CC=mpicc CXX=mpicxx F77C=mpif77 F90C=mpif90
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		make CC=mpifccpx CXX=mpiFCCpx F77C=mpifrtpx F90C=mpifrtpx
	else
		make CC=fccpx CXX=FCCpx F77C=frtpx F90C=frtpx LD=FCCpx F90FLAGS="-I$ROOTDIR/dep/mpistub/include/mpistub"
	fi
	unset TOP_DIR
	unset TYPE
	cd $ROOTDIR
fi

