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
elif [[ "$1" = *"fuji"* ]]; then
	module load FujitsuCompiler/201903
	#export PATH=$ROOTDIR/dep/mpistub/bin:$PATH	#XXX: those won't run because of missing cross-compile
	export LD_LIBRARY_PATH=$ROOTDIR/dep/mpistub/lib:$LD_LIBRARY_PATH
fi

BM="HPL"
if [ ! -f $ROOTDIR/$BM/bin/Linux_Intel64/xhpl ]; then
	mkdir -p $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/
	wget http://www.netlib.org/benchmark/hpl/hpl-2.2.tar.gz
	tar xzf ./hpl-2.2.tar.gz -C $ROOTDIR/$BM --strip-components 1
	patch -p1 < $ROOTDIR/patches/*1-${BM}*.patch
	if [ -z $1 ]; then
		sed -i -e 's/mpiicc/mpicc/' -e 's/ -L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" ./Make.Linux_Intel64
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/mpiicc/mpicc/' -e 's/-ipo -xHost/-march=native -m64 -static/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's/libmkl_intel_thread.a/libmkl_gnu_thread.a/g' -e 's/-lpthread/-lgomp -lpthread -lm/g' -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" -e 's/ -ansi-alias -i-static//g' -e 's/ -nocompchk//g' -e 's/ -mt_mpi//g' ./Make.Linux_Intel64
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		sed -i -e 's/mpiicc/mpifccpx/' -e 's/-ipo -xHost/-Kfast/g' -e 's/ -z noexecstack -z relro -z now//g' -e 's/ -Wall//g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify##g" -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" -e 's/ -ansi-alias -i-static//g' -e 's/ -nocompchk//g' -e 's/ -mt_mpi//g' -e 's/ -I$(LAinc)//g' -e "s# \$(LAlib)# -SSL2BLAMP#g" ./Make.Linux_Intel64
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"fuji"* ]]; then
		SSL2LIB=/opt/FJT/FJTMathlibs_201903/lib64			# new Fj version lacks ssl2
		ln -s $(dirname `which fccpx`)/../lib64/libfj90rt2.a $ROOTDIR/$BM/libfj90rt.a	# fix broken linker
		sed -i -e 's/mpiicc/fccpx/' -e 's/-ipo -xHost//g' -e 's/ -z noexecstack -z relro -z now//g' -e 's/ -Wall//g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -Bstatic#g" -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" -e 's/ -ansi-alias -i-static//g' -e 's/ -nocompchk//g' -e 's/ -mt_mpi//g' -e 's/ -I$(LAinc)//g' -e "s# \$(LAlib)# -L$ROOTDIR/$BM/ -L$SSL2LIB -SSL2BLAMP#g" ./Make.Linux_Intel64
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' -e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' -e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE; done
	fi
	make arch=Linux_Intel64
	cd $ROOTDIR
fi

