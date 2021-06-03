#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
if [ -z $1 ]; then
	source $ROOTDIR/conf/intel.cfg
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
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
	spack load openmpi@3.1.6%gcc@8.4.0
	export OMPI_CC=gcc
	export OMPI_CXX=g++
	export OMPI_F77=gfortran
	export OMPI_FC=gfortran
elif [[ "$1" = *"fuji"* ]]; then
	module load FujitsuCompiler/202007
	export LD_LIBRARY_PATH=$ROOTDIR/dep/mpistub/lib:$LD_LIBRARY_PATH
else
	echo 'wrong compiler'
	exit 1
fi

BM="FFVC"
VERSION="890a3f9bb3a5cf358504063a1751383b7d46f86d"
if [ ! -f $ROOTDIR/$BM/bin/ffvc_mini ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	if [ -z $1 ]; then
		rm make_setting; ln -s make_setting.intel make_setting
		sed -i -e 's/= -lifport/= -static -static-intel -qopenmp-link=static -lifport/' ./make_setting
	elif [[ "$1" = *"gnu"* ]]; then
		rm make_setting; ln -s make_setting.gcc make_setting
		#XXX: RAGE.... segfaults w/o -g, i give up, this is getting beyond stupid
		sed -i -e 's/= -lgfortran/= /g' -e 's/-O3/-O3 -g -march=native/g' ./make_setting
		sed -i -e 's/\$(LIBS).*/\$(LIBS) -static -l:libgfortran.a -l:libquadmath.a/g' ./FFV/Makefile
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		rm make_setting; ln -s make_setting.fx10 make_setting
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	elif [[ "$1" = *"fuji"* ]]; then
		ln -s $(dirname `which fccpx`)/../lib64/libfj90rt2.a $ROOTDIR/$BM/libfj90rt.a	# fix broken linker
		rm make_setting; ln -s make_setting.fx10 make_setting
		sed -i -e 's/= mpi/= /g' -e "s#^LIBS.*#LIBS = -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi\nCXXFLAGS += -I$ROOTDIR/dep/mpistub/include/mpistub\nF90FLAGS += -I$ROOTDIR/dep/mpistub/include/mpistub#g" ./make_setting
		sed -i -e "s#\$(LIBS).*#\$(LIBS) -L$ROOTDIR/$BM -Bstatic#g" ./FFV/Makefile
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' -e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' -e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE; done
		sed -i -e 's/#define message()/#define fuckthismessage()/' ./FB/mydebug.h
	fi
	make
	cd $ROOTDIR
fi

