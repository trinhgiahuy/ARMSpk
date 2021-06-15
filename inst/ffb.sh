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
	echo "DOES NOT compile in -Nclang mode"
elif [[ "$1" = *"gem5"* ]]; then
	echo "DOES NOT compile in -Nclang mode"
	#module load FujitsuCompiler/202007
	export LD_LIBRARY_PATH=$ROOTDIR/dep/mpistub/lib:$LD_LIBRARY_PATH
else
	echo 'wrong compiler'
	exit 1
fi

if [ ! -f $ROOTDIR/dep/REVOCAP_Refiner-1.1.01.tgz ]; then
	echo "ERR: Cannot find REVOCAP_Refiner-1.1.01.tgz"
	echo "Please download from: http://www.ciss.iis.u-tokyo.ac.jp/dl/index.php?pScdownload_6 and place REVOCAP_Refiner-1.1.01.tgz in ./dep subfolder"
	exit
fi

BM="FFB"
VERSION="e273244b65c7d340cc101ae596a55301359024dd"
if [ ! -f $ROOTDIR/$BM/bin/les3x.mpi ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	if [ ! -f ./metis-5.1.0/bin/graphchk ]; then
		if [ ! -f ./metis-5.1.0.tar.gz ]; then wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz; fi
		tar xzf metis-5.1.0.tar.gz
		cd ./metis-5.1.0/
		if [ -z $1 ]; then
			make config cc=icc prefix=`pwd`
		elif [[ "$1" = *"gnu"* ]]; then
			make config cc=gcc prefix=`pwd`
		elif [[ "$1" = *"fuji"* ]]; then
			export fccpx_ENV="-Kfast"
			make config cc=fccpx prefix=`pwd`
			unset fccpx_ENV
		elif [[ "$1" = *"gem5"* ]]; then
			export fccpx_ENV="-Kfast -Knolargepage"
			make config cc=fccpx prefix=`pwd`
			unset fccpx_ENV
		fi
		make install
		cd $ROOTDIR/$BM/src
	fi
	if [ ! -f ./REVOCAP_Refiner-1.1.01/lib/x86_64-linux-intel/libRcapRefiner.a ]; then
		tar xzf $ROOTDIR/dep/REVOCAP_Refiner-1.1.01.tgz
		cd ./REVOCAP_Refiner-1.1.01
		if [ -z $1 ]; then
			rm ./MakefileConfig.in; ln -s ./MakefileConfig.LinuxIntelCompiler ./MakefileConfig.in
			sed -i -e 's/O2 -w2 -wd1782/O2 -w2 -xHost -wd1782/g' ./MakefileConfig.in
			sed -i '/#include <sstream>/a #include <stdlib.h>' ./RevocapIO/kmbHecmwIO_V3.cpp
		elif [[ "$1" = *"gnu"* ]]; then
			sleep 0
		elif [[ "$1" = *"fuji"* ]]; then
			rm ./MakefileConfig.in; ln -s ./MakefileConfig.Kei ./MakefileConfig.in
		elif [[ "$1" = *"gem5"* ]]; then
			rm ./MakefileConfig.in; ln -s ./MakefileConfig.Kei ./MakefileConfig.in
			sed -i -e 's/-Kfast/-Kfast -Knolargepage/g' ./MakefileConfig.in
			sed -i -e 's/= mpi/= /g' -e "s# -lstd -lm# -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lm#g" ./MakefileConfig.in
		fi
		make
		if [[ "$1" = *"gnu"* ]]; then
			cd ./lib/; ln -s ./x86_64-linux x86_64-linux-intel; cd -
		elif [[ "$1" = *"fuji"* ]] || [[ "$1" = *"gem5"* ]]; then
			cd ./lib/; ln -s ./kei x86_64-linux-intel; cd -
		fi
		ln -s ./Refiner include
		cd $ROOTDIR/$BM/src
	fi
	sed -i -e 's/^DEFINE += -DNO_METIS/#DEFINE += -DNO_METIS/g' -e "s#\$(HOME)/opt_intel/metis5#$ROOTDIR/$BM/src/metis-5.1.0#g" ./make_setting
	if [ -z $1 ]; then
		sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt_intel/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/x86_64-linux-intel #" -e "s/-ipo -xHost -mcmodel=large -shared-intel/-xHost/g" -e 's/LIBS += -L${ADVISOR/LIBS += -static -static-intel -qopenmp-link=static -pthread -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -L${ADVISOR/' ./make_setting
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt_intel/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/x86_64-linux #" -e "s/-ipo -xHost -mcmodel=large -shared-intel/-march=native -static/g" -e 's/LIBS += -L${ADVISOR.*/LIBS += -static/' -e 's# -I${ADVISOR_2018_DIR}/include##g' ./make_setting
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"fuji"* ]]; then
		rm -f ./make_setting; cp ./make_setting.k ./make_setting
		sed -i -e "s#/opt/klocal#$ROOTDIR/$BM/src/metis-5.1.0#g" ./make_setting
		sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/kei #" -e 's/^FLAGS /FFLAGS /g' -e "s#REFINER)/include#REFINER)/Refiner#g" -e 's/-Kvisimpact,ocl/-Kvisimpact,ocl -Kfast/g' ./make_setting
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"gem5"* ]]; then
		rm -f ./make_setting; cp ./make_setting.k ./make_setting
		# crashing in some stupid yaml shit with fujitsu compilers
		sed -i -e 's/-DPROF_MAPROF/-DNO_PROF_MAPROF/g' ./make_setting
		sed -i '/CALL DDINIT/i \      LERR = 0' ./les3x.F
		sed -i -e "s/ = mpi/ = /g" ./make_setting
		sed -i -e "s#/opt/klocal#$ROOTDIR/$BM/src/metis-5.1.0#g" ./make_setting
		sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/kei #" -e 's/^FLAGS /FFLAGS /g' -e "s#REFINER)/include#REFINER)/Refiner#g" -e "s#-Kvisimpact,ocl#-Kvisimpact,ocl -Kfast -Knolargepage -I$ROOTDIR/dep/mpistub/include/mpistub#g" ./make_setting
		sed -i -e "s#^LIBS += -lRcapRefiner#LIBS += -lRcapRefiner -L$ROOTDIR/$BM -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" ./make_setting
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' -e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' -e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE; done
		sed -i -e '/use mpi/d' -e "/implicit none/a \  include 'mpif.h'" ./ma_prof/src/mod_maprof.F90
		sed -i -e '/use mpi/d' -e "/use makemesh/a \  include 'mpif.h'" ./ffb_mini_main.F90
	fi
	make
	cd $ROOTDIR
fi

