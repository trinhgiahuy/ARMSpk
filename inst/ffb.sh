#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
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

if [ ! -f $ROOTDIR/dep/REVOCAP_Refiner-1.1.01.tgz ]; then
	echo "ERR: Cannot find REVOCAP_Refiner-1.1.01.tgz"
	echo "Please download from: http://www.ciss.iis.u-tokyo.ac.jp/dl/index.php?pScdownload_6 and place REVOCAP_Refiner-1.1.01.tgz in ./dep subfolder"
	exit
fi

BM="FFB"
VERSION="e273244b65c7d340cc101ae596a55301359024dd"
if [ ! -f $ROOTDIR/$BM/bin/les3x.mpi ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	if [ ! -f ./metis-5.1.0/bin/graphchk ]; then
		wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
		tar xzf metis-5.1.0.tar.gz
		cd ./metis-5.1.0/
		make config cc=icc prefix=`pwd`
		make install
		cd $ROOTDIR/$BM/src
	fi
	if [ ! -f ./REVOCAP_Refiner-1.1.01/lib/x86_64-linux-intel/libRcapRefiner.a ]; then
		tar xzf $ROOTDIR/dep/REVOCAP_Refiner-1.1.01.tgz
		cd ./REVOCAP_Refiner-1.1.01
		rm ./MakefileConfig.in; ln -s ./MakefileConfig.LinuxIntelCompiler ./MakefileConfig.in
		sed -i -e 's/O2 -w2 -wd1782/O2 -w2 -xHost -wd1782/g' ./MakefileConfig.in
		sed -i '/#include <sstream>/a #include <stdlib.h>' ./RevocapIO/kmbHecmwIO_V3.cpp
		make
		ln -s ./Refiner include
		cd $ROOTDIR/$BM/src
	fi
	sed -i -e 's/^DEFINE += -DNO_METIS/#DEFINE += -DNO_METIS/g' -e "s#\$(HOME)/opt_intel/metis5#$ROOTDIR/$BM/src/metis-5.1.0#g" ./make_setting
	sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt_intel/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/x86_64-linux-intel #" -e "s/-ipo -xHost -mcmodel=large -shared-intel/-xHost/g" -e 's/LIBS += -L${ADVISOR/LIBS += -static -static-intel -qopenmp-link=static -pthread -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -L${ADVISOR/' ./make_setting
	make
	cd $ROOTDIR
fi

