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
	echo 'this does NOT work, because gem5 doesnt like socket API'; exit 1
	module load FujitsuCompiler/202007
	export PATH=$HOME/ompi-3.1.6/bin:$PATH
	export LD_LIBRARY_PATH=$HOME/ompi-3.1.6/lib:$LD_LIBRARY_PATH
	export OMPI_CC=fccpx
	export OMPI_CXX=FCCpx
	export OMPI_F77=frtpx
	export OMPI_FC=frtpx
else
	echo 'wrong compiler'
	exit 1
fi

BM="AMG"
VERSION="295de9693eaabf6f7330ac3a35fd9bd4ad030522"
if [ ! -f $ROOTDIR/$BM/test/amg ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	# avx512 on KNL/KNM causes errors in AMG exec
	if [ -z $1 ]; then
		if [[ $HOSTNAME = *"${IKNLHOST}"* ]] || [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
			sed -i -e 's/xHost/xCORE-AVX2/g' ./Makefile.include
		fi
		sed -i -e 's/INCLUDE_LFLAGS = /INCLUDE_LFLAGS = -static -static-intel -qopenmp-link=static /' ./Makefile.include
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/-ipo -xHost/-march=native -static/g' ./Makefile.include
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile.include
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"fuji"* ]]; then
		sed -i -e 's/-ipo -xHost/-Bstatic/g' ./Makefile.include
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile.include
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	fi
	make
	cd $ROOTDIR
fi
