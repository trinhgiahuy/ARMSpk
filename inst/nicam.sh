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
elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
	sleep 0
elif [[ "$1" = *"gem5"* ]]; then
	echo "no point, because this one needs MPI"; exit 1
fi

BM="NICAM"
VERSION="3f758000ffce6ee95a27fb6099f654ecdc5e3add"
if [ ! -f $ROOTDIR/$BM/bin/nhm_driver ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	cd $ROOTDIR/$BM/src
	if [ -z $1 ]; then
		sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./Makefile
		sed -i -e 's/mpiifort/mpif90/' -e 's/mpiicc/mpicc/' -e 's/-shared-intel/-static -static-intel -qopenmp-link=static/g' -e 's/^LFLAGS = /LFLAGS = -static -static-intel -qopenmp-link=static /' ../sysdep/Makedef.Linux64-intel-impi
		export NICAM_SYS=Linux64-intel-impi
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -static#g' ./Makefile
		sed -i -e 's/-O3/-O3 -march=native/g' -e 's/ -pedantic-errors//' -e 's/std=f2003/std=f2008/' -e 's/^LFLAGS = /LFLAGS = -static /' ../sysdep/Makedef.Linux64-gnu-openmpi
		export NICAM_SYS=Linux64-gnu-openmpi
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' ./Makefile
		export NICAM_SYS=FX10
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE; done
	fi
	make ENABLE_OPENMP=1
	cd '../test/case/jablonowski'
	make ENABLE_OPENMP=1
	unset NICAM_SYS
	cd $ROOTDIR
fi

