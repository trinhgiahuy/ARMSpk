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
else
	source `echo $INTEL_PACKAGE | cut -d'/' -f-3`/mkl/bin/mklvars.sh intel64 > /dev/null 2>&1
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
	spack load openmpi@3.1.6%gcc@8.4.0
	export OMPI_CC=gcc
	export OMPI_CXX=g++
	export OMPI_F77=gfortran
	export OMPI_FC=gfortran
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
	else
		sed -i -e 's/mpiicc/mpicc/' -e 's/-ipo -xHost/-march=native -m64 -static/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's/libmkl_intel_thread.a/libmkl_gnu_thread.a/g' -e 's/-lpthread/-lgomp -lpthread -lm/g' -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" -e 's/ -ansi-alias -i-static//g' -e 's/ -nocompchk//g' -e 's/ -mt_mpi//g' ./Make.Linux_Intel64
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	fi
	make arch=Linux_Intel64
	cd $ROOTDIR
fi

