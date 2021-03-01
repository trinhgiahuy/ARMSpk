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
else
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
	spack load openmpi@3.1.6%gcc@8.4.0
	export OMPI_CC=gcc
	export OMPI_CXX=g++
	export OMPI_F77=gfortran
	export OMPI_FC=gfortran
fi

if [ ! -f $ROOTDIR/dep/modylas-mini-1.0.0.tar.gz ]; then
	echo "ERR: Cannot find modylas-mini-1.0.0.tar.gz"
	echo "Please download from: http://hpci-aplfs.aics.riken.jp/fiber/modylas-mini.html and place modylas-mini-1.0.0.tar.gz in ./dep subfolder"
	exit
fi

BM="MODYLAS" # (req. license agreement on website)
if [ ! -f $ROOTDIR/$BM/src/modylas_mini ]; then
	mkdir -p $ROOTDIR/$BM/
	tar xzf $ROOTDIR/dep/modylas-mini-1.0.0.tar.gz -C $ROOTDIR/$BM --strip-components 1
	cd $ROOTDIR/$BM/
	patch -p1 < $ROOTDIR/patches/*1-${BM}*.patch
	cd $ROOTDIR/$BM/src
	if [ -z $1 ]; then
		rm make_setting; ln -s make_setting.intel make_setting
		sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./make_setting
	else
		rm make_setting; ln -s make_setting.gcc make_setting
		sed -i -e 's/-O3/-O3 -march=native\nLIBS += -static/g' ./make_setting
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	fi
	make
	cd $ROOTDIR
fi

