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

BM="MACSio"
VERSION1="e8bece99bfa5eab9355549bb587ee36aec9d6c67"
VERSION2="30008dc17cd5f787dedd303c51367bb5a8885271"
if [ ! -f $ROOTDIR/$BM/macsio/macsio ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION1}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	if [ ! -f $ROOTDIR/dep/json-cwx/lib/libjson-cwx.a ]; then
		cd $ROOTDIR/dep/json-cwx/
		git checkout -b precision ${VERSION2}
		cd $ROOTDIR/dep/json-cwx/json-cwx
		./autogen.sh
		if [ -z $1 ]; then
			./configure --disable-shared --enable-static --prefix=`pwd`/../ CC=icc CFLAGS="-O2 -ipo -xHost"
		else
			./configure --disable-shared --enable-static --prefix=`pwd`/../ CC=gcc CFLAGS="-O2 -march=native"
		fi
		make
		make install
		cd $ROOTDIR/$BM/
	fi
	if [ ! -f $ROOTDIR/dep/silo-4.10.2/bin/silofile ]; then
		cd $ROOTDIR/dep/
		wget https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz
		tar xzf silo-4.10.2.tar.gz
		cd silo-4.10.2/
		if [ -z $1 ]; then
			./configure --prefix=`pwd` CC=icc CFLAGS="-O2 -ipo -xHost" CXX=icpc CXXFLAGS="-O2 -ipo -xHost" FC=ifort FCFLAGS="-O2 -ipo -xHost" F77=ifort FFLAGS="-O2 -ipo -xHost"
		else
			./configure --prefix=`pwd` CC=gcc CFLAGS="-O2 -march=native" CXX=g++ CXXFLAGS="-O2 -march=native" FC=gfortran FCFLAGS="-O2 -march=native" F77=gfortran FFLAGS="-O2 -march=native"
		fi
		make install
		cd $ROOTDIR/$BM/
	fi
	rm -rf build; mkdir -p build; cd build
	if [ -z $1 ]; then
		cmake -DCMAKE_C_COMPILER=`which mpicc` -DCMAKE_C_FLAGS="-O3 -ipo -xHost -I${ADVISOR_2018_DIR}/include" -DCMAKE_CXX_COMPILER=`which mpicxx` -DCMAKE_CXX_FLAGS="-O3 -ipo -xHost -I${ADVISOR_2018_DIR}/include" -DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
		sed -i -e "s#libjson-cwx.a #libjson-cwx.a -static -static-intel -qopenmp-link=static -L${ADVISOR_2018_DIR}/lib64 -littnotify #" ./macsio/CMakeFiles/macsio.dir/link.txt
	else
		cmake -DCMAKE_C_COMPILER=`which mpicc` -DCMAKE_C_FLAGS="-O3 -march=native" -DCMAKE_CXX_COMPILER=`which mpicxx` -DCMAKE_CXX_FLAGS="-O3 -march=native" -DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
		sed -i -e "s#libjson-cwx.a #libjson-cwx.a -static #" ./macsio/CMakeFiles/macsio.dir/link.txt
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r .. | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	fi
	make
	make install
	cd $ROOTDIR
fi

