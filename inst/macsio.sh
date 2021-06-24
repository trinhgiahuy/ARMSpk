#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

if [[ "$1" = *"fuji"* ]] || [[ "$1" = *"gem5"* ]]; then
	echo "WRN: DOES NOT compile in -Nclang mode"
	. /vol0004/apps/oss/spack/share/spack/setup-env.sh
	spack load zlib@1.2.11%fj
fi

BM="MACSio"
VERSION1="v1.1" #VERSION1="e8bece99bfa5eab9355549bb587ee36aec9d6c67"
VERSION2="30008dc17cd5f787dedd303c51367bb5a8885271"
if [ ! -f $ROOTDIR/$BM/macsio/macsio ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION1}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	if [ ! -f $ROOTDIR/dep/json-cwx/lib/libjson-cwx.a ]; then
		cd $ROOTDIR/dep/json-cwx/
		if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION2}; fi
		cd $ROOTDIR/dep/json-cwx/json-cwx
		./autogen.sh
		if lscpu | grep 'aarch64' >/dev/null 2>&1 ; then
			rm -f config.guess config.sub
			wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.guess'
			wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.sub'
		fi
		if [[ "$1" = *"intel"* ]]; then
			./configure --disable-shared --enable-static --prefix=`pwd`/../ CC=icc CFLAGS="-O2 -ipo -xHost"
		elif [[ "$1" = *"gnu"* ]]; then
			./configure --disable-shared --enable-static --prefix=`pwd`/../ CC=gcc CFLAGS="-O2 -march=native"
		elif [[ "$1" = *"fujitrad"* ]] || [[ "$1" = *"fujiclang"* ]]; then
			./configure --disable-shared --enable-static --prefix=`pwd`/../ CC=fcc CFLAGS="-O2 -KA64FX,SVE -Kocl,largepage"
		elif [[ "$1" = *"gem5"* ]]; then
			./configure --disable-shared --enable-static --prefix=`pwd`/../ CC=fcc CFLAGS="-O2 -Nclang -mcpu=a64fx+sve -ffj-ocl -ffj-no-largepage -fno-lto"
		fi
		make
		make install
		cd $ROOTDIR/$BM/
	fi
	if [ ! -f $ROOTDIR/dep/silo-4.10.2/bin/silofile ]; then
		cd $ROOTDIR/dep/
		URL="https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz"; DEP=$(basename $URL)
		if [ ! -f $ROOTDIR/dep/${DEP} ]; then if ! wget ${URL} -O $ROOTDIR/dep/${DEP}; then echo "ERR: download failed for ${URL}"; exit 1; fi; fi; tar xzf $ROOTDIR/dep/${DEP}
		cd ./silo-4.10.2/
		if lscpu | grep 'aarch64' >/dev/null 2>&1 ; then
			cd config/; rm -f config.guess config.sub
			wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.guess'
			wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.sub'
			cd -
		fi
		if [[ "$1" = *"intel"* ]]; then
			./configure --prefix=`pwd` CC=icc CFLAGS="-O2 -ipo -xHost" CXX=icpc CXXFLAGS="-O2 -ipo -xHost" FC=ifort FCFLAGS="-O2 -ipo -xHost" F77=ifort FFLAGS="-O2 -ipo -xHost"
		elif [[ "$1" = *"gnu"* ]]; then
			./configure --prefix=`pwd` CC=gcc CFLAGS="-O2 -march=native" CXX=g++ CXXFLAGS="-O2 -march=native" FC=gfortran FCFLAGS="-O2 -march=native" F77=gfortran FFLAGS="-O2 -march=native"
		elif [[ "$1" = *"fujitrad"* ]] || [[ "$1" = *"fujiclang"* ]]; then
			./configure --prefix=`pwd` CC=fcc CFLAGS="-O2 -KA64FX,SVE -Kocl,largepage" \
				CXX=FCC CXXFLAGS="-O2 -KA64FX,SVE -Kocl,largepage" \
				FC=frt FCFLAGS="-O2 -KA64FX,SVE -Kocl,largepage" \
				F77=frt FFLAGS="-O2 -KA64FX,SVE -Kocl,largepage"
		elif [[ "$1" = *"gem5"* ]]; then
			./configure --prefix=`pwd` CC=fcc CFLAGS="-O2 -KA64FX,SVE -Kocl,nolargepage,nolto" \
				CXX=FCC CXXFLAGS="-O2 -KA64FX,SVE -Kocl,nolargepage,nolto" \
				FC=frt FCFLAGS="-O2 -KA64FX,SVE -Kocl,nolargepage,nolto" \
				F77=frt FFLAGS="-O2 -KA64FX,SVE -Kocl,nolargepage,nolto"
		fi
		make install
		cd $ROOTDIR/$BM/
	fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	rm -rf build; mkdir -p build; cd build
	if [[ "$1" = *"intel"* ]]; then
		cmake -DCMAKE_BUILD_TYPE=release \
			-DCMAKE_C_COMPILER=`which mpicc` -DCMAKE_C_FLAGS="-O3 -ipo -xHost -I${ADVISOR_2018_DIR}/include" \
			-DCMAKE_CXX_COMPILER=`which mpicxx` -DCMAKE_CXX_FLAGS="-O3 -ipo -xHost -I${ADVISOR_2018_DIR}/include" \
			-DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
		sed -i -e "s#libjson-cwx.a #libjson-cwx.a -static -static-intel -qopenmp-link=static -L${ADVISOR_2018_DIR}/lib64 -littnotify #" ./macsio/CMakeFiles/macsio.dir/link.txt
	elif [[ "$1" = *"gnu"* ]]; then
		cmake -DCMAKE_BUILD_TYPE=release \
			-DCMAKE_C_COMPILER=`which mpicc` -DCMAKE_C_FLAGS="-O3 -march=native" \
			-DCMAKE_CXX_COMPILER=`which mpicxx` -DCMAKE_CXX_FLAGS="-O3 -march=native" \
			-DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
		sed -i -e "s#libjson-cwx.a #libjson-cwx.a -static #" ./macsio/CMakeFiles/macsio.dir/link.txt
	elif [[ "$1" = *"fujitrad"* ]] || [[ "$1" = *"fujiclang"* ]]; then
		cmake -DCMAKE_BUILD_TYPE=release \
			-DCMAKE_C_COMPILER=`which mpifcc` -DCMAKE_C_FLAGS="-Kfast,openmp,ocl,largepage" \
			-DCMAKE_CXX_COMPILER=`which mpiFCC` -DCMAKE_CXX_FLAGS="-Kfast,openmp,ocl,largepage" \
			-DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
	elif [[ "$1" = *"gem5"* ]]; then
		cmake -DCMAKE_BUILD_TYPE=release \
			-DMPI_C_LIBRARIES="-lmpi" -DMPI_C_INCLUDE_PATH=$ROOTDIR/dep/mpistub/include/mpistub \
			-DMPI_CXX_LIBRARIES="-lmpi" -DMPI_CXX_INCLUDE_PATH=$ROOTDIR/dep/mpistub/include/mpistub \
			-DCMAKE_EXE_LINKER_FLAGS="-Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi" \
			-DCMAKE_C_COMPILER=`which fcc` -DCMAKE_C_FLAGS="-Kfast,ocl,nolargepage,nolto -I$ROOTDIR/dep/mpistub/include/mpistub" \
			-DCMAKE_CXX_COMPILER=`which FCC` -DCMAKE_CXX_FLAGS="-Kfast,ocl,nolargepage,nolto -I$ROOTDIR/dep/mpistub/include/mpistub" \
			-DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
		sed -i -e 's/ -lmpi / /g' -e 's/libsilo.a/libsilo.a -lmpi/g' ./macsio/CMakeFiles/macsio.dir/link.txt
	fi
	make VERBOSE=1
	make install
	cd $ROOTDIR
fi

