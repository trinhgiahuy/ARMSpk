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
	echo "does not compile on login node"; exit 1
elif lscpu | grep 'sve' >/dev/null 2>&1 && [[ "$1" = *"fuji"* ]]; then
	. /vol0004/apps/oss/spack/share/spack/setup-env.sh
	spack load zlib@1.2.11%fj@4.3.1
elif [[ "$1" = *"fuji"* ]]; then
	echo "ERR: Ignore this BM, because json dependency will not compile correctly"; exit 1
	module load FujitsuCompiler/202007
	export LD_LIBRARY_PATH=$ROOTDIR/dep/mpistub/lib:$LD_LIBRARY_PATH
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load zlib@1.2.11%gcc@4.8.5
else
	echo 'wrong compiler'
	exit 1
fi

BM="MACSio"
VERSION1="v1.1" #VERSION1="e8bece99bfa5eab9355549bb587ee36aec9d6c67"
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
		elif [[ "$1" = *"gnu"* ]]; then
			./configure --disable-shared --enable-static --prefix=`pwd`/../ CC=gcc CFLAGS="-O2 -march=native"
		elif lscpu | grep 'sve' >/dev/null 2>&1 && [[ "$1" = *"fuji"* ]]; then
			rm -f config.guess config.sub; wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.guess'; wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.sub'
			./configure --disable-shared --enable-static --prefix=`pwd`/../ CC=fcc CFLAGS="-O2 -march=native"
		elif [[ "$1" = *"fuji"* ]]; then
			./configure --host=aarch64-unknown-linux-gnu --build=x84_64-unknown-linux-gnu \
				--disable-shared --enable-static --prefix=`pwd`/../ CC=fccpx CFLAGS="-O2 -march=native"
		fi
		make
		make install
		cd $ROOTDIR/$BM/
	fi
	if [ ! -f $ROOTDIR/dep/silo-4.10.2/bin/silofile ]; then
		cd $ROOTDIR/dep/
		if ! [ -f silo-4.10.2.tar.gz ]; then wget https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz; fi
		tar xzf silo-4.10.2.tar.gz
		cd silo-4.10.2/
		if [ -z $1 ]; then
			./configure --prefix=`pwd` CC=icc CFLAGS="-O2 -ipo -xHost" CXX=icpc CXXFLAGS="-O2 -ipo -xHost" FC=ifort FCFLAGS="-O2 -ipo -xHost" F77=ifort FFLAGS="-O2 -ipo -xHost"
		elif [[ "$1" = *"gnu"* ]]; then
			./configure --prefix=`pwd` CC=gcc CFLAGS="-O2 -march=native" CXX=g++ CXXFLAGS="-O2 -march=native" FC=gfortran FCFLAGS="-O2 -march=native" F77=gfortran FFLAGS="-O2 -march=native"
		elif lscpu | grep 'sve' >/dev/null 2>&1 && [[ "$1" = *"fuji"* ]]; then
			cd config/; rm -f config.guess config.sub; wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.guess'; wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.sub'; cd -
			./configure --prefix=`pwd` CC=fcc CFLAGS="-O2 -march=native" CXX=FCC CXXFLAGS="-O2 -march=native" FC=frt FCFLAGS="-O2 -march=native" F77=frt FFLAGS="-O2 -march=native"
		elif [[ "$1" = *"fuji"* ]]; then
			./configure --host=aarch64-unknown-linux-gnu --build=x84_64-unknown-linux-gnu \
				--with-zlib="`spack find -p | /bin/grep zlib | cut -d' ' -f2- | head -1 | tr -d ' '`/include,`spack find -p | /bin/grep zlib | cut -d' ' -f2- | head -1 | tr -d ' '`/lib" \
				--prefix=`pwd` CC=fccpx CFLAGS="-O2 -march=native" CXX=FCCpx CXXFLAGS="-O2 -march=native" FC=frtpx FCFLAGS="-O2 -march=native" F77=frtpx FFLAGS="-O2 -march=native"
			make install	#runs after 2. call... facepalm
		fi
		make install
		cd $ROOTDIR/$BM/
	fi
	rm -rf build; mkdir -p build; cd build
	if [ -z $1 ]; then
		cmake -DCMAKE_C_COMPILER=`which mpicc` -DCMAKE_C_FLAGS="-O3 -ipo -xHost -I${ADVISOR_2018_DIR}/include" -DCMAKE_CXX_COMPILER=`which mpicxx` -DCMAKE_CXX_FLAGS="-O3 -ipo -xHost -I${ADVISOR_2018_DIR}/include" -DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
		sed -i -e "s#libjson-cwx.a #libjson-cwx.a -static -static-intel -qopenmp-link=static -L${ADVISOR_2018_DIR}/lib64 -littnotify #" ./macsio/CMakeFiles/macsio.dir/link.txt
	elif [[ "$1" = *"gnu"* ]]; then
		cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_C_COMPILER=`which mpicc` -DCMAKE_C_FLAGS="-O3 -march=native" -DCMAKE_CXX_COMPILER=`which mpicxx` -DCMAKE_CXX_FLAGS="-O3 -march=native" -DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
		sed -i -e "s#libjson-cwx.a #libjson-cwx.a -static #" ./macsio/CMakeFiles/macsio.dir/link.txt
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r .. | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif lscpu | grep 'sve' >/dev/null 2>&1 && [[ "$1" = *"fuji"* ]]; then
		cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_C_COMPILER=`which mpifcc` -DCMAKE_C_FLAGS="-O3 -march=native" -DCMAKE_CXX_COMPILER=`which mpiFCC` -DCMAKE_CXX_FLAGS="-O3 -march=native" -DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r .. | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE; done
	elif [[ "$1" = *"fuji"* ]]; then
		for FILE in `/usr/bin/grep 'include.*ittnotify' -r .. | cut -d':' -f1 | sort -u`; do sed -i -e 's/.*include.*ittnotify\.h.*/#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' -e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' -e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE; done
		cmake -DMPI_C_LIBRARIES="-lmpi" -DMPI_C_INCLUDE_PATH=$ROOTDIR/dep/mpistub/include/mpistub -DMPI_CXX_LIBRARIES="-lmpi" -DMPI_CXX_INCLUDE_PATH=$ROOTDIR/dep/mpistub/include/mpistub -DCMAKE_BUILD_TYPE=release -DCMAKE_C_COMPILER=`which fccpx` -DCMAKE_C_FLAGS="-O3 -march=native -I$ROOTDIR/dep/mpistub/include/mpistub -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -Bstatic" -DCMAKE_CXX_COMPILER=`which FCCpx` -DCMAKE_CXX_FLAGS="-O3 -march=native -I$ROOTDIR/dep/mpistub/include/mpistub -Wl,-rpath -Wl,$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -Bstatic" -DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
	fi
	make VERBOSE=1
	make install
	cd $ROOTDIR
fi

