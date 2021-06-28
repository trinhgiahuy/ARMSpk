#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

BM="Laghos"
VERSION="9a074521257434e0b9acff9e59ff10e3e881bc32"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/laghos ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	if [ ! -f ./hypre-2.10.0b/src/hypre/lib/libHYPRE.a ]; then
		URL="https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.10.0b.tar.gz"; DEP=$(basename $URL)
		if [ ! -f $ROOTDIR/dep/${DEP} ]; then if ! wget ${URL} -O $ROOTDIR/dep/${DEP}; then echo "ERR: download failed for ${URL}"; exit 1; fi; fi; tar xzf $ROOTDIR/dep/${DEP}
		cd ./hypre-2.10.0b/src
		if lscpu | grep 'aarch64' >/dev/null 2>&1 ; then
			cd config/; rm -f config.guess config.sub
			wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.guess'
			wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.sub'
			cd -
		fi
		if [[ "$1" = *"intel"* ]]; then
			./configure --disable-fortran -with-openmp \
				CC=mpicc CFLAGS="-O3 -ipo -xHost" \
				CXX=mpicxx CXXFLAGS="-O3 -ipo -xHost" \
				F77=mpif77 FFLAGS="-O3 -ipo -xHost"
		elif [[ "$1" = *"gnu"* ]]; then
			if [ -n "$FJMPI" ]; then
			./configure --disable-fortran -with-openmp \
				CC=mpifcc CFLAGS="-O3 -march=native -fno-lto" \
				CXX=mpiFCC CXXFLAGS="-O3 -march=native -fno-lto" \
				F77=mpifrt FFLAGS="-O3 -march=native -fno-lto"
			else
			./configure --disable-fortran -with-openmp \
				CC=mpicc CFLAGS="-O3 -march=native -flto" \
				CXX=mpicxx CXXFLAGS="-O3 -march=native -flto" \
				F77=mpif77 FFLAGS="-O3 -march=native -flto"
			fi
		elif [[ "$1" = *"fujitrad"* ]]; then
			./configure --disable-fortran -with-openmp \
				CC=mpifcc CFLAGS="-Kfast,openmp,ocl,largepage" \
				CXX=mpiFCC CXXFLAGS="-Kfast,openmp,ocl,largepage" \
				F77=mpifrt FFLAGS="-Kfast,openmp,ocl,largepage,lto"
		elif [[ "$1" = *"fujiclang"* ]]; then
			./configure --disable-fortran -with-openmp \
				CC=mpifcc CFLAGS="-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto" \
				CXX=mpiFCC CXXFLAGS="-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto" \
				F77=mpifrt FFLAGS="-Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto"
		elif [[ "$1" = *"gem5"* ]]; then
			./configure --disable-fortran -with-openmp --without-MPI \
				CC=fcc CFLAGS="-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto" \
				CXX=FCC CXXFLAGS="-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto" \
				F77=frt FFLAGS="-Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,nolargepage,nolto"
		elif [[ "$1" = *"llvm12"* ]]; then
			./configure --disable-fortran -with-openmp \
				CC=mpifcc CFLAGS="-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)" \
				CXX=mpiFCC CXXFLAGS="-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)" \
				F77=mpifrt FFLAGS="-mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto"
		fi
		sed -i -e 's/ -openmp/ -fopenmp/g' ./config/Makefile.config
		make -j
		cd $ROOTDIR/$BM/
	fi
	if [ ! -f ./metis-4.0.3/graphchk ]; then
		URL="http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz"; DEP=$(basename $URL)
		if [ ! -f $ROOTDIR/dep/${DEP} ]; then if ! wget ${URL} -O $ROOTDIR/dep/${DEP}; then echo "ERR: download failed for ${URL}"; exit 1; fi; fi; tar xzf $ROOTDIR/dep/${DEP}
		cd ./metis-4.0.3/
		if [[ "$1" = *"intel"* ]]; then
			sed -i -e 's/CC = cc/CC = icc/g' -e 's/OPTFLAGS = -O2\s*$/OPTFLAGS = -O2 -ipo -xHost/g' ./Makefile.in
		elif [[ "$1" = *"gnu"* ]]; then
			sed -i -e 's/CC = cc/CC = gcc/g' -e 's/OPTFLAGS = -O2\s*$/OPTFLAGS = -O2 -march=native -fno-lto/g' ./Makefile.in
		elif [[ "$1" = *"fujitrad"* ]]; then
			sed -i -e 's/CC = cc/CC = fcc/g' -e 's/OPTFLAGS = -O2\s*$/OPTFLAGS = -O2 -KA64FX,SVE -Kocl,largepage/g' -e 's/^LDOPTIONS =/LDOPTIONS = $(OPTFLAGS) /g' ./Makefile.in
		elif [[ "$1" = *"fujiclang"* ]]; then
			sed -i -e 's/CC = cc/CC = fcc/g' -e 's/OPTFLAGS = -O2\s*$/OPTFLAGS = -O2 -Nclang -mcpu=a64fx+sve -ffj-ocl -ffj-largepage -flto/g' -e 's/^LDOPTIONS =/LDOPTIONS = $(OPTFLAGS) /g' ./Makefile.in
		elif [[ "$1" = *"gem5"* ]]; then
			sed -i -e 's/CC = cc/CC = fcc/g' -e 's/OPTFLAGS = -O2\s*$/OPTFLAGS = -O2 -Nclang -mcpu=a64fx+sve -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's/^LDOPTIONS =/LDOPTIONS = $(OPTFLAGS) /g' ./Makefile.in
		elif [[ "$1" = *"llvm12"* ]]; then
			sed -i -e 's/CC = cc/CC = clang/g' -e "s#OPTFLAGS = -O2\s*$#OPTFLAGS = -O2 -mcpu=a64fx -mtune=a64fx -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./Makefile.in
		fi
		make
		cd $ROOTDIR/$BM/
	fi
	if [[ "$2" = *"rebuild"* ]]; then cd $ROOTDIR/; rm -rf dep/mfem .git/modules/dep/mfem; git submodule update --init dep/mfem; cd -; fi
	if [ ! -f $ROOTDIR/dep/mfem/libmfem.a ]; then
		cd $ROOTDIR/dep/mfem/
		if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"laghos-v1.0"* ]]; then git checkout laghos-v1.0; fi
		git apply --check $ROOTDIR/patches/*1-mfem*.patch
		if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-mfem*.patch; fi
		if [[ "$1" = *"gnu"* ]]; then
			sed -i -e 's/icpc/g++/g' -e "s/-ipo -xHost/-march=native -fno-lto/g" ./config/defaults.mk
		elif [[ "$1" = *"fujitrad"* ]]; then
			sed -i -e 's/icpc/FCC/g' -e 's/-ipo -xHost/-Kfast,openmp,ocl,largepage/g' ./config/defaults.mk
		elif [[ "$1" = *"fujiclang"* ]]; then
			sed -i -e 's/icpc/FCC/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/g' ./config/defaults.mk
		elif [[ "$1" = *"gem5"* ]]; then
			sed -i -e 's/icpc/FCC/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' ./config/defaults.mk
		elif [[ "$1" = *"llvm12"* ]]; then
			sed -i -e 's/icpc/clang++/g' -e "s#-ipo -xHost#-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./config/defaults.mk
		fi
		if [[ "$1" = *"intel"* ]]; then
			make config MFEM_USE_MPI=YES MPICXX=mpicxx MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
		elif [[ "$1" = *"gnu"* ]]; then
			if [ -n "$FJMPI" ]; then
			make config MFEM_USE_MPI=YES MPICXX=mpiFCC MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
			else
			make config MFEM_USE_MPI=YES MPICXX=mpicxx MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
			fi
		elif [[ "$1" = *"fujitrad"* ]]; then
			make config MFEM_USE_MPI=YES MPICXX=mpiFCC MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
		elif [[ "$1" = *"fujiclang"* ]]; then
			make config MFEM_USE_MPI=YES MPICXX=mpiFCC MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
		elif [[ "$1" = *"gem5"* ]]; then
			make config CMAKE_CXX_COMPILER=FCC MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
		elif [[ "$1" = *"llvm12"* ]]; then
			make config MFEM_USE_MPI=YES MPICXX="mpiFCC -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)" MFEM_USE_OPENMP=YES MFEM_THREAD_SAFE=YES MFEM_DEBUG=NO && make -j
		fi
		cd $ROOTDIR/$BM/
	fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/= -L${ADVISOR/= -static -static-intel -qopenmp-link=static -L${ADVISOR/' ./makefile
		make
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -flto ${MAYBESTATIC}#g" -e 's/-ipo -xHost/-march=native -flto/g' -e 's/ -lirc -lsvml//g' ./makefile
		make
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's/-ipo -xHost/-Kfast,openmp,ocl,largepage/g' -e 's/ -lirc -lsvml//g' ./makefile
		make
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage/g' -e 's/ -lirc -lsvml//g' ./makefile
		make
	elif [[ "$1" = *"gem5"* ]]; then
		cd serial/
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's/$(LAGHOS_LIBS) $(LDFLAGS)/$(LDFLAGS) $(LAGHOS_LIBS)/g' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's/ -lirc -lsvml//g' ./makefile
		make
		cd -
		cp serial/laghos .
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" -e "s#-ipo -xHost#-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" -e 's/ -lirc -lsvml//g' ./makefile
		make
	fi
	cd $ROOTDIR
fi

