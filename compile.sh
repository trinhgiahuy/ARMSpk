#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile AMG -> comes w/ 2 problems
if [ ! -f $ROOTDIR/AMG/test/amg ]; then
	cd $ROOTDIR/AMG/
	make
	cd $ROOTDIR
fi

# compile CANDLE -> comes w/ 7 problems
if [ ! -f $HOME/anaconda2/bin/anaconda ]; then
	cd $ROOTDIR/CANDLE/
	curl -o Anaconda2-4.3.1-Linux-x86_64.sh https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh
	chmod u+x ./Anaconda2-4.3.1-Linux-x86_64.sh
	./Anaconda2-4.3.1-Linux-x86_64.sh -b
	export PATH=$HOME/anaconda2/bin:$PATH
	conda install -y -c conda-forge tensorflow
	conda install -y -c anaconda hdf5=1.8.17
	conda install -y -c anaconda theano
	conda install -y -c conda-forge keras=2
	conda install -y -c conda-forge opencv
	conda install -y -c conda-forge tqdm
	conda update -y -c conda-forge numpy
	cd $ROOTDIR
fi

# compile CoMD
if [ ! -f $ROOTDIR/CoMD/bin/CoMD-openmp-mpi ]; then
	cd $ROOTDIR/CoMD/src-openmp/
	cp Makefile.vanilla Makefile
	make
	cd $ROOTDIR
fi

# compile Laghos
if [ ! -f $ROOTDIR/Laghos/laghos ]; then
	cd $ROOTDIR/Laghos/
	if [ ! -f ./hypre-2.10.0b/src/hypre/lib/libHYPRE.a ]; then
		wget https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/download/hypre-2.10.0b.tar.gz
		tar xzf hypre-2.10.0b.tar.gz
		cd ./hypre-2.10.0b/src
		./configure --disable-fortran
		make -j
		cd $ROOTDIR/Laghos/
	fi
	if [ ! -f ./metis-4.0.3/graphchk ]; then
		wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/metis-4.0.3.tar.gz
		tar xzf metis-4.0.3.tar.gz
		cd ./metis-4.0.3/
		make
		cd $ROOTDIR/Laghos/
	fi
	if [ ! -f $ROOTDIR/dep/mfem/libmfem.a ]; then
		cd $ROOTDIR/dep/mfem/
		git checkout laghos-v1.0
		sed -i -e 's#@MFEM_DIR@/../hypre#@MFEM_DIR@/../../Laghos/hypre#' config/defaults.mk
		sed -i -e 's#@MFEM_DIR@/../metis-4.0$#@MFEM_DIR@/../../Laghos/metis-4.0.3#' config/defaults.mk
		make parallel -j
		cd $ROOTDIR/Laghos/
	fi
	sed -i -e 's#MFEM_DIR = ../mfem$#MFEM_DIR = ../dep/mfem#' makefile
	sed -i -e 's#LAGHOS_LIBS = $(MFEM_LIBS)$#LAGHOS_LIBS = $(MFEM_LIBS) -lirc -lsvml#' makefile
	make
	cd $ROOTDIR
fi

# compile MACSio
if [ ! -f $ROOTDIR/MACSio/macsio/macsio ]; then
	cd $ROOTDIR/MACSio/
	if [ ! -f $ROOTDIR/dep/json-cwx/lib/libjson-cwx.a ]; then
		cd $ROOTDIR/dep/json-cwx/json-cwx
		./autogen.sh
		./configure --prefix=`pwd`/../
		make install
		cd $ROOTDIR/MACSio/
	fi
	if [ ! -f $ROOTDIR/dep/silo-4.10.2/bin/silofile ]; then
		cd $ROOTDIR/dep/
		wget https://wci.llnl.gov/content/assets/docs/simulation/computer-codes/silo/silo-4.10.2/silo-4.10.2.tar.gz
		tar xzf silo-4.10.2.tar.gz
		cd silo-4.10.2/
		./configure --prefix=`pwd`
		make install
		cd $ROOTDIR/MACSio/
	fi
	mkdir -p build; cd build
	cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CC_COMPILER=mpicc -DCMAKE_INSTALL_PREFIX=../ -DWITH_JSON-CWX_PREFIX=../../dep/json-cwx -DWITH_SILO_PREFIX=../../dep/silo-4.10.2 ..
	make
	make install
	cd $ROOTDIR
fi

# compile miniAMR
if [ ! -f $ROOTDIR/MiniAMR/ref/ma.x ]; then
	cd $ROOTDIR/MiniAMR/ref
	sed -i -e 's#= cc#= mpicc#' Makefile
	make
	cd $ROOTDIR
fi

# compile miniFE
if [ ! -f $ROOTDIR/MiniFE/mkl/src/miniFE.x ]; then
	cd $ROOTDIR/MiniFE/mkl/src
	make
	cd $ROOTDIR/MiniFE/openmp-opt/src
	make
	cd $ROOTDIR/MiniFE/openmp-opt-knl/src
	make
	cd $ROOTDIR
fi

# compile miniTri
if [ ! -f $ROOTDIR/MiniTri/miniTri/linearAlgebra/MPI/miniTri.exe ]; then
	cd $ROOTDIR/MiniTri/miniTri/linearAlgebra/MPI
	make
	cd $ROOTDIR/MiniTri/miniTri/linearAlgebra/openmp
	sed -i -e 's/= g++/= icpc/' Makefile
	sed -i -r '/Time to compute miniTri/ s#^(.*)$#//\1#' miniTri.cc
	make
	# get an valid input
	if [ ! -f $ROOTDIR/MiniTri/bcsstk30.mtx ]; then
		cd $ROOTDIR/MiniTri
		wget ftp://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/bcsstruc5/bcsstk30.mtx.gz
		gunzip bcsstk30.mtx.gz
	fi
	cd $ROOTDIR
fi

# compile Nekbone
if [ ! -f $ROOTDIR/Nekbone/test/nek_mgrid/nekbone ]; then
	cd $ROOTDIR/Nekbone/test/nek_mgrid
	sed -i -e 's/lp = 10)/lp = 576)/' -e 's/lelt=100/lelt=1024/' SIZE
	./makenek NotUsedCasename $ROOTDIR/Nekbone/src
	cd $ROOTDIR
fi

# compile SW4lite
if [ ! -f $ROOTDIR/SW4lite/optimize_mp_kiev/sw4lite ]; then
	cd $ROOTDIR/SW4lite
	sed -i -e "s/^HOSTNAME := /HOSTNAME := kiev #/g" Makefile
	sed -i -e "s/quadknl/kiev/g" -e "s#/opt/intel/compilers_and_libraries_2017/linux#`dirname $MKLROOT`#g"  Makefile
	sed -i -e "s/-xmic-avx512/#NOKNL-xmic-avx512/g" Makefile
	make
	sed -i -e "s/kiev/mill/g" Makefile
	sed -i -e "s/#NOKNL-xmic-avx512/-xmic-avx512/g" Makefile
	make
	cd $ROOTDIR
fi

# compile SWFFT
if [ ! -f $ROOTDIR/SWFFT/build.xeon/TestDfft ]; then
	cd $ROOTDIR/SWFFT
	if [ ! -f $ROOTDIR/fftw-3.3.4/bin/fftw-wisdom ]; then
		wget http://fftw.org/fftw-3.3.4.tar.gz
		tar xzf fftw-3.3.4.tar.gz
		cd ./fftw-3.3.4/
		./configure --prefix=`pwd`/../fftw-xmic --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=icc
		make -j CFLAGS="-O3 -xMIC-AVX512 -fp-model fast=2 -no-prec-div -qoverride-limits"
		make install
		make distclean
		./configure --prefix=`pwd`/../fftw-xeon --disable-mpi --enable-openmp --disable-fortran --enable-sse2 --enable-avx CC=icc
		make -j CFLAGS="-O3 -xCORE-AVX2 -fp-model fast=2 -no-prec-div -qoverride-limits"
		make install
		cd $ROOTDIR/SWFFT
	fi
	oldPATH=$PATH
	export PATH=$oldPATH:`pwd`/fftw-xeon/bin
	make -f GNUmakefile.openmp
	cp -r build.openmp build.xeon
	make -f GNUmakefile.openmp clean
	export PATH=$oldPATH:`pwd`/fftw-xmic/bin
	make -f GNUmakefile.openmp
	cp -r build.openmp build.xmic
	make -f GNUmakefile.openmp clean
	export PATH=$oldPATH
	cd $ROOTDIR
fi

# compile XSBench
if [ ! -f $ROOTDIR/XSBench/src/XSBench ]; then
	cd $ROOTDIR/XSBench/src
	sed -i -e 's/^COMPILER.*= gnu/COMPILER = intel/' -e 's/^MPI.* = no/MPI = yes/' -e 's/-openmp/-fopenmp/' Makefile
	make
	cd $ROOTDIR
fi

# compile CCS QCD
if [ ! -f $ROOTDIR/QCD/src/ccs_qcd_solver_bench_class1 ]; then
	cd $ROOTDIR/QCD/src
	sed -i -e 's/-openmp/-fopenmp/' make.ifort.inc
	make MAKE_INC=make.ifort.inc CLASS=1 PX=1 PY=1 PZ=1
	make MAKE_INC=make.ifort.inc CLASS=2 PX=1 PY=1 PZ=1
	make MAKE_INC=make.ifort.inc CLASS=3 PX=1 PY=1 PZ=1
	cd $ROOTDIR
fi

# compile FFVC
if [ ! -f $ROOTDIR/FFVC/bin/ffvc_mini ]; then
	cd $ROOTDIR/FFVC/src
	sed -i -e 's/-openmp/-fopenmp/' make_setting.intel
	rm make_setting; ln -s make_setting.intel make_setting
	make
	cd $ROOTDIR
fi

# compile NICAM
if [ ! -f $ROOTDIR/NICAM/bin/nhm_driver ]; then
	cd $ROOTDIR/NICAM/src
	export NICAM_SYS=Linux64-intel-impi
	make
	cd '../test/case/jablonowski'
	make
	cd $ROOTDIR
fi

# compile mVMC
if [ ! -f $ROOTDIR/MVMC/src/vmc.out ]; then
	cd $ROOTDIR/MVMC/src
	sed -i -e 's/-openmp/-fopenmp/g' -e 's/-opt-prefetch=3/-qopt-prefetch=3/g' -e "s#L/usr/local/intel/composer_xe_2013/mkl#L$MKLROOT#g"  Makefile_intel
	make intel
	cd $ROOTDIR
fi

# compile NGS Analyzer
if [ ! -f $ROOTDIR/NGSAnalyzer/bin/workflow ]; then
	cd $ROOTDIR/NGSAnalyzer
	sed -i -e 's/^N_THREADS=1/N_THREADS=$OMP_NUM_THREADS/' ./workflow/workflow_01.sh
	sed -i -e 's/CC.*=.*gcc/CC=icc/g' -e 's/CXX.*=.*g++/CXX=icpc/g' ./makefile.x86_64_gcc
	sed -i -e 's/CC.*=.*gcc/CC=icc/g' -e 's/CXX.*=.*g++/CXX=icpc/g' ./SNP_indel_caller/Makefile
	# bwa w/ intel breaks
	#sed -i -e 's/CC.*=.*gcc/CC=icc/g' -e 's/CXX.*=.*g++/CXX=icpc/g' ./bwa-0.5.9rc1_kei/bwt_gen/Makefile
	sed -i -e 's/^BWA_OPT.*= CC=$(CC) CXX=$(CXX)/BWA_OPT = CC=gcc CXX=g++/g' ./makefile.x86_64_gcc
	sed -i -e 's/CC.*=.*gcc/CC=icc/g' -e 's/CXX.*=.*g++/CXX=icpc/g' ./samtools-0.1.8_kei/examples/Makefile
	sed -i -e 's/CC.*=.*gcc/CC=icc/g' -e 's/CXX.*=.*g++/CXX=icpc/g' ./samtools-0.1.8_kei/misc/Makefile
	make -f makefile.x86_64_gcc
	# we also need to get an input data set
	if [ ! -f $ROOTDIR/NGSAnalyzer/ngsa_mini_input/reference.fa ]; then
		echo "Creating NGSAnalyzer input may take 1h or more (stay tuned) ..."
		mkdir -p ./ngsa_mini_input/work/../bwa_db; cd ./ngsa_mini_input
		echo "  (downloading and processing reference genome)"
		bash $ROOTDIR/NGSAnalyzer/bin/download_reference.sh ./work
		../bin/samtools faidx reference.fa
		cd ./work
		../../bin/bwa index -a bwtsw -p ../bwa_db/reference.fa ../reference.fa
		bash $ROOTDIR/NGSAnalyzer/bin/download_contig.sh
		echo "  (downloading and processing pseudo-genome data)"
		cd $ROOTDIR/NGSAnalyzer/ngsa_mini_input
		wget http://mt.aics.riken.jp/hpci-miniapp/ngsa-data/ngsa-dummy.tar.gz
		tar zxf ngsa-dummy.tar.gz
		echo "... done"
	fi
	cd $ROOTDIR
fi

# compile MODYLAS (req. license agreement on website)
if [ ! -f $ROOTDIR/MODYLAS/src/modylas_mini ]; then
	mkdir -p $ROOTDIR/MODYLAS
	tar xzf $ROOTDIR/dep/modylas-mini-1.0.0.tar.gz -C $ROOTDIR/MODYLAS --strip-components 1
	cd $ROOTDIR/MODYLAS/src
	sed -i -e 's/-openmp/-fopenmp/' make_setting.intel
	rm make_setting; ln -s make_setting.intel make_setting
	make
	cd $ROOTDIR
fi

# compile NTChem
if [ ! -f $ROOTDIR/NTChem/bin/rimp2.exe ]; then
	cd $ROOTDIR/NTChem
	TOP_DIR=`pwd`
	TYPE=intel
	cp platforms/config_mine.${TYPE} ./config_mine
	sed -i -e 's/-openmp/-fopenmp/g' ./config/linux64_mpif90_omp_intel_proto.makeconfig.in
	./config_mine
	mkdir $ROOTDIR/NTChem/bin
	make CC=mpicc CXX=mpicxx F77C=mpif77 F90C=mpif90
	cd $ROOTDIR
fi

# compile FFB (w/o REVOCAP_Refiner library b/c it has 'registration wall' -> fix later)
if [ ! -f $ROOTDIR/FFB/bin/les3x.mpi ]; then
	cd $ROOTDIR/FFB/src
	if [ ! -f ./metis-5.1.0/bin/graphchk ]; then
		wget http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz
		tar xzf metis-5.1.0.tar.gz
		cd ./metis-5.1.0/
		make config cc=icc prefix=`pwd`
		make install
		cd $ROOTDIR/FFB/src
	fi
	if [ ! -f ./REVOCAP_Refiner-1.1.01/lib/x86_64-linux-intel/libRcapRefiner.a ]; then
		tar xzf $ROOTDIR/dep/REVOCAP_Refiner-1.1.01.tgz
		cd ./REVOCAP_Refiner-1.1.01
		rm ./MakefileConfig.in; ln -s ./MakefileConfig.LinuxIntelCompiler ./MakefileConfig.in
		sed -i '/#include <sstream>/a #include <stdlib.h>' ./RevocapIO/kmbHecmwIO_V3.cpp
		make
		cd $ROOTDIR/FFB/src
	fi
	cp ./make_setting.intel ./make_setting
	sed -i -e 's/^DEFINE += -DNO_METIS/#DEFINE += -DNO_METIS/g' -e "s#\$(HOME)/opt_intel/metis5#$ROOTDIR/FFB/src/metis-5.1.0#g" ./make_setting
	sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt_intel/REVOCAP_Refiner#$ROOTDIR/FFB/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/x86_64-linux-intel #" ./make_setting
	make
	cd $ROOTDIR
fi

