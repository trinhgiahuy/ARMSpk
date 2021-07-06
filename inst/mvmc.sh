#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

if [[ "$1" = *"fuji"* ]] || [[ "$1" = *"gem5"* ]]; then
	echo "WRN: DOES NOT compile in -Nclang mode"
fi

BM="MVMC"
VERSION="7c58766b180ccb1035e4c220208b64ace3c49cf2"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/src/vmc.out ] || [ "x`ls -s $ROOTDIR/$BM/src/vmc.out | awk '{print $1}'`" = "x0" ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/src
	for x in `/bin/grep -r 'init_by_array(' . | cut -d':' -f1 | sort -u`; do sed -i -e 's/init_by_array(/xxxinit_by_array(/' $x; done
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/blacs_intelmpi/blacs_openmpi/' -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' -e "s#L/usr/local/intel/composer_xe_2013/mkl#L$MKLROOT#g" ./Makefile_intel
	elif [[ "$1" = *"gnu"* ]]; then
		if [ -n "$FJMPI" ]; then sed -i -e 's/ mpif90/ mpifrt/g' -e 's/ mpicc/ mpifcc/g' ./Makefile_intel; fi
		if [ -n "$MKLROOT" ]; then
			sed -i -e 's/-ipo -xHost/-march=native -flto/g' -e 's# -I${ADVISOR_2018_DIR}/include# -m64 -I${MKLROOT}/include#g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -flto ${MAYBESTATIC}#g" -e "s#L/usr/local/intel/composer_xe_2013/mkl#L$MKLROOT#g" -e 's/blacs_intelmpi/blacs_openmpi/' -e 's/mkl_intel_thread/mkl_gnu_thread/' -e 's/-lpthread/-lgomp -lpthread -lm -ldl/' -e 's/ -qopt-prefetch=3 -nofor-main//' -e 's/ -vec-report//' ./Makefile_intel
		elif [ -n "$FJBLAS" ]; then
			sed -i -e 's/-ipo -xHost/-march=native -flto/g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$(dirname `which fcc`)/../include#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -flto ${MAYBESTATIC}#g" -e "s#\$(MKL)#-L$(readlink -f $(dirname $(which mpifcc))/../lib64) -lscalapacksve -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi_usempif08 $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjhpctag.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjlang08.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjomp.o -lfj90rt2 -lssl2mtexsve -lssl2mtsve -lfj90i -lfj90fmt_sve -lfj90f -lfjsrcinfo -lfj90rt -lfjprofcore -lfjprofomp -lelf -flto ${MAYBESTATIC}#g" -e 's/ -qopt-prefetch=3 -nofor-main//' -e 's/ -vec-report//' ./Makefile_intel
			sed -i -e 's/ -DHAVE_SSE2//g' ./sfmt/Makefile_intel
		fi
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's/ ifort/ gfortran/' -e 's/-implicitnone/-fimplicit-none/' ./pfapack/Makefile_intel
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's/ icc/ gcc/' -e 's/-no-ansi-alias//' ./sfmt/Makefile_intel
	elif [[ "$1" = *"fujitrad"* ]]; then
		cp ./Makefile_kei ./Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./Makefile_intel
		sed -i -e 's/Makefile_kei/Makefile_intel/g' -e 's/-Kfast/-Kfast,openmp,ocl,largepage/g' ./Makefile_intel
		cp ./pfapack/Makefile_kei ./pfapack/Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./pfapack/Makefile_intel
		sed -i -e 's/-Kfast/-Kfast,openmp,ocl,largepage/g' ./pfapack/Makefile_intel
		cp ./sfmt/Makefile_kei ./sfmt/Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./sfmt/Makefile_intel
		sed -i -e 's/-Kfast/-Kfast,openmp,ocl,largepage/g' ./sfmt/Makefile_intel
	elif [[ "$1" = *"fujiclang"* ]]; then
		cp ./Makefile_kei ./Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./Makefile_intel
		sed -i -e 's/Makefile_kei/Makefile_intel/g' -e 's/-Kfast.*/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -fno-lto/g' -e '/^CFLAGS =/a FFLAGS = -Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto' -e 's/CFLAGS) $(LIB/FFLAGS) $(LIB/g' ./Makefile_intel
		cp ./pfapack/Makefile_kei ./pfapack/Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./pfapack/Makefile_intel
		sed -i -e 's/-Kfast/-Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' ./pfapack/Makefile_intel
		cp ./sfmt/Makefile_kei ./sfmt/Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./sfmt/Makefile_intel
		sed -i -e 's/-Kfast.*/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -fno-lto/g' ./sfmt/Makefile_intel
	elif [[ "$1" = *"gem5"* ]]; then
		# FJ's -SCALAPACK is hardwired to FJ's MPI, so we need a replacement
		URL="http://www.netlib.org/scalapack/scalapack-2.0.2.tgz"; DEP=$(basename $URL)
		if [ ! -f $ROOTDIR/dep/${DEP} ]; then if ! wget ${URL} -O $ROOTDIR/dep/${DEP}; then echo "ERR: download failed for ${URL}"; exit 1; fi; fi; tar xzf $ROOTDIR/dep/${DEP}
		cd scalapack-2.0.2
		sed -e 's/mpif90/frt/g' -e 's/mpicc/fcc/g' -e "s#-O3#-Nclang -Kfast,openmp,ocl,nolargepage -I$ROOTDIR/dep/mpistub/include/mpistub#" -e "s#^LIBS .*#LIBS = -SSL2BLAMP -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi#" SLmake.inc.example > SLmake.inc
		make lib -j $(nproc); make lib
		cd $ROOTDIR/$BM/src
		cp ./Makefile_kei ./Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./Makefile_intel
		sed -i -e 's/Makefile_kei/Makefile_intel/g' -e 's/-Kfast.*/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's/ = mpi/ = /g' -e "s#^CFLAGS = #CFLAGS = -I$ROOTDIR/dep/mpistub/include/mpistub #g" -e "s#^LIB = #LIB = $ROOTDIR/$BM/src/scalapack-2.0.2/libscalapack.a -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi #g" -e '/^CFLAGS =/a FFLAGS = -Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,nolargepage,nolto' -e 's/CFLAGS) $(LIB/FFLAGS) $(LIB/g' ./Makefile_intel
		cp ./pfapack/Makefile_kei ./pfapack/Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./pfapack/Makefile_intel
		sed -i -e 's/-Kfast/-Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,nolargepage,nolto/g' ./pfapack/Makefile_intel
		cp ./sfmt/Makefile_kei ./sfmt/Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./sfmt/Makefile_intel
		sed -i -e 's/-Kfast.*/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' ./sfmt/Makefile_intel
	elif [[ "$1" = *"llvm12"* ]]; then
		cp ./Makefile_kei ./Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./Makefile_intel
		sed -i -e 's/Makefile_kei/Makefile_intel/g' -e 's/-Kfast.*/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -fno-lto/g' -e 's/^REPORT.*/REPORT=/g' -e '/^CFLAGS =/a FFLAGS = -Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto' -e 's/CFLAGS) $(LIB/FFLAGS) $(LIB/g' ./Makefile_intel
		cp ./pfapack/Makefile_kei ./pfapack/Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./pfapack/Makefile_intel
		sed -i -e 's/-Kfast/-Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' ./pfapack/Makefile_intel
		cp ./sfmt/Makefile_kei ./sfmt/Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./sfmt/Makefile_intel
		sed -i -e 's/-Kfast.*/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=full/g' ./sfmt/Makefile_intel
	fi
	make intel
	cd $ROOTDIR
fi

