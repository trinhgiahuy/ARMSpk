#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

if [[ "$1" = *"fuji"* ]] || [[ "$1" = *"gem5"* ]]; then
	echo "WRN: DOES NOT compile in -Nclang mode"
fi

BM="MVMC"
VERSION="7c58766b180ccb1035e4c220208b64ace3c49cf2"
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
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's# -I${ADVISOR_2018_DIR}/include# -m64 -I${MKLROOT}/include#g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify# -static#g' -e "s#L/usr/local/intel/composer_xe_2013/mkl#L$MKLROOT#g" -e 's/blacs_intelmpi/blacs_openmpi/' -e 's/mkl_intel_thread/mkl_gnu_thread/' -e 's/-lpthread/-lgomp -lpthread -lm -ldl/' -e 's/ -qopt-prefetch=3 -nofor-main//' -e 's/ -vec-report//' ./Makefile_intel
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's/ ifort/ gfortran/' -e 's/-implicitnone/-fimplicit-none/' ./pfapack/Makefile_intel
		sed -i -e 's/-ipo -xHost/-march=native/g' -e 's/ icc/ gcc/' -e 's/-no-ansi-alias//' ./sfmt/Makefile_intel
	elif [[ "$1" = *"fujitrad"* ]] || [[ "$1" = *"fujiclang"* ]]; then
		cp ./Makefile_kei ./Makefile_intel
		cp ./pfapack/Makefile_kei ./pfapack/Makefile_intel
		cp ./sfmt/Makefile_kei ./sfmt/Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' -e 's/-Kfast/-Kfast,openmp,ocl,largepage/g' ./Makefile_intel
	elif [[ "$1" = *"gem5"* ]]; then
		# FJ's -SCALAPACK is hardwired to FJ's MPI, so we need a replacement
		URL="http://www.netlib.org/scalapack/scalapack-2.0.2.tgz"; DEP=$(basename $URL)
		if [ ! -f $ROOTDIR/dep/${DEP} ]; then if ! wget ${URL} -O $ROOTDIR/dep/${DEP}; then echo "ERR: download failed for ${URL}"; exit 1; fi; fi; tar xzf $ROOTDIR/dep/${DEP}
		cd scalapack-2.0.2
		sed -e 's/mpif90/frt/g' -e 's/mpicc/fcc/g' -e "s#-O3#-Kfast,openmp,ocl,nolargepage,nolto -I$ROOTDIR/dep/mpistub/include/mpistub#" -e "s#^LIBS .*#LIBS = -SSL2BLAMP -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi#" SLmake.inc.example > SLmake.inc
		make lib -j $(nproc); make lib
		cd $ROOTDIR/$BM/src
		cp ./Makefile_kei ./Makefile_intel
		cp ./pfapack/Makefile_kei ./pfapack/Makefile_intel
		cp ./sfmt/Makefile_kei ./sfmt/Makefile_intel
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' -e 's/-Kfast/-Kfast,openmp,ocl,nolargepage,nolto/g' -e 's/ = mpi/ = /g' -e "s#^CFLAGS = #CFLAGS = -I$ROOTDIR/dep/mpistub/include/mpistub #g" -e "s#^LIB = #LIB = $ROOTDIR/$BM/src/scalapack-2.0.2/libscalapack.a -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi #g" ./Makefile_intel
	fi
	make intel
	cd $ROOTDIR
fi

