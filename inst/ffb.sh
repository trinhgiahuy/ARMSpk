#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

if [[ "$1" = *"fujiclang"* ]] || [[ "$1" = *"llvm12"* ]]; then
	echo "WRN: REVOCAP_Refiner-1.1.01 DOES NOT compile in clang mode"
fi

if [ ! -f $ROOTDIR/dep/REVOCAP_Refiner-1.1.01.tgz ]; then
	echo "ERR: Cannot find REVOCAP_Refiner-1.1.01.tgz"
	echo "Please download from: http://www.ciss.iis.u-tokyo.ac.jp/dl/index.php?pScdownload_6 and place REVOCAP_Refiner-1.1.01.tgz in ./dep subfolder"
	exit
fi

BM="FFB"
VERSION="e273244b65c7d340cc101ae596a55301359024dd"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/bin/les3x.mpi ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/src
	if [ ! -f ./metis-5.1.0/bin/graphchk ]; then
		URL="http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz"; DEP=$(basename $URL)
		if [ ! -f $ROOTDIR/dep/${DEP} ]; then if ! wget ${URL} -O $ROOTDIR/dep/${DEP}; then echo "ERR: download failed for ${URL}"; exit 1; fi; fi; tar xzf $ROOTDIR/dep/${DEP}
		cd ./metis-5.1.0/
		if [[ "$1" = *"intel"* ]]; then
			make config cc=icc prefix=`pwd`
		elif [[ "$1" = *"gnu"* ]]; then
			make config cc=gcc prefix=`pwd`
		elif [[ "$1" = *"fujitrad"* ]]; then
			export fcc_ENV="-Kocl,largepage"
			make config cc=fcc prefix=`pwd`
			unset fcc_ENV
		elif [[ "$1" = *"fujiclang"* ]]; then
			export fcc_clang_ENV="-Nclang -mcpu=a64fx+sve -ffj-ocl -ffj-largepage -flto"
			make config cc=fcc prefix=`pwd`
			unset fcc_clang_ENV
		elif [[ "$1" = *"gem5"* ]]; then
			#export fcc_clang_ENV="-Nclang -mcpu=a64fx+sve -ffj-ocl -ffj-no-largepage -fno-lto"
			export fcc_ENV="-Kocl,nolargepage"
			make config cc=fcc prefix=`pwd`
			unset fcc_clang_ENV
		elif [[ "$1" = *"llvm12"* ]]; then
			make config cc=clang prefix=`pwd`
		fi
		make install
		cd $ROOTDIR/$BM/src
	fi
	if [ ! -f ./REVOCAP_Refiner-1.1.01/lib/x86_64-linux-intel/libRcapRefiner.a ]; then
		tar xzf $ROOTDIR/dep/REVOCAP_Refiner-1.1.01.tgz
		cd ./REVOCAP_Refiner-1.1.01
		if [[ "$1" = *"intel"* ]]; then
			rm ./MakefileConfig.in; ln -s ./MakefileConfig.LinuxIntelCompiler ./MakefileConfig.in
			sed -i -e 's/O2 -w2 -wd1782/O2 -w2 -xHost -wd1782/g' ./MakefileConfig.in
			sed -i '/#include <sstream>/a #include <stdlib.h>' ./RevocapIO/kmbHecmwIO_V3.cpp
		elif [[ "$1" = *"gnu"* ]]; then
			sed -i -e 's/ -Wall/ -fno-lto -Wall/g' ./MakefileConfig.in
		elif [[ "$1" = *"fujitrad"* ]]; then
			rm ./MakefileConfig.in; ln -s ./MakefileConfig.Kei ./MakefileConfig.in
			sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./MakefileConfig.in
			sed -i -e 's/-Kfast/-Kfast,openmp,ocl,largepage/g' ./MakefileConfig.in
		elif [[ "$1" = *"fujiclang"* ]]; then
			rm ./MakefileConfig.in; ln -s ./MakefileConfig.Kei ./MakefileConfig.in
			sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./MakefileConfig.in
			sed -i -e 's/-Kfast/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/g' -e 's/ -Kauto//g' ./MakefileConfig.in
		elif [[ "$1" = *"gem5"* ]]; then
			rm ./MakefileConfig.in; ln -s ./MakefileConfig.Kei ./MakefileConfig.in
			sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./MakefileConfig.in
			#sed -i -e 's/-Kfast/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's/ -Kauto//g' ./MakefileConfig.in
			sed -i -e 's/-Kfast/-Kfast,openmp,ocl,nolargepage/g' ./MakefileConfig.in
			sed -i -e 's/= mpi/= /g' -e "s# -lstd -lm# -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lm#g" ./MakefileConfig.in
		elif [[ "$1" = *"llvm12"* ]]; then
			rm ./MakefileConfig.in; ln -s ./MakefileConfig.Kei ./MakefileConfig.in
			sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./MakefileConfig.in
			sed -i -e "s#-Kfast.*#-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=full -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./MakefileConfig.in
		fi
		make
		if [[ "$1" = *"gnu"* ]]; then
			cd ./lib/; ln -s ./x86_64-linux x86_64-linux-intel; cd -
		elif [[ "$1" = *"fujitrad"* ]] || [[ "$1" = *"fujiclang"* ]] || [[ "$1" = *"gem5"* ]] || [[ "$1" = *"llvm12"* ]]; then
			cd ./lib/; ln -s ./kei x86_64-linux-intel; cd -
		fi
		ln -s ./Refiner include
		cd $ROOTDIR/$BM/src
	fi
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/^DEFINE += -DNO_METIS/#DEFINE += -DNO_METIS/g' -e "s#\$(HOME)/opt_intel/metis5#$ROOTDIR/$BM/src/metis-5.1.0#g" ./make_setting
		sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt_intel/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/x86_64-linux-intel #" -e "s/-ipo -xHost -mcmodel=large -shared-intel/-xHost/g" -e 's/LIBS += -L${ADVISOR/LIBS += -static -static-intel -qopenmp-link=static -pthread -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -L${ADVISOR/' ./make_setting
	elif [[ "$1" = *"gnu"* ]]; then
		if [ -n "$FJMPI" ]; then sed -i -e 's/^CC =.*/CC = mpifcc/g' -e 's/^FC =.*/FC = mpifrt/g' ./make_setting; fi
		sed -i -e 's/^DEFINE += -DNO_METIS/#DEFINE += -DNO_METIS/g' -e "s#\$(HOME)/opt_intel/metis5#$ROOTDIR/$BM/src/metis-5.1.0#g" ./make_setting
		sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt_intel/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/x86_64-linux #" -e "s/-ipo -xHost -mcmodel=large -shared-intel/-march=native -fallow-argument-mismatch -fno-lto ${MAYBESTATIC}/g" -e "s/LIBS += -L\${ADVISOR.*/LIBS += -fno-lto ${MAYBESTATIC}/" -e 's# -I${ADVISOR_2018_DIR}/include##g' ./make_setting
	elif [[ "$1" = *"fujitrad"* ]]; then
		rm -f ./make_setting; cp ./make_setting.k ./make_setting
		sed -i -e 's/-DPROF_MAPROF/-DNO_PROF_MAPROF/g' ./make_setting
		sed -i '/CALL DDINIT/i \      LERR = 0' ./les3x.F
		sed -i -e "s#/opt/klocal#$ROOTDIR/$BM/src/metis-5.1.0#g" ./make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/kei #" -e 's/^FLAGS /FFLAGS /g' -e "s#REFINER)/include#REFINER)/Refiner#g" -e 's/-Kvisimpact,ocl -Qt/-Kfast,openmp,ocl,largepage,lto/g' -e 's/-Kvisimpact,ocl/-Kvisimpact -Kfast,openmp,ocl,largepage/g' ./make_setting
	elif [[ "$1" = *"fujiclang"* ]]; then
		rm -f ./make_setting; cp ./make_setting.k ./make_setting
		sed -i -e 's/-DPROF_MAPROF/-DNO_PROF_MAPROF/g' ./make_setting
		sed -i '/CALL DDINIT/i \      LERR = 0' ./les3x.F
		sed -i -e "s#/opt/klocal#$ROOTDIR/$BM/src/metis-5.1.0#g" ./make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		#sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/kei #" -e 's/^FLAGS /FFLAGS /g' -e "s#REFINER)/include#REFINER)/Refiner#g" -e 's/-Kvisimpact,ocl -Qt/-Nclang -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' -e 's/-Kvisimpact,ocl/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage/g' ./make_setting
		sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/kei #" -e 's/^FLAGS /FFLAGS /g' -e "s#REFINER)/include#REFINER)/Refiner#g" -e 's/-Kvisimpact,ocl -Qt/-Kfast,openmp,ocl,nolargepage,nolto/g' -e 's/-Kvisimpact,ocl/-Kvisimpact -Kfast,openmp,ocl,nolargepage/g' ./make_setting
		#XXX: LD is frt, so clang lto does not work
	elif [[ "$1" = *"gem5"* ]]; then
		rm -f ./make_setting; cp ./make_setting.k ./make_setting
		# crashing in some stupid yaml shit with fujitsu compilers
		sed -i -e 's/-DPROF_MAPROF/-DNO_PROF_MAPROF/g' ./make_setting
		sed -i '/CALL DDINIT/i \      LERR = 0' ./les3x.F
		sed -i -e "s/ = mpi/ = /g" ./make_setting
		sed -i -e "s#/opt/klocal#$ROOTDIR/$BM/src/metis-5.1.0#g" ./make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/kei #" -e 's/^FLAGS /FFLAGS /g' -e "s#REFINER)/include#REFINER)/Refiner#g" -e "s#-Kvisimpact,ocl -Qt#-Kfast,openmp,ocl,nolargepage,nolto -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s#-Kvisimpact,ocl#-Kvisimpact -Kfast,openmp,ocl,nolargepage -I$ROOTDIR/dep/mpistub/include/mpistub#g" ./make_setting
		#XXX: does not compile in clang mode
		sed -i -e "s#^LIBS += -lRcapRefiner#LIBS += -lRcapRefiner -L$ROOTDIR/$BM -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi -lmpifort#g" ./make_setting
		sed -i -e '/use mpi/d' -e "/implicit none/a \  include 'mpif.h'" ./ma_prof/src/mod_maprof.F90
		sed -i -e '/use mpi/d' -e "/use makemesh/a \  include 'mpif.h'" ./ffb_mini_main.F90
	elif [[ "$1" = *"llvm12"* ]]; then
		rm -f ./make_setting; cp ./make_setting.k ./make_setting
		sed -i -e 's/-DPROF_MAPROF/-DNO_PROF_MAPROF/g' ./make_setting
		sed -i '/CALL DDINIT/i \      LERR = 0' ./les3x.F
		sed -i -e "s#/opt/klocal#$ROOTDIR/$BM/src/metis-5.1.0#g" ./make_setting
		sed -i -E 's/(fcc|FCC|frt)px/\1/g' ./make_setting
		sed -i -e 's/^DEFINE += -DNO_REFINER/#DEFINE += -DNO_REFINER/g' -e "s#\$(HOME)/opt/REVOCAP_Refiner#$ROOTDIR/$BM/src/REVOCAP_Refiner-1.1.01#g" -e "s#REFINER)/lib #REFINER)/lib/kei #" -e 's/^FLAGS /FFLAGS /g' -e "s#REFINER)/include#REFINER)/Refiner#g" -e 's/-Kvisimpact,ocl -Qt/-mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' -e 's/-Kvisimpact,ocl/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly/g' -e "s#^LDFLAGS =.*#LDFLAGS = -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" ./make_setting
		#XXX: LD is frt, so clang lto does not work
	fi
	make
	cd $ROOTDIR
fi

