#!/bin/bash

function do_dir {
	local OBJFILE=""
	echo "==================== $1 ===================="
	cd $1
	/bin/grep -A 4 $(basename $(pwd) | cut -d'.' -f2-) ../option.list | tail -3 > Makefile
	cat ../Makefile.inc >> Makefile

	local OBJFILES="$(find . -maxdepth 1 -iname '*.[cCsSfF]' -or -iname '*.f90' -or -iname '*.cc' | cut -c3-)"
	local MAINFILE="$(/bin/grep -l -e '[^a-zA-Z]main[[:space:]]*(' -i -e program ${OBJFILES})"
	local BIN=${MAINFILE%.*}
	local EXT=${MAINFILE##*.}

	for i in ${OBJFILES}; do
		if [[ "$i" = *"$MAINFILE"* ]]; then continue; fi
		OBJFILE="$OBJFILE ${i%.*}.o"
	done
	sed -i -e "s/OBJFILE=.*/OBJFILE= $OBJFILE/g" Makefile

	echo "${MAINFILE%.*}.o: \$(OBJFILE)" >> Makefile
	echo "${BIN}: ${MAINFILE%.*}.o \$(OBJFILE) \$(COMMON)" >> Makefile
	if /bin/grep -i program $MAINFILE; then
		echo -e '\t$(FC) $^ $(LDLIBS) -o $@' >> Makefile
	elif [[ "$MAINFILE" = *".cc" ]]; then
		echo -e '\t$(CXX) $^ $(LDLIBS) -o $@' >> Makefile
	fi

	make -B ${BIN}
	cd -
}


ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

if [ ! -f $ROOTDIR/dep/fs2020_microkernel_tmp.zip ]; then
	echo "ERR: Cannot find fs2020_microkernel_tmp.zip; if you have them, place them in ./dep subfolder"
	exit
fi

BM="fs2020"
source $ROOTDIR/conf/${BM}.sh
if [ ! -f $ROOTDIR/$BM/14.streamlike_pattern1/main ]; then
	mkdir -p $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/
	if ! [ -f $ROOTDIR/$BM/option.list ]; then
		unzip -d $ROOTDIR/$BM/ $ROOTDIR/dep/fs2020_microkernel_tmp.zip && f=($ROOTDIR/$BM/*) && mv $ROOTDIR/$BM/*/* $ROOTDIR/$BM/ && rmdir "${f[@]}"
	fi
	patch -p1 --forward < $ROOTDIR/patches/*1-${BM}*.patch
	cat <<EOF > Makefile.inc
DEF00=
DEF01=
DEF02=
DEF03=
DEF04=
DEF05=
DEF06=
DEF07=
DEF08=-DSINGLE
DEF09=-DSINGLE
DEF10=-DSINGLE
DEF11=-DSINGLE
DEF12=-DSINGLE
DEF13=-DSINGLE
DEF14=-DSINGLE
DEF15=-DSINGLE
DEF16=-DSINGLE
DEF17=
DEF18=
DEF19=
DEF20=
DEF21=
DEF22=-DRDC -DV512
DEF23=-DRDC -DVLENS=16 -DPREFETCH -D_CHECK_SIM -DTARGET_JINV
CC=fcc
CXX=FCC
FC=frt
COMMON=common/src/report.o
LDLIBS=-Kopenmp -Bstatic

FFLAGS=-Icommon/include -DDISABLE_VALIDATION \$(FOPTIMIZE) -Bstatic
CFLAGS=-Icommon/include -DDISABLE_VALIDATION \$(COPTIMIZE) -Bstatic
CXXFLAGS=-Icommon/include -DDISABLE_VALIDATION \$(COPTIMIZE) -Bstatic
OBJFILE=

%.o : %.F90
	\$(FC) -c \$(FFLAGS) $< -o \$@

%.o : %.f90
	\$(FC) -c \$(FFLAGS) $< -o \$@

EOF
	instrument_kernel "$1" $ROOTDIR/$BM/
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/CC=.*/CC=icc/g' -e 's/CXX=.*/CXX=icpc/g' -e 's/FC=.*/FC=ifort/g' -e 's#-Kopenmp -Bstatic#-qopenmp -static -static-intel -L${ADVISOR_2018_DIR}/lib64 -littnotify#g' Makefile.inc
		sed -i -e 's#$(.OPTIMIZE) -Bstatic#$(DEF00) -fpp -std=gnu99 -qopenmp -O3 -xHost -static -static-intel -I${ADVISOR_2018_DIR}/include#g' Makefile.inc
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/CC=.*/CC=gcc/g' -e 's/CXX=.*/CXX=g++/g' -e 's/FC=.*/FC=gfortran/g' -e 's/-Kopenmp -Bstatic/-fopenmp -static/g' Makefile.inc
		sed -i -e 's/$(.OPTIMIZE) -Bstatic/$(DEF00) -fopenmp -O3 -march=native -static/g' Makefile.inc
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's/-Bstatic//g' Makefile.inc
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's/$(FOPTIMIZE) -Bstatic/$(DEF00) -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -Kfast,ocl,largepage,lto/g' -e 's/$(COPTIMIZE) -Bstatic/$(DEF00) -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage -flto/g' Makefile.inc
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/$(FOPTIMIZE) -Bstatic/$(DEF00) -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -Kfast,ocl,nolargepage,nolto/g' -e 's/$(COPTIMIZE) -Bstatic/$(DEF00) -Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' Makefile.inc
	elif [[ "$1" = *"arm"* ]]; then
		sed -i -e 's/CC=.*/CC=armclang/g' -e 's/CXX=.*/CXX=armclang++/g' -e 's/FC=.*/FC=armflang/g' -e 's/-Kopenmp -Bstatic/-fopenmp/g' Makefile.inc
		sed -i -e 's/$(.OPTIMIZE) -Bstatic/$(DEF00) -Ofast -ffast-math -fopenmp -march=armv8.2-a+sve -mcpu=a64fx -mtune=a64fx -flto/g' Makefile.inc
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's/ -V / /g' -e 's/ -v / /g' $ROOTDIR/$BM/option.list		#flang doesn'f like -V
		sed -i -e 's/CC=.*/CC=clang/g' -e 's/CXX=.*/CXX=clang++/g' -e 's/FC=.*/FC=frt/g' -e 's/-Kopenmp -Bstatic/-fopenmp/g' Makefile.inc
		sed -i -e 's/$(FOPTIMIZE) -Bstatic/$(FOPTIMIZE) -Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -Kfast,ocl,largepage,lto/g' -e 's/$(COPTIMIZE) -Bstatic/$(DEF00) -fopenmp -Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -mllvm -polly -mllvm -polly-vectorizer=polly -flto/g' -e "s#-Kopenmp -Bstatic#-fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" Makefile.inc
	fi
	for BMconf in ${BINARYS}; do
		NR=`echo $BMconf | cut -d'.' -f1`; sed -i -e "s/(DEF..)/(DEF$NR)/g" Makefile.inc
		do_dir $ROOTDIR/$BM/$(dirname $BMconf)
	done
	find $ROOTDIR/$BM/ -executable -type f -exec echo {} \; -exec ldd {} \;
	cd $ROOTDIR
fi
