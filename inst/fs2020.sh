#!/bin/bash

function prep_dir {
	for FILE in `/bin/grep '__builtin_fj_prefetch' -r $1 | cut -d':' -f1 | sort -u`; do sed -i -e 's#__builtin_fj_prefetch#//__builtin_fj_prefetch#g' $FILE; done
	for FILE in `/bin/grep 'PROF_START_ALL\|PROF_STOP_ALL' -r $1 | /bin/grep -v 'define' | cut -d':' -f1 | sort -u`; do
		if [[ "$FILE" = *".c" ]] || [[ "$FILE" = *".cc" ]]; then
			sed -i -e '/#include.*"profiler.h"/a #include <time.h>\n#include <stdio.h>' -e 's/PROF_INIT/double mkrts, mkrte; struct timespec mkrtsclock;/g' -e 's/PROF_START_ALL/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/g' -e 's/PROF_STOP_ALL/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/g' -e 's/PROF_FINALIZE/fprintf(stdout,"Walltime of the main kernel: %.6lf sec\\n", mkrte - mkrts);/g' $FILE
		else
			if ! [[ "$FILE" = *"postK_nicam_divdamp_ijsplit_tune03.f90" ]]; then sed -i -n 'p; s/implicit none/real :: mkrts, mkrte/p' $FILE; fi
			sed -i -n 'p; s/PROF_START_ALL/call cpu_time(mkrts)/p' $FILE
			sed -i -n 'p; s/PROF_STOP_ALL/call cpu_time(mkrte)/p' $FILE
			if [[ "$FILE" = *".F" ]] || [[ "$FILE" = *".f" ]]; then
				sed -i -n 'p; s/PROF_FINALIZE/write(*,"(A,f11.6,A)") "Walltime of the main kernel: ", \n     +  mkrte-mkrts, " sec"/p' $FILE
			else
				sed -i -n 'p; s/PROF_FINALIZE/write(*,"(A,f11.6,A)") "Walltime of the main kernel: ", \&\n     \&  mkrte-mkrts, " sec"/p' $FILE
			fi
			if [[ "$FILE" = *"kernel_pairlist_july.f90" ]] || [[ "$FILE" = *"kernel_july.f90" ]]; then
				sed -i -n 'p; s/integer(1),.*pointer,contiguous :: exclusion_mask.*/real :: mkrts, mkrte/p' $FILE
			fi
			if [[ "$FILE" = *"postK_nicam_divdamp_ijsplit_tune03.f90" ]]; then
				sed -i -n 'p; s/integer :: i$/real :: mkrts, mkrte/p' $FILE
			fi
		fi
	done
	if [ -z $2 ] && [[ "$1" = *"18.Adventure.region1.tune01.outputcheck.160808.M24"* ]]; then
		for FILE in `/bin/grep 'implicit none' -r $1 | cut -d':' -f1 | sort -u`; do sed -i '0,/.*implicit none.*/ s///' $FILE; done
	fi
}

function do_dir {
	local OBJFILE=""
	echo "==================== $1 ===================="
	cd $1
	/bin/grep -A 4 $(basename $(pwd) | cut -d'.' -f2-) ../option.list | tail -3 > Makefile
	cat ../Makefile.inc >> Makefile

	local OBJFILES="$(find . -maxdepth 1 -iname '*.[cCsSfF]' -or -iname '*.f90' -or -iname '*.cc' | cut -c3-)"
	local MAINFILE="$(/bin/grep -l -e '[^a-zA-Z]main' -i -e program ${OBJFILES})"
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
if [ -z $1 ]; then
	source $ROOTDIR/conf/intel.cfg
	source $INTEL_PACKAGE intel64 > /dev/null 2>&1
	alias ar=`which xiar`
	alias ld=`which xild`
	export ADVISOR_2018_DIR=${ADVISOR_2019_DIR}
elif [[ "$1" = *"gnu"* ]]; then
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
elif [[ "$1" = *"fuji"* ]]; then
	module load FujitsuCompiler/201903
elif [[ "$1" = *"arm"* ]]; then
	module load /opt/arm/modulefiles/A64FX/RHEL/8/arm-linux-compiler-20.3/armpl/20.3.0
	#module load /opt/arm/modulefiles/A64FX/RHEL/8/gcc-9.3.0/armpl/20.3.0
fi

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
CC=fccpx
CXX=FCCpx
FC=frtpx
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
	if [ -z $1 ]; then
		sed -i -e 's/CC=.*/CC=icc/g' -e 's/CXX=.*/CXX=icpc/g' -e 's/FC=.*/FC=ifort/g' -e 's/-Kopenmp -Bstatic/-qopenmp -static -static-intel/g' Makefile.inc
		sed -i -e 's/$(.OPTIMIZE) -Bstatic/$(DEF00) -fpp -std=gnu99 -qopenmp -O3 -xHost -static -static-intel/g' Makefile.inc
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/CC=.*/CC=gcc/g' -e 's/CXX=.*/CXX=g++/g' -e 's/FC=.*/FC=gfortran/g' -e 's/-Kopenmp -Bstatic/-fopenmp -static/g' Makefile.inc
		sed -i -e 's/$(.OPTIMIZE) -Bstatic/$(DEF00) -fopenmp -O3 -march=native -static/g' Makefile.inc
	elif lscpu | grep 'sve' >/dev/null 2>&1 && [[ "$1" = *"fuji"* ]]; then
		sed -i -e 's/-Bstatic//g' Makefile.inc
	elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
		sed -i -e 's/-Bstatic//g' Makefile.inc
	elif [[ "$1" = *"arm"* ]]; then
		sed -i -e 's/CC=.*/CC=armclang/g' -e 's/CXX=.*/CXX=armclang++/g' -e 's/FC=.*/FC=armflang/g' -e 's/-Kopenmp -Bstatic/-fopenmp/g' Makefile.inc
		sed -i -e 's/$(.OPTIMIZE) -Bstatic/$(DEF00) -fopenmp -mcpu=a64fx -march=armv8.2-a+sve -Ofast -ffast-math -flto/g' Makefile.inc
	fi
	for BMconf in ${BINARYS}; do
		NR=`echo $BMconf | cut -d'.' -f1`; sed -i -e "s/(DEF..)/(DEF$NR)/g" Makefile.inc
		prep_dir $ROOTDIR/$BM/$(dirname $BMconf) $1
		do_dir $ROOTDIR/$BM/$(dirname $BMconf)
	done
	find $ROOTDIR/$BM/ -executable -type f -exec echo {} \; -exec ldd {} \;
	cd $ROOTDIR
fi

