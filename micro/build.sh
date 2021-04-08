#!/bin/bash
# set -e
# set -o pipefail

function do_dir {
	local OBJFILE=""
	echo "==================== $1 ===================="
	cd $1
	grep -A 4 $(basename $(pwd) | cut -d'.' -f2-) ../option.list | tail -3 > Makefile
	cat ../../../Makefile.inc >> Makefile

	local OBJFILES="$(find . -maxdepth 1 -iname '*.[cCsSfF]' -or -iname '*.f90' -or -iname '*.cc' | cut -c3-)"
	local MAINFILE="$(grep -l -e '[^a-zA-Z]main' -i -e program ${OBJFILES})"
	local BIN=${MAINFILE%.*}
	local EXT=${MAINFILE##*.}

	for i in ${OBJFILES}; do
		if [[ "$i" = *"$MAINFILE"* ]]; then continue; fi
		OBJFILE="$OBJFILE ${i%.*}.o"
	done
	sed -i -e "s/OBJFILE=.*/OBJFILE= $OBJFILE/g" Makefile

	echo "${MAINFILE%.*}.o: \$(OBJFILE)" >> Makefile
	echo "${BIN}: ${MAINFILE%.*}.o \$(OBJFILE) \$(COMMON)" >> Makefile
	if grep -i program $MAINFILE; then
		echo -e '\t$(FC) $^ $(LDLIBS) -o $@' >> Makefile
	elif [[ "$MAINFILE" = *".cc" ]]; then
		echo -e '\t$(CXX) $^ $(LDLIBS) -o $@' >> Makefile
	fi
	if lscpu | grep 'sve' >/dev/null 2>&1 && which fcc >/dev/null 2>&1 ; then
		sed -i -e 's/fccpx/fcc/g' -e 's/FCCpx/FCC/g' -e 's/frtpx/frt/g' Makefile
	fi

	make -B ${BIN}
	cd -
}

find . -mindepth 3 -maxdepth 3 -type d | while read PROG; do
	for FILE in `/bin/grep '__builtin_fj_prefetch' -r $PROG | cut -d':' -f1 | sort -u`; do sed -i -e 's#__builtin_fj_prefetch#//__builtin_fj_prefetch#g' $FILE; done
	for FILE in `/bin/grep 'PROF_START_ALL\|PROF_STOP_ALL' -r $PROG | /bin/grep -v 'define' | cut -d':' -f1 | sort -u`; do
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
	do_dir $PROG
done
