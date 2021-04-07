#!/usr/bin/bash
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
		echo '	$(FC) $^ $(LDLIBS) -o $@' >> Makefile
	fi

	make -B ${BIN}
	cd -
}

find . -mindepth 3 -maxdepth 3 -type d | while read PROG; do
	for FILE in `/usr/bin/grep '__builtin_fj_prefetch' -r $PROG | cut -d':' -f1 | sort -u`; do sed -i -e 's#__builtin_fj_prefetch#//__builtin_fj_prefetch#g' $FILE; done
	do_dir $PROG
done
