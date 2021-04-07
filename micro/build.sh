#!/usr/bin/bash
# set -e 
# set -o pipefail

function do_dir {
  local OBJFILE=""
  echo "==================== $1 ===================="
  cd $1
  grep -A 4 $(basename $(pwd)) ../option.list | tail -3 > Makefile 
  cat ../Makefile >> Makefile 
  
  local OBJFILES="$(find  -maxdepth 1 -iname '*.[cCsSfF]' -or -iname '*.f90' -or -iname '*.cc' | cut -c3-)"
  local MAINFILE="$(grep -l -e '[^a-zA-Z]main' -i -e program ${OBJFILES})"
  local BIN=${MAINFILE%.*}
  local EXT=${MAINFILE##*.}
  
  for i in ${OBJFILES}; do
    OBJFILE="$OBJFILE ${i%.*}.o"
  done
  sed -i "s/OBJFILE =.*/OBJFILE = $OBJFILE/"  Makefile 
  
  echo ${BIN}: '$(OBJFILE) $(COMMON)' >> Makefile
  if grep -i program $MAINFILE; then
    echo '	$(FC) $^ $(LDLIBS) -o $@' >> Makefile
  fi
  
  make -j -B 
  cd -
}

mv -f Adventure.region0.tune4pad-acle.armtest.20170821_resize/adventure_kernel_region0_tune0000644 Adventure.region0.tune4pad-acle.armtest.20170821_resize/adventure_kernel_region0_tune0000644.c

find -mindepth 1 -maxdepth 1 -type d | while read PROG; do
  do_dir $PROG
done
