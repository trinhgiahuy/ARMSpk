#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

if which numactl >/dev/null 2>&1; then
	if [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then PIN="numactl -l -C 24";
	else PIN="numactl -l -C 2"; fi
else PIN=""; fi

# ============================ PolyBench ====================================
source conf/polybench.sh
DEFLOG="$ROOTDIR/log/`hostname -s`/bestrun/polybench"
mkdir -p $DEFLOG
cd $APPDIR
for BEST in $BESTCONF; do
	for BMconf in $BINARYS; do
		NumMPI=1
		NumOMP=1
		BINARY="`echo ${BMconf} | tr '|' '_'`"
		BName="`basename $(echo ${BMconf} | cut -d'|' -f1)`"
		LOG="${DEFLOG}/${BName}.log"
		echo "OMP_NUM_THREADS=$NumOMP timeout --kill-after=30s $MAXTIME $PIN $BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRunsBEST`; do
			START="`date +%s.%N`"
			OMP_NUM_THREADS=$NumOMP timeout --kill-after=30s $MAXTIME $PIN $BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="`date +%s.%N`"
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		done
		BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
		echo "Best $BName run: $BEST"
	done
done
cd $ROOTDIR
