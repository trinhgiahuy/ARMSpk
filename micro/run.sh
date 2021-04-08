#!/bin/bash
ulimit -s unlimited

if [ -z $1 ]; then NTHREADS=$(nproc); else NTHREADS=$1; fi
export OMP_NUM_THREADS=$NTHREADS
export OMP_NUM_PARALELL=$NTHREADS
export OMP_PROC_BIND=close
export FLIB_FASTOMP=FALSE
export FLIB_CNTL_BARRIER_ERR=FALSE
export KMP_AFFINITY=granularity=fine,compact,1,0

find -mindepth 2 -executable  -type f | sort | while read BIN; do
	if lscpu | grep 'sve' >/dev/null 2>&1; then
		if [ $NTHREADS -le 12 ]; then
			PIN="numactl -N 4 -m 4"
		elif [ $NTHREADS -le 24 ]; then
			PIN="numactl -N 4,5 -m 4,5"
		elif [ $NTHREADS -le 36 ]; then
			PIN="numactl -N 4,5,6 -m 4,5,6"
		else
			PIN=""
		fi
	fi
	echo -n $(basename $(dirname $BIN))
	rm -f ${BIN}.log
	for x in `seq 1 10`; do
		if [[ "$BIN" = *"19.Adventure.region2.tune161011.armtest-M24.170727"* ]] && [ $NTHREADS -gt 24 ]; then continue; fi
		$PIN $BIN >> ${BIN}.log 2>&1
	done
	if /bin/grep '^Walltime' ${BIN}.log >/dev/null 2>&1; then
		/bin/grep '^Walltime' ${BIN}.log | awk -F 'kernel:' '{print $2}' | sort -g | head -1
	else
		echo ''
	fi
done
