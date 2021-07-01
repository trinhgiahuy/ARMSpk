#!/bin/bash

SELF="$(readlink -f "${BASH_SOURCE[0]}")"
ROOTDIR="$(readlink -f $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../)"
BenchID="$(basename $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
get_comp_env_name "$1"
maybe_submit_job "${COMP}" "${SELF}" "${ROOTDIR}/conf/${BenchID}.sh"
load_compiler_env "${COMP}"

if which numactl >/dev/null 2>&1; then
	if [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then PIN="numactl -l -C 24";
	else PIN="numactl -l -C 2"; fi
else PIN=""; fi

source $ROOTDIR/conf/${BenchID}.sh
LOGDIR="${ROOTDIR}/log/$(hostname -s)/bestrun/${BenchID}"
mkdir -p $LOGDIR
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for BEST in $BESTCONF; do
	for BMconf in $BINARYS; do
		NumMPI="$(echo $BEST | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
		NumOMP="$(echo $BEST | cut -d '|' -f2)"
		BINARY="$(echo ${BMconf} | tr '|' '_')"
		BName="$(basename $(echo ${BMconf} | cut -d'|' -f1))"
		LOG="${LOGDIR}/${BName}.log"
		echo "OMP_NUM_THREADS=$NumOMP timeout --kill-after=30s $MAXTIME $PIN $BINARY $INPUT" >> $LOG 2>&1
		for i in $(seq 1 $NumRunsBEST); do
			START="$(date +%s.%N)"
			OMP_NUM_THREADS=$NumOMP timeout --kill-after=30s $MAXTIME $PIN $BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="$(date +%s.%N)"
			echo "Total running time: $(echo "$ENDED - $START" | bc -l)" >> $LOG 2>&1
		done
		BEST="$(/bin/grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1)"
		echo "Best $BName run: $BEST"
	done
done
cd $ROOTDIR
