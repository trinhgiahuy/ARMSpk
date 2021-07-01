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

source $ROOTDIR/conf/${BenchID}.sh
LOG="${ROOTDIR}/log/$(hostname -s)/bestrun/${BenchID}.log"
mkdir -p $(dirname $LOG)
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for BEST in $BESTCONF; do
	NumMPI="$(echo $BEST | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
	NumOMP="$(echo $BEST | cut -d '|' -f2)"
	PX="$(echo $BEST | cut -d '|' -f3)"
	PY="$(echo $BEST | cut -d '|' -f4)"
	PZ="$(echo $BEST | cut -d '|' -f5)"
	BINARYYY="${BINARY}_${PX}${PY}${PZ}"
	echo "$(get_mpi_cmd $NumMPI $NumOMP $LOG "") $BINARYYY $INPUT" >> $LOG 2>&1
	for i in $(seq 1 $NumRunsBEST); do
		START="$(date +%s.%N)"
		timeout --kill-after=30s $MAXTIME $(get_mpi_cmd $NumMPI $NumOMP $LOG "") $BINARYYY $INPUT >> $LOG 2>&1
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
		ENDED="$(date +%s.%N)"
		echo "Total running time: $(echo "$ENDED - $START" | bc -l)" >> $LOG 2>&1
	done
done
echo "Best ${BenchID} run:"
BEST="$(/bin/grep 'BiCGStab Total FLOPS:' $LOG | awk -F 'FLOPS:' '{print $2}' | sort -r -g | head -1)"
/bin/grep "$BEST\|mpiexec" $LOG | /bin/grep -B1 "$BEST"
echo ""
cd $ROOTDIR
