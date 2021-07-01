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
source activate idp
LOG="${ROOTDIR}/log/$(hostname -s)/bestrun/${BenchID}.log"
mkdir -p $(dirname $LOG)
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for BEST in $BESTCONF; do
	for BINARY in $BINARYS; do
		NumMPI="$(echo $BEST | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
		NumOMP="$(echo $BEST | cut -d '|' -f2)"
		pushd "$(find . -name $BINARY -exec dirname {} \;)"
		make libssc.so
		# check if data is hot or must be preloaded
		python ./p1b1.py
		#echo "mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT" >> $LOG 2>&1
		echo "OMP_NUM_THREADS=$NumOMP KMP_BLOCKTIME=30 KMP_SETTINGS=1 KMP_AFFINITY='granularity=fine,compact,1,0' numactl --preferred 1 python $BINARY $INPUT" >> $LOG 2>&1
		for i in $(seq 1 $NumRunsBEST); do
			START="$(date +%s.%N)"
			#timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT >> $LOG 2>&1
			export OMP_NUM_THREADS=$NumOMP; export KMP_BLOCKTIME=30; export KMP_SETTINGS=1; export KMP_AFFINITY='granularity=fine,compact,1,0';
			timeout --kill-after=30s $MAXTIME numactl --preferred 1 python $BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="$(date +%s.%N)"
			echo "Total running time: $(echo "$ENDED - $START" | bc -l)" >> $LOG 2>&1
		done
		popd
	done
done
echo "Best ${BenchID} run:"
BEST="$(/bin//bin/grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1)"
/bin//bin/grep "$BEST\|mpiexec" $LOG | /bin//bin/grep -B1 "$BEST"
echo ""
cd $ROOTDIR
