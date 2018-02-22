#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

export PATH=$ROOTDIR/dep/sde-external-8.16.0-2018-01-30-lin:$PATH
SDE="`which sde64` -sse-sde -global_region -mix_omit_per_thread_stats -mix_omit_per_function_stats -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -iform 1 -omix oSDE/\"\$MPI_LOCALRANKID\".txt"
if [[ $HOSTNAME = *"kiev"* ]]; then
	SDE="$SDE -bdw -- "
elif [[ $HOSTNAME = *"lyon"* ]]; then
	SDE="$SDE -knl -- "
elif [[ $HOSTNAME = *"mill"* ]]; then
	SDE="$SDE -knm -- "
else
	echo "Unsupported host"
	exit
fi

# ============================ SWFFFT =========================================
source conf/swfft.sh
LOG="$ROOTDIR/log/profrun/swfft.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	mkdir -p ./oSDE
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	# test if Decomposition is valid
	INSIZE="`echo $INPUT | awk '{print $2}'`"
	`dirname $BINARY`/CheckDecomposition $INSIZE $INSIZE $INSIZE $NumMPI > /dev/null 2>&1
	if [ "x$?" = "x0" ]; then
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE $BINARY $INPUT\"" >> $LOG 2>&1
	else
		echo "INVALID Decomposition: mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE $BINARY $INPUT\"" >> $LOG 2>&1
		continue
	fi
	mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c "$SDE $BINARY $INPUT" >> $LOG 2>&1
	for P in `seq 0 $((NumMPI - 1))`; do
		echo "SDE output of MPI process $P" >> $LOG 2>&1
		cat ./oSDE/${P}.txt >> $LOG 2>&1
	done
	echo "=== SDE summary ===" >> $LOG 2>&1
	$ROOTDIR/util/analyze_sde.py `echo $LOG | sed 's#profrun#bestrun#g'` ./oSDE >> $LOG 2>&1
	rm -rf ./oSDE
done
cd $ROOTDIR
