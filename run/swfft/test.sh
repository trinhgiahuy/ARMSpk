#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ SWFFFT =========================================
source conf/swfft.sh
LOG="$ROOTDIR/log/`hostname -s`/testrun/swfft.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for TEST in $TESTCONF; do
	NumMPI="`echo $TEST | cut -d '|' -f1`"
	NumOMP="`echo $TEST | cut -d '|' -f2`"
	# test if Decomposition is valid
	INSIZE="`echo $INPUT | awk '{print $2}'`"
	`dirname $BINARY`/CheckDecomposition $INSIZE $INSIZE $INSIZE $NumMPI > /dev/null 2>&1
	if [ "x$?" = "x0" ]; then
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	else
		echo "INVALID Decomposition: mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		continue
	fi
	for i in `seq 1 $NumRunsTEST`; do
		START="`date +%s.%N`"
		timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
		ENDED="`date +%s.%N`"
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
	done
done
echo "Best SWFFT run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
