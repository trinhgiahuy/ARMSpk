#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096

# ============================ CANDLE =========================================
source conf/candle.sh
LOG="$ROOTDIR/log/`hostname -s`/testrun/candle.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for TEST in $TESTCONF; do
	for BINARY in $BINARYS; do
		NumMPI=1
		NumOMP=$TEST
		pushd "`find . -name $BINARY -exec dirname {} \;`"
		make libssc.so
		echo "mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRunsTEST`; do
			START="`date +%s.%N`"
			timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="`date +%s.%N`"
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		done
		popd
	done
done
echo "Best CANDLE run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
