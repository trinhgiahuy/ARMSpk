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
source activate idp
LOG="$ROOTDIR/log/`hostname -s`/testrun/candle.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for TEST in $TESTCONF; do
	for BINARY in $BINARYS; do
		NumMPI=1
		NumOMP=$TEST
		pushd "`find . -name $BINARY -exec dirname {} \;`"
		make libssc.so
		# check if data is hot or must be preloaded
		python ./p1b1.py
		#echo "mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT" >> $LOG 2>&1
		echo "OMP_NUM_THREADS=$NumOMP KMP_BLOCKTIME=30 KMP_SETTINGS=1 KMP_AFFINITY=\"granularity=fine,compact,1,0\" numactl --preferred 1 python $BINARY $INPUT" >> $LOG 2>&1
		NROUND=1
		for i in `seq 1 $NumRunsTEST`; do
			START="`date +%s.%N`"
			#timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT >> $LOG 2>&1
			export OMP_NUM_THREADS=$NumOMP; export KMP_BLOCKTIME=30; export KMP_SETTINGS=1; export KMP_AFFINITY="granularity=fine,compact,1,0"
			timeout --kill-after=30s $MAXTIME numactl --preferred 1 python $BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; NROUND=0; fi
			ENDED="`date +%s.%N`"
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
			if [ "x${NROUND}" = "x0" ]; then break; fi
		done
		popd
	done
done
echo "Best CANDLE run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
