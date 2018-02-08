#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"
#-genv KMP_PLACE_THREADS=1T -genv KMP_AFFINITY=compact -genv I_MPI_PIN_DOMAIN=node -genv I_MPI_PIN_ORDER=spread"

# ============================ CoMD ===========================================
source conf/comd.sh
NumRUNS=1
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/testrun/comd.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for TEST in $TESTCONF; do
	NumMPI="`echo $TEST | cut -d '|' -f1`"
	NumOMP="`echo $TEST | cut -d '|' -f2`"
	if [ "x${NumMPI}x" != "x1x" ]; then
		I="`echo $TEST | cut -d '|' -f3`"
		J="`echo $TEST | cut -d '|' -f4`"
		K="`echo $TEST | cut -d '|' -f5`"
		INPUT=$DEFINPUT
		INPUT="`echo $INPUT | sed -e \"s/i1/i$I/\"`"
		INPUT="`echo $INPUT | sed -e \"s/j1/j$J/\"`"
		INPUT="`echo $INPUT | sed -e \"s/k1/k$K/\"`"
	fi
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
	done
done
#echo "Best CoMD run:"
#BEST="`grep -A3 'export\|mpiexec\|Problem 1.*AMG-PCG Solve' $LOG | grep 'wall' | cut -d '=' -f2 | sort -n | head -1`"
#grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
#echo ""
cd $ROOTDIR
