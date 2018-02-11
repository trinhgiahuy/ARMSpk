#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"
#-genv KMP_PLACE_THREADS=1T -genv KMP_AFFINITY=compact -genv I_MPI_PIN_DOMAIN=node -genv I_MPI_PIN_ORDER=spread"

# ============================ SW4lite ========================================
source conf/sw4lite.sh
NumRUNS=20
LOG="$ROOTDIR/log/testrun/sw4lite.log"
mkdir -p `dirname $LOG`
cd $APPDIR
sed -i -e 's/corder=1/corder=0/' $INPUT
for TEST in $TESTCONF; do
	NumMPI="`echo $TEST | cut -d '|' -f1`"
	NumOMP="`echo $TEST | cut -d '|' -f2`"
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
	done
done
echo "Best SW4lite run:"
BEST="`grep 'Total running time:' $LOG | awk -F 'time:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR