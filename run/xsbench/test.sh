#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"
#-genv KMP_PLACE_THREADS=1T -genv KMP_AFFINITY=compact -genv I_MPI_PIN_DOMAIN=node -genv I_MPI_PIN_ORDER=spread"

# ============================ XSBench ========================================
source conf/xsbench.sh
NumRUNS=5
LOG="$ROOTDIR/log/testrun/xsbench.log"
mkdir -p `dirname $LOG`
cd $APPDIR
DEFINPUT=$INPUT
for TEST in $TESTCONF; do
	NumMPI="`echo $TEST | cut -d '|' -f1`"
	NumOMP="`echo $TEST | cut -d '|' -f2`"
	INPUT="`echo $DEFINPUT | sed -e \"s/OMPNT/$NumOMP/\"`"
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
	done
done
echo "Best XSBench run:"
grep 'Total Lookups' $LOG | awk -F 's:' '{print $2}' > /dev/shm/1
cat /dev/shm/1 | sed -e 's/,//g' > /dev/shm/2
BEST="`paste /dev/shm/1 /dev/shm/2 | awk '{print $2 "\t" $1}' | sort -g -r | head -1 | awk '{print $2}'`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
