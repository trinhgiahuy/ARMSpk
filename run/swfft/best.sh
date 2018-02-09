#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"
#-genv KMP_PLACE_THREADS=1T -genv KMP_AFFINITY=compact -genv I_MPI_PIN_DOMAIN=node -genv I_MPI_PIN_ORDER=spread"

# ============================ SWFFFT =========================================
source conf/swfft.sh
NumRUNS=10
LOG="$ROOTDIR/log/bestrun/swfft.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		echo "Start at " `date --iso-8601=s` >> $LOG 2>&1
		START="`date +%s.%N`"
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		ENDED="`date +%s.%N`"
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		echo "Ended at " `date --iso-8601=s` >> $LOG 2>&1
	done
done
echo "Best SWFFT run:"
BEST="`grep 'Total running time:' $LOG | awk -F 'time:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
