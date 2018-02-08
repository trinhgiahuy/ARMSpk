#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096

# ============================ miniFE =========================================
source conf/minife.sh
NumRUNS=10
LOG="$ROOTDIR/log/bestrun/minife.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	for BINARY in $BINARYS; do
		NumMPI="`echo $BEST | cut -d '|' -f1`"
		NumOMP="`echo $BEST | cut -d '|' -f2`"
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRUNS`; do
			echo "Start at " `date --iso-8601=s` >> $LOG 2>&1
			mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			echo "Ended at " `date --iso-8601=s` >> $LOG 2>&1
			cat miniFE*.yaml >> $LOG; rm miniFE*.yaml
		done
	done
done
echo "Best miniFE run:"
BEST="`grep 'Total CG Mflops' $LOG | awk -F 'Mflops:' '{print $2}' | sort -n -r | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
cd $ROOTDIR
