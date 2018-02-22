#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ miniFE =========================================
source conf/minife.sh
LOG="$ROOTDIR/log/bestrun/minife.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	for BINARY in $BBINARY; do
		NumMPI="`echo $BEST | cut -d '|' -f1`"
		NumOMP="`echo $BEST | cut -d '|' -f2`"
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRunsBEST`; do
			START="`date +%s.%N`"
			mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			ENDED="`date +%s.%N`"
			cat miniFE*.yaml >> $LOG; rm miniFE*.yaml
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		done
	done
done
echo "Best miniFE run:"
BEST="`grep 'Total CG Mflops' $LOG | awk -F 'Mflops:' '{print $2}' | sort -r -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
