#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ NTChem =========================================
source conf/ntchem.sh $ROOTDIR
LOG="$ROOTDIR/log/bestrun/ntchem.log"
mkdir -p `dirname $LOG`
cd $DATA_DIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $APPDIR/$BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRunsBEST`; do
		echo "Start at " `date --iso-8601=s` >> $LOG 2>&1
		START="`date +%s.%N`"
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $APPDIR/$BINARY $INPUT >> $LOG 2>&1
		ENDED="`date +%s.%N`"
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		echo "Ended at " `date --iso-8601=s` >> $LOG 2>&1
	done
done
echo "Best NTChem run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
