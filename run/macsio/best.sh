#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ MACSio =========================================
source conf/macsio.sh
NumRUNS=10
LOG="$ROOTDIR/log/bestrun/macsio.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	mkdir -p ./bestrun; cd ./bestrun
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		echo "Start at " `date --iso-8601=s` >> $LOG 2>&1
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		echo "Ended at " `date --iso-8601=s` >> $LOG 2>&1
	done
done
cd $ROOTDIR