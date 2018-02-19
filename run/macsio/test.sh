#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ MACSio =========================================
source conf/macsio.sh
NumRUNS=3
LOG="$ROOTDIR/log/testrun/macsio.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for TEST in $TESTCONF; do
	NumMPI="`echo $TEST | cut -d '|' -f1`"
	NumOMP="`echo $TEST | cut -d '|' -f2`"
	mkdir -p ./testrun; cd ./testrun
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT >> $LOG 2>&1
		grep 'Processor\|^Info' macsio-log.log >> $LOG 2>&1
	done
	cd ../; rm -rf ./testrun
done
echo "Best MACSio run:"
BEST="`grep 'Summed  BW:' $LOG | awk -F 'BW:' '{print $2}' | sort -r -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
cd $ROOTDIR
