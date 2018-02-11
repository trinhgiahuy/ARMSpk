#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ FFVC ===========================================
source conf/ffvc.sh
NumRUNS=5
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/testrun/ffvc.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for TEST in $TESTCONF; do
	NumMPI="`echo $TEST | cut -d '|' -f1`"
	NumOMP="`echo $TEST | cut -d '|' -f2`"
	X="`echo $TEST | cut -d '|' -f3`"
	Y="`echo $TEST | cut -d '|' -f4`"
	Z="`echo $TEST | cut -d '|' -f5`"
	INPUT=$DEFINPUT
	INPUT="`echo $INPUT | sed -e \"s/DX/$X/\"`"
	INPUT="`echo $INPUT | sed -e \"s/DY/$Y/\"`"
	INPUT="`echo $INPUT | sed -e \"s/DZ/$Z/\"`"
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
	done
done
echo "Best FFVC run:"
BEST="`grep 'Main Loop' $LOG | awk -F 'max)' '{print $2}' | awk -F '[' '{print $1}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
