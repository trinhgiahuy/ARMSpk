#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ FFB ============================================
source conf/ffb.sh
NumRUNS=5
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/testrun/ffb.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for TEST in $TESTCONF; do
	NumMPI="`echo $TEST | cut -d '|' -f1`"
	NumOMP="`echo $TEST | cut -d '|' -f2`"
	X="`echo $TEST | cut -d '|' -f3`"
	Y="`echo $TEST | cut -d '|' -f4`"
	Z="`echo $TEST | cut -d '|' -f5`"
	INPUT="`echo $DEFINPUT | sed -e \"s/PX/$X/\" -e \"s/PY/$Y/\" -e \"s/PZ/$Z/\"`"
	# try finding closest cube size per proc
	FLOAT=`echo "e((1/3)*l($MAXDCZ / $NumMPI))" | bc -l`
	DCZ=`echo "($FLOAT+0.5)/1" | bc`
	INPUT="`echo $INPUT | sed -e \"s/DCZ/$DCZ/\"`"
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		mkdir ./tmp; sleep 2; cd ./tmp
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT >> $LOG 2>&1
		cd ../; rm -rf ./tmp; sleep 2
	done
done
echo "Best FFB run:"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
echo ""
cd $ROOTDIR