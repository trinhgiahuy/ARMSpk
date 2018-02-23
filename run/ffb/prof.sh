#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

export PATH=$ROOTDIR/dep/sde-external-8.16.0-2018-01-30-lin:$PATH
SDE="`which sde64` -sse-sde -global_region -mix_omit_per_thread_stats -mix_omit_per_function_stats -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -iform 1 -omix oSDE/\"\$MPI_LOCALRANKID\".txt"
if [[ $HOSTNAME = *"kiev"* ]]; then
	SDE="$SDE -bdw -- "
elif [[ $HOSTNAME = *"lyon"* ]]; then
	SDE="$SDE -knl -- "
elif [[ $HOSTNAME = *"mill"* ]]; then
	SDE="$SDE -knm -- "
else
	echo "Unsupported host"
	exit
fi

# ============================ FFB ============================================
source conf/ffb.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/profrun/ffb.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	X="`echo $BEST | cut -d '|' -f3`"
	Y="`echo $BEST | cut -d '|' -f4`"
	Z="`echo $BEST | cut -d '|' -f5`"
	INPUT="`echo $DEFINPUT | sed -e \"s/PX/$X/\" -e \"s/PY/$Y/\" -e \"s/PZ/$Z/\"`"
	# try finding closest cube size per proc
	FLOAT=`echo "e((1/3)*l($MAXDCZ / $NumMPI))" | bc -l`
	DCZ=`echo "($FLOAT+0.5)/1" | bc`
	INPUT="`echo $INPUT | sed -e \"s/DCZ/$DCZ/\"`"
	mkdir ./tmp; sleep 1; cd ./tmp
	mkdir -p ./oSDE
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE ../$BINARY $INPUT\"" >> $LOG 2>&1
	mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c "$SDE ../$BINARY $INPUT" >> $LOG 2>&1
	for P in `seq 0 $((NumMPI - 1))`; do
		echo "SDE output of MPI process $P" >> $LOG 2>&1
		cat ./oSDE/${P}.txt >> $LOG 2>&1
	done
	echo "=== SDE summary ===" >> $LOG 2>&1
	$ROOTDIR/util/analyze_sde.py ./oSDE `echo $LOG | sed 's#profrun#bestrun#g'` >> $LOG 2>&1
	rm -rf ./oSDE
	cd ../; rm -rf ./tmp; sleep 1
done
cd $ROOTDIR
