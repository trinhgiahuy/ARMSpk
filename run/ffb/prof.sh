#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

export PATH=$ROOTDIR/dep/sde-external-8.16.0-2018-01-30-lin:$PATH
if [ ! -x "`which sde64 2>/dev/null`" ]; then echo "ERROR: SDE missing, please download from Intel sde-external-8.16.0-2018-01-30-lin.tar.bz2 and untar in ./dep folder"; exit; fi;
SDE="`which sde64` -sse-sde -global_region -mix_omit_per_thread_stats -mix_omit_per_function_stats -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -iform 1 -omix oSDE/\"\$MPI_LOCALRANKID\".txt"
if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	SDE="$SDE -bdw -- "
elif [[ $HOSTNAME = *"${IKNLHOST}"* ]]; then
	SDE="$SDE -knl -- "
elif [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
	SDE="$SDE -knm -- "
else
	echo "Unsupported host"
	exit
fi
export PATH=$ROOTDIR/dep/intel-pcm:$PATH
PCMM="pcm-memory.x 360000 -- "
PCMP="pcm-power.x 360000 -- "

# ============================ FFB ============================================
source conf/ffb.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/profrun/ffb.log"
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
	if [ "x$RUNSDE" = "xyes" ]; then
		mkdir -p ./oSDE
		echo "=== sde run ===" >> $LOG 2>&1
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE ../$BINARY $INPUT\"" >> $LOG 2>&1
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c "$SDE ../$BINARY $INPUT" >> $LOG 2>&1
		for P in `seq 0 $((NumMPI - 1))`; do
			echo "SDE output of MPI process $P" >> $LOG 2>&1
			cat ./oSDE/${P}.txt >> $LOG 2>&1
		done
		echo "=== SDE summary ===" >> $LOG 2>&1
		$ROOTDIR/util/analyze_sde.py ./oSDE `echo $LOG | sed 's#profrun#bestrun#g'` >> $LOG 2>&1
		rm -rf ./oSDE
	fi
	if [ "x$RUNPCM" = "xyes" ]; then
		echo "=== intel pcm-memory.x run ===" >> $LOG 2>&1
		echo "$PCMM mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT" >> $LOG 2>&1
		$PCMM mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT >> $LOG 2>&1
		echo "=== intel pcm-power.x run ===" >> $LOG 2>&1
		echo "$PCMP mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT" >> $LOG 2>&1
		$PCMP mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT >> $LOG 2>&1
	fi
	if [ "x$RUNVTUNE" = "xyes" ]; then
		echo "=== vtune hpc-performance ===" >> $LOG 2>&1
		echo "mpiexec -gtool "amplxe-cl -collect hpc-performance -data-limit=0 -no-auto-finalize -no-summary -trace-mpi -result-dir ./oVTP:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT" >> $LOG 2>&1
		mpiexec -gtool "amplxe-cl -collect hpc-performance -data-limit=0 -no-auto-finalize -no-summary -trace-mpi -result-dir ./oVTP:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT
		amplxe-cl -report summary -q -result-dir ./oVTP.`hostname` >> $LOG 2>&1
		rm -rf ./oVTP.`hostname`
		echo "=== vtune memory-access ===" >> $LOG 2>&1
		echo "mpiexec -gtool "amplxe-cl -collect memory-access -data-limit=0 -no-auto-finalize -no-summary -trace-mpi -result-dir ./oVTM:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT" >> $LOG 2>&1
		mpiexec -gtool "amplxe-cl -collect memory-access -data-limit=0 -no-auto-finalize -no-summary -trace-mpi -result-dir ./oVTM:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT
		amplxe-cl -report summary -q -result-dir ./oVTM.`hostname` >> $LOG 2>&1
		rm -rf ./oVTM.`hostname`
	fi
	cd ../; rm -rf ./tmp; sleep 1
done
cd $ROOTDIR
