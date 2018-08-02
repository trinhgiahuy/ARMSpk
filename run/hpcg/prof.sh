#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	MPIEXECOPT="-host `hostname` -genv I_MPI_ADJUST_ALLREDUCE=5 -genv KMP_AFFINITY=granularity=fine,compact,1,0"
else
	MPIEXECOPT="-host `hostname` -genv I_MPI_ADJUST_ALLREDUCE=5 -genv KMP_AFFINITY=compact"
fi

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
PCMB="pcm.x pcm-memory.x pcm-power.x"

# ============================ HPCG ===========================================
source conf/hpcg.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/profrun/hpcg.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	# test to identify hpcg's internal dimensions
	rm -f hpcg_log_* n*.yaml
	mpiexec $MPIEXECOPT -n $NumMPI $BINARY -n 1 > /dev/null 2>&1
	X="`grep 'npx:' n*.yaml | awk -F 'npx:' '{print $2}'`"
	Y="`grep 'npy:' n*.yaml | awk -F 'npy:' '{print $2}'`"
	Z="`grep 'npz:' n*.yaml | awk -F 'npz:' '{print $2}'`"
	rm -f hpcg_log_* n*.yaml
	X=$(($MAXXYZ / $X))
	Y=$(($MAXXYZ / $Y))
	Z=$(($MAXXYZ / $Z))
	INPUT="`echo $DEFINPUT | sed -e \"s/NX/$X/\" -e \"s/NY/$Y/\" -e \"s/NZ/$Z/\"`"
	if [ "x$RUNSDE" = "xyes" ]; then
		mkdir -p ./oSDE
		echo "=== sde run ===" >> $LOG 2>&1
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE $BINARY $INPUT\"" >> $LOG 2>&1
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c "$SDE $BINARY $INPUT" >> $LOG 2>&1
		cat hpcg_log_* >> $LOG 2>&1
		cat n*.yaml >> $LOG 2>&1
		rm -f hpcg_log_* n*.yaml
		for P in `seq 0 $((NumMPI - 1))`; do
			echo "SDE output of MPI process $P" >> $LOG 2>&1
			cat ./oSDE/${P}.txt >> $LOG 2>&1
		done
		echo "=== SDE summary ===" >> $LOG 2>&1
		$ROOTDIR/util/analyze_sde.py ./oSDE `echo $LOG | sed -e 's#profrun#bestrun#g'` >> $LOG 2>&1
		rm -rf ./oSDE
	fi
	if [ "x$RUNPCM" = "xyes" ]; then
		# reset PMU counters
		pcm.x -r -- mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=1 `lscpu | grep '^CPU(s):' | awk '{print $2}'` sleep 0.1 >> /dev/null 2>&1
		for PCM in $PCMB; do
			echo "=== intel $PCM run ===" >> $LOG 2>&1
			PCM+=" 360000 -- "
			echo "$PCM mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
			$PCM mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			rm -f hpcg_log_* n*.yaml
		done
	fi
	if [ "x$RUNVTUNE" = "xyes" ]; then
		echo "=== vtune hpc-performance ===" >> $LOG 2>&1
		echo "mpiexec -gtool \"amplxe-cl -collect hpc-performance -data-limit=0 -finalization-mode=none -no-summary -trace-mpi -result-dir ./oVTP:all\" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		mpiexec -gtool "amplxe-cl -collect hpc-performance -data-limit=0 -finalization-mode=none -no-summary -trace-mpi -result-dir ./oVTP:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		amplxe-cl -report summary -q -result-dir ./oVTP.`hostname` >> $LOG 2>&1
		rm -rf ./oVTP.`hostname`
		echo "=== vtune memory-access ===" >> $LOG 2>&1
		echo "mpiexec -gtool \"amplxe-cl -collect memory-access -data-limit=0 -finalization-mode=none -no-summary -trace-mpi -result-dir ./oVTM:all\" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		mpiexec -gtool "amplxe-cl -collect memory-access -data-limit=0 -finalization-mode=none -no-summary -trace-mpi -result-dir ./oVTM:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		amplxe-cl -report summary -q -result-dir ./oVTM.`hostname` >> $LOG 2>&1
		rm -rf ./oVTM.`hostname`
		rm -f hpcg_log_* n*.yaml
	fi
done
cd $ROOTDIR
