#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-genv I_MPI_FABRICS=shm:ofi -genv FI_PROVIDER=sockets -genv I_MPI_HBW_POLICY=hbw_preferred -host `hostname`"

export PATH=$ROOTDIR/dep/sde-external-8.35.0-2019-03-11-lin:$PATH
if [ ! -x "`which sde64 2>/dev/null`" ]; then echo "ERROR: SDE missing, please download from Intel sde-external-8.35.0-2019-03-11-lin.tar.bz2 and untar in ./dep folder"; exit; fi;
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
VTAO="hpc-performance memory-access"
VTRO="-data-limit=0 -finalization-mode=none -no-summary -trace-mpi -result-dir ./oVTP:all"

# ============================ miniAMR ========================================
source conf/miniamr.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/profrun/miniamr.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	X="`echo $BEST | cut -d '|' -f3`"
	Y="`echo $BEST | cut -d '|' -f4`"
	Z="`echo $BEST | cut -d '|' -f5`"
	INPUT="`echo $DEFINPUT | sed -e \"s/PX/$X/\" -e \"s/PY/$Y/\" -e \"s/PZ/$Z/\"`"
	if [ "x$RUNSDE" = "xyes" ]; then
		mkdir -p ./oSDE
		echo "=== sde run ===" >> $LOG 2>&1
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE $BINARY $INPUT\"" >> $LOG 2>&1
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c "$SDE $BINARY $INPUT" >> $LOG 2>&1
		for P in `seq 0 $((NumMPI - 1))`; do
			echo "SDE output of MPI process $P" >> $LOG 2>&1
			cat ./oSDE/${P}.txt >> $LOG 2>&1
		done
		echo "=== SDE summary ===" >> $LOG 2>&1
		$ROOTDIR/util/analyze_sde.py ./oSDE `echo $LOG | sed 's#profrun#bestrun#g'` >> $LOG 2>&1
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
		done
	fi
	if [ "x$RUNVTUNE" = "xyes" ]; then
		for VTO in $VTAO; do
			echo "=== vtune $VTO ===" >> $LOG 2>&1
			echo "mpiexec -gtool \"amplxe-cl -collect $VTO $VTRO\" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
			mpiexec -gtool "amplxe-cl -collect $VTO $VTRO" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			amplxe-cl -report summary -q -result-dir ./oVTP.`hostname` >> $LOG 2>&1
			rm -rf ./oVTP.`hostname`
		done
	fi
done
cd $ROOTDIR
