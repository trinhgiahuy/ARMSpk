#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096

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
VTAO="hpc-performance memory-access"
VTRO="-data-limit=0 -finalization-mode=none -no-summary -trace-mpi -result-dir ./oVTP:all"

# ============================ BabelStream ====================================
source conf/babelstream.sh
DEFLOG="$ROOTDIR/log/`hostname -s`/profrun/babelstream"
mkdir -p `dirname $LOG`
cd $APPDIR
DEFINPUT=$INPUT
for BEST in $BESTCONF; do
	for BINARY in $BINARYS; do
		NumMPI=1
		NumOMP=$BEST
		S="`echo $BINARY | cut -d '_' -f2`"
		BINARY="`echo $BINARY | cut -d '_' -f1`"
		LOG="${DEFLOG}${S}gb.log"
		S=$((S*1024*1024*1024/8))
		INPUT="`echo $DEFINPUT | sed -e \"s/SIZE/$S/\"`"
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
done
cd $ROOTDIR