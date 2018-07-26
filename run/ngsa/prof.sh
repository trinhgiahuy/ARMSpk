#!/bin/bash

function PreprocessInput {
	PROCS=$1
	INDIR=$2
	for P in `seq 0 $((PROCS - 1))`; do
		mkdir -p $INDIR/00-read-rank/${P}
	done
	for N in `seq 1 2`; do
		P=0; C=0
		for S in `seq -w 0 11`; do
			if [ "$C" = "$((12/PROCS))" ]; then
				P=$((P + 1))
				C=0
			fi
			cp $INDIR/00-read/part_${N}.${S} $INDIR/00-read-rank/${P}/
			C=$((C + 1))
		done
	done
}

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
#SDE="`which sde64` -sse-sde -global_region -mix_omit_per_thread_stats -mix_omit_per_function_stats -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -iform 1 -omix oSDE/\"\$MPI_LOCALRANKID\".txt"
SDE="`which sde64` -sse-sde -global_region -mix_omit_per_thread_stats -mix_omit_per_function_stats -iform 1 -omix oSDE/\"\$MPI_LOCALRANKID\".txt"
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
PCM="pcm-memory.x 360000 -- "

# ============================ NGSA ===========================================
source conf/ngsa.sh $ROOTDIR
LOG="$ROOTDIR/log/`hostname -s`/profrun/ngsa.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	# prep input (dep on numMPI; up to 12 supported)
	PreprocessInput $NumMPI $INPUTDIR
	if [ "x$RUNSDE" = "xyes" ]; then
		mkdir -p ./oSDE
		echo "=== sde run ===" >> $LOG 2>&1
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE $BINARY $INPUT\"" >> $LOG 2>&1
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c "$SDE $BINARY $INPUT" >> $LOG 2>&1
		# clean up
		rm -rf workflow_*
		if [ -d $INPUTDIR/00-read-rank ]; then rm -rf $INPUTDIR/00-read-rank; fi
		for P in `seq 0 $((NumMPI - 1))`; do
			echo "SDE output of MPI process $P" >> $LOG 2>&1
			cat ./oSDE/${P}.txt >> $LOG 2>&1
		done
		echo "=== SDE dump all ===" >> $LOG 2>&1
		for x in `ls ./oSDE/`; do
			echo "dump of ./oSDE/$x" >> $LOG 2>&1
			cat ./oSDE/$x >> $LOG 2>&1
		done
		echo "=== SDE summary ===" >> $LOG 2>&1
		$ROOTDIR/util/analyze_sde.py ./oSDE `echo $LOG | sed 's#profrun#bestrun#g'` >> $LOG 2>&1
		rm -rf ./oSDE
	fi
	if [ "x$RUNPCM" = "xyes" ]; then
		echo "=== intel pcm-memory.x run ===" >> $LOG 2>&1
		echo "$PCM mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		$PCM mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		rm -rf workflow_*
		if [ -d $INPUTDIR/00-read-rank ]; then rm -rf $INPUTDIR/00-read-rank; fi
	fi
	if [ "x$RUNVTUNE" = "xyes" ]; then
		echo "=== vtune hpc-performance ===" >> $LOG 2>&1
		echo "mpiexec -gtool "amplxe-cl -collect hpc-performance -data-limit=0 -no-auto-finalize -no-summary -trace-mpi -result-dir ./oVTP:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		mpiexec -gtool "amplxe-cl -collect hpc-performance -data-limit=0 -no-auto-finalize -no-summary -trace-mpi -result-dir ./oVTP:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT
		amplxe-cl -report summary -q -result-dir ./oVTP.`hostname` >> $LOG 2>&1
		rm -rf ./oVTP.`hostname`
		rm -rf workflow_*
		if [ -d $INPUTDIR/00-read-rank ]; then rm -rf $INPUTDIR/00-read-rank; fi
		echo "=== vtune memory-access ===" >> $LOG 2>&1
		echo "mpiexec -gtool "amplxe-cl -collect memory-access -data-limit=0 -no-auto-finalize -no-summary -trace-mpi -result-dir ./oVTM:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		mpiexec -gtool "amplxe-cl -collect memory-access -data-limit=0 -no-auto-finalize -no-summary -trace-mpi -result-dir ./oVTM:all" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT
		amplxe-cl -report summary -q -result-dir ./oVTM.`hostname` >> $LOG 2>&1
		rm -rf ./oVTM.`hostname`
		rm -rf workflow_*
		if [ -d $INPUTDIR/00-read-rank ]; then rm -rf $INPUTDIR/00-read-rank; fi
	fi
done
cd $ROOTDIR
