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
PCMB="pcm.x pcm-memory.x pcm-power.x"
VTAO="hpc-performance memory-access"
VTRO="-data-limit=0 -finalization-mode=none -no-summary -trace-mpi -result-dir ./oVTP:all"

# ============================ mVMC ===========================================
source conf/mvmc.sh
LOG="$ROOTDIR/log/`hostname -s`/profrun/mvmc.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	# prepare input for strong scaling (scale down a bit from the default 64 node run)
	if [ -d ./job_mpi${NumMPI} ]; then rm -rf job_mpi${NumMPI}; fi
	sed -i -e 's/^Lx = Ly = 12/Lx = Ly = 4 #12/' -e 's/^NTotalSample = 4096/NTotalSample = 512 #4096/' -e 's/^NOuterMPI = 64/NOuterMPI = 2 #64/' ./makeDef/makeDef_large.py
	./makeDef/makeDef_large.py ${NumMPI}
	cd ./job_mpi${NumMPI}
	if [ "x$RUNSDE" = "xyes" ]; then
		mkdir -p ./oSDE
		echo "=== sde run ===" >> $LOG 2>&1
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE $BINARY $INPUT\"" >> $LOG 2>&1
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c "$SDE $BINARY $INPUT" >> $LOG 2>&1
		cat Lx*Ly*/zvo_HitachiTimer.dat >> $LOG 2>&1
		rm -f Lx*Ly*/zvo_*
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
			cat Lx*Ly*/zvo_HitachiTimer.dat >> $LOG 2>&1
			rm -f Lx*Ly*/zvo_*
		done
	fi
	if [ "x$RUNVTUNE" = "xyes" ]; then
		for VTO in $VTAO; do
			echo "=== vtune $VTO ===" >> $LOG 2>&1
			echo "mpiexec -gtool \"amplxe-cl -collect $VTO $VTRO\" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
			mpiexec -gtool "amplxe-cl -collect $VTO $VTRO" $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			amplxe-cl -report summary -q -result-dir ./oVTP.`hostname` >> $LOG 2>&1
			cat Lx*Ly*/zvo_HitachiTimer.dat >> $LOG 2>&1
			rm -f Lx*Ly*/zvo_*
			rm -rf ./oVTP.`hostname`
		done
	fi
	cd ../
done
cd $ROOTDIR
