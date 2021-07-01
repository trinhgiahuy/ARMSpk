#!/bin/bash
exit 1 #ignore in this study

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096

export PATH=$ROOTDIR/dep/sde-external-8.35.0-2019-03-11-lin:$PATH
if [ ! -x "`which sde64 2>/dev/null`" ]; then echo "ERROR: SDE missing, please sde-external-8.35.0-2019-03-11-lin.tar.bz2 and untar in ./dep folder"; exit; fi;
SDE="`which sde64` -sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 -dcfg:out_base_name dcfg-out.rank-\"\$MPI_LOCALRANKID\" -align_checker_prefetch 0 -align_correct 0 -emu_fast 1"
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

# ============================ CANDLE =========================================
source conf/candle.sh
source activate idp
LOG="$ROOTDIR/log/`hostname -s`/profrun/candle.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	for BINARY in $BINARYS; do
		NumMPI=1
		NumOMP=$BEST
		pushd "`find . -name $BINARY -exec dirname {} \;`"
		make libssc.so
		# check if data is hot or must be preloaded
		python ./p1b1.py
		if [ "x$RUNSDE" = "xyes" ]; then
			echo "=== sde run ===" >> $LOG 2>&1
			echo "mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE python $BINARY $INPUT\"": >> $LOG 2>&1
			#mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c "$SDE python $BINARY $INPUT" >> $LOG 2>&1
			export OMP_NUM_THREADS=$NumOMP; export KMP_BLOCKTIME=30; export KMP_SETTINGS=1; export KMP_AFFINITY="granularity=fine,compact,1,0"
			numactl --preferred 1 $SDE python $BINARY $INPUT >> $LOG 2>&1
			mkdir -p ${LOG}_${NumMPI}_${NumOMP}_sde; mv dcfg-out.* ${LOG}_${NumMPI}_${NumOMP}_sde/
		fi
		if [ "x$RUNPCM" = "xyes" ]; then
			# reset PMU counters
			pcm.x -r -- mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=1 `lscpu | grep '^CPU(s):' | awk '{print $2}'` sleep 0.1 >> /dev/null 2>&1
			for PCM in $PCMB; do
				echo "=== intel $PCM run ===" >> $LOG 2>&1
				PCM+=" 360000 -- "
				echo "$PCM mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT" >> $LOG 2>&1
				#$PCM mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT >> $LOG 2>&1
				export OMP_NUM_THREADS=$NumOMP; export KMP_BLOCKTIME=30; export KMP_SETTINGS=1; export KMP_AFFINITY="granularity=fine,compact,1,0"
				numactl --preferred 1 python $BINARY $INPUT >> $LOG 2>&1
			done
		fi
		if [ "x$RUNVTUNE" = "xyes" ]; then
			for VTO in $VTAO; do
				echo "=== vtune $VTO ===" >> $LOG 2>&1
				echo "mpiexec -gtool \"amplxe-cl -collect $VTO $VTRO\" $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT" >> $LOG 2>&1
				#mpiexec -gtool "amplxe-cl -collect $VTO $VTRO" $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT >> $LOG 2>&1
				numactl --preferred 1 amplxe-cl -collect $VTO $VTRO python $BINARY $INPUT >> $LOG 2>&1
				amplxe-cl -report summary -q -result-dir ./oVTP.`hostname` >> $LOG 2>&1
				rm -rf ./oVTP.`hostname`
			done
		fi
		popd
	done
done
cd $ROOTDIR
