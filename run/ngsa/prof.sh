#!/bin/bash
exit 1 #ignore in this study

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
ulimit -n 8192; ulimit -u 8192
source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load openmpi@3.1.6%intel@19.0.1.144
MPIEXECOPT="--mca btl ^openib,tcp --oversubscribe --host `hostname`"
NumCORES=$((`lscpu | /bin/grep ^Socket | cut -d ':' -f2` * `lscpu | /bin/grep ^Core | cut -d ':' -f2`))

export PATH=$ROOTDIR/dep/sde-external-8.35.0-2019-03-11-lin:$PATH
if [ ! -x "`which sde64 2>/dev/null`" ]; then echo "ERROR: SDE missing, please download from Intel sde-external-8.35.0-2019-03-11-lin.tar.bz2 and untar in ./dep folder"; exit; fi;
SDE="`which sde64` -sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 -align_checker_prefetch 0 -align_correct 0 -emu_fast 1"
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

# ============================ NGSA ===========================================
source conf/ngsa.sh $ROOTDIR
LOG="$ROOTDIR/log/`hostname -s`/profrun/ngsa.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF $SCALCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	TenRanks="`python3 -c \"from random import sample, seed; seed(0); need=10; print(0, ' '.join(['%s' % x for x in sorted(sample(range(1, $NumMPI), k=need-1))])) if $NumMPI > need else print(' '.join(['%s' % y for y in range($NumMPI)]))\"`"
	ProcElem=$(((NumCORES / NumMPI) + (NumCORES < NumMPI)))
	if [ "x$RUNSDE" = "xyes" ]; then
		# prep input (dep on numMPI; up to 12 supported)
		echo "=== sde run ===" >> $LOG 2>&1
		echo "mpiexec $MPIEXECOPT --map-by slot:pe=$ProcElem -x SDE="$SDE" -x SDE_TEST_RANK=$R -x OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE $BINARY $INPUT\"" >> $LOG 2>&1
		for R in $TenRanks; do
			PreprocessInput $NumMPI $INPUTDIR
			mpiexec $MPIEXECOPT --map-by slot:pe=$ProcElem -x SDE="$SDE" -x SDE_TEST_RANK=$R -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			mkdir -p ${LOG}_${NumMPI}_${NumOMP}_sde; mv dcfg-out.* ${LOG}_${NumMPI}_${NumOMP}_sde/
			# clean up
			rm -rf workflow_*
			if [ -d $INPUTDIR/00-read-rank ]; then rm -rf $INPUTDIR/00-read-rank; fi
		done
	fi
	if [ "x$RUNPCM" = "xyes" ]; then
		# reset PMU counters
		pcm.x -r -- mpiexec $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=1 `lscpu | grep '^CPU(s):' | awk '{print $2}'` sleep 0.1 >> /dev/null 2>&1
		for PCM in $PCMB; do
			# prep input (dep on numMPI; up to 12 supported)
			PreprocessInput $NumMPI $INPUTDIR
			echo "=== intel $PCM run ===" >> $LOG 2>&1
			PCM+=" 360000 -- "
			echo "$PCM mpiexec $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
			$PCM mpiexec $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			rm -rf workflow_*
			if [ -d $INPUTDIR/00-read-rank ]; then rm -rf $INPUTDIR/00-read-rank; fi
		done
	fi
	if [ "x$RUNVTUNE" = "xyes" ]; then
		for VTO in $VTAO; do
			# prep input (dep on numMPI; up to 12 supported)
			PreprocessInput $NumMPI $INPUTDIR
			echo "=== vtune $VTO ===" >> $LOG 2>&1
			echo "mpiexec -gtool \"amplxe-cl -collect $VTO $VTRO\" $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
			mpiexec -gtool "amplxe-cl -collect $VTO $VTRO" $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			amplxe-cl -report summary -q -result-dir ./oVTP.`hostname` >> $LOG 2>&1
			rm -rf ./oVTP.`hostname`
			rm -rf workflow_*
			if [ -d $INPUTDIR/00-read-rank ]; then rm -rf $INPUTDIR/00-read-rank; fi
		done
	fi
done
cd $ROOTDIR
