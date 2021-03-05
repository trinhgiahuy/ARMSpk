#!/bin/bash

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
if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	MPIEXECOPT+=" -x KMP_AFFINITY=granularity=fine,compact,1,0"
else
	MPIEXECOPT=" -x KMP_AFFINITY=compact"
fi

export PATH=$ROOTDIR/dep/sde-external-8.35.0-2019-03-11-lin:$PATH
if [ ! -x "`which sde64 2>/dev/null`" ]; then echo "ERROR: SDE missing, please download from Intel sde-external-8.35.0-2019-03-11-lin.tar.bz2 and untar in ./dep folder"; exit; fi;
SDE="`which sde64` -sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 -dcfg:out_base_name dcfg-out.rank-\"\$PMIX_RANK\" -align_checker_prefetch 0 -align_correct 0 -emu_fast 1 -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat"
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

# ============================ HPCG ===========================================
source conf/hpcg.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/profrun/hpcg.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF $SCALCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	TenRanks="`python3 -c \"from random import sample, seed; seed(0); need=10; print(0, ' '.join(['%s' % x for x in sorted(sample(range(1, $NumMPI), k=need-1))])) if $NumMPI > need else print(' '.join(['%s' % y for y in range($NumMPI)]))\"`"
	ProcElem=$(((NumCORES / NumMPI) + (NumCORES < NumMPI)))
	# test to identify hpcg's internal dimensions
	rm -f hpcg_log_* n*.yaml
	rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
	mpiexec $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=1 -n $NumMPI $BINARY -n 1 > /dev/null 2>&1
	if [ ! "x$?" = "x0" ]; then continue; fi
	if [ -f n*.yaml ]; then
		X=$(($MAXXYZ / `grep 'npx:' n*.yaml | awk -F 'npx:' '{print $2}'`))
		Y=$(($MAXXYZ / `grep 'npy:' n*.yaml | awk -F 'npy:' '{print $2}'`))
		Z=$(($MAXXYZ / `grep 'npz:' n*.yaml | awk -F 'npz:' '{print $2}'`))
	elif [ -f HPCG-Benchmark_3*.txt ]; then
		#non-intel version needs to be div8 https://github.com/hpcg-benchmark/hpcg/issues/47
		X=$((($MAXXYZ / `grep 'npx=' HPCG-Benchmark_3*.txt | awk -F 'npx=' '{print $2}'` / 8) * 8))
		Y=$((($MAXXYZ / `grep 'npy=' HPCG-Benchmark_3*.txt | awk -F 'npy=' '{print $2}'` / 8) * 8))
		Z=$((($MAXXYZ / `grep 'npz=' HPCG-Benchmark_3*.txt | awk -F 'npz=' '{print $2}'` / 8) * 8))
	else continue; fi
	rm -f hpcg_log_* n*.yaml
	rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
	INPUT="`echo $DEFINPUT | sed -e \"s/NX/$X/\" -e \"s/NY/$Y/\" -e \"s/NZ/$Z/\"`"
	if [ "x$RUNSDE" = "xyes" ]; then
		echo "=== sde run ===" >> $LOG 2>&1
		echo "mpiexec $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c \"$SDE $BINARY $INPUT\"" >> $LOG 2>&1
		for R in $TenRanks; do
			mpiexec $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=$NumOMP -n $NumMPI bash -c "if [[ \$PMIX_RANK = ${R} ]]; then $SDE $BINARY $INPUT; else $BINARY $INPUT; fi" >> $LOG 2>&1
			cat hpcg_log_* >> $LOG 2>&1
			cat n*.yaml >> $LOG 2>&1
			rm -f hpcg_log_* n*.yaml
			cat hpcg20*.txt >> $LOG 2>&1
			cat HPCG-Benchmark_3*.txt >> $LOG 2>&1
			rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
		done
		mkdir -p ${LOG}_${NumMPI}_${NumOMP}_sde; mv dcfg-out.* ${LOG}_${NumMPI}_${NumOMP}_sde/
	fi
	if [ "x$RUNPCM" = "xyes" ]; then
		# reset PMU counters
		pcm.x -r -- mpiexec $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=1 `lscpu | grep '^CPU(s):' | awk '{print $2}'` sleep 0.1 >> /dev/null 2>&1
		for PCM in $PCMB; do
			echo "=== intel $PCM run ===" >> $LOG 2>&1
			PCM+=" 360000 -- "
			echo "$PCM mpiexec $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
			$PCM mpiexec $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			rm -f hpcg_log_* n*.yaml
			rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
		done
	fi
	if [ "x$RUNVTUNE" = "xyes" ]; then
		for VTO in $VTAO; do
			echo "=== vtune $VTO ===" >> $LOG 2>&1
			echo "mpiexec -gtool \"amplxe-cl -collect $VTO $VTRO\" $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
			mpiexec -gtool "amplxe-cl -collect $VTO $VTRO" $MPIEXECOPT --map-by slot:pe=$ProcElem -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			amplxe-cl -report summary -q -result-dir ./oVTP.`hostname` >> $LOG 2>&1
			rm -rf ./oVTP.`hostname`
			rm -f hpcg_log_* n*.yaml
			rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
		done
	fi
done
cd $ROOTDIR
