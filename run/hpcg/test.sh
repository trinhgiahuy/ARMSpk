#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then
	MPIEXECOPT="-genv I_MPI_FABRICS=shm:ofi -genv FI_PROVIDER=sockets -genv I_MPI_HBW_POLICY=hbw_preferred -host `hostname` -genv I_MPI_ADJUST_ALLREDUCE=5 -genv KMP_AFFINITY=granularity=fine,compact,1,0"
else
	MPIEXECOPT="-genv I_MPI_FABRICS=shm:ofi -genv FI_PROVIDER=sockets -genv I_MPI_HBW_POLICY=hbw_preferred -host `hostname` -genv I_MPI_ADJUST_ALLREDUCE=5 -genv KMP_AFFINITY=compact"
fi

# ============================ HPCG ===========================================
source conf/hpcg.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/testrun/hpcg.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for TEST in $TESTCONF; do
	NumMPI="`echo $TEST | cut -d '|' -f1`"
	NumOMP="`echo $TEST | cut -d '|' -f2`"
	# test to identify hpcg's internal dimensions
	rm -f hpcg_log_* n*.yaml
	rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
	mpiexec $MPIEXECOPT -n $NumMPI $BINARY -n 1 > /dev/null 2>&1
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
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRunsTEST`; do
		START="`date +%s.%N`"
		timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
		ENDED="`date +%s.%N`"
		cat hpcg_log_* >> $LOG 2>&1
		cat n*.yaml >> $LOG 2>&1
		rm -f hpcg_log_* n*.yaml
		cat hpcg20*.txt >> $LOG 2>&1
		cat HPCG-Benchmark_3*.txt >> $LOG 2>&1
		rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
	done
done
echo "Best HPCG run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
