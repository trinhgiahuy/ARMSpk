#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname` -genv I_MPI_ADJUST_ALLREDUCE=5"

# ============================ HPCG ===========================================
source conf/hpcg.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/bestrun/hpcg.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	if [[ $HOSTNAME = *"kiev"* ]]; then
		MPIEXECOPT="$MPIEXECOPT -genv KMP_AFFINITY=granularity=fine,compact,1,0 -genv OMP_NUM_THREADS=$NumOMP"
	else
		MPIEXECOPT="$MPIEXECOPT –u OMP_NUM_THREADS –u KMP_AFFINITY -genv MIC_OMP_NUM_THREADS=$NumOMP"
	fi
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
	echo "mpiexec $MPIEXECOPT -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRunsBEST`; do
		START="`date +%s.%N`"
		mpiexec $MPIEXECOPT -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		ENDED="`date +%s.%N`"
		cat hpcg_log_* >> $LOG 2>&1
		cat n*.yaml >> $LOG 2>&1
		rm -f hpcg_log_* n*.yaml
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
	done
done
echo "Best HPCG run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
