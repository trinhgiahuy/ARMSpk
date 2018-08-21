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
export PATH=$ROOTDIR/dep/likwid/bin:$PATH
export LD_LIBRARY_PATH=$ROOTDIR/dep/likwid/lib:$LD_LIBRARY_PATH
if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then UNCORE="--umin 2.7 --umax 2.7"; fi

# ============================ HPCG ===========================================
source conf/hpcg.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/freqrun/hpcg.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for FREQ in $FREQR; do
	echo "likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE" >> $LOG 2>&1
	likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE
	likwid-setFrequencies -c 0 -p | grep 'CPU 0' >> $LOG 2>&1
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
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRunsBEST`; do
			START="`date +%s.%N`"
			timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="`date +%s.%N`"
			cat hpcg_log_* >> $LOG 2>&1
			cat n*.yaml >> $LOG 2>&1
			rm -f hpcg_log_* n*.yaml
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		done
	done
done
cd $ROOTDIR
