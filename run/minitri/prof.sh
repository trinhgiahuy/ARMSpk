#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

export PATH=$ROOTDIR/dep/sde-external-8.16.0-2018-01-30-lin:$PATH
SDE="`which sde64` -sse-sde -global_region -mix_omit_per_thread_stats -mix_omit_per_function_stats -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -iform 1 -omix oSDE/\"\$MPI_LOCALRANKID\".txt"
if [[ $HOSTNAME = *"kiev"* ]]; then
	SDE="$SDE -bdw -- "
elif [[ $HOSTNAME = *"lyon"* ]]; then
	SDE="$SDE -knl -- "
elif [[ $HOSTNAME = *"mill"* ]]; then
	SDE="$SDE -knm -- "
else
	echo "Unsupported host"
	exit
fi

# ============================ miniTri ========================================
source conf/minitri.sh
LOG="$ROOTDIR/log/bestrun/minitri.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	if [ "x$NumMPI" != "x1" ]; then
		BINARY=$BINARYMPI
		INPUT=$INPUTMPI
	else
		BINARY=$BINARYOMP
		INPUT="`echo $INPUTOMP | sed -e \"s/OMPNT/$NumOMP/\"`"
	fi
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRunsBEST`; do
		START="`date +%s.%N`"
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		ENDED="`date +%s.%N`"
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
	done
done
echo "Best miniTri run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
