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

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ NGSA ===========================================
source conf/ngsa.sh $ROOTDIR
NumRUNS=10
LOG="$ROOTDIR/log/bestrun/ngsa.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		# prep input (dep on numMPI; up to 12 supported)
		PreprocessInput $NumMPI $INPUTDIR
		START="`date +%s.%N`"
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		ENDED="`date +%s.%N`"
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		# clean up
		rm -rf workflow_*
		if [ -d $INPUTDIR/00-read-rank ]; then
			rm -rf $INPUTDIR/00-read-rank
		fi
	done
done
echo "Best NGS Analyzer run:"
BEST="`grep 'Total running time:' $LOG | awk -F 'time:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR