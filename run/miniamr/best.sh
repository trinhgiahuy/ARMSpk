#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"
#-genv KMP_PLACE_THREADS=1T -genv KMP_AFFINITY=compact -genv I_MPI_PIN_DOMAIN=node -genv I_MPI_PIN_ORDER=spread"

# ============================ miniAMR ========================================
source conf/miniamr.sh
NumRUNS=10
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/bestrun/miniamr.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	if [ "x${NumMPI}x" != "x1x" ]; then
		X="`echo $BEST | cut -d '|' -f3`"
		Y="`echo $BEST | cut -d '|' -f4`"
		Z="`echo $BEST | cut -d '|' -f5`"
		INPUT=$DEFINPUT
		INPUT="`echo $INPUT | sed -e \"s/npx 1/npx $X/\"`"
		INPUT="`echo $INPUT | sed -e \"s/npy 1/npy $Y/\"`"
		INPUT="`echo $INPUT | sed -e \"s/npz 1/npz $Z/\"`"
	fi
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		echo "Start at " `date --iso-8601=s` >> $LOG 2>&1
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		echo "Ended at " `date --iso-8601=s` >> $LOG 2>&1
	done
done
echo "Best miniAMR run:"
BEST="`grep 'total GFLOPS' $LOG | awk -F 'GFLOPS:' '{print $2}' | sort -r -n | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR