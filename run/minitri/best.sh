#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ miniTri ========================================
source conf/minitri.sh
NumRUNS=10
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
	for i in `seq 1 $NumRUNS`; do
		echo "Start at " `date --iso-8601=s` >> $LOG 2>&1
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		echo "Ended at " `date --iso-8601=s` >> $LOG 2>&1
	done
done
echo "Best miniTri run:"
grep 'compute.*B:' $LOG | awk -F 'B:' '{print $2}' > /dev/shm/1
grep 'compute ver' $LOG | awk -F 'vertex degrees:' '{print $2}' > /dev/shm/2
grep 'compute edg' $LOG | awk -F 'edge degrees:' '{print $2}' > /dev/shm/3
BEST="`paste /dev/shm/1 /dev/shm/2 /dev/shm/3 | awk '{print $1+$2+$3 "\t" $1}' | sort -g | head -1 | awk '{print $2}'`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR