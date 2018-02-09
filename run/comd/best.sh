#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ CoMD ===========================================
source conf/comd.sh
NumRUNS=10
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/bestrun/comd.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	if [ "x${NumMPI}x" != "x1x" ]; then
		I="`echo $BEST | cut -d '|' -f3`"
		J="`echo $BEST | cut -d '|' -f4`"
		K="`echo $BEST | cut -d '|' -f5`"
		INPUT=$DEFINPUT
		INPUT="`echo $INPUT | sed -e \"s/i1/i$I/\"`"
		INPUT="`echo $INPUT | sed -e \"s/j1/j$J/\"`"
		INPUT="`echo $INPUT | sed -e \"s/k1/k$K/\"`"
	fi
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		echo "Start at " `date --iso-8601=s` >> $LOG 2>&1
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		echo "Ended at " `date --iso-8601=s` >> $LOG 2>&1
	done
done
rm *.yaml
echo "Best CoMD run:"
grep 'Starting simulation' $LOG | awk -F '=' '{print $2}' > /dev/shm/1
grep 'Ending simulation' $LOG | awk -F '=' '{print $2}' > /dev/shm/2
BEST="`paste /dev/shm/1 /dev/shm/2 | awk '{print $2-$1 "\t" $1}' | sort -g | head -1 | awk '{print $2}'`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
