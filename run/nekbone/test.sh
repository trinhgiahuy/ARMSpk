#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ Nekbone ========================================
source conf/nekbone.sh
NumRUNS=5
LOG="$ROOTDIR/log/testrun/nekbone.log"
mkdir -p `dirname $LOG`
cd $APPDIR
if [ ! -f ./data.rea.bak ]; then cp ./data.rea ./data.rea.bak; fi
for TEST in $TESTCONF; do
	NumMPI="`echo $TEST | cut -d '|' -f1`"
	NumOMP="`echo $TEST | cut -d '|' -f2`"
	# prep input for strong scaling test
	NEPP=$(($ielN / $NumMPI))
	sed -e "s/1   50  1 = iel0/$NEPP  $NEPP  1 = iel0/" -e 's/8   10  2 = nx0/8    8  2 = nx0/' ./data.rea.bak > ./data.rea
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
	done
done
echo "Best Nekbone run:"
BEST="`grep '^Avg MFlops' $LOG | awk -F 'MFlops =' '{print $2}' | sort -g -r | head -1`"
grep "^Avg.*$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
