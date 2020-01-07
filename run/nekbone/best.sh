#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-genv I_MPI_FABRICS=shm:ofi -genv FI_PROVIDER=sockets -genv I_MPI_HBW_POLICY=hbw_preferred -host `hostname`"

# ============================ Nekbone ========================================
source conf/nekbone.sh
LOG="$ROOTDIR/log/`hostname -s`/bestrun/nekbone.log"
mkdir -p `dirname $LOG`
cd $APPDIR
if [ ! -f ./data.rea.bak ]; then cp ./data.rea ./data.rea.bak; fi
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	# prep input for strong scaling test
	NEPP=$(($ielN / $NumMPI))
	sed -e "s/1   50  1 = iel0/$NEPP  $NEPP  1 = iel0/" -e 's/8   10  2 = nx0/8    8  2 = nx0/' ./data.rea.bak > ./data.rea
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRunsBEST`; do
		START="`date +%s.%N`"
		timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
		ENDED="`date +%s.%N`"
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
	done
done
echo "Best Nekbone run:"
BEST="`grep '^Avg MFlops' $LOG | awk -F 'MFlops =' '{print $2}' | sort -g -r | head -1`"
grep "^Avg.*$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
