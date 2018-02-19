#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"

# ============================ mVMC ===========================================
source conf/mvmc.sh
NumRUNS=10
LOG="$ROOTDIR/log/bestrun/mvmc.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	# prepare input for strong scaling (scale down a bit from the default 64 node run)
	if [ -d ./job_mpi${NumMPI} ]; then rm -rf job_mpi${NumMPI}; fi
	sed -i -e 's/^Lx = Ly = 12/Lx = Ly = 4 #12/' -e 's/^NTotalSample = 4096/NTotalSample = 512 #4096/' -e 's/^NOuterMPI = 64/NOuterMPI = 2 #64/' ./makeDef/makeDef_large.py
	./makeDef/makeDef_large.py ${NumMPI}
	cd ./job_mpi${NumMPI}
	echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRUNS`; do
		echo "Start at " `date --iso-8601=s` >> $LOG 2>&1
		mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
		cat Lx*Ly*/zvo_HitachiTimer.dat >> $LOG 2>&1
		rm -f Lx*Ly*/zvo_*
		echo "Ended at " `date --iso-8601=s` >> $LOG 2>&1
	done
	cd ../
done
echo "Best mVMC run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
