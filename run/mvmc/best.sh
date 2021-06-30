#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"
if [ -n "${FUJIHOST}" ]; then module load Python2-CN; export FLIB_CNTL_BARRIER_ERR=FALSE; fi

# ============================ mVMC ===========================================
source conf/mvmc.sh
LOG="$ROOTDIR/log/`hostname -s`/bestrun/mvmc.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	# prepare input for strong scaling (scale down a bit from the default 64 node run)
	if [ -d ./job_mpi${NumMPI} ]; then rm -rf job_mpi${NumMPI}; fi
	sed -i -e 's/^Lx = Ly = 12/Lx = Ly = 4 #12/' -e 's/^NTotalSample = 4096/NTotalSample = 512 #4096/' -e 's/^NOuterMPI = 64/NOuterMPI = 2 #64/' ./makeDef/makeDef_large.py
	python2 ./makeDef/makeDef_large.py ${NumMPI}
	cd ./job_mpi${NumMPI}
	echo "$(get_mpi_cmd $NumMPI $NumOMP $LOG "") $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRunsBEST`; do
		START="`date +%s.%N`"
		timeout --kill-after=30s $MAXTIME $(get_mpi_cmd $NumMPI $NumOMP $LOG "") $BINARY $INPUT >> $LOG 2>&1
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
		ENDED="`date +%s.%N`"
		cat Lx*Ly*/zvo_HitachiTimer.dat >> $LOG 2>&1
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		rm -f Lx*Ly*/zvo_*
	done
	cd ../
done
echo "Best mVMC run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
