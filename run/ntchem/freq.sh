#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname`"
export PATH=$ROOTDIR/dep/likwid/bin:$PATH
export LD_LIBRARY_PATH=$ROOTDIR/dep/likwid/lib:$LD_LIBRARY_PATH
if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then UNCORE="--umin 2.7 --umax 2.7"; fi

# ============================ NTChem =========================================
source conf/ntchem.sh $ROOTDIR
LOG="$ROOTDIR/log/`hostname -s`/freqrun/ntchem.log"
mkdir -p `dirname $LOG`
cd $DATA_DIR
for FREQ in $FREQR; do
	echo "likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE" >> $LOG 2>&1
	likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE
	likwid-setFrequencies -c 0 -p | grep 'CPU 0' >> $LOG 2>&1
	for BEST in $BESTCONF; do
		NumMPI="`echo $BEST | cut -d '|' -f1`"
		NumOMP="`echo $BEST | cut -d '|' -f2`"
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $APPDIR/$BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRunsBEST`; do
			START="`date +%s.%N`"
			timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $APPDIR/$BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="`date +%s.%N`"
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		done
	done
done
cd $ROOTDIR
