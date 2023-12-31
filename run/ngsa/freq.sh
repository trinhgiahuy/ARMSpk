#!/bin/bash
exit 1 #ignore in this study

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

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-genv I_MPI_FABRICS=shm:ofi -genv FI_PROVIDER=sockets -genv I_MPI_HBW_POLICY=hbw_preferred -host `hostname`"
export PATH=$ROOTDIR/dep/likwid/bin:$PATH
export LD_LIBRARY_PATH=$ROOTDIR/dep/likwid/lib:$LD_LIBRARY_PATH
if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then UNCORE="--umin 2.7 --umax 2.7"; fi

# ============================ NGSA ===========================================
source conf/ngsa.sh $ROOTDIR
LOG="$ROOTDIR/log/`hostname -s`/freqrun/ngsa.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for FREQ in $FREQR; do
	echo "likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE" >> $LOG 2>&1
	likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE
	likwid-setFrequencies -c 0 -p | grep 'CPU 0' >> $LOG 2>&1
	for BEST in $BESTCONF; do
		NumMPI="`echo $BEST | cut -d '|' -f1`"
		NumOMP="`echo $BEST | cut -d '|' -f2`"
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRunsBEST`; do
			# prep input (dep on numMPI; up to 12 supported)
			PreprocessInput $NumMPI $INPUTDIR
			START="`date +%s.%N`"
			timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="`date +%s.%N`"
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
			# clean up
			rm -rf workflow_*
			if [ -d $INPUTDIR/00-read-rank ]; then
				rm -rf $INPUTDIR/00-read-rank
			fi
		done
	done
done
cd $ROOTDIR
