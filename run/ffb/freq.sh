#!/bin/bash

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

# ============================ FFB ============================================
source conf/ffb.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/freqrun/ffb.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for FREQ in $FREQR; do
	echo "likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE" >> $LOG 2>&1
	likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE
	likwid-setFrequencies -c 0 -p | grep 'CPU 0' >> $LOG 2>&1
	for BEST in $BESTCONF; do
		NumMPI="`echo $BEST | cut -d '|' -f1`"
		NumOMP="`echo $BEST | cut -d '|' -f2`"
		X="`echo $BEST | cut -d '|' -f3`"
		Y="`echo $BEST | cut -d '|' -f4`"
		Z="`echo $BEST | cut -d '|' -f5`"
		INPUT="`echo $DEFINPUT | sed -e \"s/PX/$X/\" -e \"s/PY/$Y/\" -e \"s/PZ/$Z/\"`"
		# try finding closest cube size per proc
		FLOAT=`echo "e((1/3)*l($MAXDCZ / $NumMPI))" | bc -l`
		DCZ=`echo "($FLOAT+0.5)/1" | bc`
		INPUT="`echo $INPUT | sed -e \"s/DCZ/$DCZ/\"`"
		echo "mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRunsBEST`; do
			mkdir ./tmp; sleep 1; cd ./tmp
			START="`date +%s.%N`"
			timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI ../$BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="`date +%s.%N`"
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
			cd ../; rm -rf ./tmp; sleep 1
		done
	done
done
cd $ROOTDIR
