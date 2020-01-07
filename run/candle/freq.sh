#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
export PATH=$ROOTDIR/dep/likwid/bin:$PATH
export LD_LIBRARY_PATH=$ROOTDIR/dep/likwid/lib:$LD_LIBRARY_PATH
if [[ $HOSTNAME = *"${XEONHOST}"* ]]; then UNCORE="--umin 2.7 --umax 2.7"; fi

# ============================ CANDLE =========================================
source conf/candle.sh
source activate idp
LOG="$ROOTDIR/log/`hostname -s`/freqrun/candle.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for FREQ in $FREQR; do
	echo "likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE" >> $LOG 2>&1
	likwid-setFrequencies -g performance --freq $FREQ --turbo 0 $UNCORE
	likwid-setFrequencies -c 0 -p | grep 'CPU 0' >> $LOG 2>&1
	for BEST in $BESTCONF; do
		for BINARY in $BINARYS; do
			NumMPI=1
			NumOMP=$BEST
			pushd "`find . -name $BINARY -exec dirname {} \;`"
			make libssc.so
			# check if data is hot or must be preloaded
			python ./p1b1.py
			#echo "mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT" >> $LOG 2>&1
			echo "OMP_NUM_THREADS=$NumOMP KMP_BLOCKTIME=30 KMP_SETTINGS=1 KMP_AFFINITY=\"granularity=fine,compact,1,0\" numactl --preferred 1 python $BINARY $INPUT" >> $LOG 2>&1
			for i in `seq 1 $NumRunsBEST`; do
				START="`date +%s.%N`"
				#timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT -genvall -genv OMP_NUM_THREADS=$NumOMP -n $NumMPI python $BINARY $INPUT >> $LOG 2>&1
				export OMP_NUM_THREADS=$NumOMP; export KMP_BLOCKTIME=30; export KMP_SETTINGS=1; export KMP_AFFINITY="granularity=fine,compact,1,0"
				timeout --kill-after=30s $MAXTIME numactl --preferred 1 python $BINARY $INPUT >> $LOG 2>&1
				if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
				ENDED="`date +%s.%N`"
				echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
			done
			popd
		done
	done
done
cd $ROOTDIR
