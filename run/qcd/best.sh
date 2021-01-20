#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load openmpi@3.1.6%intel@19.0.1.144
MPIEXECOPT="--mca btl ^openib,tcp --oversubscribe --host `hostname`"
NumCORES=$((`lscpu | /bin/grep ^Socket | cut -d ':' -f2` * `lscpu | /bin/grep ^Core | cut -d ':' -f2`))

# ============================ CCS QCD ========================================
source conf/qcd.sh
LOG="$ROOTDIR/log/`hostname -s`/bestrun/qcd.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	PX="`echo $BEST | cut -d '|' -f3`"
	PY="`echo $BEST | cut -d '|' -f4`"
	PZ="`echo $BEST | cut -d '|' -f5`"
	BINARYYY="${BINARY}_${PX}${PY}${PZ}"
	echo "mpiexec $MPIEXECOPT --map-by slot:pe=$(((NumCORES / NumMPI) + (NumCORES < NumMPI))) -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARYYY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRunsBEST`; do
		START="`date +%s.%N`"
		timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT --map-by slot:pe=$(((NumCORES / NumMPI) + (NumCORES < NumMPI))) -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARYYY $INPUT >> $LOG 2>&1
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
		ENDED="`date +%s.%N`"
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
	done
done
echo "Best QCD run:"
BEST="`grep 'BiCGStab Total FLOPS:' $LOG | awk -F 'FLOPS:' '{print $2}' | sort -r -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
