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

# ============================ miniFE =========================================
source conf/minife.sh
LOG="$ROOTDIR/log/`hostname -s`/bestrun/minife.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	for BINARY in $BBINARY; do
		NumMPI="`echo $BEST | cut -d '|' -f1`"
		NumOMP="`echo $BEST | cut -d '|' -f2`"
		echo "mpiexec $MPIEXECOPT --map-by slot:pe=$(((NumCORES / NumMPI) + (NumCORES < NumMPI))) -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRunsBEST`; do
			START="`date +%s.%N`"
			timeout --kill-after=30s $MAXTIME mpiexec $MPIEXECOPT --map-by slot:pe=$(((NumCORES / NumMPI) + (NumCORES < NumMPI))) -x OMP_NUM_THREADS=$NumOMP -n $NumMPI $BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="`date +%s.%N`"
			cat ./miniFE.*.yaml >> $LOG 2>&1
			rm -f ./miniFE.*.yaml
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		done
	done
done
echo "Best miniFE run:"
BEST="`grep 'Total CG Mflops' $LOG | awk -F 'Mflops:' '{print $2}' | sort -r -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
