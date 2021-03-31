#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
#source $ROOTDIR/conf/intel.cfg
#source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
#source $ROOTDIR/dep/spack/share/spack/setup-env.sh
#spack load gcc@8.4.0
### all static, no need for compilers and envs

if which numactl >/dev/null 2>&1; then PIN="numactl -l -C 1"; else PIN=""; fi

# ============================ PolyBench ====================================
source conf/polybench.sh
DEFLOG="$ROOTDIR/log/`hostname -s`/bestrun/polybench"
mkdir -p $DEFLOG
cd $APPDIR
for BEST in $BESTCONF; do
	for BMconf in $BINARYS; do
		NumMPI=1
		NumOMP=1
		BINARY="`echo ${BMconf} | cut -d '|' -f1`"
		BName="`basename $BINARY`"
		LOG="${DEFLOG}/${BName}/conf${1}.log"
		mkdir -p `dirname $LOG`
		echo "$PIN OMP_NUM_THREADS=$NumOMP $BINARY $INPUT" >> $LOG 2>&1
		for i in `seq 1 $NumRunsBEST`; do
			START="`date +%s.%N`"
			timeout --kill-after=30s $MAXTIME $PIN OMP_NUM_THREADS=$NumOMP $BINARY $INPUT >> $LOG 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
			ENDED="`date +%s.%N`"
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		done
		BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
		echo "Best $BName run: $BEST"
	done
done
cd $ROOTDIR
