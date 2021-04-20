#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load gcc@8.4.0
### all static, no need for compilers and envs -> XXX: not true for 2 intel version...

SPECCMD="runspec --config=nedo.cfg --nobuild --action=run --noreportable"

# ============================ SPEC OMP =======================================
source conf/spec_omp.sh
LOGDIR="$ROOTDIR/log/`hostname -s`/bestrun/spec_omp"
mkdir -p $LOGDIR
cd $APPDIR
for BENCH in $BINARY; do
	BM="`echo $BENCH | cut -d '|' -f1`"
	if [ -z $1 ]; then  COMP="`echo $BENCH | cut -d '|' -f2`"; else COMP=$1; fi
	SIZE="`echo $BENCH | cut -d '|' -f3`"
	LOG="$LOGDIR/${BM}.log"
	for NumOMP in $BESTCONF; do
		echo -e "=== runing $BENCH ===\nsource ./shrc; $SPECCMD --iterations=$NumRunsBEST --size=$SIZE --threads=$NumOMP --define COMP=$COMP --define RESDIR=0 $BM" >> $LOG 2>&1
		START="`date +%s.%N`"
		bash -c "source ./shrc; timeout --kill-after=30s $MAXTIME $SPECCMD --iterations=$NumRunsBEST --size=$SIZE --threads=$NumOMP --define COMP=$COMP --define RESDIR=0 $BM" >> $LOG 2>&1
		ENDED="`date +%s.%N`"
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
		REPORT="`find $ROOTDIR/$APPDIR/result -type f -name '*.log' | sort -g | tail -1`"
		cat $REPORT >> $LOG 2>&1
		/bin/grep 'Run Reported' $REPORT >> $LOG 2>&1
		echo "Best $BENCH run: " "`/bin/grep 'Reported' $REPORT | awk -F'Reported:[[:blank:]]+' '{print $2}' | tr -s ' ' | cut -d ' ' -f3 | sort -g | head -1`"
	done
done
cd $ROOTDIR
