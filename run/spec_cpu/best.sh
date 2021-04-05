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

SPECCMD="runcpu --config=nedo.cfg --nobuild --action=run --noreportable --use_submit_for_speed"

# ============================ SPEC CPU =======================================
source conf/spec_cpu.sh
LOGDIR="$ROOTDIR/log/`hostname -s`/bestrun/spec_cpu"
mkdir -p $LOGDIR
cd $APPDIR
for BENCH in $BINARY; do
	BM="`echo $BENCH | cut -d '|' -f1`"
	COMP="`echo $BENCH | cut -d '|' -f2`"
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
		echo "Best $BENCH run: " "`/bin/grep 'Run Reported' $REPORT | cut -d'(' -f2 | cut -d')' -f1 | tr -s ' ' | cut -d ' ' -f3 | sort -g | head -1`"
	done
done
cd $ROOTDIR
