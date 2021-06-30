#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

# ============================ MACSio =========================================
source conf/macsio.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/bestrun/macsio.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	INPUT="`echo $DEFINPUT | sed -e \"s/NDPP/$(($MAXNDPP / $NumMPI))/\"`"
	mkdir -p ./bestrun; cd ./bestrun
	echo "$(get_mpi_cmd $NumMPI $NumOMP $LOG "") ../$BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRunsBEST`; do
		START="`date +%s.%N`"
		timeout --kill-after=30s $MAXTIME $(get_mpi_cmd $NumMPI $NumOMP $LOG "") ../$BINARY $INPUT >> $LOG 2>&1
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
		ENDED="`date +%s.%N`"
		grep 'Processor\|^Info' macsio-log.log >> $LOG 2>&1
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
	done
	cd ../; rm -rf ./bestrun
done
echo "Best MACSio run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
