#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

# ============================ miniTri ========================================
source conf/minitri.sh
LOG="$ROOTDIR/log/fugaku/bestrun/minitri.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI="`echo $BEST | cut -d '|' -f1`"
	NumOMP="`echo $BEST | cut -d '|' -f2`"
	if [ "x$NumMPI" != "x1" ]; then
		BINARY=$BINARYMPI
		INPUT=$INPUTMPI
	else
		BINARY=$BINARYOMP
		INPUT="`echo $INPUTOMP | sed -e \"s/OMPNT/$NumOMP/\"`"
	fi
	echo "$(get_mpi_cmd $NumMPI $NumOMP $LOG "") $BINARY $INPUT" >> $LOG 2>&1
	for i in `seq 1 $NumRunsBEST`; do
		START="`date +%s.%N`"
		timeout --kill-after=30s $MAXTIME $(get_mpi_cmd $NumMPI $NumOMP $LOG "") $BINARY $INPUT >> $LOG 2>&1
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding $MAXTIME timeout" >> $LOG 2>&1; fi
		ENDED="`date +%s.%N`"
		echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
	done
done
echo "Best miniTri run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpirun" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
