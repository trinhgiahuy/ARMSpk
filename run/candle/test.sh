#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096

# ============================ CANDLE =========================================
source conf/candle.sh
LOG="$ROOTDIR/log/testrun/candle.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for TEST in $TESTCONF; do
	for BINARY in $BINARYS; do
		pushd "`find . -name $BINARY -exec dirname {} \;`"
		for i in `seq 1 $NumRunsTEST`; do
			python $BINARY >> $LOG 2>&1
		done
		popd
	done
done
echo "Best CANDLE run:"
BEST="`grep '^Walltime' $LOG | awk -F 'kernel:' '{print $2}' | sort -g | head -1`"
grep "$BEST\|mpiexec" $LOG | grep -B1 "$BEST"
echo ""
cd $ROOTDIR
